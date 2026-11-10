/* Symmetric thickness mode of a planar gas film between two deep liquids. */
#include "grid/multigrid1D.h"
#include <stdint.h>
#include <string.h>

scalar gap[], pressure[];

static const double gamma_film = 1., mu_gas = .02, mu_liquid = 1.;
static double pi_value, dx, attraction, temperature, step, neighbor_kernel;
static int cells, radix;
static double rotation_cos[4096], rotation_sin[4096];
static uint64_t rng;

/* Radix-two transform adapted from bm_fft1 in SLE_TCW. */
static void fft_radix2 (double * re, double * im, int n, int inverse)
{
  for (int i = 1, j = 0; i < n; i++) {
    int bit = n >> 1;
    for (; j & bit; bit >>= 1)
      j ^= bit;
    j ^= bit;
    if (i < j) {
      double tmp = re[i]; re[i] = re[j]; re[j] = tmp;
      tmp = im[i]; im[i] = im[j]; im[j] = tmp;
    }
  }
  for (int len = 2; len <= n; len *= 2) {
    double angle = (inverse ? 2. : -2.)*pi_value/len;
    double cr = cos(angle), ci = sin(angle);
    for (int i = 0; i < n; i += len) {
      double wr = 1., wi = 0.;
      for (int j = 0; j < len/2; j++) {
        int a = i + j, b = a + len/2;
        double vr = re[b]*wr - im[b]*wi, vi = re[b]*wi + im[b]*wr;
        re[b] = re[a] - vr; im[b] = im[a] - vi;
        re[a] += vr; im[a] += vi;
        double next = wr*cr - wi*ci;
        wi = wi*cr + wr*ci; wr = next;
      }
    }
  }
  if (inverse)
    for (int i = 0; i < n; i++) {
      re[i] /= n;
      im[i] /= n;
    }
}

static void fft (double * re, double * im, int inverse)
{
  if (radix == 1) {
    fft_radix2(re, im, cells, inverse);
    return;
  }
  // Split a 3*2^m transform into three radix-two transforms of stride three.
  int n = cells/3;
  double r[cells], s[cells];
  for (int phase = 0; phase < 3; phase++) {
    for (int j = 0; j < n; j++) {
      r[phase*n + j] = re[3*j + phase];
      s[phase*n + j] = im[3*j + phase];
    }
    fft_radix2(r + phase*n, s + phase*n, n, inverse);
  }
  double sign = inverse ? 1. : -1.;
  double scale = inverse ? 1./3. : 1.;
  for (int j = 0; j < cells; j++) {
    int a = j % n, b = a + n, c = b + n, k = 2*j % cells;
    double w1 = sign*rotation_sin[j], w2 = sign*rotation_sin[k];
    re[j] = scale*(r[a] + r[b]*rotation_cos[j] - s[b]*w1 +
                   r[c]*rotation_cos[k] - s[c]*w2);
    im[j] = scale*(s[a] + r[b]*w1 + s[b]*rotation_cos[j] +
                   r[c]*w2 + s[c]*rotation_cos[k]);
  }
}

static double uniform (void)
{
  uint64_t bits = (rng += UINT64_C(0x9e3779b97f4a7c15));
  bits = (bits ^ (bits >> 30))*UINT64_C(0xbf58476d1ce4e5b9);
  bits = (bits ^ (bits >> 27))*UINT64_C(0x94d049bb133111eb);
  bits ^= bits >> 31;
  return ((bits >> 12) + .5)*0x1p-52;
}

static double normal (void)
{
  double radius = sqrt(-2.*log(uniform()));
  return radius*cos(2.*pi_value*uniform());
}

static int index_of (double coordinate)
{
  return ((int) floor(coordinate/dx + 1e-8) + cells) % cells;
}

static double wave_number (int j)
{
  return 2.*pi_value*(j <= cells/2 ? j : cells - j)/L0;
}

static void apply_symbol (double * values, double * imag, const double * symbol)
{
  memset(imag, 0, cells*sizeof(double));
  fft(values, imag, 0);
  for (int j = 0; j < cells; j++) {
    values[j] *= symbol[j];
    imag[j] *= symbol[j];
  }
  fft(values, imag, 1);
}

int main (int argc, char ** argv)
{
  if (argc < 11 || argc > 13 || (argc == 13 && strcmp(argv[12], "stdin"))) {
    fprintf(stderr, "usage: gas_film N L mode amplitude A theta dt end stride seed [output_modes [stdin]]\n");
    return 2;
  }
  cells = atoi(argv[1]);
  radix = cells % 3 == 0 ? 3 : 1;
  int block = cells/radix;
  double length = atof(argv[2]);
  int mode = atoi(argv[3]);
  double amplitude = atof(argv[4]);
  attraction = atof(argv[5]);
  temperature = atof(argv[6]);  // k_B T divided by the physical transverse width.
  step = atof(argv[7]);
  double end = atof(argv[8]);
  int stride = atoi(argv[9]);
  rng = strtoull(argv[10], NULL, 10);
  if (cells < 16 || cells > 4096 || (block & (block - 1)) ||
      mode < 1 || mode >= cells/2 || stride < 1 ||
      !(length > 0. && fabs(amplitude) < .2 && attraction >= 0. &&
        temperature >= 0. && step > 0. && end >= step && end/step <= 1e8))
    return 2;
  long steps = lround(end/step);
  if (fabs(steps*step - end) > 1e-9*end)
    return 2;

  pi_value = acos(-1.);
  size(length);
  dimensions(radix);
  periodic(right);
  init_grid(cells);
  if (N != cells)
    return 2;
  dx = length/cells;
  double * work = calloc(17*cells, sizeof(double));
  if (!work)
    return 2;
  double * surface = work, * root = surface + cells, * rate = root + cells;
  double * propagator = rate + cells, * weight = propagator + cells;
  double * random_weight = weight + cells, * height = random_weight + cells;
  double * height_imag = height + cells, * rhs = height_imag + cells;
  double * rhs_imag = rhs + cells, * stress = rhs_imag + cells;
  double * shear_noise = stress + cells, * imag = shear_noise + cells;
  double * kick = imag + cells, * kick_imag = kick + cells;
  double * flow = kick_imag + cells, * random_flow = flow + cells;

  for (int j = 0; j < cells; j++) {
    rotation_cos[j] = cos(2.*pi_value*j/cells);
    rotation_sin[j] = sin(2.*pi_value*j/cells);
    double q = wave_number(j), qh = 2.*sin(q*dx/2.)/dx;
    if (j)
      surface[j] = 1./(4.*mu_liquid*q);
    root[j] = sqrt(surface[j]);
    neighbor_kernel += surface[j]*cos(2.*pi_value*j/cells)/cells;
    /* Exact propagation of the reference capillary operator only. */
    rate[j] = -.5*gamma_film*sq(sq(qh))*(1./(12.*mu_gas) + surface[j]);
    double z = rate[j]*step;
    propagator[j] = exp(z);
    weight[j] = z == 0. ? step : expm1(z)/rate[j];
    random_weight[j] = z == 0. ? 1. : sqrt(expm1(2.*z)/(2.*z));
  }

  if (argc == 13)
    for (int j = 0; j < cells; j++)
      if (scanf("%lf", &height[j]) != 1) {
        fprintf(stderr, "expected N initial thicknesses on stdin\n");
        free(work);
        return 2;
      }
  foreach()
    gap[] = argc == 13 ? height[index_of(x)] : 1. + amplitude*cos(2.*pi_value*mode*x/L0);
  boundary({gap});
  printf("t,mean,hmin,hmax");
  int observed = argc >= 12 ? atoi(argv[11]) : (cells/2 - 1 < 12 ? cells/2 - 1 : 12);
  if (observed < 1 || observed >= cells/2) {
    fprintf(stderr, "output_modes must lie between 1 and N/2 - 1\n");
    free(work);
    return 2;
  }
  for (int j = 1; j <= observed; j++)
    printf(",cos%d,sin%d", j, j);
  putchar('\n');

  for (long iteration = 0; iteration <= steps; iteration++) {
    double hmin = HUGE_VAL, hmax = 0.;
    int invalid = 0;
    foreach() {
      int j = index_of(x);
      height[j] = gap[] - 1.;
      height_imag[j] = 0.;
      invalid |= !isfinite(gap[]) || gap[] <= .4 || gap[] >= 1.6;
      hmin = fmin(hmin, gap[]);
      hmax = fmax(hmax, gap[]);
    }
    if (invalid) {
      fprintf(stderr, "left the pre-contact validation interval at t=%.17g\n", iteration*step);
      free(work);
      return 3;
    }
    fft(height, height_imag, 0);
    if (iteration % stride == 0 || iteration == steps) {
      printf("%.17g,%.17g,%.17g,%.17g", iteration*step, 1. + height[0]/cells, hmin, hmax);
      for (int j = 1; j <= observed; j++) {
        double angle = pi_value*j/cells;
        printf(",%.17g,%.17g", 2.*(height[j]*cos(angle) + height_imag[j]*sin(angle))/cells,
               2.*(height[j]*sin(angle) - height_imag[j]*cos(angle))/cells);
      }
      putchar('\n');
    }
    if (iteration == steps)
      break;

    foreach()
      pressure[] = -.5*gamma_film*(gap[1] - 2.*gap[] + gap[-1])/sq(Delta)
                   + attraction/cube(gap[]);
    boundary({pressure});
    foreach_face(x) {
      int j = index_of(x);
      double hface = .5*(gap[] + gap[-1]);
      kick[j] = hface;
      stress[j] = hface*(pressure[] - pressure[-1])/Delta;
    }
    apply_symbol(stress, imag, surface);
    if (temperature > 0.) {
      for (int j = 0; j < cells; j++)
        shear_noise[j] = normal();
      apply_symbol(shear_noise, imag, root);
    }
    foreach_face(x) {
      int j = index_of(x);
      double hface = .5*(gap[] + gap[-1]);
      double local = cube(hface)/(12.*mu_gas);
      flow[j] = -local*(pressure[] - pressure[-1])/Delta - hface*stress[j];
      if (temperature > 0.) {
        /* Ito drift of G^T diag(h_face) S diag(h_face) G on the periodic mesh. */
        double neighbors = kick[(j + cells - 1) % cells] - kick[(j + 1) % cells];
        flow[j] += temperature/dx*neighbor_kernel*neighbors/(2.*Delta);
        random_flow[j] = sqrt(local)*normal() + hface*shear_noise[j];
      }
    }
    foreach() {
      int j = index_of(x);
      // Use the same periodic face sample for both cells sharing the face.
      int right_face = (j + 1) % cells;
      rhs[j] = (flow[j] - flow[right_face])/Delta;
      rhs_imag[j] = 0.;
      if (temperature > 0.) {
        kick[j] = sqrt(2.*temperature*step/dx)*(random_flow[j] - random_flow[right_face])/Delta;
        kick_imag[j] = 0.;
      }
    }
    fft(rhs, rhs_imag, 0);
    if (temperature > 0.)
      fft(kick, kick_imag, 0);
    if (fabs(rhs[0]) > 1e-10 || (temperature > 0. && fabs(kick[0]) > 1e-10)) {
      fprintf(stderr, "periodic divergence produced a nonzero mean at t=%.17g\n", iteration*step);
      free(work);
      return 3;
    }
    for (int j = 1; j < cells; j++) {
      height[j] = propagator[j]*height[j] + weight[j]*(rhs[j] - rate[j]*height[j]);
      height_imag[j] = propagator[j]*height_imag[j] + weight[j]*(rhs_imag[j] - rate[j]*height_imag[j]);
      if (temperature > 0.) {
        height[j] += random_weight[j]*kick[j];
        height_imag[j] += random_weight[j]*kick_imag[j];
      }
    }
    fft(height, height_imag, 1);
    foreach()
      gap[] = 1. + height[index_of(x)];
    boundary({gap});
  }
  free(work);
  return 0;
}
