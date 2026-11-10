/* Frozen-mobility 2D gas-gap model with fully implicit midpoint pressure.
 * Binary input: one initial gap field, then one normal field per noisy step.
 * Binary output: gap fields at t=0, first step, and every stride thereafter.
 */
#include "grid/multigrid.h"
#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif

scalar gap[];
static int n;
static double pi_value, dx;
#include "gas_fft_2d.h"

static int cell_index (double xcoord, double ycoord)
{
  return (int) floor(xcoord/dx)*n + (int) floor(ycoord/dx);
}

int main (int argc, char **argv)
{
  if (argc != 9) {
    fprintf(stderr, "usage: gas_film_midpoint_2d N L a theta dt steps stride noisy_steps\n");
    return 2;
  }
  n = atoi(argv[1]);
  double length = atof(argv[2]), a = atof(argv[3]), theta = atof(argv[4]);
  double dt = atof(argv[5]);
  int steps = atoi(argv[6]), stride = atoi(argv[7]), noisy = atoi(argv[8]);
  if (n < 16 || n > 512 || (n & (n-1)) || steps < 1 || stride < 1 ||
      noisy < 0 || noisy > steps || !(length > 0. && a >= 0. && theta >= 0. && dt > 0.))
    return 2;
#ifdef _WIN32
  _setmode(_fileno(stdin), _O_BINARY);
  _setmode(_fileno(stdout), _O_BINARY);
#endif
  pi_value = acos(-1.);
  dx = length/n;
  size(length);
  periodic(right);
  periodic(top);
  init_grid(n);
  int count = n*n;
  double * work = calloc(10*count, sizeof(double));
  if (!work)
    return 2;
  double * old = work, * re = old + count, * im = re + count;
  double * base_re = im + count, * base_im = base_re + count;
  double * p_re = base_im + count, * p_im = p_re + count;
  double * gain = p_im + count, * pressure = gain + count, * noise = pressure + count;
  int status = 0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      int k = i*n + j;
      int mx = i <= n/2 ? i : i-n, my = j <= n/2 ? j : j-n;
      double q = 2.*pi_value*hypot(mx,my)/length;
      double qx = sin(pi_value*mx/n), qy = sin(pi_value*my/n);
      double qd2 = 4.*(qx*qx + qy*qy)/(dx*dx);
      double mobility = q > 0. ? qd2*(1./.24 + 1./(4.*q)) : 0.;
      double half = .5*dt*mobility*(3.*a - .5*qd2);
      double denominator = 1. - half;
      if (!(denominator > 0.)) {
        status = 2;
        goto done;
      }
      gain[k] = (1. + half)/denominator;
      pressure[k] = -dt*mobility/denominator;
      noise[k] = sqrt(2.*theta*mobility*dt)/(dx*denominator);
    }
  if (fread(old, sizeof(double), count, stdin) != (size_t) count) {
    status = 4;
    goto done;
  }
  foreach()
    gap[] = old[cell_index(x,y)];

  for (int step = 0; step <= steps; step++) {
    int invalid = 0;
    foreach() {
      int k = cell_index(x,y);
      double h = gap[];
      old[k] = h;
      re[k] = h - 1.;
      im[k] = 0.;
      invalid |= !isfinite(h) || fabs(h-1.) > .08;
    }
    if (invalid) {
      fprintf(stderr, "left the near-flat interval at step %d\n", step);
      status = 3;
      goto done;
    }
    if ((step % stride == 0 || step == 1 || step == steps) &&
        fwrite(old, sizeof(double), count, stdout) != (size_t) count) {
      status = 4;
      goto done;
    }
    if (step == steps)
      break;
    int thermal = step < noisy;
    if (thermal) {
      if (fread(p_re, sizeof(double), count, stdin) != (size_t) count) {
        status = 4;
        goto done;
      }
      memset(p_im, 0, count*sizeof(double));
      fft2(p_re, p_im, 0);
    }
    fft2(re, im, 0);
    for (int k = 0; k < count; k++) {
      base_re[k] = gain[k]*re[k] + (thermal ? noise[k]*p_re[k] : 0.);
      base_im[k] = gain[k]*im[k] + (thermal ? noise[k]*p_im[k] : 0.);
      re[k] = base_re[k];
      im[k] = base_im[k];
    }
    fft2(re, im, 1);
    if (a > 0.) {
      double change = HUGE_VAL;
      for (int iter = 0; iter < 16 && change > 2e-13; iter++) {
        for (int k = 0; k < count; k++) {
          double u = .5*(old[k] - 1. + re[k]);
          double h = 1. + u;
          p_re[k] = a*u*u*(6. + u*(8. + 3.*u))/(h*h*h);
          p_im[k] = 0.;
        }
        fft2(p_re, p_im, 0);
        for (int k = 0; k < count; k++) {
          p_re[k] = base_re[k] + pressure[k]*p_re[k];
          p_im[k] = base_im[k] + pressure[k]*p_im[k];
        }
        fft2(p_re, p_im, 1);
        change = 0.;
        for (int k = 0; k < count; k++) {
          change = fmax(change, fabs(p_re[k]-re[k]));
          re[k] = p_re[k];
        }
      }
      if (change > 2e-13) {
        fprintf(stderr, "midpoint residual did not converge at step %d\n", step);
        status = 3;
        goto done;
      }
    }
    foreach()
      gap[] = 1. + re[cell_index(x,y)];
  }
done:
  free(work);
  return status;
}
