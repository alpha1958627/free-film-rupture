/* Radix-two transform adapted from bm_fft1 in SLE_TCW, with a column stride. */
static void fft_line (double * re, double * im, int stride, int inverse)
{
  for (int i = 1, j = 0; i < n; i++) {
    int bit = n >> 1;
    for (; j & bit; bit >>= 1)
      j ^= bit;
    j ^= bit;
    if (i < j) {
      int a = i*stride, b = j*stride;
      double tmp = re[a]; re[a] = re[b]; re[b] = tmp;
      tmp = im[a]; im[a] = im[b]; im[b] = tmp;
    }
  }
  for (int len = 2; len <= n; len *= 2) {
    double angle = (inverse ? 2. : -2.)*pi_value/len;
    double cr = cos(angle), ci = sin(angle);
    for (int i = 0; i < n; i += len) {
      double wr = 1., wi = 0.;
      for (int j = 0; j < len/2; j++) {
        int a = (i + j)*stride, b = (i + j + len/2)*stride;
        double vr = re[b]*wr - im[b]*wi, vi = re[b]*wi + im[b]*wr;
        re[b] = re[a] - vr; im[b] = im[a] - vi;
        re[a] += vr; im[a] += vi;
        double next = wr*cr - wi*ci;
        wi = wi*cr + wr*ci; wr = next;
      }
    }
  }
}

static void fft2 (double * re, double * im, int inverse)
{
  for (int j = 0; j < n; j++)
    fft_line(re + j*n, im + j*n, 1, inverse);
  for (int j = 0; j < n; j++)
    fft_line(re + j, im + j, n, inverse);
  if (inverse)
    for (int j = 0; j < n*n; j++) {
      re[j] /= n*n;
      im[j] /= n*n;
    }
}
