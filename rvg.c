#include "rvg.h"
#include <math.h>
#include <stdlib.h>

enum rvg_error rvg_unif(int (*rng)(), double *points, int n, int a, int b)
{
    for (int i = 0; i < n; i++)
        points[i] = (double) rng() / RAND_MAX;
   return RVG_ERR_OK; 
}

enum rvg_error rvg_exp(int (*rng)(), double *points, int n, double lambda)
{
    for (int i = 0; i < n; i++) {
        u = rng() / RAND_MAX;
        points[i] = -(1.0 / lambda) * log(u);
    }
    return RVG_ERR_OK;
}

enum rvg_error rvg_weibull(int (*rng)(), double *points, int n, double shape, double scale)
{
    for (int i = 0; i < n; i++) {
        u = rng() / RAND_MAX;
        points[i] = -scale * pow(log(u), 1.0 / shape;
    }
    return RVG_ERR_OK;
}

static double rvg_std_normal(int (*rng)())
{
    c0 = 2.515517;
    c1 = 0.802853;
    c2 = 0.010328;
    d1 = 1.432788;
    d2 = 0.189269;
    d3 = 0.001308;
    u = rng() / RAND_MAX
}

int main(int argc, char* argv[])
{
    return 0;
}
	
