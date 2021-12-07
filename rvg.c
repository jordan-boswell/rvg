#include "rvg.h"
#include <math.h>
#include <stdlib.h>


enum rvg_error rvg_exp(int (*rng)(), double *points, int n, double rate)
{
    double u;
    for (int i = 0; i < n; i++) {
        u = (double) rng() / RAND_MAX;
        points[i] = -(1.0 / rate) * log(u);
    }
    return RVG_ERR_OK;
}


enum rvg_error rvg_normal(int (*rng)(), double *points, int n, double mean, double sd)
{
    const double c0 = 2.515517;
    const double c1 = 0.802853;
    const double c2 = 0.010328;
    const double d1 = 1.432788;
    const double d2 = 0.189269;
    const double d3 = 0.001308;
    double t;
    double sign;
    double u;
    double z;
    for (int i = 0; i < n; i++) {
        u = (double) rng() / RAND_MAX
        if (u < 1 - u)
            t = sqrt(-log(u * u))
        else
            t = sqrt(-log((1 - u) * (1 - u)))
        if (u - 0.5 < 0)
            sign = -1;
        else if (u - 0.5 > 0)
            sign = 1;
        else
            sign = 0;
        z = sign * (t - (c0 + c1 * t + c2 * t * t) / (1 + d1 * t + d2 * t * t + d3 * t * t * t));
        points[i] = mean + sd * z;
    }
    return RVG_ERR_OK;
}


enum rvg_error rvg_triangular(int (*rng)(), double *points, double a, double b, double c)
{
    double u;
    for (int i = 0; i < n; i++) {
        u = (double) rng() / RAND_MAX;
        if (u <= (c - a) / (b - a))
            points[i] = a + sqrt((b - a) * (c - a) * u);
        else
            points[i] = b - sqrt((b - a) * (b - c) * (1 - u));
    }
    return RVG_ERR_OK;
}


enum rvg_error rvg_unif(int (*rng)(), double *points, int n, int a, int b)
{
    double u;
    for (int i = 0; i < n; i++) {
        u = (double) rng() / RAND_MAX;
        points[i] = a + (b - a) * u);
    }
   return RVG_ERR_OK; 
}


enum rvg_error rvg_weibull(int (*rng)(), double *points, int n, double scale, double shape)
{
    double u;
    for (int i = 0; i < n; i++) {
        u = rng() / RAND_MAX;
        points[i] = scale * pow(-log(u), 1 / shape);
    }
    return RVG_ERR_OK;
}


int main(int argc, char* argv[])
{
    return 0;
}
	
