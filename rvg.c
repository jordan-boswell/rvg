#include "rvg.h"
#include <math.h>
#include <stdlib.h>


#ifndef M_PI
#define M_PI 3.141592653589
#endif /* M_PI */


enum rvg_error rvg_bernoulli(int (*rng)(), int *points, int len, double p)
{
    int i;
    for (i = 0; i < len; i++) {
        double u = (double) rng() / RAND_MAX;
        if (u < p)
            points[i] = 1;
        else
            points[i] = 0;
    }
    return RVG_ERR_OK;
}


enum rvg_error rvg_binomial(int (*rng)(), int *points, int len, int n, double p)
{
    int i;
    for (i = 0; i < len; i++) {
        int count = 0;
        int j;
        for (j = 0; j < n; j++) {
            double u = (double) rng() / RAND_MAX;
            if (u < p)
                count++;
        }
        points[i] = count;
    }
    return RVG_ERR_OK;
}


enum rvg_error rvg_chi_squared(int (*rng)(), double *points, int len, int df)
{
    double u;
    int k = df / 2;
    double accum;
    if (df % 2 == 0) {
        int i;
        for (i = 0; i < len; i++) {
            int j;
            accum = 0;
            for (j = 0; j < k; j++) {
                u = (double) rng() / RAND_MAX;
                accum += -2 * log(u);
            }
            points[i] = accum;
        }
    } else {
        int i;
        for (i = 0; i < len; i++) {
            int j;
            accum = 0;
            for (j = 0; j < k; j++) {
                u = (double) rng() / RAND_MAX;
                accum += -2 * log(u);
            }
            points[i] = accum + log(rng()) * pow(sin(2 * M_PI * rng()), 2);
        }
    }
    return RVG_ERR_OK;
}


enum rvg_error rvg_exp(int (*rng)(), double *points, int len, double rate)
{
    int i;
    for (i = 0; i < len; i++) {
        double u = (double) rng() / RAND_MAX;
        points[i] = -(1.0 / rate) * log(u);
    }
    return RVG_ERR_OK;
}


enum rvg_error rvg_f(int (*rng)(), double *points, int len, int df1, int df2)
{
    double u;
    int k = df2 / 2;
    double accum;
    rvg_chi_squared(rng, points, len, df1);
    if (df2 % 2 == 0) {
        int i;
        for (i = 0; i < len; i++) {
            int j;
            accum = 0;
            for (j = 0; j < k; j++) {
                u = (double) rng() / RAND_MAX;
                accum += -2 * log(u);
            }
            points[i] *= df2 / (df1 * accum);
        }
    } else {
        int i;
        for (i = 0; i < len; i++) {
            int j;
            accum = 0;
            for (j = 0; j < k; j++) {
                u = (double) rng() / RAND_MAX;
                accum += -2 * log(u);
            }
            points[i] *= df2 / (df1 * (accum + log(rng()) * pow(sin(2 * M_PI * rng()), 2)));
        }
    }
    return RVG_ERR_OK;
}


enum rvg_error rvg_gamma(int (*rng)(), double *points, int len, double shape, double rate)
{
    double u1;
    double u2;
    if (shape < 1) {
        double b = (exp(1) + shape) / exp(1);
        int i;
        for (i = 0; i < len; i++) {
            while (1) {
                double y;
                double w;
                u1 = (double) rng() / RAND_MAX;
                w = b * u1;
                if (w < 1) {
                    y = pow(w, 1.0 / shape);
                    u2 = (double) rng() / RAND_MAX;
                    if (u2 <= exp(-y)) {
                        points[i] = y / rate;
                        break;
                    }
                } else {
                    y = -log((b - w) / shape); 
                    u2 = (double) rng() / RAND_MAX;
                    if (u2 <= pow(y, shape - 1)) {
                        points[i] = y / rate;
                        break;
                    }
                }
            }    
        }
    } else {
        double a = 1.0 / sqrt(2 * shape - 1);
        double b = shape - log(4);
        double c = shape + 1.0 / a;
        double d = 1 + log(4.5);
        int i;
        for (i = 0; i < len; i++) {
            while (1) {
                double v;
                double y;
                double z;
                double w;
                u1 = (double) rng() / RAND_MAX;
                u2 = (double) rng() / RAND_MAX;
                v = a * log(u1 / (1 - u2));
                y = shape * exp(v);
                z = u1 * u1 * u2;
                w = b + c * v - y;
                if ( w + d - 4.5 * z >= 0) {
                    points[i] = y / rate;
                    break;
                } else if (w >= log(z)) {
                    points[i] = y / rate;
                    break;
                }
            }
        }
    }
    return RVG_ERR_OK;
}


enum rvg_error rvg_geometric(int (*rng)(), int *points, int len, double p)
{
    int i;
    for (i = 0; i < len; i++) {
        double u = (double) rng() / RAND_MAX;
        points[i] = ceil(log(u) / log(1 - p));
    }
    return RVG_ERR_OK;
}


enum rvg_error rvg_negative_binomial(int (*rng)(), int *points, int len, int n, double p)
{
    int i;
    for (i = 0; i < len; i++) {
        double sum = 0.0;
        int j;
        for (j = 0; j < n; j++) {
            double u = (double) rng() / RAND_MAX;
            sum += ceil(log(u) / log(1 - p));
        }
        points[i] = sum;
    }
    return RVG_ERR_OK;
}


enum rvg_error rvg_normal(int (*rng)(), double *points, int len, double mean, double sd)
{
    const double c0 = 2.515517;
    const double c1 = 0.802853;
    const double c2 = 0.010328;
    const double d1 = 1.432788;
    const double d2 = 0.189269;
    const double d3 = 0.001308;
    double t;
    double sign;
    double z;
    int i;
    for (i = 0; i < len; i++) {
        double u = (double) rng() / RAND_MAX;
        if (u < 1 - u)
            t = sqrt(-log(u * u));
        else
            t = sqrt(-log((1 - u) * (1 - u)));
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


enum rvg_error rvg_poisson(int (*rng)(), int *points, int len, double rate)
{
    int i;
    for (i = 0; i < len; i++) {
        double a, p, x, u;
        a = exp(-rate);
        p = 1;
        x = -1;
        while (a <= p) {
            u = (double) rng() / RAND_MAX;
            p = p * u;
            x = x + 1;
        }
        points[i] = x;
    }
    return RVG_ERR_OK;
}


enum rvg_error rvg_t(int (*rng)(), double *points, int len, int df)
{
    double u;
    int k = df / 2;
    double accum;
    rvg_normal(rng, points, len, 0, 1);
    if (df % 2 == 0) {
        int i;
        for (i = 0; i < len; i++) {
            int j;
            accum = 0;
            for (j = 0; j < k; j++) {
                u = (double) rng() / RAND_MAX;
                accum += -2 * log(u);
            }
            points[i] /= sqrt(accum / df);
        }
    } else {
        int i;
        for (i = 0; i < len; i++) {
            int j;
            accum = 0;
            for (j = 0; j < k; j++) {
                u = (double) rng() / RAND_MAX;
                accum += -2 * log(u);
            }
            points[i] /= sqrt((accum + log(rng()) * pow(sin(2 * M_PI * rng()), 2)) / df);
        }
    }
    return RVG_ERR_OK;
}


enum rvg_error rvg_triangular(int (*rng)(), double *points, int len, double a, double b, double c)
{
    int i;
    for (i = 0; i < len; i++) {
        double u = (double) rng() / RAND_MAX;
        if (u <= (c - a) / (b - a))
            points[i] = a + sqrt((b - a) * (c - a) * u);
        else
            points[i] = b - sqrt((b - a) * (b - c) * (1 - u));
    }
    return RVG_ERR_OK;
}


enum rvg_error rvg_unif(int (*rng)(), double *points, int len, int a, int b)
{
    int i;
    for (i = 0; i < len; i++) {
        double u = (double) rng() / RAND_MAX;
        points[i] = a + (b - a) * u;
    }
   return RVG_ERR_OK; 
}


enum rvg_error rvg_weibull(int (*rng)(), double *points, int len, double scale, double shape)
{
    int i;
    for (i = 0; i < len; i++) {
        double u = rng() / RAND_MAX;
        points[i] = scale * pow(-log(u), 1 / shape);
    }
    return RVG_ERR_OK;
}