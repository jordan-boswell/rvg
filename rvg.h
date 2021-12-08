#ifdef RVG_H
#define RVG_H


enum rvg_error {
    RVG_ERR_OK
};

enum rvg_error rvg_bernoulli(int (*rng)(), int *points, int len, double p);

enum rvg_error rvg_binomial(int (*rng)(), int *points, int len, int n, double p);

enum rvg_error rvg_chi_squared(int (*rng)(), double *points, int len, int df);

enum rvg_error rvg_exp(int (*rng)(), double *points, int len, double rate);

enum rvg_error rvg_f(int (*rng)(), double *points, int len, int df1, int df2);

enum rvg_error rvg_gamma(int (*rng)(), double *points, int len, double shape, double rate);

enum rvg_error rvg_geometric(int (*rng)(), int *points, int len, double p);

enum rvg_error rvg_negative_binomial(int (*rng)(), int *points, int len, int n, double p);

enum rvg_error rvg_normal(int (*rng)(), double *points, int len, double mean, double sd);

enum rvg_error rvg_poisson(int (*rng)(), int *points, int len, double rate);

enum rvg_error rvg_t(int (*rng)(), double *points, int len, int df);

enum rvg_error rvg_triangular(int (*rng)(), double *points, int len, double a, double b, double c);

enum rvg_error rvg_unif(int (*rng)(), double *points, int len, int a, int b);

enum rvg_error rvg_weibull(int (*rng)(), double *points, int len, double scale, double shape);


#endif /* RVG_H */