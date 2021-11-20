#ifdef RVG_H
#define RVG_H

struct rvg_rng {
    int (*rand)();
    int max;
};

enum rvg_error {
    RVG_ERR_OK
};

enum rvg_error rvg_unif_c(int (*rng)(), double *points, int n, int a, int b);

#endif /* RVG_H */
