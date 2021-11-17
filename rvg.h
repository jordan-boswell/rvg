#ifdef RVG_H
#define RVG_H

int rvg_unif(int (*rng)(), double *points, int n, int a, int b)
{
    for (int i = 0; i < n; i++)
        points[i] = (double) rng() / RAND_MAX;
   return 0; 
}

#endif /* RVG_H */
