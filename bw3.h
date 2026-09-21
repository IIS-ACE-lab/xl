#ifndef _BW3_h_
#define _BW3_h_

#include "matrix_array.h"
#include "matrix_polynomial.h"

#include "matrix.h"
#include "sp_matrix.h"

#include "monomial.h"

class BW3_opt
{
    template <unsigned N, unsigned M, unsigned width_sys>
    static bool sol_check(gfv<N> &sol, matrix<width_sys, M> &sys)
    {
       boost::scoped_ptr<matrix<M, width_sys> > _sys(new matrix<M, width_sys>);

        _sys->set_zero();

        for (unsigned i = 0; i < width_sys; i++)
            for (unsigned j = 0; j < M; j++)
                _sys->L[i].set(j, sys.L[j][i]);


        monomial<N> mon;

        gfv<M> tmp_v, sum_v;

        sum_v.set_zero();

        for (unsigned i = 0; i < width_sys; i++)
        {
            tmp_v = _sys->L[ i ];

            for (unsigned j = 0; j < N; j++)
                for (unsigned d = 0; d < mon[j]; d++)
                    tmp_v *= sol[ j ];

            mon.step();

            sum_v += tmp_v;
        }

        if (sum_v.is_zero())
           return true;
        else
           return false;
    }

 public:

   template <unsigned WIDTH_M, unsigned m, unsigned n,
             unsigned num_var, unsigned nsol, class orig_sys,
             unsigned deg_poly, unsigned deg_ai>
   static double bw_3(matrix<num_var, nsol> &sol, matrix_polynomial<n, m+n, deg_poly> &min_poly,
         matrix<n, m> &z, matrix_array<m, n, deg_ai> &ai, orig_sys &sys)
   {
       ECHO("BW_3 OPT\n");
   
   #ifdef MPI
       // no work for other nodes
       if (mpi_rank > 0)
         return 0;
   #endif
   
       unsigned rowDeg[ m+n ];
   
       min_poly.getRowDeg(rowDeg);
   
   
       unsigned maxDeg = 0;
   
       for (unsigned i = 0; i < n; i++)
           if (min_poly.nom_deg[i] > maxDeg)
               maxDeg = min_poly.nom_deg[i];
   
       unsigned minDeg = maxDeg;
   
       for (unsigned i = 0; i < n; i++)
           if (rowDeg[i] < minDeg)
               minDeg = rowDeg[i];

       ECHO("\n");
       for (unsigned i = 0; i < n; i++)
         ECHO("deg %3i: %4i - %4i\n", i, rowDeg[i], min_poly.nom_deg[i]);
       ECHO("\n");
   
   
   ////////////////////////////////////////////
   
       matrix<m, n+m> tmp;
   
       matrix<m, n> Bz0;
       Bz0.set_zero();
   
       // copy top part of sparse z to Bz0 *transposed*
       for (unsigned i = 0; i < n; i++)
          for (unsigned j = 0; j < m; j++)
             Bz0.L[i].set(j, z.L[j][i]);
   
   
       ECHO("min_poly.deg: %i\n", min_poly.deg);
       ECHO("min deg: %i\n", minDeg);
       ECHO("max deg: %i\n", maxDeg);
   
   
       unsigned num_sol = 0;
   
       double t = -get_ms_time();
   
       ECHO_NL();
   
       for (unsigned deg = minDeg; deg <= maxDeg; deg++)
       {
          matrix_prod(tmp, min_poly[deg], Bz0, 0, n);
   
          for (unsigned i = 0; i < deg; i++)
             matrix_mad(tmp, min_poly[deg - 1 - i], ai[i]);
   
          for (unsigned i = 0; i < n; i++)
             if ((rowDeg[i] <= deg) and (deg <= min_poly.nom_deg[i]))
             {
                gfv<num_var> sol_vec; sol_vec.set_zero();
   
                if (tmp.L[i][0] != 0)
                {
                   gf gf_tmp = tmp.L[i][0].inv();
   
                   for (uint64_t j = num_var; j >= 1; j--)
                      sol_vec.set(num_var - j, tmp.L[i][j] * gf_tmp);
   
                   if (sol_check(sol_vec, sys))
                   {
                      for (uint64_t j = num_var; j >= 1; j--)
                         printf(GF_FMT " ", (tmp.L[i][j] * gf_tmp).v);
   
                      printf(" is sol\n");
   
   
                      sol.L[num_sol++] = sol_vec;
   
                      if (num_sol >= nsol)
                         goto end;
                   }
   
                }
             }
       }
   
   end:
       t += get_ms_time();
   
       if (num_sol == 0)
          printf("\n!!! NO sol !!!\n");
       else if (num_sol < nsol)
          printf("\ndid not find %i solutions, only got %i...\n", nsol, num_sol);
   
       ECHO_NL();
   
       ECHO("BW_3 time: %.2f ms\n", t);
   
       return 0;
   }
 
};

#endif // ifndef _BW3_h_

