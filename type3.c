/*
 *  author: 
 *    Christiane Peters, 
 *    c.p.peters@mat.dtu.dk
 *
 * 
 *  date: 12-sep-2008
 *
 * 
 *  acknowledgements:
 *    thanks to Dan Bernstein for 
 *    helpful comments,
 *    bug-fixing and in particular for recommending
 *    using the MPFI library in order to achieve
 *    high precision
 *
 * 
 *  type-3 Markov analysis as in:    
 *    Daniel J. Bernstein, Tanja Lange, Christiane Peters. 
 *    Attacking and defending the McEliece cryptosystem. 
 *    To appear in: Proceedings of PQCrypto 2008, 
 *    Lecture Notes in Computer Science, Springer, 2008. 
 *    
 *    http://www2.mat.dtu.dk/people/C.Peters/publications/2008-mceliece.pdf
 *
 *
 *  this program uses the MPFI library
 *    See http://perso.ens-lyon.fr/nathalie.revol/mpfi.html
 *    The MPFI library is built on top of the MPFR library,
 *    which is built on top of the GMP library.
 *
 *  compile: 
 *    gcc type3.c -o type3 -lm -lgmp -lmpfr -lmpfi
 * 
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "mpfi.h"
#include "mpfi_io.h"


int n = 1024;
int k = 525;
int x = 262;
int w = 50;

mpfi_t **P;
mpfi_t **Tmp;
mpfi_t **R;

mpfi_t beta;
mpfi_t tmp0, tmp1, tmp2, Anum, Bnum;
mpfi_t num;
mpfr_t numleft;
mpfr_t numright;


#define NMAX 100000

mpfi_t *factorial = 0;

void factorial_init(void)
{
  int i;
  
  factorial = (mpfi_t *) malloc((NMAX + 1) * sizeof(mpfi_t));
  
  if (!factorial) abort();
  
  for (i = 0; i <= NMAX; i++)
  {  
    mpfi_init(factorial[i]);
  }
  
  mpfi_set_ui(factorial[0], 1);

  for (i = 1; i <= NMAX; i++) 
  {
    mpfi_mul_ui(factorial[i], factorial[i - 1], i);
  }

}


void factorial_free(void)
{
  int i;
  
  for(i = 0; i <= NMAX; i++)
  {
    mpfi_clear(factorial[i]);
  }
  
  free(factorial);
}


void C(mpfi_t result, int a, int b)
{
  if (b < 0 || b > a)
  {
    mpfi_set_ui(result, 0);
    return;
  }
  if (a > NMAX) abort();
  if (!factorial) abort();
  
  mpfi_div(result, factorial[a], factorial[b]);
  mpfi_div(result, result, factorial[a - b]);
}


/* matrix ops */
void id_mat(mpfi_t **matrix, int length)
{
  int i, j;
  
  for(i = 0; i < length; i++)
  {
    for(j = 0; j < length; j++)
    {
      if (i==j)  mpfi_set_ui(matrix[i][j], 1);
      else    mpfi_set_ui(matrix[i][j], 0); 
    }
  }
}


void copy_mat(mpfi_t **matA, mpfi_t **matB, int length)
{
  int i, j;
  
  for(i = 0; i < length; i++)
  {
    for(j = 0; j < length; j++)
    {
      mpfi_set(matA[i][j], matB[i][j]);
    }
  }
}


void clear_mat(mpfi_t **matrix, int length)
{
  int i, j;

  for(i = 0; i < length; i++ )
  {
    for(j = 0; j < length; j++)
    {
      mpfi_set_ui(matrix[i][j], 0);
    }
  }
}


void print_mat(mpfi_t **matrix, int length)
{
  int i, j;
  
  for(i = 0; i < length; i = i+1)
  {
    for(j = 0; j < length; j = j+1)
    {
      mpfi_out_str(stdout, 10, 10, matrix[i][j]);
      printf("  ");
//      printf("%lf\t", mpfi_get_d(matrix[i][j]));
    }
    printf("\n");
  }

  printf("\n");
}


void prod_mat(mpfi_t **matC, mpfi_t **matA, mpfi_t **matB, int length)
{
  int i, j, k;

  for(i = 0; i < length; i++)
  {
    for(j = 0; j < length; j++)
    {
      mpfi_set_ui(matC[i][j], 0);
      
      for(k = 0; k < length; k++)
      {
        mpfi_mul(tmp0, matA[i][k], matB[k][j]);
        mpfi_add(matC[i][j], matC[i][j], tmp0);
      }
    }
  }
}


mpfi_t** init_mat(int length)
{
  int i, j;

  mpfi_t ** matrix = (mpfi_t **) malloc((length)*sizeof(mpfi_t *));
  
  if(matrix != NULL)
  {
    for(i = 0; i < length; i++)
    {
      matrix[i] = (mpfi_t*) malloc((length)* sizeof(mpfi_t));

      if(matrix[i] == NULL)
      {
        printf("%i allocation error\n", i);
        abort();
      }
    }
  }
  else
  {
    printf("allocation error\n");
    abort();
  }

  for(i = 0; i < length; i++)
  {

    for(j = 0; j < length; j++)
    {
      mpfi_init(matrix[i][j]);
    }
  }

  return matrix;
}


void free_mat(mpfi_t **matrix, int length)
{
   int i, j;

   for(i = 0; i < length; i++)
   {
     for(j = 0; j < length; j++)
     {
       mpfi_clear(matrix[i][j]);
     }
   }
  
   for(i = 0; i < length; i++)
     free(matrix[i]);
   free(matrix);
}


void swaprow_mat(mpfi_t **matA, int i, int k, int length)
{
  int j;

  for(j = 0; j < length; j++)
  {
     mpfi_set(tmp0, matA[i][j]);
     mpfi_set(matA[i][j], matA[k][j]);
     mpfi_set(matA[k][j], tmp0);
  }

}


void multline_mat(mpfi_t **matA, mpfi_t v, int j, int length)
{
  int k;

  for(k = 0; k < length; k++)
  {
    mpfi_mul(matA[j][k], matA[j][k], v);
  }
}


void addmultline_mat(mpfi_t **matA, mpfi_t v, int upper, int lower, int length)
{
  int k;
  /* tmp0, tmp1, tmp2 are in use */
  mpfi_t tmp3;
  mpfi_init(tmp3);

  for(k = 0; k < length; k++)
  {
    mpfi_mul(tmp3, matA[upper][k], v);
    mpfi_add(matA[lower][k], matA[lower][k], tmp3);
  }
      
  mpfi_clear(tmp3);
}


void gaussjord_mat(mpfi_t **matB, mpfi_t **matA, int length)
{
  int i, j, r_val, j_max;
  mpfi_t max_piv;
  mpfi_init(max_piv);
  
  id_mat(matB, length);

  for(i = 0; i < length; i++)
  {  
    mpfi_set(max_piv, matA[i][i]);
    j_max = i;
    
    for(j = i+1; j < length; j++)
    {
      r_val = mpfi_abs(tmp0, max_piv);
      r_val = mpfi_abs(tmp1, matA[j][i]);
      
      if(mpfi_cmp(tmp0, tmp1) < 0)
      {
        mpfi_set(max_piv, matA[j][i]);
        j_max = j;
      }
    }

    if(mpfi_is_zero(max_piv) <= 0)
    {
      if(i != j_max)
      {
        swaprow_mat(matA, i, j_max, length);
        swaprow_mat(matB, i, j_max, length);
      }

      mpfi_neg(tmp1, matA[i][i]);
      mpfi_ui_div(tmp2, 1, tmp1); /* (-1/matA[i][i]) */

      for(j = i+1; j < length; j++)
      {
        mpfi_set(tmp0, matA[j][i]);

        multline_mat(matA, tmp1, j, length);
        addmultline_mat(matA, tmp0, i, j, length);
        multline_mat(matA, tmp2, j, length);

        multline_mat(matB, tmp1, j, length);
        addmultline_mat(matB, tmp0, i, j, length);
        multline_mat(matB, tmp2, j, length);
      }
      
    }
  }
  
  for(i = length-1; i >= 0; i--)
  {
    if(mpfi_is_zero(matA[i][i]) <= 0)
    {
      mpfi_ui_div(tmp0, 1, matA[i][i]);

    multline_mat(matA, tmp0, i, length);
      multline_mat(matB, tmp0, i, length);

      for(j = i-1; j >= 0; j--)
      {
        mpfi_set(tmp0, matA[j][i]);
        mpfi_set_si(tmp1, -1);

        multline_mat(matA, tmp1, j, length);
        addmultline_mat(matA, tmp0, i, j, length);
        multline_mat(matB, tmp1, j, length);
        addmultline_mat(matB, tmp0, i, j, length);

      }
    }
  }

  mpfi_clear(max_piv);
}


void si_inc(mpfi_t rop, int c, int u, int d, int i)
{
  C(tmp2, w-u, i);        mpfi_set(rop, tmp2);  
  C(tmp2, n-k-w+u, c-i);  mpfi_mul(rop, rop, tmp2);
  C(tmp2, u, d+i);        mpfi_mul(rop, rop, tmp2);
  C(tmp2, k-u, c-d-i);    mpfi_mul(rop, rop, tmp2);
  C(tmp2, n-k, c);        mpfi_div(rop, rop, tmp2);
  C(tmp2, k, c);          mpfi_div(rop, rop, tmp2);
}


void si_dec(mpfi_t rop, int c, int u, int d, int i)
{
  C(tmp2, w-u, d+i);         mpfi_set(rop, tmp2);
  C(tmp2, n-k-w+u, c-d-i);   mpfi_mul(rop, rop, tmp2);
  C(tmp2, u, i);             mpfi_mul(rop, rop, tmp2);
  C(tmp2, k-u, c-i);         mpfi_mul(rop, rop, tmp2);
  C(tmp2, n-k, c);           mpfi_div(rop, rop, tmp2);
  C(tmp2, k, c);             mpfi_div(rop, rop, tmp2);
}


void P_compute(int p, int c)
{
  int u, d, i;

  /* all entries are set to zero */
  clear_mat(P, w+2);

  for(u = 0; u < w+1; u++)
  {
    for(d = c; d >= 0; d--)
    {
      /* increasing errors */
      if((u-d >= 0) && (u-d < w+1))
      {
        mpfi_set_ui(tmp0, 0);

        for(i = 0; i<=(c-d); i++)
        {
          if((w-u >= i) && (u >= d+i))
          {
            si_inc(tmp1, c, u, d, i);
            mpfi_add(tmp0, tmp0, tmp1);
          }
        }
        mpfi_set(P[u][u-d], tmp0);
      }
      
      /* decreasing errors */
      if((u+d > 0) && (u+d < w+1))
      {
        mpfi_set_ui(tmp0, 0);

        for(i = 0; i<=(c-d); i++)
        {
          if((w-u >= d+i) && (u >= i))
          {
            si_dec(tmp1, c, u, d, i);
            mpfi_add(tmp0, tmp0, tmp1);
          }
        }
        mpfi_set(P[u][u+d], tmp0);
      }
    }

  }
  
  /* transition: 2p+1 -> (2p)_S */
  mpfi_set(P[2*p+1][w+1], P[2*p+1][2*p]); 
  
  /* transition: 2p-1 -> (2p)_S */
  mpfi_set(P[2*p-1][w+1], P[2*p-1][2*p]); 

  /* transition: (2p)_F -> (2p)_S */
  mpfi_set(P[2*p][w+1], P[2*p][2*p]); 

  /* transition: (2p)_S -> (2p)_S */
  mpfi_ui_div(P[w+1][w+1], 1, beta);
}


void R_compute(int p)
{
  int i;

  id_mat(Tmp, w+2);
  mpfi_ui_sub(Tmp[2*p][2*p], 1, beta);
  mpfi_set(Tmp[w+1][w+1], beta);
  
  prod_mat(R, P, Tmp, w+2);
  
  /* -(Id-R) = R-Id */
  /* subtract -1 from all diagonal elements of R */
  for(i = 0; i < w+1; i++)
  {
    mpfi_sub_ui(R[i][i], R[i][i], 1);
  }
  
  /* inverse */
  copy_mat(Tmp, R, w+1);
  gaussjord_mat(R, Tmp, w+1);
}


void pi(mpfi_t pi0, int u, int p)
{
  /* tmp0, tmp1 are in use --> tmp2 */

  if(u < 0 || u > w)  abort();

  if(u == 2*p)
  {
    mpfi_ui_sub(pi0, 1, beta);
  }
  else
  {
    mpfi_set_ui(pi0, 1);
  }

  C(tmp2, w, u); mpfi_mul(pi0, pi0, tmp2);
  C(tmp2, n-w, k-u); mpfi_mul(pi0, pi0, tmp2);
  C(tmp2, n, k); mpfi_div(pi0, pi0, tmp2);

}


void it_count(int p)
{
  int u, v;
  
  mpfi_set_ui(num, 0);
  
  for(u = 0; u <= w; u++)
  {
    mpfi_set_ui(tmp0, 0);
    
    /* we have -R, since we considered (Id-Q)^(-1) */
    for(v = 0; v <= w; v++)
    {
      mpfi_sub(tmp0, tmp0, R[u][v]);
    }

    pi(tmp1, u, p);
    mpfi_mul(tmp0, tmp0, tmp1);
    mpfi_add(num, num, tmp0);

  }

}


int main(int argc, char **argv)
{
  int prec = 150;
  mpfr_set_default_prec(prec);

  int p = 2;
  int m = 2;
  int l = 20;
  int c = 7;
  int r = 0;
  int bestr = r;
  int i;
  int sign;
  double opsperiteration;
  double ops;
  double iterations;
  double bestops = 0.;
  double twor;
  double twol = exp(l * log(2.0));


  factorial_init();
  mpfi_init(Anum);
  mpfi_init(Bnum);
  mpfi_init(beta);
  mpfi_init(tmp0);
  mpfi_init(tmp1);
  mpfi_init(tmp2);
  mpfi_init(num);
  mpfr_init(numleft);
  mpfr_init(numright);

//  for (w = 4; w<=50; w++)
  {
    P = init_mat(w+2);
    Tmp = init_mat(w+2);
    R = init_mat(w+2);

    C(Anum, x, p);
    C(Bnum, k-x, p);

    /* set success probability beta */
    mpfi_set_ui(tmp1, 0);
    sign = -1;

    for(i = 1; i<=m; i++)
    {
      sign *= -1;
      C(tmp0, m, i);  mpfi_mul_si(tmp0, tmp0, sign);
      C(tmp2, n-k-w+2*p, i*l);  mpfi_mul(tmp0, tmp0, tmp2);
      C(tmp2, n-k, i*l);  mpfi_div(tmp0, tmp0, tmp2);
      mpfi_add(tmp1, tmp1, tmp0);
    }

    mpfi_mul(beta, Anum, Bnum);
    mpfi_mul(beta, beta, tmp1);
    C(tmp2, k, 2*p);  mpfi_div(beta, beta, tmp2);

    /* initialize Markov matrix P */
    P_compute(p, c);

    R_compute(p);
    it_count(p);

    if (!mpfi_bounded_p(num)) abort();
    mpfi_get_left(numleft, num);
    mpfi_get_right(numright, num);
    if (mpfr_get_d(numleft, GMP_RNDD) <= 0) abort();
    if (mpfr_get_d(numright, GMP_RNDU)/mpfr_get_d(numleft, GMP_RNDD) > 1.001) abort();

    while(r < c)
    {
      ++r;
      twor = exp(r * log(2.0));
      opsperiteration = (n-1) * ((n-k-1)*(1-1/twor)+(twor-r)) * ceil(c * 1.0 / r);
      opsperiteration += m * l * (mpfi_get_d(Anum) + mpfi_get_d(Bnum));
      opsperiteration += m * 2.0 * (w - 2 * p) * (p + p - 1) * mpfi_get_d(Anum) * mpfi_get_d(Bnum) / twol;

      iterations = mpfi_get_d(num) + 1;
      ops = opsperiteration * iterations;
    
      if (!bestops || ops < bestops)
      {
        bestr = r;
        bestops = ops;
      }
    }
       
    r = bestr;
    ops = bestops;

    printf("type-3 iteration count:\nn=%d k=%d x=%d w=%d p=%d m=%d l=%d c=%d r=%d: \n", n, k, x, w, p, m, l, c, r);
    printf("ops \t%lf\n", ops);
    printf("ops/it \t%lf\n", ops / iterations);
    printf("it \t%lf\n", iterations);
    fflush(stdout);

    free_mat(P, w+2);
    free_mat(Tmp, w+2);
    free_mat(R, w+2);

  }

  mpfi_clear(Anum);
  mpfi_clear(Bnum);
  mpfi_clear(beta);
  mpfi_clear(tmp0);
  mpfi_clear(tmp1);
  mpfi_clear(tmp2);
  mpfi_clear(num);
  mpfr_clear(numleft);
  mpfr_clear(numright);

  factorial_free();
  
  return 0;

}
