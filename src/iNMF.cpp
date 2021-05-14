// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <time.h>

// Record running time
clock_t start_t, end_t;

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
List iNMF_BCD(List count_list,
              NumericVector Label,
              int K_l,
              int K,
              double lda1,
              double lda2,
              double eps=1e-5,
              int thre_niter=5000,
              double thre_delta=0.00001,
              int seed = 1,
              bool loop=true
)
{
  /*
   iNMF version 8

   BCD for iNMF algorithm

   Args:
   count_list: A list of cell by gene matrices(length L).
   Label: Label of each Data
   K_l: Column dimensions of W_1(length J)
   K: Row dimension of V
   lda1: lambda1
   lda2: lambda2
   eps: a number to ensure numerical stability
   thre_niter: maximum number of iterations
   thre_delta: if the values of objective funciton in two
   iterations differ less than thre_delta, stop the iteration
   loop: whether cycle at the same time

   Returns:
   W_1, W_2, H, V: the factory matrices of iNMF output

   */

  // setseed
  Rcpp::Environment baseEnv("package:base");
  Rcpp::Function setSeed = baseEnv["set.seed"];
  setSeed(seed);

  // Initialize output
  int L = count_list.size();
  int J = Label.size();

  // Rcpp::Rcout<<"Start iNMF"<<std::endl;
  std::vector<arma::mat> X(L), W_1(L), W_2(L), H(J);
  for(int i=0;i<L;i++)
  {
    mat temp = count_list[i];
    X[i] = temp;
  }

  int N = X[0].n_cols;


  mat V = randu<mat>(K,N);

  NumericVector M_l(L);
  for(int i=0; i<L; i++)
  {
    M_l[i] = X[i].n_rows;
  }


  for(int i=0;i<L;i++)
  {
    W_2[i] = randu<mat>(M_l[i], K);
    W_1[i] = randu<mat>(M_l[i], K_l);
  }

  for(int i=0; i<J; i++)
  {
    H[i] = randu<mat>(K_l, N);
    for(int k=0; k<K_l; k++)
    {
      double n_H = sqrt(sum(H[i].row(k)%H[i].row(k)));
      H[i].row(k) /= n_H;
      for(int m=0; m<L; m++)
      {
        if(Label[m]==i)
        {
          W_1[m].col(k) *= n_H;
        }
      }

    }
  }

  // Loss and Iteration termination criterion
  double Loss_1 = 0.0;
  double Loss_2 = 1000000.0;

  // Iterations
  start_t = clock(); // Record the start time of iterations
  int iter; // Record the number of iterations

  for(iter=0; iter<thre_niter; iter++)
  {
    if(iter % 20 == 0)
    {
      Rcpp::Rcout<<iter<<std::endl;
    }

    std::vector<arma::mat> Res_v(L);
    for(int i=0; i<L; i++)
    {
      Res_v[i] = X[i] - W_1[i]*H[Label[i]];

    }

    for(int k=0; k<K; k++)
    {
      rowvec V_a = zeros<rowvec>(N);
      double V_b = 0;
      for(int i=0;i<L;i++)
      {
        V_a += W_2[i].col(k).t()*(Res_v[i]-W_2[i]*V)/M_l[i];
        V_b += sum(W_2[i].col(k)%W_2[i].col(k))/M_l[i];
      }
      V.row(k) = max(V.row(k) + V_a/V_b, eps*ones<rowvec>(N));
    }

    if(loop) // Whether to solve W_2 and W_1 at the same time
    {
      for(int i=0; i<L; i++)
      {
        mat Res_2 = X[i] - W_1[i]*H[Label[i]];
        for(int k=0; k<K; k++)
        {
          W_2[i].col(k) = max(W_2[i].col(k) + (Res_2 - W_2[i]*V)*(V.row(k).t())/sum(V.row(k)%V.row(k)), eps*ones<colvec>(M_l[i]));
        }

        mat Res_h = X[i] - W_2[i]*V;

        for(int k=0; k<K_l; k++)
        {
          W_1[i].col(k) = max(W_1[i].col(k) + ((Res_h - (1 + lda1)*W_1[i]*H[Label[i]])*(H[Label[i]].row(k).t()))/sum((1+lda1)*H[Label[i]].row(k)%H[Label[i]].row(k)),eps*ones<colvec>(M_l[i]));
        }

      }
    }
    else
    {
      for(int i=0; i<L; i++)
      {
        mat Res_2 = X[i] - W_1[i]*H[Label[i]];
        for(int k=0; k<K; k++)
        {
          W_2[i].col(k) = max(W_2[i].col(k) + ((Res_2-W_2[i]*V)*(V.row(k).t()))/sum(V.row(k)%V.row(k)),eps*ones<colvec>(M_l[i]));
        }
      }

      for(int i=0;i<L;i++)
      {
        mat Res_h = X[i]-W_2[i]*V;
        for(int k=0;k<K_l;k++)
        {
          W_1[i].col(k) = max(W_1[i].col(k)+(Res_h-(1+lda1)*W_1[i]*H[Label[i]])*H[Label[i]].row(k).t()/sum((1+lda1)*H[Label[i]].row(k)%H[Label[i]].row(k)),eps*ones<colvec>(M_l[i]));
        }

      }
    }

    // update H
    for(int i=0; i<J; i++)
    {
      for(int k=0; k<K_l; k++)
      {
        rowvec temp_H = zeros<rowvec>(N);
        double temp_W = 0.0;
        for(int m=0; m<L; m++)
        {
          if(Label[m]==i)
          {
            temp_H += (W_1[m].col(k).t()*(X[m]-W_2[m]*V-(1+lda1)*W_1[m]*H[i]))/M_l[m];
            temp_W += (1+lda1)*sum(W_1[m].col(k)%W_1[m].col(k))/M_l[m];
          }
        }


        rowvec s_H = zeros<rowvec>(N);
        for(int j=0; j<J; j++)
        {
          if(j==i)
          {
            continue;
          }
          for(int t=0; t<K_l; t++)
          {
            s_H += H[j].row(t);
          }
        }

        H[i].row(k) = max(H[i].row(k) + (temp_H - (lda2/4)*s_H)/temp_W, eps*ones<rowvec>(N));

        double n_H = sqrt(sum(H[i].row(k)%H[i].row(k)));
        H[i].row(k) /= n_H;
        for(int m=0; m<L; m++)
        {
          if(Label[m]==i)
          {
            W_1[m].col(k) *= n_H;
          }
        }

      }
    }


    // Calculate the value of Objective function
    Loss_1 = Loss_2;
    Loss_2 = 0;

    for(int i=0;i<L;i++)
    {
      Loss_2 += pow(norm(X[i]-W_1[i]*H[Label[i]]-W_2[i]*V,"fro"),2)/M_l[i];
      Loss_2 += lda1*pow(norm(W_1[i]*H[Label[i]],"fro"),2)/M_l[i];
    }


    for(int i=0; i<J; i++)
    {
      for(int j=0;j<J; j++)
      {
        if(j==i)
        {
          continue;
        }
        Loss_2 += lda2/2*accu(H[i]*(H[j].t()));
      }
    }


    if((Loss_1-Loss_2)<thre_delta)
    {
      break;
    }


  }

  end_t = clock(); // Record the end time of iterations

  // Convert values that less than eps to zero
  V.elem(find(V<=eps)).zeros();

  for(int l=0;l<L;l++)
  {
    W_1[l].elem(find(W_1[l]<=eps)).zeros();
    W_2[l].elem(find(W_2[l]<=eps)).zeros();
  }

  for(int l=0;l<J;l++)
  {
    H[l].elem(find(H[l]<=eps)).zeros();
  }


  // Normalize V
  for(int i=0;i<K;i++)
  {
    double n_V = sqrt(sum(V.row(i)%V.row(i)));
    V.row(i) /= n_V;
    for(int l=0;l<L;l++)
    {
      W_2[l].col(i) *= n_V;
    }
  }


  // Convert the results to several Lists (Can not directly assign the array of mat into List)
  // It means that it is illegal to write down 'Named("W_1") = W_1', W_1 needs to be converted to a List W_1r
  List W_1r(L),W_2r(L),H_r(J);
  for(int i=0;i<L;i++)
  {
    W_1r[i] = W_1[i];
    W_2r[i] = W_2[i];
  }

  for(int i=0;i<J;i++)
  {
    H_r[i] = H[i];
  }


  // Return the result as a List
  return List::create(Named("W_1") = W_1r,
                      Named("W_2") = W_2r,
                      Named("H")  =  H_r,
                      Named("V")  =  V,
                      Named("iters") = iter+1,
                      Named("loss")  =  Loss_2,
                      Named("times") = (double)(end_t-start_t)/CLOCKS_PER_SEC);



  // // Return the result as a List
  // return List::create(Named("W_1") = W_1,
  //                     Named("W_2") = W_2,
  //                     Named("H")  =  H,
  //                     Named("V")  =  V,
  //                     Named("iters") = iter+1,
  //                     Named("loss")  =  Loss_2,
  //                     Named("times") = (double)(end_t-start_t)/CLOCKS_PER_SEC);
}
