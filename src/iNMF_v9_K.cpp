// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <time.h>

// Record running time
clock_t start, end;

using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
List iNMF_BCD_Decrease(List count_list,
                       NumericVector Label,
                       int K_l,
                       IntegerVector K,
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
   iNMF version 9

   BCD for iNMF algorithm

   Args:
   count_list: A list of cell by gene matrices(length L).
   Label: Label of each Data
   K_l: Column dimensions of W_1
   K: Row dimensions of V
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
  int J = max(Label)+1;


  int n_K = K.size();
  K.sort(true);

  int K_1 = K[0];

  std::vector<arma::mat> X(L), W_2(L), W_1(L), H(J);
  for(int i=0;i<L;i++)
  {
    mat temp = count_list[i];
    X[i] = temp;
  }

  int N = X[0].n_cols;

  mat V = randu<mat>(K_1, N);

  NumericVector M_l(L);
  for(int i=0; i<L; i++)
  {
    M_l[i] = X[i].n_rows;
  }
  for(int i=0; i<L; i++)
  {
    W_2[i] = randu<mat>(M_l[i], K_1);
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
  List output(n_K);
  int num;
  for(num=0;num<n_K;num++)
  {
    int K_2 = K[num];
    Rcpp::Rcout<<"K="<<K_2<<std::endl;

    if(num>0)
    {
      NumericVector gamma(K_1);
      for(int k=0;k<K_1;k++)
      {
        gamma[k] = 0;
        for(int i=0;i<L;i++)
        {
          gamma[k] += sum(W_2[i].col(k)%W_2[i].col(k))/M_l[i];
        }
        gamma[k] *= sum(V.row(k)%V.row(k));
      }

      NumericVector ind(K_1);

      // sort and get index
      for(int i=0;i<K_1;i++)
      {
        ind[i] = i;
      }
      for(int i=0;i<K_1;i++)
      {
        for(int j=0;j<K_1-i-1;j++)
        {
          if(gamma[j]>gamma[j + 1])
          {
            double temp = gamma[j];
            gamma[j] = gamma[j+1];
            gamma[j+1] = temp;

            int ind_temp = ind[j];
            ind[j] = ind[j+1];
            ind[j+1] = ind_temp;
          }
        }
      }
      // finish sort

      uvec index(K_2);
      for(int k=0;k<K_2;k++)
      {
        index[k] = ind[k];
      }

      V = V.rows(index);
      for(int i=0;i<L;i++)
      {
        W_2[i] = W_2[i].cols(index);
      }
      K_1 = K_2;
    }


    // Iterations
    double Loss_1 = 0.0;
    double Loss_2 = 1000000.0;

    start = clock(); // Record the start time of iterations
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
      for(int k=0; k<K_2; k++)
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
          for(int k=0; k<K_2; k++)
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
          for(int k=0; k<K_2; k++)
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

    end = clock(); // Record the end time of iterations

    // temporary variables
    std::vector<arma::mat> oW_2(L), oW_1(L), oH(J);
    mat oV = V;

    for(int i=0; i<L; i++)
    {
      oW_2[i] = W_2[i];
      oW_1[i] = W_1[i];
    }

    for(int j=0; j<J; j++)
    {
      oH[j] = H[j];
    }


    // Convert values that less than eps to zero
    oV.elem(find(oV<=eps)).zeros();

    for(int l=0;l<L;l++)
    {
      oW_1[l].elem(find(oW_1[l]<=eps)).zeros();
      oW_2[l].elem(find(oW_2[l]<=eps)).zeros();
    }

    for(int l=0;l<J;l++)
    {
      oH[l].elem(find(oH[l]<=eps)).zeros();
    }

    // Normalize V
    for(int i=0;i<K_2;i++)
    {
      double n_V = sqrt(sum(oV.row(i)%oV.row(i)));
      oV.row(i) /= n_V;
      for(int l=0;l<L;l++)
      {
        oW_2[l].col(i) *= n_V;
      }
    }

    // Convert the results to several Lists (Can not directly assign the array of mat into List)
    // It means that it is illegal to write down 'Named("W_1") = W_1', W_1 needs to be converted to a List W_1r
    List W_1r(L),W_2r(L),H_r(J);
    for(int i=0;i<L;i++)
    {
      W_1r[i] = oW_1[i];
      W_2r[i] = oW_2[i];
    }

    for(int j=0;j<J;j++)
    {
      H_r[j] = oH[j];
    }



    output[num] = List::create(Named("W_1") = W_1r,
                               Named("W_2") = W_2r,
                               Named("H")  =  H_r,
                               Named("V")  =  oV,
                               Named("iters") = iter,
                               Named("loss")  =  Loss_2,
                               Named("times") = (double)(end-start)/CLOCKS_PER_SEC);
  }
  // Return the result as a List
  return output;
}

// [[Rcpp::export]]
List iNMF_BCD_Increase(List count_list,
                       NumericVector Label,
                       int K_l,
                       IntegerVector K,
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
   iNMF version 9

   BCD for iNMF algorithm

   Args:
   count_list: A list of cell by gene matrices(length L).
   Label: Label of each Data
   K_l: Column dimensions of W_1
   K: Row dimensions of V
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
  int J = max(Label)+1;


  int n_K = K.size();
  K.sort();
  int K_1 = K[0];

  std::vector<arma::mat> X(L), W_2(L), W_1(L), H(J);
  for(int i=0;i<L;i++)
  {
    mat temp = count_list[i];
    X[i] = temp;
  }

  int N = X[0].n_cols;


  mat V = randu<mat>(K_1,N);

  NumericVector M_l(L);
  for(int i=0; i<L; i++)
  {
    M_l[i] = X[i].n_rows;
  }


  for(int i=0; i<L; i++)
  {
    W_2[i] = randu<mat>(M_l[i], K_1);
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
  List output(n_K);
  int num, K_add;
  std::vector<arma::mat> E(L), W_add(L);
  mat V_add;

  for(num=0;num<n_K;num++)
  {
    int K_2 = K[num];
    Rcpp::Rcout<<"K="<<K_2<<std::endl;

    if(num>0)
    {
      K_add = K_2 - K_1;

      // initiate V_add & W_add
      V_add = randu<mat>(K_add, N);
      for(int i=0;i<L;i++)
      {
        E[i] = X[i] - W_1[i]*H[Label[i]] - W_2[i]*V;
        W_add[i] = randu<mat>(M_l[i], K_add);
      }

      // iterations for V_add & W_add
      int iter;
      int add_thre_iter = 300;

      double add_Loss_1 = 0.0;
      double add_Loss_2 = 1000000.0;
      double add_thre_delta = 0.0001;

      // Rcpp::Rcout<<"start initializing V_add & W_add"<<std::endl;

      for(iter=0; iter<add_thre_iter; iter++)
      {
        for(int k=0; k<K_add; k++)
        {
          rowvec V_a = zeros<rowvec>(N);
          double V_b = 0;
          for(int i=0;i<L;i++)
          {
            V_a += W_add[i].col(k).t()*(E[i] - W_add[i]*V_add)/M_l[i];
            V_b += sum(W_add[i].col(k)%W_add[i].col(k))/M_l[i];
          }
          V_add.row(k) = max(V_add.row(k) + V_a/V_b, eps*ones<rowvec>(N));
        }

        for(int i=0; i<L; i++)
        {
          for(int k=0; k<K_add; k++)
          {
            W_add[i].col(k) = max(W_add[i].col(k) + (E[i] - W_add[i]*V_add)*(V_add.row(k).t())/sum(V_add.row(k)%V_add.row(k)), eps*ones<colvec>(M_l[i]));
          }
        }

        add_Loss_1 = add_Loss_2;
        add_Loss_2 = 0;

        for(int i=0;i<L;i++)
        {
          add_Loss_2 += pow(norm(E[i]-W_add[i]*V_add,"fro"),2)/M_l[i];
        }

        // Rcpp::Rcout<<"init_loss:"<<add_Loss_1-add_Loss_2<<std::endl;

        if((add_Loss_1-add_Loss_2)<add_thre_delta)
        {
          break;
        }

      }

      // Rcpp::Rcout<<"Finish initializing V_add & W_add with time: "<<(double)(end-start)/CLOCKS_PER_SEC<<"s"<<std::endl;

      for(int k=0; k<K_add; k++)
      {
        double n_V = sqrt(sum(V_add.row(k)%V_add.row(k)));
        V_add.row(k) /= n_V;
        for(int i=0; i<L; i++)
        {
          W_add[i].col(k) *= n_V;
        }
      }

      mat temp_V(K_2, N);
      std::vector<arma::mat> temp_W(L);
      for(int i=0;i<L;i++)
      {
        temp_W[i] = zeros<mat>(M_l[i],K_2);
      }

      for(int k=0;k<K_2;k++)
      {
        if(k<K_1)
        {
          temp_V.row(k) = V.row(k);
          for(int i=0;i<L;i++)
          {
            temp_W[i].col(k) = W_2[i].col(k);
          }
        }
        else
        {
          temp_V.row(k) = V_add.row(k-K_1);
          for(int i=0;i<L;i++)
          {
            temp_W[i].col(k) = W_add[i].col(k-K_1);
          }
        }
      }

      V = temp_V;
      for(int i=0;i<L;i++)
      {
        W_2[i] = temp_W[i];
      }

      K_1 = K_2;
    }


    // Iterations
    double Loss_1 = 0.0;
    double Loss_2 = 1000000.0;

    start = clock(); // Record the start time of iterations
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

      for(int k=0; k<K_2; k++)
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
          for(int k=0; k<K_2; k++)
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
          for(int k=0; k<K_2; k++)
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

    end = clock(); // Record the end time of iterations

    // temporary variables
    std::vector<arma::mat> oW_2(L), oW_1(L), oH(J);
    mat oV = V;


    for(int i=0; i<L; i++)
    {
      oW_2[i] = W_2[i];
      oW_1[i] = W_1[i];
    }


    for(int j=0; j<J; j++)
    {
      oH[j] = H[j];
    }

    // Convert values that less than eps to zero
    oV.elem(find(oV<=eps)).zeros();

    for(int l=0;l<L;l++)
    {
      oW_1[l].elem(find(oW_1[l]<=eps)).zeros();
      oW_2[l].elem(find(oW_2[l]<=eps)).zeros();
    }

    for(int l=0;l<J;l++)
    {
      oH[l].elem(find(oH[l]<=eps)).zeros();
    }


    // Normalize V
    for(int i=0;i<K_2;i++)
    {
      double n_V = sqrt(sum(oV.row(i)%oV.row(i)));
      oV.row(i) /= n_V;
      for(int l=0;l<L;l++)
      {
        oW_2[l].col(i) *= n_V;
      }
    }



    // Convert the results to several Lists (Can not directly assign the array of mat into List)
    // It means that it is illegal to write down 'Named("W_1") = W_1', W_1 needs to be converted to a List W_1r
    List W_1r(L),W_2r(L),H_r(J);
    for(int i=0;i<L;i++)
    {
      W_1r[i] = oW_1[i];
      W_2r[i] = oW_2[i];
    }

    for(int j=0;j<J;j++)
    {
      H_r[j] = oH[j];
    }


    output[num] = List::create(Named("W_1") = W_1r,
                               Named("W_2") = W_2r,
                               Named("H")  =  H_r,
                               Named("V")  =  oV,
                               Named("iters") = iter,
                               Named("loss")  =  Loss_2,
                               Named("times") = (double)(end-start)/CLOCKS_PER_SEC);
  }
  // Return the result as a List
  return output;
}
