#' The scINSIGHT Class
#'
#' The scINSIGHT object is created from two or more single cell datasets. To construct a
#' scINSIGHT object, the user needs to provide at least two normalized expression (or another
#' single-cell modality) matrices and the condition vector.
#'
#' The key slots used in the scINSIGHT object are described below.
#'
#' @slot norm.data List of normalized expression matrices (genes by cells). Each matrix should have the same number and name of genes.
#' @slot condition Vector specifying each sample's condition name.
#' @slot W_1 List of \eqn{W_{\ell1}} estimated by scINSIGHT, names correspond to sample names.
#' @slot W_2 List of \eqn{W_{\ell2}} estimated by scINSIGHT, names correspond to sample names.
#' @slot H List of \eqn{H} estimated by scINSIGHT, names correspond to condition names.
#' @slot V Matrix \eqn{V} estimated by scINSIGHT.
#' @slot norm.W_2 List of \eqn{W_{\ell2}} after normalization. Recommended for downstream analysis.
#' @slot parameters List of selected parameters, including \eqn{K} and \eqn{\lambda}.
#'
#' @exportClass scINSIGHT
#' @name scINSIGHT-class
#' @rdname scINSIGHT-class

scINSIGHT <- methods::setClass(
  "scINSIGHT",
  slots = c(
    norm.data = "list",
    condition = "character",
    W_1 = "list",
    W_2 = "list",
    H = "list",
    V = "matrix",
    norm.W_2 = "list",
    parameters = "list"
  )
)

#' Create an scINSIGHT object.
#'
#' This function initializes an scINSIGHT object with normalized data passed in.
#'
#' @param norm.data List of normalized expression matrices (genes by cells). Gene names should be the same in all matrices.
#' @param condition Vector specifying sample conditions.

#' @return \code{scINSIGHT} object with norm.data slot set.

#' @importFrom methods new is as
#'
#' @export
#' @examples
#' # Demonstration using matrices with randomly generated numbers
#' S1 <- matrix(runif(50000,0,2), 500,100)
#' S2 <- matrix(runif(60000,0,2), 500,120)
#' S3 <- matrix(runif(80000,0,2), 500,160)
#' S4 <- matrix(runif(75000,0,2), 500,150)
#' data = list(S1, S2, S3, S4)
#' sample = c("sample1", "sample2", "sample3", "sample4")
#' condition = c("control", "activation", "control", "activation")
#' names(data) = sample
#' names(condition) = sample
#' scINSIGHTx <- create_scINSIGHT(data, condition)

create_scINSIGHT <- function(norm.data, condition) {

  L = length(norm.data)

  if(L<=1){
    stop("Please provide at least two expression matrices.")
  }


  for(i in 1:(L-1)){
    if(!all(rownames(norm.data[[i]]) == rownames(norm.data[[i+1]]))){
      stop("Please ensure that each sample's number and name of gene are consistent.")
    }
  }

  if(length(condition) != L){
    stop("Please ensure that the number of 'condition' is the same as the number of 'norm.data'.")
  }

  if(!all(names(norm.data) == names(condition))){
    stop("Please ensure that sample and condition share the same names.")
  }

  object <- methods::new(Class = "scINSIGHT", norm.data = norm.data, condition = condition)


  return(object)
}


#' Perform scINSIGHT on normalized datasets
#'
#' @description
#' Perform INterpreting single cell gene expresSIon bioloGically Heterogeneous daTa (scINSIGHT) to return factorized \eqn{W_{\ell1}}, \eqn{W_{\ell2}}, \eqn{H} and \eqn{V} matrices.
#'
#' This factorization produces a \eqn{W_{\ell1}} matrix (cells by \eqn{K_j}), a \eqn{W_{\ell2}} matrix (cells by \eqn{K}), a shared \eqn{V} matrix (\eqn{K} by genes)
#' for each sample, and a \eqn{H} (\eqn{K_j} by genes) matrix for each condition. \eqn{W_{\ell2}} are the expression matrices of \eqn{K} common gene pathways for all samples,
#' \eqn{V} is the membership matrix of \eqn{K} common gene pathways, and it's shared by all samples.
#' \eqn{W_{\ell1}} are the expression matrices of \eqn{K_j} condition-specific gene pathways for all samples,
#' and \eqn{H} are the membership matrices of \eqn{K_j} condition-specific gene pathways for all conditions.
#'
#' @param object \code{scINSIGHT} object.
#' @param K Number of common gene pathways. (default \code{c(5, 7, 9, 11, 13, 15)})
#' @param K_j Number of dataset-specific gene pathways. (default 2)
#' @param LDA Regularization parameters. (default \code{c(0.001, 0.01, 0.1, 1, 10)})
#' @param num.cores Number of cores used for optimizing factorizations in parallel (default 1).
#' @param B Number of repeats with random seed from 1 to B. (default 5)
#' @param out.dir Output directory of scINSIGHT results. (default NULL)
#' @param method Method of updating the factorization (default "increase"). If provide multiple \eqn{K}, user can choose method between "increase" and "decrease".
#'
#' For "increase", the algorithm will first perform factorization with the least \eqn{K=K_1}. Then initialize \eqn{K_2-K_1} facotrs,
#'  where \eqn{K_2} is the \eqn{K} sightly larger than \eqn{K_1}, and perform facotrization with these new facotrs. Continue this process until the largest \eqn{K}.
#'
#' For "increase", the algorithm will first perform factorization with the largest \eqn{K=K_1}. Then choose \eqn{K_2} facotrs,
#'  where \eqn{K_2} is the \eqn{K} sightly less than \eqn{K_1}, and perform facotrization with these new facotrs. Continue this process until the least \eqn{K}.
#' @param thre.niter Maximum number of block coordinate descent iterations to perform. (default 500)
#' @param thre.delta Stop iteration when the reduction of objective function is less than the threshold. (default 0.01)
#
#' @return \code{scINSIGHT} object with \eqn{W_1}, \eqn{W_2}, \eqn{H}, \eqn{V} and parameters slots set.
#'
#' @importFrom stringr str_sub
#' @importFrom parallel mclapply
#' @importFrom Rcpp evalCpp
#' @useDynLib scINSIGHT
#'
#' @export


run_scINSIGHT<- function(
  object,
  K = seq(5,15,2),
  K_j = 2,
  LDA = c(0.001, 0.01, 0.1, 1, 10),
  thre.niter = 500,
  thre.delta = 0.01,
  num.cores = 1,
  B = 5,
  out.dir = NULL,
  method = "increase"
){
  if(!(method %in% c("increase", "decrease"))){
    stop("Please input right method")
  }

  object@parameters[["K_j"]] = K_j

  LDA = sort(LDA)
  index_lda = which.min(abs(LDA-0.01))
  lda0 = LDA[[index_lda]]
  object@parameters[["lda"]] = lda0

  n_K = length(K)
  K = sort(K)
  cnt_list = object@norm.data
  cnt_list = lapply(cnt_list, function(x){t(x)})
  labels = object@condition

  if(!is.null(out.dir)){
    if(str_sub(out.dir, -1)!="/"){
      out.dir = paste0(out.dir,"/")
    }
    dir.create(out.dir, recursive = TRUE)
  }

  seeds = 1:B

  message(paste("Run scINSIGHT with method:", method, "K, and with lda =",lda0))

  # Run iNMF BCD
  if(method == "increase")
  {
    res_all = mclapply(seeds, function(seed){
      res = iNMF_BCD_Increase(cnt_list, labels, K_j, K = K, lda1 = lda0, lda2 = lda0, eps = 1e-5,
                              thre.niter, thre.delta, seed, TRUE)
      gc()
      return(res)
    }, mc.cores = num.cores)

    res_parallel = list()
    for(i in 1:n_K){
      res_name = paste0("K_",as.character(K[[i]]))
      res_parallel[[res_name]] = lapply(seeds, function(seed){res_all[[seed]][[i]]})
      names(res_parallel[[res_name]]) = sapply(seeds, function(seed){paste0("seed_",seed)})
      if(!is.null(out.dir)){
        saveRDS(res_parallel[[res_name]], file = paste0(out.dir, "res-k", K[[i]], "-lda", lda0, ".rds"))
      }

      object@parameters[["stability"]][i] = get_stability_StrictAndWeight_k(res_parallel[[res_name]], nk = 20)

    }

    names(object@parameters[["stability"]]) = sapply(1:n_K,function(i){paste0("K_",as.character(K[[i]]))})
    print(object@parameters[["stability"]])

    message("run increaseK completed!")

  }

  if(method == "decrease")
  {
    res_all = mclapply(seeds, function(seed){
      res = iNMF_BCD_Decrease(cnt_list, labels, K_j, K = K, lda1 = lda0, lda2 = lda0, eps = 1e-5,
                              thre.niter, thre.delta, seed, TRUE)
      gc()
      return(res)
    }, mc.cores = num.cores)

    res_parallel = list()
    for(i in 1:n_K){
      res_name = paste0("K_",as.character(K[[i]]))
      res_parallel[[res_name]] = lapply(seeds, function(seed){res_all[[seed]][[i]]})
      names(res_parallel[[res_name]]) = sapply(seeds, function(seed){paste0("seed_",seed)})
      if(!is.null(out.dir)){
        saveRDS(res_parallel[[res_name]], file = paste0(out.dir, "res-k", K[[i]], "-lda", lda0, ".rds"))
      }

      object@parameters[["stability"]][[i]] = get_stability_StrictAndWeight_k(res_parallel[[res_name]], nk = 20)

    }

    names(object@parameters[["stability"]]) = sapply(1:n_K, function(i){paste0("K_",as.character(K[[i]]))})
    message("run decreaseK completed!")
  }


  ### Select K
  if(length(K)>1){
    message("Select optimal K...")
  }
  bestK = K[which.max(object@parameters[["stability"]])]
  object@parameters[["K"]] = bestK
  if(length(K)>1){
    message(paste("Best K is selected as K =", bestK))
  }

  ### Select lda
  if(length(LDA)>1){
    message("Select optimal lda...")
  }
  res_parallel = res_parallel[[paste0("K_",as.character(bestK))]]

  if((!is.null(out.dir)) && (length(LDA)>1)){
    out.dir = paste0(out.dir, "select_lda/")
    object = select_LDA(object, res_parallel, LDA = LDA, out.dir = out.dir, ncores = num.cores,
                        thre.delta = thre.delta, thre.niter = thre.niter, topn = 100)
  }
  else{
    object = select_LDA(object, res_parallel, LDA = LDA, ncores = num.cores,
                        thre.delta = thre.delta, thre.niter = thre.niter, topn = 100)
  }

  if(length(LDA)>1){
    message(paste("Best lambda is selected as lambda =", object@parameters$lda))
  }

  if(length(K)>1 || length(LDA)>1){
    message("Finish selecting parametesr...")
  }

  return(object)

}

