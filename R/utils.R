# Utility functions for scINSIGHT methods. Some published, some not.

#' Normalize \eqn{W_2} and cluster cells.
#' @description
#' Quantile normalization and clustering for the list of \eqn{W_2} (expression matrices of common gene modules). Use weights in Louvain.
#'
#' @param W_2 List of \eqn{W_2}.
#' @param Knn The maximum number of nearest neighbors to search (default 20).
#'
#' @return List of normalized \eqn{W_2}.
#'
#' @importFrom RANN nn2
#' @import igraph
#' @importFrom stats approxfun
#'
#' @examples
#' \dontrun{
#' # Requires scINSIGHT object
#' # Get factorization using K from 5 to 15 and LDA from 0.001 to 10, repeat 5 times
#' # (default setting, can be adjusted for ideal results)
#' scINSIGHTex <- run_scINSIGHT(scINSIGHTex, K = seq(5,15,2),
#'                              LDA = c(0.001, 0.01, 0.1, 1, 10), B = 5)
#' }

norm_clust_strict_weighted = function(W2, Knn = 20, ...){
  ### cosine normalization
  W2norm = lapply(W2, function(x){
    l2 = sqrt(rowSums(x^2))
    l2[l2 < 1e-10] = 1e-10
    x = sweep(x, 1, 1/l2, FUN = "*")
    return(x)
  })
  ### cell number
  nc = sapply(W2, nrow)
  nc = c(0,nc)
  ### find knns
  L = length(W2norm)
  find_knns = lapply(1:L, function(l1){
    tp = lapply(1:L, function(l2){
      nns = RANN::nn2(data = W2norm[[l2]], query = W2norm[[l1]],
                      min(Knn, nrow(W2norm[[l2]])),
                      treetype = "kd", searchtype = "standard")
      from = rep((sum(nc[1:l1])+1):sum(nc[1:(l1+1)]), each = min(Knn, nrow(W2norm[[l2]])))
      to = as.numeric(t(nns$nn.idx)) + sum(nc[1:l2])
      weight = as.numeric(t(nns$nn.dists))
      da = data.frame(from = from, to = to, weight = weight)
      return(da)
    })
    return(Reduce(rbind, tp))
  })
  adjmat = Reduce(rbind, find_knns)
  ig = graph_from_data_frame(adjmat, directed = TRUE)
  # ig = subgraph.edges(ig, eids = E(ig)[which_mutual(ig)], delete.vertices = F)
  ig = as.undirected(ig, mode = "mutual", edge.attr.comb="first")
  ig = simplify(ig)
  E(ig)$weight = max(E(ig)$weight) - E(ig)$weight

  cl_louvain = cluster_louvain(ig)$membership
  clusters = split(cl_louvain, rep(1:L, time = nc[-1]))
  Klv = length(unique(cl_louvain))
  ### normalization
  quantiles = 10
  for(j in 1:Klv){
    ns = sapply(clusters, function(cl) sum(cl==j))
    if(sum(ns > 5) < 2) {next}
    refid = min(which(ns > 5))

    cells1 = which(clusters[[refid]] == j)
    num_cells1 = length(cells1)
    #if (num_cells1 < 2) {next}
    for(l in 1:L){
      if(l == refid) next
      cells2 = which(clusters[[l]] == j)
      num_cells2 = length(cells2)
      if (num_cells2 < 2) {next}
      for(i in 1:ncol(W2[[refid]])){
        q2 = quantile(W2[[l]][cells2, i], seq(0, 1, by = 1 / quantiles))
        q1 = quantile(W2[[refid]][cells1, i], seq(0, 1, by = 1 / quantiles))
        if (sum(q1) == 0 | sum(q2) == 0 | length(unique(q1)) <
            2 | length(unique(q2)) < 2) {
          new_vals = rep(0, num_cells2)
        }
        else {
          warp_func = approxfun(q2, q1, rule = 2)
          new_vals = warp_func(W2[[l]][cells2, i])
        }
        W2[[l]][cells2, i] = new_vals
      }
    }

  }
  return(list(W2=W2, clusters = clusters))
}



clust2connect = function(assignment){
  n = length(assignment)
  mat = matrix(0, nrow = n, ncol = n)
  for(l in unique(assignment)){
    id = which(assignment == l)
    mat[id,id] = 1
  }
  return(mat)
}



#' Calculate stability score.
#' @description
#' Calculate stability of result for a certain \eqn{K}.
#'
#' @param res_parallel List of result calculated by scINSIGHT (default NULL).
#' @param res.dir Path to the result (default NULL).
#' @param K Number of \eqn{K}, needed if res.dir is provided (default NULL).
#' @param nk The maximum number of nearest neighbors to search (default 20).
#'
#' @importFrom stats cophenetic hclust as.dist cor quantile
#'
#' @return stability of result for a certain \eqn{K}
#'
#' @examples
#' \dontrun{
#' # If provide result: res_parallel
#' stability <- get_stability_StrictAndWeight_k(res_parallel, nk = 20)
#' # Or provide K and the directory of result: res_parallel
#' stability <- get_stability_StrictAndWeight_k(K = 5, c = 'some_path', nk = 20)
#' }

get_stability_StrictAndWeight_k <- function(
  res_parallel = NULL,
  res.dir = NULL,
  K = NULL,
  nk = 20){

  if(is.null(res_parallel)){
    if(is.null(res.dir)){
      stop("Please input the Path of result file")
    }
    res_parallel = readRDS(file = paste0(res.dir, "res-k", K, ".rds"))
  }

  ### stability -----------------------
  n_cell = sum(sapply(res_parallel[[1]][['W_2']], nrow))
  n_res = length(res_parallel)
  
  x = res_parallel[[1]]
  W2 = x[["W_2"]]
  clust = norm_clust_strict_weighted(W2, Knn = nk)
  W2 = Reduce(rbind, W2)
  assign = unlist(clust$clusters)
  if(n_cell > 20000){
    n_sample = min(2^16-1, floor(n_cell*0.2))
    index = sample(1:n_cell, n_sample)
    assign = assign[index]
  }
  consmat = clust2connect(assign)/n_res
  
  if(n_res > 1){
    
    for(i in c(2:n_res)){
      x = res_parallel[[i]]
      W2 = x[["W_2"]]
      clust = norm_clust_strict_weighted(W2, Knn = nk)
      W2 = Reduce(rbind, W2)
      assign = unlist(clust$clusters)
      if(n_cell > 20000){
        n_sample = min(2^16-1, floor(n_cell*0.2))
        index = sample(1:n_cell, n_sample)
        assign = assign[index]
      }
      consmat = consmat + clust2connect(assign)/n_res
    }
  }
  
  gc()
  
  obs_hclust = hclust(as.dist(1-consmat), method = "average")
  cop = cophenetic(obs_hclust)
  d1 = (1-consmat)[upper.tri(consmat, diag = FALSE)]
  d2 = as.matrix(cop)[upper.tri(consmat, diag = FALSE)]
  stability = cor(d1, d2)

  return(stability)
}



#' Calculate specificity score and select the \eqn{lambda}.
#' @description
#' Calculate specificity of each lambda and select the \eqn{lambda} that leads to the largest specificity score.
#'
#' @param object \code{scINSIGHT} object.
#' @param res_parallel List of result calculated by scINSIGHT.
#' @param LDA Regularization parameters to select (default \code{c(0.001, 0.01, 0.1, 1, 10)}).
#' @param ncores Number of cores to use for optimizing factorizations in parallel (default NULL).
#' @param out.dir Output directory of intermediate results (default NULL).
#' @param thre.niter Maximum number of block coordinate descent iterations to perform (default 500).
#' @param num.cores Number of cores to use for optimizing factorizations in parallel (default 1).
#' @param thre.delta Stop iteration when the reduction of objective function is less than delta (default 0.05).
#' @param topn Number of genes that have the largest loadings on the module (default 100).
#'
#' @useDynLib scINSIGHT
#'
#' @return \code{scINSIGHT} object with \eqn{W_1}, \eqn{W_2}, \eqn{H}, \eqn{V} and parameters slots set.
#'
#' @examples
#' \dontrun{
#' # Requires scINSIGHT object and results calculated by scINSIGHT
#' scINSIGHTex = select_LDA(scINSIGHTex, res_parallel)
#' }



select_LDA = function(
  object,
  res_parallel,
  LDA = c(0.001, 0.01, 0.1, 1, 10),
  ncores = NULL,
  out.dir = NULL,
  thre.niter = 500,
  thre.delta = 0.01,
  num.cores = 1,
  topn = 100){


  LDA = sort(LDA)

  cnt_list = object@norm.data
  cnt_list = lapply(cnt_list, function(x){t(x)})
  samples = names(cnt_list)
  L = length(cnt_list)
  uLabels = unique(object@condition)
  labels = sapply(object@condition, function(condition){which(uLabels==condition)-1})

  J = length(uLabels)

  bestK = object@parameters[["K"]]
  K_j = object@parameters[["K_j"]]
  lda0 = object@parameters[["lda"]]

  bseed_index = which.min(sapply(res_parallel, function(x) x$loss))
  bseed = as.numeric(str_sub(names(res_parallel)[[bseed_index]], 6))
  res0 = res_parallel[[bseed_index]]

  res_parallel = mclapply(LDA, function(lda){
    if(lda == lda0){
      res = res0
    }else{
      res = iNMF_BCD(cnt_list, labels, K_j, bestK, lda1 = lda, lda2 = lda, eps = 1e-5,
                     thre.niter, thre.delta, bseed, loop = TRUE)
      gc()
    }
    return(res)
  }, mc.cores = ncores)

  names(res_parallel) = sapply(LDA, function(lda){paste0("lda_",lda)})

  if(!is.null(out.dir)){
    if(str_sub(out.dir, -1)!="/"){
      out.dir = paste0(out.dir,"/")
    }
    dir.create(out.dir, recursive = TRUE)

    for(lda in LDA){
      res_name = paste0("lda_",lda)
      saveRDS(res_parallel[[res_name]], paste0(out.dir, "res-bestk", bestK, "-lda", lda, ".rds"))
    }
  }

  idx = match(object@condition, uLabels)
  specificity = sapply(res_parallel, function(res){
    H = res$H
    score = sapply(1:J, function(j){
      hmat = H[[j]]
      within = lapply(1:nrow(hmat), function(k){
        t = sort(hmat[k, ], decreasing = T)[topn]
        sapply(cnt_list[idx == j], function(x){
          mean(x[, hmat[k,]>t])
        })
      })
      within = unlist(within)
      between = lapply(1:nrow(hmat), function(k){
        t = sort(hmat[k, ], decreasing = T)[topn]
        sapply(cnt_list[idx != j], function(x){
          mean(x[, hmat[k,]>t])
        })
      })
      between = unlist(between)
      return(mean(within)/mean(between))
    })
    return(mean(score))
  })

  names(specificity) = sapply(LDA, function(lda){paste0("lda_",lda)})
  print(specificity)

  lda_index = which.max(specificity)
  object@parameters[["lda"]] = LDA[lda_index]
  object@parameters[["specificity"]] = specificity

  object@W_2 = res_parallel[[lda_index]]$W_2
  object@W_1 = res_parallel[[lda_index]]$W_1
  object@H = res_parallel[[lda_index]]$H
  object@V = res_parallel[[lda_index]]$V
  
  names(object@W_2) = samples
  names(object@W_1) = samples
  names(object@H) = uLabels

  clust = norm_clust_strict_weighted(object@W_2 , Knn=20)

  object@norm.W_2 = clust$W2
  object@clusters = clust$clusters
  
  names(object@clusters) = samples
  
  
  return(object)
}

