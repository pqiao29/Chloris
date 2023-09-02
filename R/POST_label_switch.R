#' Label Switching correction
#'
#' Wrapper of the Equivalence Classes Representatives (ECR) algorithm in \href{https://cran.r-project.org/web/packages/label.switching/index.html}{\code{label.switching}}
#'
#' @param cluster_record T by N matrix of MCMC clustering samples,
#' where T is the total number of Gibbs iteration and N is the number of cells being clustered.
#' @param loglik_record likelihood of each Gibbs iteration.
#' @return The order of cluster labels in each MCMC iteration.
#' @keywords internal

ECR <- function(cluster_record, loglik_record) {
  cluster_record <- apply(cluster_record, 1, function(x) {
    as.numeric(factor(x, levels = unique(x), labels = 1:length(unique(x))))
  })
  cluster_record <- t(cluster_record)

  flag <- which.max(loglik_record)
  flag_I <- cluster_record[flag, ]

  invisible(capture.output(
    ls_ECR <- label.switching::label.switching(method = "ECR", zpivot = flag_I, z = cluster_record, K = max(cluster_record))
  ))

  return(ls_ECR)
}

#' Reorder cluster labels after correction
#'
#' Align the cluster labels as well as all cluster-specific parameters
#' according to the label order returned by \code{ECR()}
#' for all MCMC samples.
#'
#' @param res output from \code{ECR()}.
#' @param cluster_record same input as in \code{ECR()}.
#' @param state_record A T by K by U array of clonal copy number states, where T is the number of MCMC iterations, K is the number of clusters and U is the number of genes.
#' @param Q_record A list of T MCMC samples of transition matrices, each sample contains K S by S matrices, where K is the number of clusters and S is the number of states.
#' @param Pi_record State probability in the same format as \code{Q_record}.
#' @keywords internal
align_MC_samples <- function(res, cluster_record, state_record, Q_record = NULL, Pi_record = NULL) {
  Gibbs_tol <- dim(state_record)[1]
  state_record_aligned <- array(NA, dim = dim(state_record))
  if (!is.null(Q_record)) Q_record_aligned <- list()
  if (!is.null(Pi_record)) Pi_record_aligned <- list()

  for (iter in 1:Gibbs_tol) {
    if (!is.null(Q_record)) tmp_Q_iter <- list()
    if (!is.null(Pi_record)) tmp_Pi_iter <- list()

    for (k in sort(unique(c(res$clusters)))) {
      I_iter <- cluster_record[iter, ]

      tmp_cluster_overlap <- table(I_iter[res$clusters == k])
      k_overlap <- names(tmp_cluster_overlap)
      prob_overlap <- tmp_cluster_overlap / table(I_iter)[k_overlap]

      if (any(prob_overlap == 1)) k_selected <- c(names(prob_overlap == 1)) else k_selected <- names(prob_overlap)[which.max(prob_overlap)]

      k_selected <- as.numeric(k_selected)
      if (length(k_selected) == 1) {
        state_record_aligned[iter, k, ] <- state_record[iter, k_selected, ]
        if (!is.null(Q_record)) tmp_Q_iter[[k]] <- Q_record[[iter]][[k_selected]]
        if (!is.null(Pi_record)) tmp_Pi_iter[[k]] <- Pi_record[[iter]][[k_selected]]
      } else {
        states_concensus_k <- apply(state_record[iter, k_selected, ], 2, function(x) {
          tmp <- table(x)
          tmp_idx <- which(tmp == max(tmp))
          if (length(tmp_idx) == 1) ret <- tmp_idx else ret <- sample(tmp_idx, 1)
          as.integer(names(ret))
        })
        state_record_aligned[iter, k, ] <- states_concensus_k
        if (!is.null(Q_record)) tmp_Q_iter[[k]] <- add_matrix(Q_record[[iter]][k_selected]) / length(k_selected)
        if (!is.null(Pi_record)) tmp_Pi_iter[[k]] <- add_matrix(Pi_record[[iter]][k_selected]) / length(k_selected)
      }
    }

    if (!is.null(Q_record)) Q_record_aligned[[iter]] <- tmp_Q_iter
    if (!is.null(Pi_record)) Pi_record_aligned[[iter]] <- tmp_Pi_iter
  }

  res_K <- length(unique(c(res$clusters)))
  res_cluster <- as.numeric(factor(res$clusters, levels = sort(unique(c(res$clusters))), labels = 1:res_K))
  state_record_aligned <- state_record_aligned[, sort(unique(c(res$clusters))), ]

  U <- dim(state_record)[3]
  res_H_aligned <- matrix(NA, res_K, U)
  for (k in 1:res_K) {
    for (u in 1:U) {
      tmp <- table(state_record_aligned[, k, u])
      res_H_aligned[k, u] <- as.numeric(names(tmp)[which.max(tmp)])
    }
  }

  if (!is.null(Q_record)) {
    res_Q_aligned <- replicate(res_K, Q_record[[1]][[1]], simplify = FALSE)
    for (k in 1:res_K) {
      res_Q_aligned[[k]] <- add_matrix(sapply(Q_record_aligned, "[[", sort(unique(c(res$clusters)))[k], simplify = FALSE)) / Gibbs_tol
    }
  }

  if (!is.null(Pi_record)) {
    res_Pi_aligned <- replicate(res_K, Pi_record[[1]][[1]], simplify = FALSE)
    for (k in 1:res_K) {
      res_Pi_aligned[[k]] <- add_matrix(sapply(Pi_record_aligned, "[[", sort(unique(c(res$clusters)))[k], simplify = FALSE)) / Gibbs_tol
    }
  }


  ret <- list("cluster_est" = res_cluster, "state_est" = res_H_aligned)
  if (!is.null(Q_record)) ret$"Q_est" <- res_Q_aligned
  if (!is.null(Pi_record)) ret$"Pi_est" <- res_Pi_aligned
  return(ret)
}
