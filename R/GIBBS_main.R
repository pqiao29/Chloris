Gibbs_main <- function(signal_RDR, signal_BAF, RDR, A, D,
                       priors, init, break_idx = NULL,
                       burnin_tol, Gibbs_tol, cluster_shrink_tol, min_cluster_size) {
  S <- priors$S
  K <- priors$K

  N <- ifelse(signal_RDR, ncol(RDR), ncol(A))
  U <- ifelse(signal_RDR, nrow(RDR), nrow(A))

  est <- Gibbs_init(signal_RDR, signal_BAF, RDR, A, D, priors, break_idx = break_idx, U, clust_method = init)

  if (is.null(cluster_shrink_tol)) {
    cat("Shrinkage disabled.\n")
    pb <- txtProgressBar(min = 0, max = burnin_tol, style = 3)
    
    for (iter in 1:burnin_tol) {
      setTxtProgressBar(pb, iter)

      ### main update
      est <- update_states(signal_RDR, signal_BAF, RDR, A, D, est, U, break_idx = break_idx)
      est <- update_clustering(signal_RDR, signal_BAF, RDR, A, D, est, priors, save_prob = FALSE)
      est <- update_transition(est, priors)
      est <- update_emission(signal_RDR, signal_BAF, RDR, A, D, est, priors)
    }
  } else {
    cat("Shrinkage enabled.\n")
    pb <- txtProgressBar(min = 0, max = burnin_tol, style = 3)
    
    empty_cluster_tracker <- NULL
    iter <- 1
    while (iter <= burnin_tol) {
      setTxtProgressBar(pb, iter)
      
      ### Main update
      est <- update_states(signal_RDR, signal_BAF, RDR, A, D, est, U, break_idx = break_idx)
      est <- update_clustering(signal_RDR, signal_BAF, RDR, A, D, est, priors, save_prob = FALSE)
      est <- update_transition(est, priors)
      est <- update_emission(signal_RDR, signal_BAF, RDR, A, D, est, priors)

      cluster_size <- apply(est$I, 2, sum)

      if (length(empty_cluster_tracker) >= cluster_shrink_tol) {
        ### Shrink K in prior
        cat("\nK shrinks from ", K)
        K <- K - min(empty_cluster_tracker)
        priors$I <- rep(1 / K, K)
        priors$K <- K
        ### restart
        cat(" to ", K, "\n")
        iter <- 0
        empty_cluster_tracker <- NULL
        pb <- txtProgressBar(min = 0, max = burnin_tol, style = 3)
        
        est <- Gibbs_init(signal_RDR, signal_BAF, RDR, A, D, priors, break_idx = break_idx, U, clust_method = init)
      } else {
        if (min(cluster_size) <= min_cluster_size) {
          empty_cluster_tracker <- c(empty_cluster_tracker, sum(cluster_size <= min_cluster_size))
        }
      }

      iter <- iter + 1
    }
    cat("\n")
  }

  N <- ifelse(signal_RDR, ncol(RDR), ncol(A))
  cluster_record <- matrix(NA, Gibbs_tol, N)
  state_record <- array(NA, dim = c(Gibbs_tol, est$K, U))
  theta_record <- array(NA, dim = c(Gibbs_tol, S, signal_RDR * 2 + signal_BAF))
  loglik_record <- rep(NA, Gibbs_tol)
  # steven_input_cluster <- array(NA, dim = c(Gibbs_tol, N, est$K))
  # all_par_record_cluster <- array(NA, dim = c(Gibbs_tol, est$K, N + U))
  Q_record <- list()
  Pi_record <- list()

  pb <- txtProgressBar(min = 0, max = Gibbs_tol, style = 3)
  for (Gibbs_iter in 1:Gibbs_tol) {
    cat("\033[0;31m")
    setTxtProgressBar(pb, Gibbs_iter)
    cat("\033[0m")
    
    est <- update_states(signal_RDR, signal_BAF, RDR, A, D, est, U, break_idx = break_idx)
    est <- update_clustering(signal_RDR, signal_BAF, RDR, A, D, est, priors, save_prob = FALSE)
    est <- update_transition(est, priors)
    est <- update_emission(signal_RDR, signal_BAF, RDR, A, D, est, priors)

    if (!is.null(break_idx)) {
      H_est <- est$H_array
    } else {
      H_est <- est$H
    }

    ################################  Simple label switch using mu ##################################
    if (signal_RDR) { ## full model or RDR model
      tmp_par_order <- sort(est$mu, index.return = T)$ix
      if (!identical(tmp_par_order, 1:S)) { ## update H_est
        est$mu <- est$mu[tmp_par_order]
        est$sigma <- est$sigma[tmp_par_order]
        if (signal_BAF) est$theta <- est$theta[tmp_par_order]
        est$Q <- lapply(est$Q, function(x) x[tmp_par_order, ][, tmp_par_order])
        tmp_H <- H_est
        for (s in 1:S) {
          tmp_H[H_est == tmp_par_order[s]] <- s
        }
        H_est <- tmp_H
      }
    } else { # BAF model
      tmp_par_order <- sort(est$theta, index.return = T)$ix
      if (!identical(tmp_par_order, 1:S)) { ## update H_est
        est$theta <- est$theta[tmp_par_order]
        est$Q <- lapply(est$Q, function(x) x[tmp_par_order, ][, tmp_par_order])
        tmp_H <- H_est
        for (s in 1:S) {
          tmp_H[H_est == tmp_par_order[s]] <- s
        }
        H_est <- tmp_H
      }
    }

    ## update est$H and est$H_array
    if (!is.null(break_idx)) {
      est$H_array <- H_est
      updated_H <- vector("list", est$K)
      full_break_idx <- c(1, break_idx, U + 1)
      for (k in 1:est$K) {
        for (chr_idx in 1:(length(full_break_idx) - 1)) {
          gene_idx <- full_break_idx[chr_idx]:(full_break_idx[chr_idx + 1] - 1)
          updated_H[[k]][[chr_idx]] <- est$H_array[k, gene_idx]
        }
      }
      est$H <- updated_H
    } else {
      est$H <- H_est
    }
    #####################################  End of sweep wrap up #####################################
    cluster_record[Gibbs_iter, ] <- est$I_label
    state_record[Gibbs_iter, , ] <- H_est
    theta_record[Gibbs_iter, , ] <- cbind(est$mu, est$sigma, est$theta)
    loglik_record[Gibbs_iter] <- lik(RDR, est$I, H_est, est$mu, est$sigma, A, D, est$theta) + lik_H(H_est, est$Q)
    Q_record[[Gibbs_iter]] <- est$Q
    Pi_record[[Gibbs_iter]] <- est$Pi
  }

  list(
    "cluster_record" = cluster_record, "state_record" = state_record, "par_record" = theta_record, "est" = est,
    "loglik_record" = loglik_record, "Q_record" = Q_record, "Pi_record" = Pi_record, "K" = K
  )
}
