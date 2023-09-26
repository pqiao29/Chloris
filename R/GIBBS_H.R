update_H_k_chr2 <- function(signal_RDR, signal_BAF, RDR_k_chr, A_k_chr, D_k_chr, S,
                            Pi_k, Q_k, mu_est = NULL, sigma_est = NULL, theta_est = NULL) {
  U <- ifelse(signal_RDR, nrow(RDR_k_chr), nrow(A_k_chr))
  ####################################################################################################################
  #############################################   Forward    #########################################################
  ####################################################################################################################
  fwd <- array(NA, dim = c(U, S, S))
  log_liklihood <- 0

  for (u in 1:U) {
    tmp_pi <- Pi_k[u, ]
    Pt <- matrix(NA, S, S)
    ## observations in all cells of clones k at position u
    if (signal_RDR) RDR_k_u <- RDR_k_chr[u, , drop = F]
    if (signal_BAF) {
      keep_snp <- (D_k_chr[u, , drop = F] > 1)
      if (sum(keep_snp) > 0) {
        A_k_u <- A_k_chr[u, keep_snp, drop = F]
        D_k_u <- D_k_chr[u, keep_snp, drop = F]
      }
    }

    for (s in 1:S) {
      ### log emission density
      ### scaled s.t not affected by sample size
      Ps <- 0
      if (signal_RDR) Ps <- Ps + mean(dnorm(RDR_k_u, mu_est[s], sqrt(sigma_est[s]), log = TRUE))
      if (signal_BAF && (sum(keep_snp) > 0)) Ps <- Ps + mean(dbinom(A_k_u, D_k_u, prob = theta_est[s], log = TRUE))

      for (r in 1:S) {
        tmp_Q <- log(Q_k[r, s]) ## transition probability
        tmp_Pi <- log(tmp_pi[r])
        Pt[r, s] <- Ps + tmp_Q + tmp_Pi
      }
    }

    Mt <- max(Pt)
    Pt <- exp(Pt - Mt)
    tmp_nc <- sum(Pt)

    log_liklihood <- log_liklihood + log(tmp_nc) + Mt
    Pt <- Pt / tmp_nc

    Pi_k[u + 1, ] <- colSums(Pt)
    fwd[u, , ] <- Pt
  }

  ## Output: fwd, Pi_k, log_liklihood

  ####################################################################################################################
  #############################################   Backward    ########################################################
  ####################################################################################################################
  H_k <- rep(NA, U)
  H_k[U] <- which(rmultinom(1, 1, prob = Pi_k[U + 1, ]) == 1)

  for (u in (U - 1):1) {
    tmp_H_next <- H_k[u + 1]
    tmp_P <- fwd[u + 1, , tmp_H_next]
    tmp_P <- tmp_P / sum(tmp_P)
    H_k[u] <- which(rmultinom(1, 1, prob = tmp_P) == 1)
  }
  ## Output: H_k

  list(
    "log_liklihood" = log_liklihood,
    "Pi" = Pi_k, "H" = H_k
  )
}



update_states <- function(signal_RDR, signal_BAF, RDR, A, D, est, U, break_idx) {
  S <- est$S
  Pi_est <- est$Pi
  mu_est <- est$mu
  sigma_est <- est$sigma
  Q_est <- est$Q
  I_label <- est$I_label
  theta_est <- est$theta

  K <- est$K
  log_liklihood <- rep(0, K)

  for (k in 1:K) {
    cell_idx <- which(I_label == k)

    if (length(cell_idx) > 0) {
      if (signal_RDR) {
        RDR_k <- RDR[, cell_idx, drop = F]
      } else {
        RDR_k <- NULL
      }

      if (signal_BAF) {
        A_k <- A[, cell_idx, drop = F]
        D_k <- D[, cell_idx, drop = F]
      } else {
        A_k <- D_k <- NULL
      }

      if (is.null(break_idx)) {
        update_k <- update_H_k_chr2(signal_RDR, signal_BAF, RDR_k,
          A_k_chr = A_k, D_k_chr = D_k, S,
          Pi_k = Pi_est[[k]], Q_k = Q_est[[k]], mu_est, sigma_est, theta_est
        )
        log_liklihood[k] <- update_k$log_liklihood
        est$Pi[[k]] <- update_k$Pi
        est$H[k, ] <- update_k$H
      } else {
        tmp_break_idx <- c(1, break_idx, U + 1)
        logliklihood_k <- 0

        for (chr_idx in 1:(length(tmp_break_idx) - 1)) {
          gene_idx <- tmp_break_idx[chr_idx]:(tmp_break_idx[chr_idx + 1] - 1)

          if (signal_RDR) {
            RDR_k_chr <- RDR_k[gene_idx, , drop = F]
          } else {
            RDR_k_chr <- NULL
          }

          if (signal_BAF) {
            A_k_chr <- A_k[gene_idx, , drop = F]
            D_k_chr <- D_k[gene_idx, , drop = F]
          } else {
            A_k_chr <- D_k_chr <- NULL
          }

          update_chr <- update_H_k_chr2(signal_RDR, signal_BAF, RDR_k_chr, A_k_chr, D_k_chr, S,
            Pi_k = Pi_est[[k]][[chr_idx]], Q_k = Q_est[[k]], mu_est, sigma_est, theta_est
          )

          logliklihood_k <- logliklihood_k + update_chr$log_liklihood
          est$H[[k]][[chr_idx]] <- update_chr$H
          est$Pi[[k]][[chr_idx]] <- update_chr$Pi
        }

        log_liklihood[k] <- logliklihood_k
      }
    }
  }

  est$log_liklihood <- sum(log_liklihood)

  if (!is.null(break_idx)) {
    H_array <- matrix(NA, K, U)
    for (k in 1:K) {
      H_array[k, ] <- unlist(est$H[[k]])
    }
    est$H_array <- H_array
  }

  return(est)
}
