## Assumes dimensions of A, D, RDR have been checked
## Assumes the relationship of between signal and input data has been checked
Gibbs_init <- function(signal_RDR, signal_BAF, RDR, A, D,
                       priors, break_idx, U, clust_method = "hclust") {
  ret <- list()
  N <- ifelse(signal_RDR, ncol(RDR), ncol(A))
  S <- priors$S
  K <- priors$K
  neutral_idx <- priors$neutral_idx

  if (signal_RDR) {
    if (clust_method == "hclust") {
      dist_RDR <- dist(t(RDR))
      if (sum(is.na(dist_RDR)) == 0) {
        hc <- hclust(dist_RDR, method = "ward.D2")
        I_label_init <- cutree(hc, k = K)
      } else {
        warning("NA exists in distance of RDR, starting with random initialization! \n")
        clust_method <- "random"
      }
    }
    if (clust_method == "kmeans") hc <- kmeans(t(RDR), centers = K)
  } else {
    if (clust_method == "hclust") {
      dist_BAF <- dist(t(A / D))
      if (sum(is.na(dist_BAF)) == 0) {
        hc <- hclust(dist_BAF, method = "ward.D2")
        I_label_init <- cutree(hc, k = K)
      } else {
        warning("NA exists in distance of BAF, starting with random initialization! \n")
        clust_method <- "random"
      }
    }
    if (clust_method == "kmeans") hc <- kmeans(t(A / D), centers = K)
  }

  if (clust_method == "random") I_label_init <- apply(rmultinom(N, 1, rep(1 / K, K)), 2, which.max)


  tmp_I <- matrix(0, N, K)
  for (i in 1:N) {
    tmp_I[i, I_label_init[i]] <- 1
  }
  ret$I <- tmp_I
  ret$I_label <- I_label_init

  #### Marginal distribution parameter: mu, sigma, theta
  if (signal_RDR) {
    ret$mu <- priors$mu[1, ]
    ret$sigma <- rep(0.1, S)
  }
  if (signal_BAF) ret$theta <- apply(priors$theta, 2, function(x) rbeta(1, x[1], x[2]))

  #### Transition matrix: Q
  ret$Q <- replicate(K, priors$Q, simplify = FALSE)

  #### State: Pi and H
  init_pi <- rep(NA, S) ## most likely to start at neutral
  init_pi[neutral_idx] <- 0.9
  init_pi[-neutral_idx] <- 0.1 / (S - 1)

  if (is.null(break_idx)) {
    H_est <- rep(NA, U)
    Pi_k <- matrix(0, U + 1, S)
    Pi_k[1, ] <- init_pi
    ret$H <- matrix(NA, K, U)
  } else {
    tmp_break_idx <- c(1, break_idx, U + 1)
    H_k <- list()
    Pi_k <- list()
    for (chr_idx in 1:(length(tmp_break_idx) - 1)) {
      chr_length <- tmp_break_idx[chr_idx + 1] - tmp_break_idx[chr_idx]
      H_k[[chr_idx]] <- rep(NA, chr_length)

      tmp_Pi <- matrix(0, chr_length + 1, S)
      tmp_Pi[1, ] <- init_pi
      Pi_k[[chr_idx]] <- tmp_Pi
    }
    ret$H <- rep(list(H_k), K)
  }

  ret$Pi <- rep(list(Pi_k), K)

  ### index size
  ret$K <- K
  ret$S <- S

  ### label switching signal
  ret$cluster_switch_signal <- FALSE
  ret$states_switch_signal <- FALSE

  ret
}
