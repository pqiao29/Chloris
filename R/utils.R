get_cell_level_states <- function(I, states) {
  if (ncol(I) != nrow(states)) {
    return(NULL)
  }
  N <- nrow(I)
  K <- ncol(I)
  U <- ncol(states)

  ret <- matrix(NA, U, N)
  for (k in 1:K) {
    cell_idx <- (I[, k] == 1)
    cluster_states <- states[k, ]
    ret[, cell_idx] <- cluster_states
  }
  ret
}

my_v_match <- function(needle, haystack) {
  a <- haystack[-length(haystack)] == needle[1]
  b <- haystack[-1] == needle[2]
  sum(a & b)
}

add_matrix <- function(x) Reduce("+", x)

lik <- function(RDR, I, H, mu, sigma, A, D, theta) {
  ret <- 0

  cell_level_states <- get_cell_level_states(I, H)
  S <- sort(unique(c(cell_level_states)))

  for (s in S) {
    if (!is.null(RDR)) ret <- ret + sum(dnorm(RDR[cell_level_states == s], mu[s], sqrt(sigma[s]), log = TRUE))
    if (!is.null(A)) ret <- ret + sum(dbinom(A[cell_level_states == s], D[cell_level_states == s], theta[s], log = TRUE))
  }

  ret
}

lik_H <- function(H, Q) {
  ret <- 0
  K <- nrow(H)
  S <- nrow(Q[[1]])

  for (k in 1:nrow(H)) {
    for (s1 in 1:S) {
      for (s2 in 1:S) {
        ret <- ret + log(Q[[k]][s1, s2]) * my_v_match(c(s1, s2), H[k, ])
      }
    }
  }

  return(ret)
}

#' Basic preprocess
#'
#' @param counts Raw gene count matrix, gene*cell
#' @param ref_cells A vector of reference or normal cell index.
#' @param window_smooth window size for window smoothing, NULL then no window smoothing applied
#' @param recentering Recenter RDR for each cell independently. 
#' @param psudo_count Psudo count in log2 transform
#' @export
preprocess <- function(counts, ref_cells, window_smooth = NULL, recentering = FALSE, psudo_count = 0) {
  ## Equalize library size
  counts <- t(t(counts) / colSums(counts))

  ## Truncate zero's before taking log
  if (psudo_count == 0) counts[counts == 0] <- min(counts[counts > 0])

  ## log ratio
  ref <- rowMeans(counts[, ref_cells])
  RDR <- log2(counts + psudo_count) - log2(ref + psudo_count)

  ## recentering
  if (recentering) {
    col_median <- apply(RDR, 2, function(x) {
      median(x, na.rm = TRUE)
    })
    RDR <- t(apply(RDR, 1, "-", col_median))
  }

  return(RDR)
}
