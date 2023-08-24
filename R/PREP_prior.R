#' Generate prior parameters
#'
#' @param K The maximum number of clusters.
#' @param S The number of states.
#' @param signal_RDR TRUE if prior for mu and sigma is required.
#' @param signal_BAF TRUE if prior for beta is required.
#' @param mu The prior distribution of mean of RDR for each state is N(\code{mu}, 0.1).
#' @param sigma_shape The prior distribution of the stand error of RDR for each state is Inv-Gamma(\code{sigma_shape}, 1).
#' @param beta_shape1 The prior distribution of BAF related parameter is Beta(\code{beta_shape1}, \code{beta_shape2}).
#' @param beta_shape2 See \code{beta_shape1}
#' @param Q_diag The prior of transition probability is a \code{S}\eqn{\times}\code{S} matrix, with off-diagonal entries 1 and diagonal elements \code{Q_diag}.
#' @param neutral_idx The index of the neutral state.
#' @param to_neutral Entries in the prior of transition probability that represent P(neutral state | non-neutral state).
#' @return A list of prior parameters for MCMC sampling.
#' @keywords internal

get_priors <- function(K, S, signal_RDR, signal_BAF,
                       mu = log2(c(0.5, 1, 1.5, 2)), sigma_shape = 10, beta_shape1 = 1, beta_shape2 = 5,
                       Q_diag = 10, neutral_idx = 2, to_neutral = 1.5){


    ### RDR
    if(signal_RDR){

        if(length(mu) != S) stop("Exception from prior setup:\nThe length of mu should be the same as S.")

        par_mu_0 <- matrix(.1, 2, S) ## sd in prior of mu: 0.1
        par_mu_0[1, ] <- mu
        par_sigma_0 <- c(sigma_shape, 1)     ## Inv-Gamma parameters
    }

    ### BAF
    if(signal_BAF){  ### Beta parameters
        par_theta_0 <- matrix(beta_shape1, 2, S)
        par_theta_0[2, -neutral_idx] <- beta_shape1*beta_shape2
    }

    ### Transition
    par_Q_gamma_0 <- matrix(1, S, S)
    par_Q_gamma_0[, neutral_idx] <- to_neutral
    diag(par_Q_gamma_0) <- Q_diag
    ### clustering
    par_I_phi_0 <- rep(1/K, K)

    ret <- list("Q" = par_Q_gamma_0, "I" = par_I_phi_0, "S" = S, "K" = K, "neutral_idx" = neutral_idx)
    if(signal_BAF) ret$theta <- par_theta_0
    if(signal_RDR){
        ret$sigma <- par_sigma_0
        ret$mu <- par_mu_0
    }

    ret
}
