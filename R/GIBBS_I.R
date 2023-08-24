update_clustering <- function(signal_RDR, signal_BAF, RDR, A, D, est, priors, save_prob){

    N <- ifelse(signal_RDR, ncol(RDR), ncol(A))
    K <- est$K
    S <- est$S

    if(signal_RDR){
        mu_est <- est$mu
        sigma_est <- est$sigma
    }
    if(signal_BAF) theta_est <- est$theta

    if(is.list(est$H)){
        H_est <- est$H_array
    }else{
        H_est <- est$H
    }

    I_est <- matrix(0, N, K)
    if(save_prob) phi_record <- matrix(NA, N, K)

    for(i in 1:N){

        I_i <- rep(0, K)

        for(k in 1:K){

            H_est_k <- H_est[k, ]

            for(s in 1:S){

                pos_idx <- which(H_est_k == s) ## Extract u s.t. X_{ku} = s

                if(signal_RDR){

                    tmp_mu_est <- mu_est[s]
                    tmp_sigma_est <- sigma_est[s]

                    sub_RDR <- RDR[pos_idx, i]
                    I_i[k] <- I_i[k] + sum(dnorm(sub_RDR, mean = tmp_mu_est, sd = sqrt(tmp_sigma_est), log = TRUE))
                }

                if(signal_BAF){

                    tmp_theta_est <- theta_est[s]

                    keep_snp <- (D[pos_idx, i] > 1)
                    if(sum(keep_snp) > 0){
                        sub_A <- A[pos_idx, i][keep_snp]
                        sub_D <- D[pos_idx, i][keep_snp]
                        I_i[k] <- I_i[k] + sum(dbinom(sub_A, sub_D, prob = tmp_theta_est, log = TRUE))
                    }
                }

            }
        }

        Mi <- max(I_i)
        I_i <- I_i - Mi
        I_i <- exp(I_i)/sum(exp(I_i))
        I_est[i, ] <- t(rmultinom(1, 1, I_i))
       # I_est[i, which.max(I_i)] <- 1
        if(save_prob) phi_record[i, ] <- I_i
    }

    est$I_label <- apply(I_est, 1, which.max)
    est$I <- I_est
    if(save_prob) est$I_prob <- phi_record
    return(est)
}

