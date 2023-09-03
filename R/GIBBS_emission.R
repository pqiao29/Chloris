#' @import invgamma
#' @keywords internal

update_emission <- function(signal_RDR, signal_BAF, RDR, A, D, est, priors){

    K <- est$K
    S <- est$S
    I_est <- est$I
    I_label <- est$I_label

    if(is.list(est$H)){
        H_est <- est$H_array
    }else{
        H_est <- est$H
    }

    ####################################################################################################################
    ###########################################      RDR       #########################################################
    ####################################################################################################################
    if(signal_RDR){

        # ======================================== mean ========================================
        mu_0 <- priors$mu
        mu_est <- est$mu
        sigma_0 <- priors$sigma
        sigma_est <- est$sigma

        for(s in 1:S){

            tmp_g <- 0
            tmp_Ns <- 0
            for(k in 1:K){

                idx_cells <- which(I_label == k)
                H_est_k <- H_est[k, ]
                idx_pos <- which(H_est_k == s)

                RDR_H <- RDR[idx_pos, idx_cells]

                tmp_g <- tmp_g + sum(RDR_H)
                tmp_Ns <- tmp_Ns + length(RDR_H) #length(idx_cells)*length(idx_pos)
            }

            tmp_eps_0 <- mu_0[2, s]
            tmp_lmd_0 <- mu_0[1, s]
            tmp_sigma_est <- sigma_est[s]

            tmp_lmd_est <- tmp_eps_0*tmp_g + tmp_sigma_est*tmp_lmd_0
            tmp_lmd_est <- tmp_lmd_est/(tmp_eps_0*tmp_Ns + tmp_sigma_est)

            tmp_eps_est <- (tmp_Ns/tmp_sigma_est) + (1/tmp_eps_0)
            tmp_eps_est <- 1/tmp_eps_est

            mu_est[s] <- rnorm(1, tmp_lmd_est, sqrt(tmp_eps_est))
        }## Output: latent_mu_est
        est$mu <- mu_est

        # ======================================== var ========================================
        tmp_nv1_0 <- sigma_0[1]
        tmp_nv2_0 <- sigma_0[2]
        for(s in 1:S){

            tmp_mu_est <- mu_est[s]
            tmp_g <- 0
            tmp_Ns <- 0
            for(k in 1:K){
                idx_cells <- which(I_label == k)
                H_est_k <- H_est[k, ]
                idx_pos <- which(H_est_k == s)

                RDR_H <- RDR[idx_pos, idx_cells]

                tmp_g <- tmp_g + sum((RDR_H - tmp_mu_est)^2)
                tmp_Ns <- tmp_Ns + length(RDR_H) #length(idx_cells)*length(idx_pos)
            }

            tmp_nv1_est <- tmp_nv1_0 + tmp_Ns/2
            tmp_nv2_est <- tmp_nv2_0 + tmp_g/2

            sigma_est[s] <- invgamma::rinvgamma(1, shape = tmp_nv1_est, rate = tmp_nv2_est)
        }## Output: latent_sigma_est
        est$sigma <- sigma_est

    }


    ####################################################################################################################
    ###############################################   BAF    ###########################################################
    ####################################################################################################################
    if(signal_BAF){

        theta_est <- est$theta
        theta_0 <- priors$theta

        for(s in 1:S){

            tmp_a_0 <- theta_0[1, s]
            tmp_b_0 <- theta_0[2, s]

            tmp_g1 <- 0
            tmp_g2 <- 0
            for(k in 1:K){
                idx_cells <- which(I_label == k)
                H_est_k <- H_est[k, ]
                idx_pos <- which(H_est_k == s)

                A_H <- A[idx_pos, idx_cells]
                D_H <- D[idx_pos, idx_cells]

                tmp_g1 <- tmp_g1 + sum(A_H)
                tmp_g2 <- tmp_g2 + sum(D_H - A_H)
            }

            tmp_a_est <- tmp_a_0 + tmp_g1
            tmp_b_est <- tmp_b_0 + tmp_g2

            theta_est[s] <- rbeta(1, tmp_a_est, tmp_b_est)
        }## Output: theta_est
        est$theta <- theta_est
    }

    return(est)
}

