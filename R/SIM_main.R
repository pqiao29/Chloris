#' Clonal CNA simulation
#'
#' Create scRNA-seq count data, RDR or BAF for testing \code{Chloris}
#'
#' @param K The number of clusters.
#' @param N The number of cells.
#' @param U The number of genes.
#' @param expr TRUE if gene count is to be generated.
#' @param RDR TRUE if RDR is to be generated. If expr == TRUE, RDR will be obtained from expr, otherwise generated from Normal distribution.
#' @param BAF TRUE if BAF is to be generated.
#' @param RDR_var The variance of RDR if generated from Normal distribution, needed iff expr == FALSE and RDR == TRUE.
#' @param norm_clone_prob The proportion of normal cluster among all cells.
#' @param CNV_overlap TRUE if each cluster has an identical interval with another cluster, which increases the similarity between different clusters.
#' @param RDR_outlier_cnt An integer. The number of outlier cells with random CN profiles.
#' @param BAF_missing_percent The percentage of NA in BAF data, effective only if BAF == TRUE.
#' @param BAF_Dmax The largest number of total count in BAF. The smaller it is the weaker the BAF signal is.
#' @return Simulated data as well as cell cluster labels and CN state profile for each cluster
#'
#' @export
#'
#' @examples 
#' \dontrun{
#' sims <- get_sim_data(K = 4, N = 100, U = 200)
#' names(sims)
#' plot_inout(sims$RDR, list(sims$cluster_true), sims$states_true, state_mean = c(-1, 0, 0.5, 1.5))
#' }
get_sim_data <- function(K, N, U, S = 4, 
                         expr = TRUE, RDR = TRUE, BAF = TRUE,
                         RDR_var = NULL,
                         norm_clone_prob = 0.2, CNV_overlap = TRUE,
                         RDR_outlier_cnt = 0,
                         BAF_missing_percent = 0.9, BAF_Dmax = 10, 
                         Q_to_neutral = 2){

    ### =================================== fixed parameters =============================================
    RDR_levels <- c(0.5, 1, 1.5, 2)[1:S]
    theta_pars <- list("1" = c(10, 200), "2" = c(50, 50), "3" = c(20, 100), "4" = c(20, 150))[1:S]
    Q_diag <- seq(100, 25, length.out = K)
    clone_prob <- c(norm_clone_prob, rep((1 - norm_clone_prob)/(K - 1), K - 1))

    ### =================================== true latents =============================================
    Q_true <- get_true_transition(K = K, S = S, diag = Q_diag, to_neutral = Q_to_neutral, neutral_idx = 2)
    states_true <- generate_cluster_states(K = K, U, Q_true, neutral_idx = 2, overlap = CNV_overlap)

    cluster_true <- sim_clusters(N, K, clone_prob)
    I_true <- apply(cluster_true, 1, which.max)
    cell_level_states <- get_cell_level_states(cluster_true, states_true)

    ### =================================== generate data =============================================
    if(expr){
        if(RDR_outlier_cnt >= 1){
            tmp_cell_level_states <- cbind(cell_level_states, matrix(sample(1:S, U*RDR_outlier_cnt, replace = T), U, RDR_outlier_cnt))
            RDR_data_full <- sim_RDR(tmp_cell_level_states, log2(RDR_levels), from_Splat = TRUE, var_true = NULL, RDR_outlier_cnt = RDR_outlier_cnt)
        }else{
            RDR_data_full <- sim_RDR(cell_level_states, log2(RDR_levels), from_Splat = TRUE, var_true = NULL, RDR_outlier_cnt = RDR_outlier_cnt)
        }

        # keep_gene <- RDR_data_full$keep_gene
        # if(sum(keep_gene) < U){
        #     U <- sum(keep_gene)
        #     states_true <- states_true[, keep_gene]
        #     cell_level_states <- RDR_data_full$cell_level_states
        # }
    }else{
        if(RDR){
            if(is.null(RDR_var)) stop("If (expr == FALSE && RDR == TRUE), RDR_var must be provided.")
            var_invgamma = 1e-06
            var_true <- matrix(NA, U, N)
            for(s in 1:S){
                var_true[cell_level_states == s] <- invgamma::rinvgamma(sum(cell_level_states == s), 
                                                                        shape = (RDR_var^2)/var_invgamma + 2, 
                                                                        rate = RDR_var*(1 + (RDR_var^2)/var_invgamma))
            } 
            RDR_data <- sim_RDR(cell_level_states, log2(RDR_levels), from_Splat = FALSE, var_true = var_true)$RDR
            
            ## add outlier
            if(RDR_outlier_cnt >= 1){
                for(i in 1:RDR_outlier_cnt){
                    outlier_H <- sample(1:S, U, replace = T)
                    outlier_Y <- rep(NA, U)
                    for(s in 1:S){
                        outlier_Y[outlier_H == s] <- rnorm(sum(outlier_H == s), RDR_levels[s], sd = 4*RDR_var)
                    }
                    RDR_data <- cbind(RDR_data, outlier_Y)
                }
            }
            N_all <- N + RDR_outlier_cnt
        }
    }

    if(BAF){
        theta_true <- matrix(NA, U, N)
        for(s in 1:S){
            theta_true[cell_level_states == s] <- rbeta(sum(cell_level_states == s), theta_pars[[s]][1], theta_pars[[s]][2])
        }
        BAF_data <- sim_BAF(cell_level_states, theta_true, Md = BAF_Dmax, BAF_missing_percent)

        if(RDR_outlier_cnt >= 1){
            BAF_data$A <- cbind(BAF_data$A, matrix(0, U, RDR_outlier_cnt))
            BAF_data$D <- cbind(BAF_data$D, matrix(0, U, RDR_outlier_cnt))
        }

    }

    ### =================================== return =============================================
    ret <- list("N" = N, "U" = U,
                "cluster_true" = I_true, "states_true" = states_true, "cell_level_states" = cell_level_states)

    if(BAF){
        ret$A <- BAF_data$A
        ret$D <- BAF_data$D
        ret$BAF_missing_percent <- BAF_missing_percent
        ret$theta_true <- theta_true
    }

    if(expr){
        ret$expr <- RDR_data_full#$raw
        ret$RDR_outlier_cnt <- RDR_outlier_cnt
        #if(RDR) ret$RDR <- RDR_data_full$RDR
    }else{
        if(RDR){
            ret$RDR <- RDR_data
            ret$RDR_var <- RDR_var
            ret$RDR_outlier_cnt <- RDR_outlier_cnt
        }
    }


    return(ret)
}



