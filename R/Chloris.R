#' R package 'Chloris'
#' @name Chloris-description
#' @importFrom stats cutree dbinom dist dnorm hclust kmeans median rbeta rbinom rchisq rlnorm rmultinom rnbinom rnorm sd
#' @importFrom utils capture.output
NULL


#' Main function
#'
#' 1. Cluster cells into clones,
#' 2. Infer copy number state profile for each clone,
#' 3. Identify the optimal number of clones.
#'
#' @param RDR A matrix of Relative Depth Ratio (RDR), where each row and column represents a gene and a cell respectively.
#' @param A A matrix of count from the alternative allele in a heterozygous SNP, same dimension as \code{RDR}.
#' @param D A matrix of total count from both alleles in a SNP, same dimension as \code{RDR}.
#' @param break_idx The index of the first gene in each chromosome, so that genes from different chromosomes are treated as independent.
#'                  If \code{break_idx = NULL}, it is assumed that there is spatial dependency between all genes with adjacent indices.
#' @param K An integer that is likely to be larger than the true number of clusters.
#'          The resulting optimal number of clusters is <= K.
#' @param min_cluster_size Cluster size constraint. All resulting clusters need to contain at least \code{min_cluster_size} cells.
#' @param cluster_shrink_tol If a cluster contains less than \code{min_cluster_size} cells for \code{cluster_shrink_tol} consecutive MCMC iterations,
#'                           this cluster will be removed and \code{K} will be decreased accordingly.
#'                           If \code{cluster_shrink_tol == NULL}, then no cluster size constraint is imposed so the function always return \code{K} cluster.
#' @param burnin_tol The number of Gibbs iterations for burn-in.
#' @param Gibbs_tol The number of Gibbs iterations after burn-in.
#'
#' @export
#' @examples
#' \dontrun{
#' sims <- get_sim_data(K = 5, N = 100, U = 200)
#' res <- Chloris(sims$RDR)
#' plot_inout(sims$RDR, list(res$cluster_est, sims$cluster_true), res$state_est) ## model result
#' plot_inout(sims$RDR, list(sims$cluster_true, res$cluster_est), sims$states_true) ## simulation truth
#' }


Chloris <- function(RDR = NULL, A = NULL, D = NULL, break_idx = NULL, init = "hclust",
                    K = 10, min_cluster_size = 1, cluster_shrink_tol = 20, S = 4,
                    burnin_tol = 300, Gibbs_tol = 300){

    #### cross checks =================================================================================
    if(is.null(A) != is.null(D)) stop("A and D need to be both provided or both NULL.")

    signal_RDR <- !is.null(RDR)
    signal_BAF <- !is.null(A)

    if(signal_RDR && signal_BAF){
        if(!(identical(dim(RDR), dim(A)) && identical(dim(D), dim(A)))) stop("Input RDR, A and D must have the same dimensions!")
    }else{
        if(!(signal_RDR || signal_BAF)) stop("Seems like RDR, A and D are all empty? At least give me something!")
    }

    #### prior ========================================================================================
    neutral_idx = 2
    prior_mu = log2(c(0.5, 1, 1.5, 2))
    prior_Q_diag = 20
    prior_beta_shape1 = 100
    prior_beta_shape2 = c(50, 5, 10)[1:(S - 1)]

    priors <- get_priors(K, S, signal_RDR, signal_BAF,
                         mu = prior_mu, beta_shape1 = prior_beta_shape1, beta_shape2 = prior_beta_shape2,
                         Q_diag = prior_Q_diag, neutral_idx = neutral_idx, to_neutral = 2)


    #### Gibbs samplling ===============================================================================
    res <- Gibbs_main(signal_RDR, signal_BAF, RDR, A, D,
                      priors, init, break_idx, burnin_tol, Gibbs_tol, cluster_shrink_tol, min_cluster_size)
    K <- res$K


    #### post process (label switching) =====================================================================
    if(!is.null(break_idx)){
        Pi_record <- list()
        for(iter in 1:Gibbs_tol){
            Pi_record[[iter]] <- list()
            for(k in 1:K){
                tmp_Pi_record <- lapply(res$Pi_record[[iter]][[k]], function(x) x[-1, ])
                Pi_record[[iter]][[k]] <- rbind(res$Pi_record[[iter]][[k]][[1]][1, ], abind::abind(tmp_Pi_record, along=1) )
            }
        }
    }else{
        Pi_record <- res$Pi_record
    } ## Pi_record
    res_post <- ECR(res$cluster_record, res$loglik_record)
    res_aligned <- align_MC_samples(res_post, res$cluster_record, res$state_record, res$Q_record, Pi_record)

    ## Final check for identity clones
    tmp <-  do.call(paste, as.data.frame(res$state_est))
    rep_cnt <- match(tmp, unique(tmp))

    if(any(table(rep_cnt) > 1)){
        cluster_est <- res_aligned$cluster_est
        rep_clones <- names(table(rep_cnt)[table(rep_cnt) > 1])
        for(k in rep_clones){
            clone_idx <- which(rep_cnt == k)
            cluster_est[cluster_est %in% clone_idx[-1]] <- k
        }
        res_aligned$cluster_est <- cluster_est
        res_aligned$state_est <- res_aligned$state_est[unique(rep_cnt), ]
    }


    return(list("cluster_est" = res_aligned$cluster_est, "state_est" = res_aligned$state_est,
                "par_record" = res$par_record, "Q_record" = res$Q_record, "K_hat" = K))

}
