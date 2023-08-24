#' Simulate transition matrix
#' @keywords internal
get_true_transition <- function(K, S, diag, to_neutral, neutral_idx, off_diag_variance = 0.1){

    if(length(diag) != K && length(diag) != 1) stop("The length of diag has to be either 1 or equal to K")
    if(length(diag) == 1) diag <- rep(diag, K)

    Q <- array(NA, dim = c(S, S, K))

    for(k in 1:K){
        Q_k <- matrix(rnorm(S^2, 1, off_diag_variance), S, S)
        if(!is.null(to_neutral)) Q_k[-neutral_idx, neutral_idx] <- rnorm(1, to_neutral, off_diag_variance)
        diag(Q_k) <- diag[k]
        Q[, , k] <- Q_k/rowSums(Q_k)
    }

    return(Q)
}


#' Simulate state profiles for a clone.
#' @keywords internal
generate_states <- function(Q, U, neutral_idx){

    if(nrow(Q) != ncol(Q)) stop("Q must be a square matrix!")
    S <- nrow(Q)

    Pi_init <- rep(1, S)
    Pi_init[neutral_idx] <- 10
    Pi_init <- Pi_init/sum(Pi_init)

    states <- vector(length = U)
    states[1] <- sample(1:S, 1, prob = Pi_init)
    for(t in 2:U){
        states[t] <- sample(1:S, 1, prob = Q[states[t - 1], ])
    }

    states
}

#' Simulate cell clustering and state profile for each cluster
#' @keywords internal
generate_cluster_states <- function(K, U, Q, neutral_idx, overlap = TRUE){

    states <- matrix(neutral_idx, K, U)
    for(k in 2:K){
        while(sum(states[k, ] != neutral_idx) == 0){
            states[k, ] <- generate_states(Q[, , k], U, neutral_idx)
        }
    }

    if(overlap && K >= 3){
        copy_from_cluster <- c(3:K, 2)
        for(k in 2:K){
            next_cluster <- copy_from_cluster[k - 1]
            tmp <-  table(states[next_cluster, ])
            bd <- sample(as.numeric(names(tmp[tmp > 2])), 1) ## find a state
            bd_idx <- sort(sample(which(states[next_cluster, ] == bd), 2))
            states[k, bd_idx[1]:bd_idx[2]] <- states[next_cluster, bd_idx[1]:bd_idx[2]]
        }
    }

    return(states)
}

#' Simulate cell clustering
#' @keywords internal
sim_clusters <- function(N, K, clone_prob = NULL){
    if(is.null(clone_prob)) clone_prob <- rep(1/K, K)
    t(rmultinom(N, 1, clone_prob))
}

#' Simulate BAF
#' @keywords internal
sim_BAF <- function(cell_level_states, theta_true, Md, missing_percent){

    N <- ncol(cell_level_states)
    U <- nrow(cell_level_states)

    A <- matrix(NA, U, N)
    D <- matrix(sample(2:Md, U*N, replace = T), U, N)
    if(missing_percent > 0) D[sample(1:length(D), length(D)*missing_percent, replace = F)] <- 0

    for(i in 1:N){
        for(u in 1:U){
            A[u, i] <- rbinom(n = 1, size = D[u, i], prob = theta_true[u, i])
        }
    }

    return(list("A" = A, "D" = D))
}

#' Simulate RDR
#' @keywords internal
sim_RDR <- function(cell_level_states, CN_levels, from_Splat, RDR_sd = NULL, RDR_outlier_cnt = 0){

    S <- length(CN_levels)
    N <- ncol(cell_level_states)
    U <- nrow(cell_level_states)

    if(!from_Splat){
        if(is.null(RDR_sd)) stop("If from_Splat == FALSE, RDR_sd must be provided!")
        RDR <- matrix(NA, nrow(cell_level_states), ncol(cell_level_states))
        for(s in 1:S){
            RDR[cell_level_states == s] <- rnorm(sum(cell_level_states == s), CN_levels[s], RDR_sd)
        }
        ret <- list("RDR" = RDR)
    }else{
        raw_count <- sim_expr_Splat(cell_level_states, 2^CN_levels, outlier_cnt = RDR_outlier_cnt)
        processed_data <- preprocess_sims(raw_count)
        ret <- list("RDR" = processed_data$RDR, "RDR_summary" = processed_data$summary,
                    "raw" = list("tumor" = raw_count$tumor, "ref" = raw_count$ref),
                    "keep_gene" = raw_count$keep_gene, "cell_level_states" = raw_count$cell_level_RDR_state)
    }

    return(ret)
}


#' Simulate gene means from Gamma distribution for gene count
#' @keywords internal
Splat_GeneMeans <- function(U){
    gamma_shape <- 0.5
    gamma_rate <- 0.15
    target_mean <- gamma_shape/gamma_rate
    target_var <- gamma_shape/(gamma_rate)^2
    logvar_est <- log((target_var/(target_mean^2)) + 1)
    logmean_est <- log((target_mean^2)/sqrt(target_var + target_mean^2))
    GeneMean <- rlnorm(U, logmean_est, sqrt(logvar_est))
    #GeneMean <- rgamma(U, shape = gamma_shape, rate = gamma_rate)
    return(GeneMean)
}

#' Simulate mean-variance trend for gene count
#' @keywords internal
Splat_meanvar <- function(BaseMeans, outlier_cnt = 0){
    N <- ncol(BaseMeans)
    U <- nrow(BaseMeans)
    bcv.common = 0.1
    bcv.df = 16
    BCV <- (bcv.common + (1 / rowMeans(BaseMeans))) * sqrt(bcv.df / rchisq(U, df = bcv.df))
    BCV <- matrix(BCV, U, N)
    #if(outlier_cnt >= 1) BCV[, N - (outlier_cnt - 1):0] <- 0.5
    Counts <- matrix(rnbinom(U*N, size = 1/(BCV ^ 2), mu = BaseMeans), nrow = U, ncol = N)
    return(Counts)
}

#' Simulate gene count
#' @keywords internal
sim_expr_Splat <- function(cell_level_states, CN_levels, outlier_cnt = 0, verbose = FALSE){

    N <- ncol(cell_level_states)
    intended_U <- nrow(cell_level_states)

    ## Base mean: tumor and ref share the same prob_gene
    ExpLibSize <- rlnorm(N, 11, 0.7)
    GeneMean <- Splat_GeneMeans(intended_U)
    prob_gene <- GeneMean/sum(GeneMean)
    BaseMeans_tumor <- prob_gene %*% t(ExpLibSize)
    BaseMeans_ref <- prob_gene * median(ExpLibSize)

    ## Add CNV to tumor
    copy_number <- as.numeric(as.character(factor(cell_level_states, levels = 1:length(CN_levels), labels = CN_levels)))
    copy_number <- matrix(copy_number, intended_U, N)
    BaseMeans_tumor <- BaseMeans_tumor * copy_number

    ## NB with common dispersion
    counts_tumor <- Splat_meanvar(BaseMeans_tumor, outlier_cnt)
    counts_ref <- Splat_meanvar(matrix(BaseMeans_ref, ncol = 1))

    ## zero's in bulk ref considered uninterpretable
    keep_gene <- (counts_ref != 0)
    U <- sum(keep_gene)
    if(verbose) cat("Simulated ", U, "genes with non-zero bulk reference\n")
    counts_ref <- counts_ref[keep_gene]
    counts_tumor <- counts_tumor[keep_gene, ]

    return(list("tumor" = counts_tumor, "ref" = counts_ref,
                "keep_gene" = keep_gene, "cell_level_RDR_states" = copy_number[keep_gene, ]))

}

#' Preprocess simulated gene expression data to RDR and summary(RDR|state)
#' @keywords internal
preprocess_sims <- function(count_data){

    for(x in names(count_data)) assign(x, count_data[[x]])

    ## Equalize library size
    tumor <- t(t(tumor)/colSums(tumor))
    ref <- ref/sum(ref)

    ## Truncate zero's before taking log
    tumor[tumor == 0] <- min(tumor[tumor > 0])

    ## log ratio
    RDR <- log2(tumor) - log2(ref)

    ## summarize
    CNs <- sort(unique(c(cell_level_RDR_states)))
    print_summary <- NULL
    for(CN in CNs){
        tmp_RDR <- RDR[cell_level_RDR_states == CN]
        print_summary <- rbind(print_summary,
                               c(log2(CN),
                                 mean(tmp_RDR), median(tmp_RDR), sd(tmp_RDR),
                                 sum(cell_level_RDR_states == CN)*100/prod(dim(RDR))))
    }
    colnames(print_summary) <- c("Expected", "mean", "median", "SD", "proportion (%)")
    rownames(print_summary) <- CNs

    list("RDR" = RDR, "summary" = print_summary)
}
