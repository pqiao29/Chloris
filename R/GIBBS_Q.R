update_transition <- function(est, priors){

    K <- est$K
    S <- est$S
    Q_gamma_0 <- priors$Q
    Q_est <- est$Q
    if(is.list(est$H)){
        H_est <- est$H_array
    }else{
        H_est <- est$H
    }

    for(k in 1:K){

        H_est_k <- H_est[k, ]

        for(r in 1:S){
            tmp_gamma <- rep(NA, S)
            for(s in 1:S){
                tmp_gamma[s] <- Q_gamma_0[r, s] + my_v_match(c(r, s), H_est_k)
            }

            Q_est[[k]][r, ] <- DirichletReg::rdirichlet(1, tmp_gamma)
        }
    }

    est$Q <- Q_est
    est
}
