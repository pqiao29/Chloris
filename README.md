# Chloris
A Bayesian framework for analysing CNA intra-tumoral heterogeneity from scRNA data. 

The model specification is detailed in the following publication:  
[Qiao, P., Kwok, C. F., Qian, G., & McCarthy, D. J. (2023). Bayesian inference for copy number intra-tumoral heterogeneity from single-cell RNA-sequencing data. bioRxiv, 2023-10.](https://www.biorxiv.org/content/10.1101/2023.10.22.563455v1)

## Installation
```
require(devtools)
devtools::install_github('pqiao29/Chloris')
```
## Quick Start
```
sims <- get_sim_data(K = 5, N = 100, U = 200)
res <- Chloris(sims$RDR)
plot_inout(sims$RDR, list(res$cluster_est, sims$cluster_true), res$state_est) ## model result
plot_inout(sims$RDR, list(sims$cluster_true, res$cluster_est), sims$states_true) ## simulation truth
```
