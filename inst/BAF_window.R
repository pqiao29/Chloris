BAF_window <- function(A, D, window_size){
  
  U <- nrow(D)

  for(cell in colnames(A)){
    
    snp_cell <- which(D[, cell] > 0)
    snp_distance <- snp_cell[-1] - snp_cell[-length(snp_cell)] - 1
    
    ### Left window
    snp_cell_idx <- snp_cell[1]
    D[snp_cell_idx - min(window_size, snp_cell_idx - 1):1, cell] <- D[snp_cell_idx, cell]
    A[snp_cell_idx - min(window_size, snp_cell_idx - 1):1, cell] <- A[snp_cell_idx, cell]
    
    if(length(snp_cell) > 1){
      for(i in 2:length(snp_cell)){
        snp_cell_idx <- snp_cell[i]
        D[snp_cell_idx - min(window_size, round((snp_distance[i - 1] - 1)/2)):1, cell] <- D[snp_cell_idx, cell]
        A[snp_cell_idx - min(window_size, round((snp_distance[i - 1] - 1)/2)):1, cell] <- A[snp_cell_idx, cell]
      }
    }
    
    ### right window
    if(length(snp_cell) > 1){
      for(i in 1:(length(snp_cell) - 1)){
        snp_cell_idx <- snp_cell[i]
        D[snp_cell_idx + min(window_size, round((snp_distance[i] - 1)/2)):1, cell] <- D[snp_cell_idx, cell]
        A[snp_cell_idx + min(window_size, round((snp_distance[i] - 1)/2)):1, cell] <- A[snp_cell_idx, cell]
      }
    }
    
    snp_cell_idx <- snp_cell[length(snp_cell)]
    D[snp_cell_idx + min(window_size, U - snp_cell_idx):1, cell] <- D[snp_cell_idx, cell]
    A[snp_cell_idx + min(window_size, U - snp_cell_idx):1, cell] <- A[snp_cell_idx, cell]
    
  }
  
  return(list(A = A, D = D))
  
}