### omitted dropout from original "get_simulated_cell_matrix_using_meanvar_trend"
sim_meanvar <- function(ref_gexp, gene_means, num_cells){
    
    mean_var_table = data.frame(m = rowMeans(ref_gexp), v = apply(ref_gexp, 1, var))
    
    #sim_cell_matrix <- .get_simulated_cell_matrix_using_meanvar_trend_helper(gene_means, mean_var_table, num_cells, dropout_logistic_params)
    # ----------------------------------------------------------
    ngenes = length(gene_means)
    logm = log(mean_var_table$m + 1)
    logv = log(mean_var_table$v + 1)
    mean_var_spline = smooth.spline(logv ~ logm)
    spike_cell_names = paste0('sim_cell_', seq_len(num_cells))
    
    sim_cell_matrix = matrix(rep(0, ngenes*num_cells), nrow = ngenes)
    rownames(sim_cell_matrix) = names(gene_means)
    colnames(sim_cell_matrix) = spike_cell_names
    
    sim_expr_vals <- function(gene_idx) {
        m = gene_means[gene_idx]
        val = 0
        if (m > 0) {
            logm = log(m+1)
            pred_log_var = predict(mean_var_spline, logm)$y
            var = max(exp(pred_log_var)-1, 0)
            val = round(max(rnorm(n=1, mean=m, sd=sqrt(var)), 0))
        }
        return(val)
    }
    
    for (i in seq_len(num_cells)) {
        newvals = sapply(seq_len(ngenes), FUN=sim_expr_vals)
        sim_cell_matrix[,i] = newvals
    }
    # ----------------------------------------------------------
    
    return(sim_cell_matrix)
}

smooth_helper <- function(obs_data, window_length) {
    # strip NAs out and replace after smoothing
    orig_obs_data = obs_data
    
    nas = is.na(obs_data)
    
    obs_data = obs_data[!nas]
    
    obs_length <- length(obs_data)
    end_data <- obs_data
    
    tail_length = (window_length - 1)/2
    if (obs_length >= window_length) {
        end_data <- .smooth_center_helper(obs_data, window_length)
    }
    
    # end_data will have the end positions replaced with mean values, smoothing just at the ends.
    
    obs_count <- length(obs_data)
    
    numerator_counts_vector = c(seq_len(tail_length), tail_length + 1, c(tail_length:1))
    
    # defining the iteration range in cases where the window size is larger than the number of genes. In that case we only iterate to the half since the process is applied from both ends.
    iteration_range = ifelse(obs_count > window_length, tail_length, ceiling(obs_count/2))
    
    for (tail_end in seq_len(iteration_range)) {
        end_tail = obs_count - tail_end + 1
        
        d_left = tail_end - 1
        d_right = obs_count - tail_end
        d_right = ifelse(d_right > tail_length, tail_length, d_right)
        
        r_left = tail_length - d_left
        r_right = tail_length - d_right
        
        denominator = (((window_length - 1)/2)^2 + window_length) - ((r_left * (r_left + 1))/2) - ((r_right * (r_right + 1))/2)
        
        left_input_vector_chunk = obs_data[seq_len(tail_end + d_right)]
        right_input_vector_chunk = obs_data[(end_tail - d_right):obs_length]
        
        numerator_range = numerator_counts_vector[(tail_length + 1 - d_left):(tail_length + 1 + d_right)]
        
        end_data[tail_end] = sum(left_input_vector_chunk * numerator_range)/denominator
        end_data[end_tail] = sum(right_input_vector_chunk * rev(numerator_range))/denominator
    }
    
    orig_obs_data[! nas] = end_data  # replace original data with end-smoothed data
    
    return(orig_obs_data)
}


.smooth_center_helper <- function(obs_data, window_length){
    
    nas = is.na(obs_data)
    vals = obs_data[! nas]
    
    custom_filter_denominator = ((window_length-1)/2)^2 + window_length
    custom_filter_numerator = c(seq_len((window_length-1)/2), ((window_length-1)/2)+1, c(((window_length-1)/2):1))
    
    custom_filter = custom_filter_numerator/rep(custom_filter_denominator, window_length)
    
    smoothed = stats::filter(vals, custom_filter, sides=2)
    
    ind = which(! is.na(smoothed))
    vals[ind] = smoothed[ind]
    
    obs_data[! nas] = vals
    
    return(obs_data)
}