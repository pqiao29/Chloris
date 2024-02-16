#' Main plot function
#'
#' Plot input RDR and output clustering and clonal copy number profiles.
#'
#' @param input Matrix of RDR or BAF, where each row is a gene and each column represents a cell.
#' @param type "RDR" or "BAF".
#' @param cluster_labels A list of one or more cell partitions, each partition is represented as
#' an N dimensional vector of cluster labels, N is the number of cells.
#' The columns in \code{input} will be ordered according to the first element in \code{cluster_labels}.
#' If \code{cluster_labels == NULL}, the order of columns in \code{input} will be maintained.
#' @param CN_states A list of K by U matrices of copy number states in K clusters with U genes.
#' @param state_mean
#' Only needed if \code{!is.null(CN_states)} for adjusting state colour representation in clonal profiles.
#' The \code{CN_states} will be plotted together with \code{input} as one heatmap, where each state is
#' represented by the according value in \code{state_mean}. 
#' If \code{state_mean == "standard"}, state labels take the most saturated colours.  
#' If \code{state_mean == "input"}, state colours reflect observed RDR assigned to each state.
#' @param lim The imposed range of \code{input}.
#' @param break_idx The index of the first gene in each chromosome. This will generate a horizontal line for each index.
#' @param cluster_colour A vector of colours for each cluster.
#' @return gg object
#'
#' @import ggplot2
#' @export
#'
#' @examples
#' \dontrun{
#' sims <- get_sim_data(K = 4, N = 100, U = 200)
#' plot_inout(sims$RDR)
#' plot_inout(sims$RDR, list(sims$cluster_true), sims$states_true, state_mean = c(-1, 0, 0.5, 1.5))
#' }
## cluster_labels: list of cell cluster labels. The cells or columns in input will be ordered according to the first element
plot_inout <- function(input, type = "RDR", cluster_labels = NULL,
                       CN_states = NULL, state_mean = "standard",
                       lim = NULL, break_idx = NULL,
                       cluster_colour = NULL) {
  U <- nrow(input)
  N <- ncol(input)

  if (type == "RDR") {
    heatmap_colour <- list("low" = "steelblue", "mid" = "white", "high" = "red1", "na" = "gray")
    midpoint <- 0
    background <- 0
  } else {
    heatmap_colour <- list("mid" = "turquoise", "low" = "red", "high" = "red", "na" = "white")
    midpoint <- 0.5
    background <- NA
  }


  ### axis span
  if (is.null(cluster_labels)) {
    if (!is.null(CN_states)) stop("cluster_labels needs to be provided for plotting CN_states.")
    df <- expand.grid(X = 1:N, Y = 1:U)
    df$input <- c(t(input))
  } else {
    gap_length_y <- max(round(U / 100), 1)
    cluster_length_y <- round(U / 30)

    X_span <- N
    Y_span <- U + (gap_length_y + cluster_length_y) * length(cluster_labels)

    I <- cluster_labels[[1]]
    
    ## Cluster count of all groups
    if(length(cluster_labels) >= length(CN_states)){
        Ks <- lapply(cluster_labels, function(x) length(unique(x))) 
    }else{
        Ks <- lapply(CN_states, function(x) nrow(x)) 
    }
    Ks <- unlist(Ks)
    
    if (N != length(I)) stop("The number of columns of RDR needs to be the same as the length of I_est!")

    if (is.null(cluster_colour)) {
      cluster_colour <- list(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"))
      cluster_colour <- rep(cluster_colour, length(cluster_labels))
    }

    if (!is.null(CN_states)) {
      fixed_gap_x <- max(round(N / 25), 1)
      gap_length_x <- max(round(N / 80), 1)
      CNV_length_x <- max(round(N / 10), 1)
      for(group in 1:length(CN_states)) X_span <- X_span + (gap_length_x + CNV_length_x) * Ks[group] + fixed_gap_x

      ### Obtain label colors for states in clonal profile 
      S <- max(unlist(CN_states))
      state_means <- rep(0, S)
      if (state_mean == "input") {
          tmp_I <- matrix(0, N, K1)
          for (k in 1:K) tmp_I[I == k, k] <- 1
          input_state <- get_cell_level_states(tmp_I, CN_states[[1]])
          for (s in 1:S) state_means[s] <- mean(input[input_state == s], na.rm = T)
      }
      if(state_mean == "standard"){ ## Only takes S = 3 or 4 for now
          if(S == 4) probs <- c(0.1, 0.8, 0.98) else probs <- c(0.1, 0.8) 
          state_means[-2] <- quantile(c(input), probs = probs)
      }
      
      CN_states <- lapply(CN_states, function(x){
          for(s in 1:S) x[x == s] <- state_means[s]
          return(x)
      })
    }

    df <- expand.grid(X = 1:X_span, Y = 1:Y_span)

    ### Order cells according to clustering
    cell_order <- sort(I, index.return = TRUE)$ix
    input_sorted <- input[, cell_order]
    cluster_labels_sorted <- lapply(cluster_labels, function(x) x[cell_order])

    ### Initialize data.frame for ggplot
    inout_matrix <- rbind(input_sorted, matrix(background, Y_span - U, N))

    ### Add states
    if (!is.null(CN_states)) {
        for(group in 1:length(CN_states)){
            inout_matrix <- cbind(inout_matrix, matrix(background, nrow(inout_matrix), fixed_gap_x))
            for (k in 1:Ks[group]) {
                states_k <- rbind(
                    matrix(CN_states[[group]][k, ], U, CNV_length_x),
                    matrix(background, Y_span - U, CNV_length_x)
                )
                inout_matrix <- cbind(inout_matrix, matrix(background, Y_span, gap_length_x), states_k)
            }
        }
    }
    df$input <- c(t(inout_matrix))
  }

  ### truncate max and min
  if (!is.null(lim)) {
    df$input[df$input > lim[2]] <- lim[2]
    df$input[df$input < lim[1]] <- lim[1]
  }

  ### Plot heatmap
  gg_ret <- ggplot(df, aes(X, Y)) +
    geom_tile(aes(fill = input)) +
    scale_fill_gradient2(low = heatmap_colour$low, high = heatmap_colour$high, mid = heatmap_colour$mid, midpoint = midpoint, na.value = heatmap_colour$na, limits = lim) +
    theme(
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none",
      panel.background = element_blank(),
      axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()
    )

  ### Add clustering labels and CN states
  if (!is.null(cluster_labels)) {
    ### First set of cluster labels
    y_low <- U + gap_length_y
    y_high <- y_low + cluster_length_y
    I_sorted <- cluster_labels_sorted[[1]]
    for (k in 1:Ks[1]) {
      cluster_idx <- which(I_sorted == k)
      gg_ret <- gg_ret + annotate("rect",
        xmin = cluster_idx[1] - 0.5, xmax = cluster_idx[length(cluster_idx)] + 0.5,
        ymin = y_low, ymax = y_high,
        fill = cluster_colour[[1]][k],
        alpha = .8
      )
    }

    ### Other cluster labels for comparison
    if (length(cluster_labels) > 1) {
      for (t in 2:length(cluster_labels)) {
        y_low <- U + gap_length_y + (t - 1) * (gap_length_y + cluster_length_y)
        y_high <- y_low + cluster_length_y

        for (i in 1:N) {
          k <- cluster_labels_sorted[[t]][i]
          gg_ret <- gg_ret + annotate("rect",
            xmin = i - 0.5, xmax = i + 0.5,
            ymin = y_low, ymax = y_high,
            fill = cluster_colour[[t]][k],
            alpha = .8
          )
        }
      }
    }

    ### Add CN states
    if (!is.null(CN_states)) {
      y_low <- U + gap_length_y
      y_high <- U + gap_length_y + cluster_length_y
      x_state_left <- N + fixed_gap_x + gap_length_x
      x_state_step <- CNV_length_x + gap_length_x
      
      xmin <- x_state_left + 1.5
      xmax <- xmin + CNV_length_x + 1.5
      
      for(group in 1:length(CN_states)){
          for (k in 1:Ks[group]) {
              cluster_idx <- which(I_sorted == k)
              gg_ret <- gg_ret + annotate("rect",
                                          xmin = xmin,
                                          xmax = xmax,
                                          ymin = y_low, ymax = y_high,
                                          fill = cluster_colour[[1]][k], alpha = .8)
              xmin <- xmin + x_state_step
              xmax <- xmin + CNV_length_x + 1
          }
          xmin <- xmin + fixed_gap_x
          xmax <- xmax + fixed_gap_x
      }
      
    }
  }

  if (!is.null(break_idx)) gg_ret <- gg_ret + geom_hline(yintercept = break_idx)

  return(gg_ret)
}
