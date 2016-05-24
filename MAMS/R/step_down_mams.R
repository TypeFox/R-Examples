step_down_mams <- function(nMat = matrix(c(10, 20), nrow=2, ncol=4), alpha_star = c(0.01, 0.025), lb = 0, selection = "all_promising"){

    # checking input parameters
    if (!all(diff(nMat) >= 0)) {stop("total sample size per arm cannot decrease between stages.")}
    J <- dim(nMat)[1]
    K <- dim(nMat)[2] - 1
    if ((J != 2) && (J != 3)) {stop("number of stages must be 2 or 3")}
    if (K < 2) {stop("must have at least two experimental treatments")}
    if (length(alpha_star) != J) {stop("length of error spending vector must be same as number of stages")}
    if (!all(diff(alpha_star) >= 0)) {stop("cumulative familywise error must increase.")}
    if (length(lb) != J - 1) {stop("lower boundary must be specified at all analysis points except the last")}
    if ((selection != "all_promising") && (selection != "select_best")) {stop("invalid selection method")}
    
  get.hyp <- function(n){ # find the nth intersection hypothesis (positions of 1s in binary n)
    indlength = ceiling(log(n)/log(2)+.0000001)
    ind = rep(0,indlength)
    newn=n
        
    for (h in seq(1,indlength)){
      ind[h] = (newn/(2^(h-1))) %% 2
      newn = newn - ind[h]*2^(h-1)
    }
    seq(1,indlength)[ind==1]
  }

  create_block <- function(control_ratios = 1:2, active_ratios = matrix(1:2, 2, 3)){ # for argument c(i,j) this gives covariance between statistics in stage i with statistics in stage j
    K <- dim(active_ratios)[2]
    block <- matrix(NA, K, K)
    for(i in 1:K){
      block[i, i] <- sqrt(active_ratios[1, i] * control_ratios[1] * (active_ratios[2, i] + control_ratios[2]) / (active_ratios[1, i] + control_ratios[1]) / active_ratios[2, i] / control_ratios[2])
    }
    for (i in 2:K){
      for (j in 1:(i - 1)){
        block[i, j] <- sqrt(active_ratios[1, i] * control_ratios[1] * active_ratios[2, j] / (active_ratios[1, i] + control_ratios[1]) / (active_ratios[2, j] + control_ratios[2]) / control_ratios[2])
        block[j, i] <- sqrt(active_ratios[1, j] * control_ratios[1] * active_ratios[2, i] / (active_ratios[1, j] + control_ratios[1]) / (active_ratios[2, i] + control_ratios[2]) / control_ratios[2])
      }
    }
    block
  }

  create_cov_matrix <- function(control_ratios = 1:2, active_ratios = matrix(1:2, 2, 3)){ # create variance-covariance matrix of the test statistics

    J <- dim(active_ratios)[1]
    K <- dim(active_ratios)[2]

    cov_matrix <- matrix(NA, J * K, J * K)
    for (i in 1:J){
      for (j in i:J){
        cov_matrix[((i - 1) * K + 1):(i * K), ((j - 1) * K + 1):(j * K)] <- create_block(control_ratios[c(i, j)], active_ratios[c(i, j), ])
        cov_matrix[((j - 1) * K + 1):(j * K), ((i - 1) * K + 1):(i * K)] <- t(cov_matrix[((i - 1) * K + 1):(i * K), ((j - 1) * K + 1):(j * K)])
      }
    }
    cov_matrix
  }

  get_path_prob <- function(surviving_subset1, surviving_subset2 = NULL, cut_off, treatments, cov_matrix, lb, upper_boundary, K, stage){ # find the probability that no test statistic crosses the upper boundary + only treatments in surviving_subsetj reach the jth stage
    treatments2 <- treatments[surviving_subset1]
    if (stage == 2){
      lower <- c(rep(-Inf, length(treatments)), rep(-Inf, length(treatments2)))
      lower[surviving_subset1] <- lb[1]

      upper <- c(rep(lb[1], length(treatments)), rep(cut_off, length(treatments2)))
      upper[surviving_subset1] <- upper_boundary[1]

      return(pmvnorm(lower = lower, upper = upper, sigma = cov_matrix[c(treatments, K + treatments2), c(treatments, K + treatments2)])[1])
    }
    treatments3 <- treatments2[surviving_subset2]

    lower <- c(rep(-Inf, length(treatments)), rep(-Inf, length(treatments2)), rep(-Inf, length(treatments3)))
    lower[surviving_subset1] <- lb[1]
    lower[length(treatments) + surviving_subset2] <- lb[2]

    upper <- c(rep(lb[1], length(treatments)), rep(lb[2], length(treatments2)), rep(cut_off, length(treatments3)))
    upper[surviving_subset1] <- upper_boundary[1]
    upper[length(treatments) + surviving_subset2] <- upper_boundary[2]

    pmvnorm(lower = lower, upper = upper, sigma = cov_matrix[c(treatments, K + treatments2, 2 * K + treatments3), c(treatments, K + treatments2, 2 * K + treatments3)])[1]
  }


  rejection_paths <- function(selected_treatment, cut_off, treatments, cov_matrix, lb, upper_boundary, K, stage){ # for the "select_best" method, find the probability that "select_treatment" is selected and subsequently crosses the upper boundary

    contrast <- diag(-1, K + stage - 1)
    contrast[1:K, selected_treatment] <- 1
    for (i in 1:(stage - 1)) contrast[K + i, K + i] <- 1
  
    bar_cov_matrix <- contrast %*% cov_matrix[c(1:K, 1:(stage - 1) * K + selected_treatment), c(1:K, 1:(stage - 1) * K + selected_treatment)] %*% t(contrast)

    lower <- c(rep(0, length(treatments)), cut_off)
    if (stage > 2) lower <- c(rep(0, length(treatments)), lb[2:(stage - 1)], cut_off)
    lower[which(treatments == selected_treatment)] <- lb[1]

    upper <- c(rep(Inf, length(treatments)), Inf)
    if (stage > 2) upper <- c(rep(Inf, length(treatments)), upper_boundary[2:(stage - 1)], Inf)
    upper[which(treatments == selected_treatment)] <- upper_boundary[1]

    pmvnorm(lower = lower, upper = upper, sigma = bar_cov_matrix[c(treatments, K + 1:(stage - 1)), c(treatments, K + 1:(stage - 1))])[1]
    
  }

  excess_alpha <- function(cut_off, alpha_star, treatments, cov_matrix, lb, upper_boundary, selection, K, stage){ # for "all_promising" rule, this gives the cumulative typeI error for 'stage' stages

      # for "select_best" rule, this gives the Type I error spent at the 'stage'th stage

      if (stage == 1) return(1 - alpha_star[1] - pmvnorm(lower = rep(-Inf, length(treatments)), upper = rep(cut_off, length(treatments)), sigma = cov_matrix[treatments, treatments])[1])
      if (selection == "select_best") return(alpha_star[stage] - alpha_star[stage - 1] - sum(unlist(lapply(treatments, rejection_paths, cut_off = cut_off, treatments = treatments, cov_matrix = cov_matrix, lb = lb, upper_boundary = upper_boundary, K = K, stage = stage)))) # any of 'treatments' could be selected, so we add all these probabilities
      if (stage == 2){
          surviving_subsets <- c(list(numeric(0)), lapply(as.list(1:(2 ^ length(treatments) - 1)), get.hyp)) # list all possible subsets of surviving treatments after the first stage
          return(1 - alpha_star[2] - sum(unlist(lapply(surviving_subsets, get_path_prob, cut_off = cut_off, treatments = treatments, cov_matrix = cov_matrix, lb = lb, upper_boundary = upper_boundary, K = K, stage = stage))))
      }
      surviving_subsets1 <- c(list(numeric(0)), lapply(as.list(1:(2 ^ length(treatments) - 1)), get.hyp)) # all possible subsets of surviving treatments after the first stage
      surviving_subsets2 <- c(list(list(numeric(0))), lapply(surviving_subsets1[-1], function(x) c(list(numeric(0)), lapply(as.list(1:(2 ^ length(x) - 1)), get.hyp)))) # for each possible subset of survivng subsets after stage 1, list the possible subsets still surviving after stage 2
      1 - alpha_star[3] - sum(unlist(Map(function(x, y) sum(unlist(lapply(y, get_path_prob, surviving_subset1 = x, cut_off = cut_off, treatments = treatments, cov_matrix = cov_matrix, lb = lb, upper_boundary = upper_boundary, K = K, stage = stage))), surviving_subsets1, surviving_subsets2)))
  }


  # get sample size ratios
  R <- nMat[, -1] / nMat[1, 1]
  r0 <- nMat[, 1] / nMat[1, 1]

  cov_matrix <- create_cov_matrix(r0, R)
  
  l <- u <- as.list(1:(2 ^ K - 1))

  alpha_star <- rep(list(alpha_star), 2 ^ K - 1)
        
  for (i in 1:(2 ^ K - 1)){

    names(u)[i] <- paste("U_{",paste(get.hyp(i), collapse = " "),"}",sep="")
    names(l)[i] <- paste("L_{",paste(get.hyp(i), collapse = " "),"}",sep="")
    names(alpha_star)[i] <- paste("alpha_star_{",paste(get.hyp(i), collapse = " "),"}",sep="")

    for (j in 1:J){

      try(new_u <- uniroot(excess_alpha, c(0, 10), alpha_star = alpha_star[[i]], treatments = get.hyp(i), cov_matrix = cov_matrix, lb = lb, upper_boundary = u[[i]], selection = selection, K = K, stage = j)$root, silent = TRUE)
      if (is.null(new_u)) {stop("upper boundary not between 0 and 10")}
      u[[i]][j] <- round(new_u, 2)
                
    }

    l[[i]] <- c(lb, u[[i]][J])
  }

  res <- NULL
  res$l <- l
  res$u <- u
  res$sample_sizes <- nMat
  res$K <- K
  res$J <- J
  res$alpha_star <- alpha_star
  res$selection <- selection
  res$z_scores <- NULL
  res$selected_trts <- list(1:K)
  class(res) <- "MAMS.step_down"

  return(res)
}
