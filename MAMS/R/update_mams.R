

update_mams <- function(current_mams = step_down_mams(), n_obs = NULL, z_scores = NULL, selected_trts = NULL, n_future = NULL) {

    
    z_scores <- c(current_mams$z_scores, list(z_scores))

    if (!is.null(selected_trts)) selected_trts <- c(current_mams$selected_trts, list(selected_trts))
    
    # checking input parameters
    if (class(current_mams) != "MAMS.step_down") {stop("current_mams must be a 'MAMS.step_down' object")}
    if (length(n_obs) != current_mams$K + 1) {stop("must provide observed cumulative sample size for each treatment")}

    completed_stages <- length(z_scores)
    for (i in 1:completed_stages) {
        if (length(z_scores[[i]]) != current_mams$K) {stop("vector of statistics is wrong length")}
    }

    if (is.null(selected_trts)){
        if (current_mams$J > completed_stages) {stop("must specify treatments selected for next stage")}
    }
    
    for (i in seq_along(selected_trts)){
        if (length(setdiff(selected_trts[[i]], 1:current_mams$K) > 0)) {stop("inappropriate treatment selection")}
    }
    if (is.matrix(n_future)){
        if (dim(n_future)[1] != current_mams$J - completed_stages) {stop("must provide future sample sizes for all remaining stages")}
        if (dim(n_future)[2] != current_mams$K + 1) {stop("must provide future sample sizes for all treatment arms")}
    }

    # load all necessary functions
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
    create_block <- function(control_ratios = 1:2, active_ratios = matrix(1:2, 2, 3)){  # for argument c(i,j) this gives covariance between statistics in stage i with statistics in stage j

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
  create_cond_cov_matrix <- function(cov_matrix, K, completed_stages){ # find the conditional covariance of future test statistics given data so far

    sigma_1_1 <- cov_matrix[((completed_stages - 1) * K + 1):(completed_stages * K), ((completed_stages - 1) * K + 1):(completed_stages * K)]
    sigma_1_2 <- cov_matrix[((completed_stages - 1) * K + 1):(completed_stages * K), -(1:(completed_stages * K))]
    sigma_2_1 <- t(sigma_1_2)
    sigma_2_2 <- cov_matrix[-(1:(completed_stages * K)), -(1:(completed_stages * K))]
    sigma_2_2 - sigma_2_1 %*% solve(sigma_1_1) %*% sigma_1_2

  }
  create_cond_mean <- function(cov_matrix, K, completed_stages, z_scores){ # find the conditional mean of future test statistics given data so far

    sigma_1_1 <- cov_matrix[((completed_stages - 1) * K + 1):(completed_stages * K), ((completed_stages - 1) * K + 1):(completed_stages * K)]
    sigma_1_2 <- cov_matrix[((completed_stages - 1) * K + 1):(completed_stages * K), -(1:(completed_stages * K))]
    sigma_2_1 <- t(sigma_1_2)
    sigma_2_1 %*% solve(sigma_1_1) %*% z_scores

  } 
  get_path_prob <- function(surviving_subset1, surviving_subset2 = NULL, cut_off, treatments, cov_matrix, lower_boundary, upper_boundary, K, stage, z_means){ # find the probability that no test statistic crosses the upper boundary + only treatments in surviving_subsetj reach the jth stage

    treatments2 <- treatments[surviving_subset1]
    if (stage == 2){
        lower <- c(rep(-Inf, length(treatments)), rep(-Inf, length(treatments2)))
        lower[surviving_subset1] <- lower_boundary[1]

        upper <- c(rep(lower_boundary[1], length(treatments)), rep(cut_off, length(treatments2)))
        upper[surviving_subset1] <- upper_boundary[1]

        return(pmvnorm(lower = lower, upper = upper, mean = z_means[c(treatments, K + treatments2)], sigma = cov_matrix[c(treatments, K + treatments2), c(treatments, K + treatments2)])[1])
    }
    treatments3 <- treatments2[surviving_subset2]

    lower <- c(rep(-Inf, length(treatments)), rep(-Inf, length(treatments2)), rep(-Inf, length(treatments3)))
    lower[surviving_subset1] <- lower_boundary[1]
    lower[length(treatments) + surviving_subset2] <- lower_boundary[2]

    upper <- c(rep(lower_boundary[1], length(treatments)), rep(lower_boundary[2], length(treatments2)), rep(cut_off, length(treatments3)))
    upper[surviving_subset1] <- upper_boundary[1]
    upper[length(treatments) + surviving_subset2] <- upper_boundary[2]

    pmvnorm(lower = lower, upper = upper, mean = z_means[c(treatments, K + treatments2, 2 * K + treatments3)], sigma = cov_matrix[c(treatments, K + treatments2, 2 * K + treatments3), c(treatments, K + treatments2, 2 * K + treatments3)])[1]
  }
  rejection_paths <- function(selected_treatment, cut_off, treatments, cov_matrix, lower_boundary, upper_boundary, K, stage, z_means){  # for the "select_best" method, find the probability that "select_treatment" is selected and subsequently crosses the upper boundary

    contrast <- diag(-1, K + stage - 1)
    contrast[1:K, selected_treatment] <- 1
    for (i in 1:(stage - 1)) contrast[K + i, K + i] <- 1

    bar_mean <- contrast %*% z_means[c(1:K, 1:(stage - 1) * K + selected_treatment)]
  
    bar_cov_matrix <- contrast %*% cov_matrix[c(1:K, 1:(stage - 1) * K + selected_treatment), c(1:K, 1:(stage - 1) * K + selected_treatment)] %*% t(contrast)
    
    lower <- c(rep(0, length(treatments)), cut_off)
    if (stage > 2) lower <- c(rep(0, length(treatments)), lower_boundary[2:(stage - 1)], cut_off)
    lower[which(treatments == selected_treatment)] <- lower_boundary[1]

    upper <- c(rep(Inf, length(treatments)), Inf)
    if (stage > 2) upper <- c(rep(Inf, length(treatments)), upper_boundary[2:(stage - 1)], Inf)
    upper[which(treatments == selected_treatment)] <- upper_boundary[1]
    
    pmvnorm(lower = lower, upper = upper, mean = bar_mean[c(treatments, K + 1:(stage - 1))], sigma = bar_cov_matrix[c(treatments, K + 1:(stage - 1)), c(treatments, K + 1:(stage - 1))])[1]
  }
  excess_alpha <- function(cut_off, alpha_star, treatments, cov_matrix, lower_boundary, upper_boundary, selection_method, K, stage, z_means){# for "all_promising" rule, this gives the cumulative typeI error for 'stage' stages

      # for "select_best" rule, this gives the Type I error spent at the 'stage'th stage
      if (stage == 1) return(1 - alpha_star[1] - pmvnorm(lower = rep(-Inf, length(treatments)), upper = rep(cut_off, length(treatments)), mean = z_means[treatments], sigma = cov_matrix[treatments, treatments])[1])
      if (selection_method == "select_best") return(sum(unlist(lapply(treatments, rejection_paths, cut_off = cut_off, treatments = treatments, cov_matrix = cov_matrix, lower_boundary = lower_boundary, upper_boundary = upper_boundary, K = K, stage = stage, z_means = z_means))) - (alpha_star[stage] - alpha_star[stage - 1])) # any of 'treatments' could be selected, so we add all these probabilities
      if (stage == 2){
          surviving_subsets <- c(list(numeric(0)), lapply(as.list(1:(2 ^ length(treatments) - 1)), get.hyp))  # list all possible subsets of surviving treatments after the first stage
          return(1 - alpha_star[2] - sum(unlist(lapply(surviving_subsets, get_path_prob, cut_off = cut_off, treatments = treatments, cov_matrix = cov_matrix, lower_boundary = lower_boundary, upper_boundary = upper_boundary, K = K, stage = stage, z_means = z_means))))
      }
      surviving_subsets1 <- c(list(numeric(0)), lapply(as.list(1:(2 ^ length(treatments) - 1)), get.hyp)) # all possible subsets of surviving treatments after the first stage
      surviving_subsets2 <- c(list(list(numeric(0))), lapply(surviving_subsets1[-1], function(x) c(list(numeric(0)), lapply(as.list(1:(2 ^ length(x) - 1)), get.hyp)))) # for each possible subset of survivng subsets after stage 1, list the possible subsets still surviving after stage 2
      1 - alpha_star[3] - sum(unlist(Map(function(x, y) sum(unlist(lapply(y, get_path_prob, surviving_subset1 = x, cut_off = cut_off, treatments = treatments, cov_matrix = cov_matrix, lower_boundary = lower_boundary, upper_boundary = upper_boundary, K = K, stage = stage, z_means = z_means))), surviving_subsets1, surviving_subsets2)))
  }

    # give everything the correct name
    alpha_star <- current_mams$alpha_star
    l <- current_mams$l
    u <- current_mams$u
    selection_method <- current_mams$selection
    sample_sizes <- current_mams$sample_sizes
    sample_sizes[completed_stages, ] <- n_obs  # Update given the sample sizes actually observed
    if (!all(diff(sample_sizes) >= 0)) {stop("total sample size per arm cannot decrease between stages.")}
    J <- dim(sample_sizes)[1]
    K <- dim(sample_sizes)[2] - 1
    R <- sample_sizes[, -1] / sample_sizes[1, 1]
    r0 <- sample_sizes[, 1] / sample_sizes[1, 1]

    # get conditional distributions BEFORE seeing the new z scores
    cond_cov_matrix <- cov_matrix <- create_cov_matrix(r0, R)
    cond_mean <- rep(0, K * J)
    if (completed_stages > 1){
        cond_cov_matrix <- create_cond_cov_matrix(cov_matrix, K, completed_stages - 1)
        cond_mean <- create_cond_mean(cov_matrix, K, completed_stages - 1, z_scores = z_scores[[completed_stages - 1]])
    }
    

    # adjust upper boundaries in light of observed sample sizes:
    for (i in 1:(2 ^ K - 1)){ 
        treatments <- intersect(selected_trts[[completed_stages]], get.hyp(i))
        if ((length(treatments > 0)) && (alpha_star[[i]][J] > 0) && (alpha_star[[i]][J] < 1)){
            for (j in completed_stages:J){

                try(new_u <- uniroot(excess_alpha, c(0, 10), alpha_star = alpha_star[[i]][completed_stages:J], treatments = treatments, cov_matrix = cond_cov_matrix, lower_boundary = l[[i]][completed_stages:J], upper_boundary = u[[i]][completed_stages:J], selection_method = selection_method, K = K, stage = j - completed_stages + 1, z_means = cond_mean)$root, silent = TRUE)
                if (is.null(new_u)) {stop("upper boundary not between 0 and 10")}
                u[[i]][j] <- round(new_u, 2)

            }
            l[[i]][J] <- u[[i]][J]
        }
    }
    if (J > completed_stages) {
        cond_cov_matrix <- create_cond_cov_matrix(cov_matrix, K, completed_stages)
        cond_mean <- create_cond_mean(cov_matrix, K, completed_stages, z_scores[[completed_stages]])
    }
    for (i in 1:(2 ^ K - 1)) { # get conditional errors
        treatments <- intersect(selected_trts[[completed_stages]], get.hyp(i))
        if ((length(treatments > 0)) && (alpha_star[[i]][J] > 0) && (alpha_star[[i]][J] < 1)){
            max_z <- max(z_scores[[completed_stages]][treatments])
            best_treatment <- treatments[which.max(z_scores[[completed_stages]][treatments])]
            if (max_z <= u[[i]][completed_stages]) alpha_star[[i]][completed_stages] <- 0
            if (max_z > u[[i]][completed_stages]) {
                alpha[[i]][completed_stages:J] <- 1
                if (J > completed_stages) {
                    l[[i]][(completed_stages + 1):J] <- u[[i]][(completed_stages + 1):J] <- -Inf
                }
            }
            else if (max_z <= l[[i]][completed_stages]){
                alpha_star[[i]][completed_stages:J] <- 0
                if (J > completed_stages) {
                    l[[i]][(completed_stages + 1):J] <- u[[i]][(completed_stages + 1):J] <- Inf
                }
            }
            else if (selection_method == "select_best") {
                for (j in (completed_stages + 1):J){
                    alpha_star[[i]][j] <- excess_alpha(cut_off = u[[i]][j], alpha_star = rep(0, J - completed_stages), treatments = best_treatment, cov_matrix = cond_cov_matrix, lower_boundary = l[[i]][(completed_stages + 1):J], upper_boundary = u[[i]][(completed_stages + 1):J], selection_method = selection_method, K = K, stage = j - completed_stages, z_means = cond_mean) + alpha_star[[i]][j - 1]
                }
            }
            else {
                for (j in (completed_stages + 1):J){
                    alpha_star[[i]][j] <- excess_alpha(cut_off = u[[i]][j], alpha_star = rep(0, J - completed_stages), treatments = treatments, cov_matrix = cond_cov_matrix, lower_boundary = l[[i]][(completed_stages + 1):J], upper_boundary = u[[i]][(completed_stages + 1):J], selection_method = selection_method, K = K, stage = j - completed_stages, z_means = cond_mean) 
                }
            }
        }
    }
    if (is.matrix(n_future)){
        sample_sizes[(completed_stages + 1):J, ] <- n_future
        if (!all(diff(sample_sizes) >= 0)) {stop("total sample size per arm cannot decrease between stages.")}
        R <- sample_sizes[, -1] / sample_sizes[1, 1]
        r0 <- sample_sizes[, 1] / sample_sizes[1, 1]
        cov_matrix <- create_cov_matrix(r0, R)
        cond_cov_matrix <- create_cond_cov_matrix(cov_matrix, K, completed_stages)
        cond_mean <- create_cond_mean(cov_matrix, K, completed_stages, z_scores = z_scores[[completed_stages]])
    }
    if (J > completed_stages){
        for (i in 1:(2 ^ K - 1)){ 
            treatments <- intersect(selected_trts[[completed_stages + 1]], get.hyp(i))
            if ((length(treatments > 0)) && (alpha_star[[i]][J] > 0) && (alpha_star[[i]][J] < 1)){
                for (j in (completed_stages + 1):J){
                    try(new_u <- uniroot(excess_alpha, c(0, 10), alpha_star = alpha_star[[i]][(completed_stages + 1):J], treatments = treatments, cov_matrix = cond_cov_matrix, lower_boundary = l[[i]][(completed_stages + 1):J], upper_boundary = u[[i]][(completed_stages + 1):J], selection_method = selection_method, K = K, stage = j - completed_stages, z_means = cond_mean)$root, silent = TRUE)
                    if (is.null(new_u)) {stop("upper boundary not between 0 and 10")}
                    u[[i]][j] <- round(new_u, 2)
                    
                }
                l[[i]][J] <- u[[i]][J]
            }
        }
    }


    res <- NULL
    res$l <- l
    res$u <- u
    res$sample_sizes <- sample_sizes
    res$K <- K
    res$J <- J
    res$alpha_star <- alpha_star
    res$selection <- selection_method
    res$z_scores <- z_scores
    res$selected_trts <- selected_trts
    class(res) <- "MAMS.step_down"

    return(res)
      
}









