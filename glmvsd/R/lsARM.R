lsARM <- function(x, y, candidate_models, n_train, no_rep, psi, prior = TRUE) {
  p <- NCOL(x)
  n <- length(y)
  n_mo <- NROW(candidate_models)
  sk <- rowSums(candidate_models)
  
  wt_calc <- function(rep_id) {
	lw_num <- rep(NA, n_mo)
    tindex <- sample(n, n_train, replace = FALSE)
    if (any(candidate_models[1, ] == 1)) {
      for (j in seq(n_mo)) {
        varindex <- (candidate_models[j, ] == 1)
        glmfit <- lm(y[tindex] ~ x[tindex, varindex])
        sigma_k <- summary(glmfit)$sigma
        if(any(is.na(glmfit$coef))) {
          lw_num[j] <- -Inf
        } else {
          dk <- sum((y[-tindex] - cbind(1, x[-tindex, varindex]) %*% glmfit$coef)^2)
		  lw_num[j] <- (-n/2) * log(sigma_k) - ((sigma_k)^(-2)) * dk/2
        }
      }
    } else {
      dk <- sum((y[-tindex] - mean(y[tindex]))^2)
      sigma_k <- sd(y[tindex])
	  lw_num[1] <- (-n/2) * log(sigma_k) - ((sigma_k)^(-2)) * dk/2
      for (j in seq(2, n_mo)) {
        varindex <- (candidate_models[j, ] == 1)
        glmfit <- lm(y[tindex] ~ x[tindex, varindex])
        sigma_k <- summary(glmfit)$sigma  	
        if (any(is.na(glmfit$coef))) {
          lw_num[j] <- -Inf
        } else {
          dk <- sum((y[-tindex] - cbind(1, x[-tindex, varindex]) %*% glmfit$coef)^2)
		  lw_num[j] <- (-n/2) * log(sigma_k) - ((sigma_k)^(-2)) * dk/2
        }
      }
    }
    return(lw_num)
  }
  lw_num <- matrix(unlist(mclapply(seq(no_rep), wt_calc)), nrow = no_rep, ncol = n_mo, byrow = TRUE)
  if (prior == TRUE) {
    ck <- ck_compute(n_mo, sk, p)
    lw_num <- sweep(lw_num, MARGIN = 2, psi * ck, "-")
  }
  lw_num <- sweep(lw_num, MARGIN = 1, apply(lw_num, 1, max), "-")
  w_num <- exp(lw_num)
  weight <- colMeans(w_num/rowSums(w_num))
  list(weight = weight)
}