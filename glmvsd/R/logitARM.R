logitARM <- function(x, y, candidate_models, n_train, no_rep, psi, prior = TRUE, reduce_bias = FALSE) {
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
        glmfit <- if(reduce_bias==TRUE) brglm(y[tindex] ~ x[tindex, varindex], family = binomial) else glm(y[tindex] ~ x[tindex, varindex], family = binomial)
        if (any(is.na(glmfit$coef))) {
          lw_num[j] <- -Inf  
        } else {
          dk <- as.vector(cbind(1, x[-tindex, varindex]) %*% glmfit$coef)
          fk <- ifelse(dk < 0, log(1 + exp(dk)), dk + log(1 + exp(-dk)))
          lw_num[j] <- sum(y[-tindex] * dk) - sum(fk)
        }
      }
    } else {
      if (mean(y[tindex])==0 | mean(y[tindex]==1)) {
        lw_num[1] <- -Inf
      } else {
        lw_num[1] <- sum(log(mean(y[tindex])) * (y[-tindex]) + log(1 - mean(y[tindex])) * (1 - y[-tindex]))
      }
      for (j in seq(2, n_mo)) {
        varindex <- (candidate_models[j, ] == 1)
        glmfit <- if(reduce_bias==TRUE) brglm(y[tindex] ~ x[tindex, varindex], family = binomial) else glm(y[tindex] ~ x[tindex, varindex], family = binomial) 
        if(any(is.na(glmfit$coef))) {
          lw_num[j] <- -Inf  
        } else {
          dk <- as.vector(cbind(1, x[-tindex, varindex]) %*% glmfit$coef)
          fk <- ifelse(dk < 0, log(1 + exp(dk)), dk + log(1 + exp(-dk)))
          lw_num[j] <- sum(y[-tindex] * dk) - sum(fk)
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
