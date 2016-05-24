logitIC <- function(x, y, candidate_models, psi,
                    type = c("BIC", "AIC"), prior = TRUE, reduce_bias = FALSE) {
  p <- NCOL(x)
  n <- length(y)
  type <- match.arg(type)
  n_mo <- NROW(candidate_models)
  sk <- rowSums(candidate_models)
  ik <- rep(NA, n_mo)
  if (any(candidate_models[1, ] == 1)) {
    for (i in seq(n_mo)) {
      glmfit <- if(reduce_bias==TRUE) brglm(y ~ x[, candidate_models[i, ] == 1], family = binomial) else glm(y ~ x[, candidate_models[i, ] == 1], family = binomial)
      ik[i] <- if (type == "BIC") extractAIC(glmfit, k=log(n))[2] else extractAIC(glmfit)[2]
    }   
  } else {
    glmfit <- if(reduce_bias==TRUE) brglm(y ~ 1, family=binomial) else glm(y ~ 1, family=binomial)
    ik[1] <- if (type == "BIC") extractAIC(glmfit, k=log(n))[2] else extractAIC(glmfit)[2]
    for (i in seq(2, n_mo)) {
      glmfit <- if(reduce_bias==TRUE) brglm(y ~ x[, candidate_models[i, ] == 1], family = binomial) else glm(y ~ x[, candidate_models[i, ] == 1], family = binomial)
      ik[i] <- if (type == "BIC") extractAIC(glmfit, k=log(n))[2] else extractAIC(glmfit)[2]
    }
  }
  if (prior == TRUE) {
    ck <- ck_compute(n_mo, sk, p)
    ik <- ik + 2 * psi * ck
  }
  ik <- ik - min(ik)
  weight <- exp(-ik/2)/sum(exp(-ik/2))
  list(weight = weight)
}
