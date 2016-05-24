lsIC <- function(x, y, candidate_models, psi, 
	type = c("BIC","AIC"), prior = TRUE) {
    p <- NCOL(x)
    n <- length(y)
    n_mo <- NROW(candidate_models)
    sk <- rowSums(candidate_models)
    ik <- rep(NA, n_mo)
    if (any(candidate_models[1, ] == 1)) {
        for (i in 1:n_mo) {
            LSL <- lm(y ~ x[, candidate_models[i, ] == 1])
            rss <- crossprod(summary(LSL)$res, summary(LSL)$res)
			if(type == "BIC") ik[i] <- n * log(rss/n) + sk[i] * log(n)
            else ik[i] <- n * log(rss/n) + sk[i] * 2
        }
    } else {
        rss <- sum((y - mean(y))^2)
        ik[1] <- n * log(rss/n) + sk[1] * log(n)
        for (i in 2:n_mo) {
            LSL <- lm(y ~ x[, candidate_models[i, ] == 1])
            rss <- crossprod(summary(LSL)$res, summary(LSL)$res)
            if(type == "BIC") ik[i] <- n * log(rss/n) + sk[i] * log(n)
            else ik[i] <- n * log(rss/n) + sk[i] * 2
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