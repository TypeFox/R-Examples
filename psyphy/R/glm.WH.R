glm.WH <- function(formula, data, NumAlt = 2, lambda.init = 0.01,
					interval = c(0, 0.05),
					trace = FALSE, tol = 1e-6, ...) 
{
	p.l <- function(lambda = lambda, m = NumAlt) {
		probit.lambda(m, lambda)
		}
	if (missing(data)) 
        data <- environment(formula)
    est.glm <- glm(formula, 
    		family = binomial(probit.lambda(NumAlt, lambda.init)),
    		data = data, ...)
    ll <- function(lam, NumAlt, x) {
    	rr <- x$model[[1]]
    	-sum(rr[, 1] * 
    		log(p.l(lam, NumAlt)$linkinv(x$linear.predictors)) +
    		 rr[, 2] * log(1 - 
    			p.l(lam, NumAlt)$linkinv(x$linear.predictors)))
    	}
    dlogcur <- dd <- -as.vector(logLik(est.glm))
    new.glm <- est.glm
    n <- 0
    while (dd > tol) {
    	n <- n + 1
    	lam.c <- optimize(ll, interval = interval, 
    			NumAlt = NumAlt,
    			x = new.glm)$minimum
    	new.glm <- glm(formula, 
    		family = binomial(probit.lambda(NumAlt, lam.c)),
    		data = data, ...)
    	dd <- abs(-as.vector(logLik(new.glm)) - dlogcur)
    	dlogcur <- -as.vector(logLik(new.glm))
    	if (trace) print(data.frame(n = n, logLik = dd, lambda = lam.c))
    	}
    new.glm$lambda <- lam.c
    cat("lambda = \t", lam.c, "\n")
    new.glm$df.residual <- new.glm$df.residual - 1
    class(new.glm) <- c("lambda", "glm", "lm")
    new.glm
}

