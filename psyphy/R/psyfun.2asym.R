psyfun.2asym <- function(formula, data, link = logit.2asym,
	init.g = 0.01, init.lam = 0.01,
	trace = FALSE, tol = 1e-6, mxNumAlt = 50, ...) {
	p.l <- function(p) {
        link(p[1], p[2])
        }
	if (missing(data)) 
        data <- environment(formula)
    est.glm <- glm(formula, 
    	family = binomial(link(g = init.g, lam = init.lam)), 
    	data = data, ...)
    ll <- function(p, x) {
    	p <- plogis(p)
        rr <- x$model[[1]]
        -sum(rr[, 1] * log(p.l(p)$linkinv(x$linear.predictors)) + 
            rr[, 2] * log(1 - p.l(p)$linkinv(x$linear.predictors)))
    }
    dlogcur <- dd <- -as.vector(logLik(est.glm))
    new.glm <- est.glm
    n <- 0
    p <- c(init.g, init.lam)
	while (dd > tol) {
        n <- n + 1
        p <- qlogis(p)
        p.opt <- optim(p, ll, x = new.glm, hessian = TRUE)
        p <- plogis(p.opt$par)
        new.glm <- glm(formula, 
        	family = binomial(link(g = p[1], lam = p[2])), 
        	data = data, ...)
        dd <- abs(-as.vector(logLik(new.glm)) - dlogcur)/dlogcur
        dlogcur <- -as.vector(logLik(new.glm))
        if (trace) 
            print(data.frame(n = n, logLik = dd, lambda = p[2], gamma = p[1]))
        if (n > mxNumAlt) {
        	print("Number of iterations exceeded without finding best fit. \n")
        	break}    
    }
    p.svd <- svd(p.opt$hessian)
    SEp <- sqrt(diag(with(p.svd, v %*% diag(1/ifelse(zapsmall(d), d, Inf)) %*% t(u))))
#    X <<- SEp
    SEp[zapsmall(SEp) == 0] <- NA
	new.glm$lambda <- p[2]
	new.glm$SElambda <- SEp[2]
	new.glm$gam <- p[1]
	new.glm$SEgam <- SEp[1]
	cat("lambda = \t", p[2], "\t", "gamma = ", p[1], "\n")
	cat("+/-SE(lambda) = \t(", plogis(qlogis(p[2]) + c(-SEp[2], SEp[2])), ")\n")
	cat("+/-SE(gamma) = \t(", plogis(qlogis(p[1]) + c(-SEp[1], SEp[1])), ")\n")
	new.glm$df.residual <- new.glm$df.residual - 2
	new.glm$call[[3]][[2]][[1]] <- as.name(substitute(link))
	class(new.glm) <- c("lambda", "glm", "lm")
	new.glm
}