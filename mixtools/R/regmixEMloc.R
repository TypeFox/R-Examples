regmixEM.loc=function (y, x, lambda = NULL, beta = NULL, sigma = NULL, k = 2, 
    addintercept = TRUE, kern.l = c("Gaussian", "Beta", "Triangle", 
        "Cosinus", "Optcosinus"), epsilon = 1e-08, maxit = 10000, 
    kernl.g = 0, kernl.h = 1, verb = FALSE) 
{
    diff <- 1
    iter <- 0
    x.1 <- cbind(1, x)
    n <- length(y)
    kern.l <- match.arg(kern.l)
    out.EM <- regmixEM.lambda(y, x, lambda = lambda, beta = beta, 
        sigma = sigma, k = k, epsilon = epsilon, maxit = maxit)
    ll = out.EM$loglik
    restarts <- 0
    while (diff > epsilon && iter < maxit) {
        old.l.x = out.EM$lambda
        old.loglik = out.EM$loglik
        l.x = matrix(nrow = n, ncol = (k - 1))
        for (i in 1:(k - 1)) {
            l.x[, i] <- c((cbind(rep(1,n),0)*t(apply(matrix(x, ncol = 1), 1, lambda, 
                z = out.EM$post[, i], xi = x, kernel = kern.l, 
                g = kernl.g, h = kernl.h)))[,1])
	    l.x[,i]=apply(cbind(l.x[,i]),1,max,0)
	    l.x[,i]=apply(cbind(l.x[,i]),1,min,1)
        }
        l.x <- cbind(l.x, 1 - apply(l.x, 1, sum))
        out.EM.loc <- regmixEM.lambda(y, x, beta = out.EM$beta, 
            sigma = out.EM$sigma, lambda = l.x, k = k)
        loglik.loc <- out.EM.loc$loglik
        out.EM.old <- regmixEM.lambda(y, x, beta = out.EM$beta, 
            sigma = out.EM$sigma, lambda = old.l.x, k = k)
        loglik.old <- out.EM.old$loglik
        if (loglik.loc > old.loglik) {
            out.EM <- out.EM.loc
        }
        else out.EM <- out.EM.old
        loglik.chosen <- out.EM$loglik
        ll <- c(ll, loglik.chosen)
        diff <- loglik.chosen - old.loglik
        if (diff < 0) {
            cat("Generating new starting values...", "\n")
            out.EM <- regmixEM.lambda(y, x, lambda = lambda, 
                beta = beta, sigma = sigma, k = k, epsilon = epsilon, 
                maxit = maxit)
            restarts <- restarts + 1
            if (restarts > 15) 
                stop("Too many tries!")
            iter <- 0
            diff <- 1
        }
        else {
            iter <- iter + 1
            if (verb) {
                cat("iteration=", iter, "diff=", diff, "log-likelihood", 
                  loglik.chosen, "\n")
            }
        }
    }
    if (iter == maxit) {
        cat("WARNING! NOT CONVERGENT!", "\n")
    }
    cat("number of overall iterations=", iter, "\n")
    a = list(x = x, y = y, lambda.x = out.EM$lambda, beta = out.EM$beta, 
        sigma = out.EM$sigma, loglik = loglik.chosen, posterior = out.EM$post, 
        all.loglik = ll, restarts = restarts, ft = "regmixEM.loc")
    class(a) = "mixEM"
    a
}

