regmix.init<-function (y, x, lambda = NULL, beta = NULL, s = NULL, k = 2, 
    addintercept = TRUE, arbmean = TRUE, arbvar = TRUE) 
{
    x <- as.matrix(x)
    n <- length(y)
    p <- ncol(x)
    if (addintercept) {
        x = x[, -1]
    }
    else x = x
    if (is.null(lambda)) {
        lambda = runif(k)
        lambda = lambda/sum(lambda)
    }
    else k = length(lambda)
    A = round(lambda * n)
	while(min(A)<4){
	       lambda = runif(k)
             lambda = lambda/sum(lambda)
		 A = round(lambda * n)
	}
    A[k] = n - sum(A[1:(k - 1)])
    B = c(0, cumsum(A))
    w = cbind(y, x)
    w = w[order(w[, 1]), ]
    w.bin = list()
    for (j in 1:k) {
        w.bin[[j]] <- w[(B[j] + 1):B[j + 1], ]
    }
    if (addintercept) {
        lm.out <- lapply(1:k, function(i) lm(w.bin[[i]][, 1] ~ 
            w.bin[[i]][, 2:p]))
    }
    else lm.out <- lapply(1:k, function(i) lm(w.bin[[i]][, 1] ~ 
        w.bin[[i]][, 2:(p + 1)] - 1))

        s.hyp = lapply(lm.out, anova)
        s.hyp = as.vector(sqrt(sapply(1:k, function(i) s.hyp[[i]]$Mean[length(s.hyp[[i]]$Mean)])))
        s.hyp[(s.hyp<=0)|(is.na(s.hyp)==1)]=1
    if (is.null(s)) {
        if (arbvar) {
            s = 1/rexp(k, rate = s.hyp)
        }
        else {
            s.hyp = mean(s.hyp)
            s = 1/rexp(1, rate = s.hyp)
        }
    }
    if (is.null(s) == FALSE && arbvar == TRUE) {
        k = length(s)
    }
    if (is.null(beta)) {
	  x.x <- solve(t(x)%*%x)
        beta.hyp = matrix(sapply(lm.out, coef), ncol = k)
        beta = matrix(0, nrow = p, ncol = k)
        for (j in 1:k) {
            beta[, j] = rnorm(p, mean = as.vector(beta.hyp[, 
#                j]), sd = (s[arbvar * j + (1 - arbvar)]*sqrt(diag(x.x))) ) 
                j]), sd = (s.hyp[arbvar * j + (1 - arbvar)]*sqrt(diag(x.x))) ) 
        }
        if (arbmean == FALSE) {
            beta = apply(beta, 1, mean)
        }
    }
    if (is.null(beta) == FALSE && arbmean == TRUE) {
        k = ncol(beta)
    }
    list(lambda = lambda, beta = beta, s = s, k = k)
}