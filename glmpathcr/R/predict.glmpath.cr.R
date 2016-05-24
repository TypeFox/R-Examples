predict.glmpath.cr<-
function (object, newx=NULL, which = "BIC", type = "class", ...) 
{
    summary.glmpath <- summary(object)
	newy<-object$y
	if (is.null(newx)) newx<-object$x
	method <- object$method
    k <- length(unique(newy))
	if (class(newx)=="numeric") newx<-matrix(newx,ncol=dim(object$x)[2])
    n <- dim(newx)[1]
    p <- dim(newx)[2]
    if (c("BIC", "AIC")[charmatch(which, c("BIC", "AIC"))] == 
        "BIC") {
        min <- which.min(summary.glmpath[, which])
    }
    if (c("BIC", "AIC")[charmatch(which, c("BIC", "AIC"))] == 
        "AIC") {
        min <- which.min(summary.glmpath[, which])
    }
    best <- as.numeric(gsub("Step ", "", rownames(summary.glmpath)[min]))
    beta <- coef(object, s = best, type = "coefficients", mode = "step")
    logit <- matrix(rep(0, n * (k - 1)), ncol = k - 1)
    for (j in 1:(k - 1)) {
        cp <- paste("cp", j, sep = "")
        logit[, j] <- beta[dimnames(beta)[[2]] == "Intercept"] + 
		beta[dimnames(beta)[[2]] == cp] + beta[2:(p + 1)] %*% 
		t(as.matrix(newx))
    }
    delta <- matrix(rep(0, n * (k - 1)), ncol = k - 1)
    for (j in 1:(k - 1)) {
        delta[, j] <- exp(logit[, j])/(1 + exp(logit[, j]))
    }
    pi <- matrix(rep(0, n * k), ncol = k)
    minus.delta <- 1 - delta
    if (c("backward", "forward")[charmatch(method, c("backward", 
													 "forward"))] == "backward") {
        for (j in k:2) {
            if (j == k) {
                pi[, j] <- delta[, k - 1]
            }
            else if (class(minus.delta[, j:(k - 1)]) == "numeric") {
                pi[, j] <- delta[, j - 1] * minus.delta[, j]
            }
            else if (dim(minus.delta[, j:(k - 1)])[2] >= 2) {
                pi[, j] <- delta[, j - 1] * apply(minus.delta[, 
												  j:(k - 1)], 1, prod)
            }
        }
        if (n==1) pi[, 1] <- 1 - sum(pi[, 2:k]) else pi[, 1] <- 1 - apply(pi[, 2:k], 1, sum)
    }
    if (c("backward", "forward")[charmatch(method, c("backward", 
													 "forward"))] == "forward") {
        for (j in 1:(k - 1)) {
            if (j == 1) {
                pi[, j] <- delta[, j]
            }
            else if (j == 2) {
                pi[, j] <- delta[, j] * minus.delta[, j - 1]
            }
            else if (j > 2 && j < k) {
                pi[, j] <- delta[, j] * apply(minus.delta[, 1:(j - 
															   1)], 1, prod)
            }
        }
        if (n==1) pi[, k] <- 1 - sum(pi[, 1:(k - 1)]) else pi[, k] <- 1 - apply(pi[, 1:(k - 1)], 1, sum)
    }
	p.class <- apply(pi, 1, which.max)
	if (sum(match(class(object$y), "factor"),na.rm=TRUE)>0) class<-levels(object$y)[p.class] else class<-p.class
    if (sum(match(class(object$y), "factor"),na.rm=TRUE)>0) dimnames(pi)[[2]] <-levels(object$y) else dimnames(pi)[[2]] <- unique(object$y)
    if (c("class", "probs")[charmatch(type, c("class", "probs"))] == 
        "class") {
        object <- class
    }
    if (c("class", "probs")[charmatch(type, c("class", "probs"))] == 
        "probs") {
        object <- pi
    }
    object
}