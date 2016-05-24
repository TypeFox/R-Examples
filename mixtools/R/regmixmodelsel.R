regmixmodel.sel = function (x, y, w=NULL, k = 2, type = c("fixed", "random", "mixed"), ...) 
{
    aic <- NULL
    bic <- NULL
    caic <- NULL
    icl <- NULL
    type <- match.arg(type)
        AIC <- function(emout, p) {
            emout$loglik - (p - 1)
        }
        BIC <- function(emout, p, n) {
            emout$loglik - log(n) * (p - 1)/2
        }
        CAIC <- function(emout, p, n) {
            emout$loglik - (log(n) + 1) * (p - 1)/2
        }
        ICL <- function(emout, p, n) {
            BIC(emout, p, n) - sum(emout$lambda * log(emout$lambda))
        }

    if (type == "fixed") {
	p <- function(emout) {length(emout$beta) + length(emout$sigma) + (length(emout$scale)-1)*is.null(emout$scale) +length(emout$lambda)}
	n <- length(y)
        for (i in 1:k) {
            if (i == 1) {
                a <- glm(y ~ x)
                beta <- matrix(a$coef, ncol = 1)
                loglik <- log(prod(dnorm(y, mean = beta[1, ] + 
                  as.matrix(x) %*% beta[2:nrow(beta), ], sd = sd(a$res))))
                emout <- list(beta = beta, sigma = sd(a$res), 
                  lambda = 1, loglik = loglik)
            }
            else emout <- regmixEM(y, x, k = i, ...)
		P = p(emout)
            aic[i] <- AIC(emout, p=P)
            bic[i] <- BIC(emout, p=P, n=n)
            caic[i] <- CAIC(emout, p=P, n=n)
            icl[i] <- ICL(emout, p=P, n=n)
        }
    }
    else if (type == "random") {
	p <- function(emout) {
		if(is.list(emout$R)){ 
			m<-nrow(emout$R[[1]])
			v<-length(emout$R)
		} else {
		m<-nrow(emout$R)
		v<-1}
length(emout$mu) +(v*m*(m+1)/2)+ length(emout$sigma)+  length(emout$lambda)}
	n <- sum(sapply(y,length))
        for (j in 1:k) {
            if (j == 1) {
                a <- lapply(1:length(y), function(i) glm(y[[i]] ~ 
                  x[[i]]))
                aa <- sapply(1:length(a), function(i) as.vector(a[[i]]$res))
                b <- sapply(1:length(a), function(i) as.matrix(a[[i]]$coef))
                mu <- apply(b, 1, mean)
                R <- cov(t(b))
                a.res <- NULL
                for (i in 1:length(aa)) {
                  a.res <- c(a.res, aa[[i]])
                }
                sd.a <- sd(a.res)
                loglik <- sapply(1:length(y), function(i) dmvnorm(as.vector(y[[i]]), 
                  b[1, i] + as.matrix(x[[i]]) %*% b[2:nrow(b), 
                    i], sd.a^2 * diag(length(y[[i]]))))
                loglik <- log(prod(loglik))
                emout <- list(mu = mu, R = R, sigma = sd.a, lambda = 1, 
                  loglik = loglik)
            }
            else emout <- regmixEM.mixed(y, x, k = j, ...)
		P = p(emout)
            aic[j] <- AIC(emout, p=P)
            bic[j] <- BIC(emout, p=P, n=n)
            caic[j] <- CAIC(emout, p=P, n=n)
            icl[j] <- ICL(emout, p=P, n=n)
        }
    }
    else if (type == "mixed") {
	p <- function(emout) {
		if(is.list(emout$R)){ 
			m<-nrow(emout$R[[1]])
			v<-length(emout$R)
		} else {
		m<-nrow(emout$R)
		v<-1}
length(emout$alpha)+ length(emout$mu) +(v*m*(m+1)/2)+ length(emout$sigma)+  length(emout$lambda)}
	n <- sum(sapply(y,length))
        for (i in 1:k) {
            if (i == 1) {
                a <- lapply(1:length(y), function(i) glm(y[[i]] ~ 
                  x[[i]] + w[[i]]))
                aa <- sapply(1:length(a), function(i) as.vector(a[[i]]$res))
                alpha <- sapply(1:length(a), function(i) matrix(a[[i]]$coef[(ncol(x[[i]]) + 
                  2):(ncol(x[[i]]) + 1 + ncol(w[[i]]))]))
                alpha <- apply(as.matrix(alpha), 1, mean)
                b <- sapply(1:length(a), function(i) matrix(a[[i]]$coef[1:(ncol(x[[i]]) + 
                  1)]))
                mu <- apply(b, 1, mean)
                R <- cov(t(b))
                a.res <- NULL
                for (i in 1:length(aa)) {
                  a.res <- c(a.res, aa[[i]])
                }
                sd.a <- sd(a.res)
                loglik <- sapply(1:length(y), function(i) dmvnorm(as.vector(y[[i]]), 
                  b[1, i] + as.matrix(x[[i]]) %*% b[2:nrow(b), 
                    i], sd.a^2 * diag(length(y[[i]]))))
                loglik <- log(prod(loglik))
                emout <- list(mu = mu, R = R, alpha = alpha, 
                  sigma = sd.a, lambda = 1, loglik = loglik)
            }
            else emout <- regmixEM.mixed(y, x, k = i, ...)
		P = p(emout)
            aic[i] <- AIC(emout, p=P)
            bic[i] <- BIC(emout, p=P, n=n)
            caic[i] <- CAIC(emout, p=P, n=n)
            icl[i] <- ICL(emout, p=P, n=n)
        }
    }
    out = rbind(aic, bic, caic, icl)
    Winner = apply(out, 1, function(x) (1:length(x))[x == max(x)])
    colnames(out) = 1:k
    rownames(out) = c("AIC", "BIC", "CAIC", "ICL")
    cbind(out, Winner)
}