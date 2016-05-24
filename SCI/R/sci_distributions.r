#########################################################
### Generalized Logistic
#########################################################

dgenlog <- function(x, shape, scale, location) {
    e <- location
    a <- scale
    k <- shape
    x.correct <- (1-(k*(x-e))/a)>0  & !is.na(x)
    if(k==0) {
        y <- (x-e)/a
    } else {
        y <- -k^(-1)*log(1-(k*(x[x.correct]-e))/a)
    }
    top <- a^(-1)*exp(-(1-k)*y)
    bottom <- (1+exp(-y))^2
    pdf <- vector()
    pdf[x.correct] <- top/bottom
    pdf[x.correct==FALSE] <- .Machine$double.eps
    pdf[is.na(x)] <- NA
    return(pdf)
}

pgenlog <- function(q, shape, scale, location) {
    e <- location
    a <- scale
    k <- shape
    q.correct <- (1-(k*(q-e))/a)>0 & !is.na(q)
    if(k==0) {
        y <- (q-e)/a
    } else {
        y <- -k^(-1)*log(1-(k*(q[q.correct]-e))/a)
    } 
    cdf <- vector()
    cdf[q.correct] <- 1/(1+exp(-y))
    if(k >0) {
	cdf[q.correct==FALSE] <- 1-.Machine$double.eps
    } else {
	cdf[q.correct==FALSE] <- .Machine$double.eps
    }
    cdf[is.na(q)] <- NA
    return(cdf)
}

qgenlog <- function(p, shape, scale, location) {
    e <- location
    a <- scale
    k <- shape
    p.correct <- !is.na(p)
    x <- vector()
    if(k==0) {
        x[p.correct] <- e-a*log((1-p[p.correct])/p[p.correct])
    } else {
        x[p.correct] <- e+(a/k)*(1-((1-p[p.correct])/p[p.correct])^k)
    }
    x[is.na(p)] <- NA
    return(x)
}


#########################################################
### Pearson Type III
#########################################################

dpe3 <- function(x, shape, scale, location) {
    e <- location
    a <- scale
    k <- shape
    x.correct <- !is.na(x)
    if (k == 0) {
        pdf.correct = dnorm((x[x.correct] - e)/a)/a
    } else {
        alpha <- 4/k^2
        tmp <- gamma(alpha)
	beta <- 0.5 * a * abs(k)
	xi <- e - 2 * a/k
	if (k > 0) {
            Y <- x[x.correct] - xi
            pdf.correct = (dgamma((Y)/beta, alpha))/beta
	}  else {
            Y <- xi - x[x.correct]
            pdf.correct = (dgamma((Y)/beta, alpha))/beta
	}
    }
    pdf <- vector()
    pdf.correct[pdf.correct <= 0] <- .Machine$double.eps
    pdf[x.correct] <- pdf.correct
    pdf[is.na(pdf)] <- .Machine$double.eps
    pdf[is.na(x)] <- NA
    return(pdf)
}

ppe3 <- function(q, shape, scale, location) {
    e <- location
    a <- scale
    k <- shape 
    q.correct <- !is.na(q)
    if (k == 0 | 4/k^2 > 170) {
	cdf.correct <- pnorm((q[q.correct] - e)/a)
    } else {
        alpha <- 4/k^2
        tmp <- gamma(alpha)
	beta <- 0.5 * a * abs(k)
	xi <- e - 2 * a/k
	if (k > 0) {
            Y <- q[q.correct] - xi
            cdf.correct <- pgamma((Y)/beta, alpha)
	}  else {
            Y <- xi - q[q.correct]
            cdf.correct <- 1-pgamma((Y)/beta, alpha)
	} }
    cdf.correct[cdf.correct == 0] <- .Machine$double.eps
    cdf.correct[cdf.correct == 1] <- 1-.Machine$double.eps
    cdf <- vector()
    cdf[q.correct] <- cdf.correct
    cdf[is.na(q)] <- NA
    return(cdf)
}


qpe3 <- function(p, shape, scale, location) {
    e <- location
    a <- scale
    k <- shape
    p.correct <- !is.na(p)
    x <- vector()
    if (4/k^2 > 170) {
	x[p.correct] <- e + a * qnorm(p[p.correct])
    } else {
        alpha <- 4/k^2
	beta <- abs(0.5 * a * k)
	if (k > 0) {
            x[p.correct] <- e - alpha*beta + qgamma(p[p.correct], alpha, scale=beta)
	}  else {
            x[p.correct] <- e - alpha*beta - qgamma(1 - p[p.correct], alpha, scale=beta)
	} }
    x[is.na(p)] <- NA
    return(x)
}

