Mstep.glm <- function(x, cond, pm, pn, family, link){
    m <- ncol(cond$u)
    beta0 <- rep(NA, m)
    beta1 <- rep(NA, m)
    sigma <- rep(NA, m)
    if (family=="binomial") x <- cbind(x, pn$size-x)
    for (k in 1:m){
        w <- cond$u[,k]
        fmly <- eval(parse(text=paste(family,
                     "(link=\"", link, "\")", sep="")))
        y <- glm(x ~ pn$x1, weights=w, family=fmly,
                        start=c(pm$beta0[k], pm$beta1[k]))
        fitval <- fitted(y)
        if (family=="poisson") var <- fitval
        if (family=="gaussian") var <- 1
        if (family=="Gamma") var <- fitval^2
        if (family=="binomial") var <- fitval*(1-fitval)/pn$size
        if (family!="binomial") sigma[k] <- sqrt(sum(w*(x-fitval)^2/var)/sum(w))
        else sigma[k] <- sqrt(sum(w*(x[,1]/pn$size-fitval)^2/var)/sum(w))
        beta0[k] <- as.vector(coefficients(y))[1]
        beta1[k] <- as.vector(coefficients(y))[2]
    }
    return(list(beta0=beta0, beta1=beta1, sigma=sigma))
}

