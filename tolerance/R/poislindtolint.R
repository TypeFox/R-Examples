poislindtol.int <- function(x, m = NULL, alpha=0.05, P=0.99, side = 1, ...){
    if (side != 1 && side != 2) {
        stop(paste("Must specify a one-sided or two-sided procedure!", 
            "\n"))
    }
	n <- ifelse(class(x)=="table",sum(x),length(x))
    if (side == 2) {
        alpha <- alpha/2
        P <- (P + 1)/2
    }
	if(is.null(m)) m <- n
	out <- poislind.ll(x, ...)
	theta <- as.numeric(stats4::coef(out))
	CI <- pmax(theta+c(-1,1)*qnorm(1-alpha)*sqrt(stats4::vcov(out)[1])*sqrt(n/m),0)
	lower <- max(qpoislind(1 - P, theta = CI[2]), 0)
    upper <- qpoislind(P, theta = CI[1])
    if (side == 2) {
        alpha <- 2 * alpha
        P <- (2 * P) - 1
    }
    temp <- data.frame(cbind(alpha, P, theta, lower, upper))
    if (side == 2) {
        colnames(temp) <- c("alpha", "P", "theta", "2-sided.lower", "2-sided.upper")
    }
    else {
        colnames(temp) <- c("alpha", "P", "theta", "1-sided.lower", "1-sided.upper")
    }
    temp	
}

