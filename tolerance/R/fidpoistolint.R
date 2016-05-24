fidpoistol.int <- function(x1, x2, n1, n2, m1 = NULL, m2 = NULL, FUN, alpha = 0.05, P = 0.99, side = 1, K = 1000, B = 1000)
{
	FUN <- match.fun(FUN)
    if (side != 1 && side != 2) {
        stop(paste("Must specify a one-sided or two-sided procedure!", 
            "\n"))
    }
    if (side == 2) {
        alpha <- alpha/2
        P <- (P + 1)/2
    }
	if(is.null(m1)) m1 <- n1
	if(is.null(m2)) m2 <- n2
	est <- FUN(x1/n1,x2/n2)
	Q.l1 <- rchisq(K,2*x1+1)/(2*n1)
	Q.l2 <- rchisq(K,2*x2+1)/(2*n2)
	TEMP <- sapply(1:1000, function(i) quantile(FUN(rpois(B,m1*Q.l1[i])/m1,rpois(B,m2*Q.l2[i])/m2),c(1-P,P),na.rm=TRUE ))   
	lower <- quantile(TEMP[1,],alpha)
	upper <- quantile(TEMP[2,],1-alpha)
    if (side == 2) {
        alpha <- 2 * alpha
        P <- (2 * P) - 1
    }
	FUN.name <- rFUN(FUN,r1="lambda.1",r2="lambda.2")
    temp <- data.frame(cbind(alpha, P, est, lower, upper))
    if (side == 2) {
        colnames(temp) <- c("alpha", "P", "fn.est", "2-sided.lower", "2-sided.upper")
    }
    else {
        colnames(temp) <- c("alpha", "P", "fn.est", "1-sided.lower", "1-sided.upper")
    }
    rownames(temp) <- 1
    temp <- list(tol.limits=temp,fn=noquote(FUN.name))
    temp
}