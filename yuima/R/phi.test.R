
phi.test <- function(yuima, H0, H1, phi, print=FALSE,...){
    
    phiname <- deparse(substitute(phi))
    if(missing(phi)){
        phi <- function(x) 1-x+x*log(x)
        phiname <- "1-x+x*log(x)"
    }
	d.size <- yuima@model@equation.number
	n <- length(yuima)[1]
	
	env <- new.env()
	assign("X",  as.matrix(onezoo(yuima)), envir=env)
	assign("deltaX",  matrix(0, n-1, d.size), envir=env)
	for(t in 1:(n-1))
	env$deltaX[t,] <- env$X[t+1,] - env$X[t,]
	
	assign("h", deltat(yuima@data@zoo.data[[1]]), envir=env)
	assign("time", as.numeric(index(yuima@data@zoo.data[[1]])), envir=env) 
    est <- FALSE
    if(missing(H1)){
     cat("\nestimating parameters via QMLE...\n")
     H1 <- coef(qmle(yuima, ...))
        est <- TRUE
    }
	g0 <- quasiloglvec(yuima=yuima, param=H0, print=print, env)
	g1 <- quasiloglvec(yuima=yuima, param=H1, print=print, env)
    y <- exp(g1-g0)
    div <- mean(phi(y), na.rm=TRUE)
    stat <- 2*sum(phi(y), na.rm=TRUE)
    df <- length(H0)
    val <- list(div=div, stat=stat, H0=H0, H1=H1, phi=phiname, phi=phi, pvalue=1-pchisq(stat, df=df), df=df,est=est)
    attr(val, "class") <- "phitest"
    return( val )
}



print.phitest <- function(x,...){
    
    symnum(x$pvalue, corr = FALSE, na = FALSE, 
    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
    symbols = c("***", "**", "*", ".", " "))->Signif
    cat(sprintf("Phi-Divergence test statistic based on phi = '%s'\n",x$phi) )
    nm <- names(x$H0)
    cat("H0: ")
    cat(sprintf("%s = %s", nm, format(x$H0,digits=3,nsmall=3)))
    cat("\nversus\n")
    cat("H1: ")
    cat(sprintf("%s = %s", nm, format(x$H1[nm],digits=3,nsmall=3)))
    if(x$est)
    cat("\nH1 parameters estimated using QMLE")
    cat(sprintf("\n\nTest statistic = %s, df = %d, p-value = %s %s\n", format(x$stat,digits=2,nsmall=3), x$df, format(x$pvalue), Signif))
    cat("---\nSignif. codes: ", attr(Signif, "legend"), "\n")
}
