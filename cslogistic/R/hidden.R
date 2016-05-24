
"print.BayesCslogistic" <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", deparse(x$call), "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Posterior Inference of Coefficients:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No coefficients\n")
    cat("\nAcceptation Rate for the Metropolis Algorihtm = ",x$arate,"\n")    
    cat("\n\n")
    invisible(x)
}



"print.MleCslogistic" <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\n",x$modelname,"\n\nCall:\n", deparse(x$call), "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(x$coefficients, digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No coefficients\n")
    cat("\nlog-likelihood",x$loglike,"\n")    
    cat("\n\n")
    invisible(x)
}


"summary.BayesCslogistic" <- function(object, ...) 
{
    dimen <- object$dimen
    coef.p <- object$coefficients
    
    coef.sd <- rep(0,dimen)
    coef.se <- rep(0,dimen)
    coef.l <- rep(0,dimen)
    coef.u <- rep(0,dimen)
    coef.m <- rep(0,dimen)
    
    names(coef.sd) <- object$pnames
    names(coef.l) <- object$pnames
    names(coef.u) <- object$pnames
    
    alpha <- 0.05
    
    for(i in 1:dimen)
	{
        alow <- rep(0,2)
        aupp <- rep(0,2)
        coef.sd[i] <- sqrt(var(object$mat[,i]))
        coef.m[i] <- median(object$mat[,i])
        vec <- object$mat[,i]
        n <- length(vec)
        
        a<-.Fortran("hpd",n=as.integer(n),alpha=as.double(alpha),x=as.double(vec),
                     alow=as.double(alow),aupp=as.double(aupp),PACKAGE="cslogistic")
        coef.l[i] <- a$alow[1]            
        coef.u[i] <- a$aupp[1]            
    }

    coef.se <- coef.sd/sqrt(n)

    coef.table <- cbind(coef.p, coef.m, coef.sd, coef.se , coef.l , coef.u)
    dimnames(coef.table) <- list(names(coef.p), c("Mean", "Median", "Std. Dev.", "Naive Std.Error",
                "95%HPD-Lower","95%HPD-Upper"))

    ans <- c(object[c("call", "modelname","arate")])
    ans$coefficients <- coef.table
    class(ans) <- "BayesCslogistic"
    return(ans)
}



"summary.MleCslogistic" <- function(object, ...) 
{
    p <- object$dimx
    coef.p <- object$coefficients
    s.err <- object$se
    tvalue <- object$tvalue
    pvalue <- 2 * pnorm(-abs(tvalue))
    or <- exp(coef.p)
    ic.l <- exp(coef.p-qnorm(0.975)*s.err)
    ic.u <- exp(coef.p+qnorm(0.975)*s.err)
    coef.table <- cbind(coef.p, s.err, or,ic.l,ic.u, pvalue)
    dimnames(coef.table) <- list(names(coef.p), c("Estimate", "Std. Error", 
                " OR ", "Lower","Upper" ,"Pr(>|z|)"))

    ans <- c(object[c("call", "loglike","modelname")])
    ans$coefficients <- coef.table
    class(ans) <- "MleCslogistic"
    return(ans)
}



"plot.BayesCslogistic" <- function(x, ...) 
{

	if(is(x, "BayesCslogistic"))
	{
	
       dimen <- x$dimen
       par(mfrow=c(3,2))
       for(i in 1:dimen)
	   {
		   nx <- round(i/4,0)
           mx <- nx*4
		   if(i==mx)
		   {
              par(mfrow=c(3,2))
           }
           title1 <- paste("Trace of",x$pnames[i],sep=" ")
           title2 <- paste("Density of",x$pnames[i],sep=" ")
           plot(x$mat[,i],type='l',main=title1,xlab="Iterations",ylab=" ")
           plot(density(x$mat[,i]),type='l',main=title2,xlab="Values", ylab="Density")
         }
	}	
}
