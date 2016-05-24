#' Print method for \code{idm} objects
#' 
#' Print a summary of a fitted illness-death model 
#' 
#' 
#' @param x Class \code{idm} object, i.e. the result of a call to the
#' \code{\link{idm}} function with \code{intensities}="Weib".
#' @param conf.int Confidence level.
#' @param digits Number of digits to print.
#' @param pvalDigits Number of digits to print for p-values.
#' @param eps Passed to \code{format.pval}.
#' @param \dots Not used.
#' @author Celia Touraine <Celia.Touraine@@isped.u-bordeaux2.fr>, Thomas A. Gerds <tag@@biostat.ku.dk> 
#' @seealso \code{\link{summary.idm}}, \code{\link{plot.idm}}
#' @keywords methods
#' @examples
#' 
#' \dontrun{
#' data(Paq1000)
#' library(prodlim)
#' fit.splines <-  idm(formula02=Hist(time=t,event=death,entry=t0)~certif,
#' 		formula01=Hist(time=list(l,r),event=dementia)~certif,
#'                 formula12=~1,
#'                 method="Splines",
#' 		data=Paq1000)
#' print(fit.splines)
#' 
#' }
#' @S3method print idm
print.idm <- function(x,conf.int=.95,digits=4,pvalDigits=4,eps=0.0001,...){
    # {{{  call
    cl <- x$call
    cat("Call:\n")
    dput(cl)
    cat("\n")
    # }}}
    # {{{ number of subjects etc
    cat("Illness-death regression model using",
        ifelse(x$method=="Splines"," M-spline approximations","Weibull parametrization"),
        "\nto estimate the baseline transition intensities.\n")
    cat("\n")
    cat("number of subjects: ", x$N,"\n")
    cat("number of events '0 -> 1': ", x$events1,"\n")
    cat("number of events '0 -> 2' or '0 -> 1 -> 2': ", x$events2,"\n")
    cat("number of covariates: ", x$NC,"\n")
    # }}}
    # {{{ convergence
    # FIXME: what is the difference between maximum number of iterations reached and
    #        model did not converge?
    if( ((x$converged[1]!=1)&(x$converged[2]!=1)) ){
        warning("The model did not converge.","\n")
        switch(as.character(x$converged[1]),
               "2"={ warning("Maximum number of iterations reached.",call.=FALSE)},
               "3"={ warning("Fisher information matrix non-positive definite.",call.=FALSE)})	
    }else{
        if((sum(x$NC)>0)&(x$converged[1]==1)&(x$converged[2]!=1) ){
            warning("The model did converge without covariates but did not converge with covariates","\n")
            switch(as.character(x$converged[2]),
                   "2"={ warning("Maximum number of iterations reached.",call.=FALSE)},
                   "3"={ warning("Model did not converge.",call.=FALSE)})
            cat("Log-likelihood without covariates: ",x$loglik[1], "\n")
        }else{
            if((sum(x$NC)>0)&(x$converged[1]!=1)&(x$converged[2]==1) ){
                cat("The model did converge with covariates but did not converge without covariates","\n")
                switch(as.character(x$converged[1]),
                       "2"={ warning("Maximum number of iterations reached.",call.=FALSE)},
                       "3"={ warning("Model did not converge.",call.=FALSE)})}}}
    if ((x$converged[1]== 1)&(x$converged[2] %in% c(0,1))  ){
        cat("----\nModel converged.\n")
        cat("number of iterations: ", x$niter,"\n")
        cat("convergence criteria: parameters=", signif(x$cv[1],2), "\n")
        cat("                    : likelihood=", signif(x$cv[2],2), "\n") 
        cat("                    : second derivatives=", signif(x$cv[3],2), "\n")
        loglik <- matrix(x$loglik,nrow=1)
        dimnames(loglik) <- list(paste(ifelse(x$method=="Splines","Penalized",""),
                                       "log-likelihood"),
                                 c("Without covariates", "With covariates"))
        cat("\n")
        prmatrix(loglik)
        cat("\n")
    }
    # }}}
    # {{{ Spline: baseline parameters
    if (x$method=="Splines"){
        splinepars <- data.frame("transition01"=c(x$nknots01,x$kappa[1]),
                                 "transition02"=c(x$nknots02,x$kappa[2]),
                                 "transition12"=c(x$nknots12,x$kappa[3]))
        rownames(splinepars) <- c("knots","kappa")
        cat("\n")
        if(x$CV){
            cat("Smoothing parameters estimated by cross validation:\n")
            print(splinepars,row.names=T)
            cat("Cross validation criterion:",x$CVcrit,"\n")
            cat("DoF: ", formatC(-x$DoF, format="f",digits=2),"\n")
        }else{
            cat("Smoothing parameters:\n")
            print(splinepars,row.names=T)
        }
        # }}}
        # {{{ Weibull: baseline parameters
    }else{
        cat("Parameters of the Weibull distributions: 'S(t) = exp(-(b*t)^a)'\n")
        wpars <- matrix(x$modelPar,nrow=2)
        dimnames(wpars) <- list(c("shape (a)","scale (b)"),
                                c("transition 0 -> 1",
                                  "transition 0 -> 2",
                                  "transition 1 -> 2"))
        prmatrix(wpars)
    }
    # }}}
    # {{{  Regression coefficients
    if(sum(x$NC)>0){
        wald <- (x$coef/x$se)**2
        z <- abs(qnorm((1 + conf.int)/2))
        coefmat <- data.frame("coef"=format(round(x$coef,digits)),
                              "SE coef"=format(round(x$se,digits)),
                              "exp(coef)"=format(round(x$HR,digits)),
                              "CI"=paste("[",
                                  format(round(exp(x$coef - z * x$se),2)),
                                  ";",
                                  format(round(exp(x$coef + z * x$se),2)),
                                  "]",
                                  sep=""),
                              ## "Wald"=format(wald,digits),
                              "p-value"=format.pval(1 - pchisq(wald, 1),digits=pvalDigits,eps=eps),
                              check.names=FALSE)
        coefmat <- cbind(Factor=names(x$coef),coefmat)
        coeflist <- split(coefmat,rep(c("transition 0 -> 1","transition 0 -> 2","transition 1 -> 2"),x$NC))
        cat("\n\nRegression coefficients:\n\n")
        print(coeflist)
    }
}

