#' Print method for \code{shrSplines} objects
#' 
#' Print a summary of a fitted illness-death model using the penalized
#' likelihood approach.
#' 
#' 
#' @param x a \code{shr} object, i.e., the result of a call to the
#' \code{\link{shr}} function with \code{hazard}="Splines".
#' @param conf.int confiance level.
#' @param digits number of digits to print.
#' @param pvalDigits number of digits to print for p-values.
#' @param eps convergence criterion used for p-values.
#' @param \dots other unusued arguments.
#' @author R: Celia Touraine <Celia.Touraine@@isped.u-bordeaux2.fr> Fortran:
#' Pierre Joly <Pierre.Joly@@isped.u-bordeaux2.fr>
#' @seealso \code{\link{summary.shr}}, \code{\link{plot.shr}}
#' @keywords methods
#' @examples
#' 
#' \dontrun{
#' # a penalized survival model 
#' library(prodlim)
#' data(testdata)
#' fit.su <- shr(Hist(time=list(l,r),id)~cov,data=testdata,method="Splines") 
#' print(fit.su) 
#' }
#' @S3method print shr
print.shr <- function(x,conf.int=.95,digits=4,pvalDigits=4,eps=0.0001,...){
    if (!inherits(x,"shr")) stop("Object must be of class 'shr'")
    cl <- x$call
    cat("Call:\n")
    dput(cl)
    cat("\n")
    if (x$converged == 1){
        if(x$NC >0){
            wald <- (x$coef/x$se)**2
            z <- abs(qnorm((1 + conf.int)/2))
            tmp <- data.frame("coef"=format(round(x$coef,digits)),
                              "SE coef"=format(round(x$se,digits)),
                              "HR"=format(round(x$HR,digits)),
                              "CI"=paste("[",format(round(exp(x$coef - z * x$se),2)),";",format(round(exp(x$coef + z * x$se),2)),"]",sep=""),
                              "Wald"=format(wald,digits),
                              "P-value"=format.pval(1 - pchisq(wald, 1),digits=pvalDigits,eps=eps))
            rownames(tmp) <- names(x$coef)
        }
        tmp1 <- matrix(x$loglik,nrow=1)
        dimnames(tmp1) <- list("Log likelihood", c("Without covariates", "With covariates"))
        if (x$method=="Splines"){
            cat("Survival model using a penalized likelihood approach \nwith M-splines approximation of the baseline hazard function.\n")
            cat("number of nodes: ", x$nknots,"\n")
            if(x$CV){
                cat("Smoothing parameters estimated by Cross validation: ",x$kappa,"\n")
                cat("Cross validation criterion:",x$CVcrit,"\n")
                cat("DoF: ", formatC(-x$DoF, format="f",digits=2),"\n")
            }else{
                cat("Smoothing parameters: ",x$kappa,"\n")
            }
        } else{
            cat("Parameters of the Weibull distribution: 'S(t) = exp(-(b*t)^a)'\n")
            cat("                   a = ",x$modelPar[1]," "," b = ",x$modelPar[2],"\n")
        }
        print(x$modelResponse)
        if(length(x$na.action))cat("number of deleted observations due to missing: ",length(x$na.action),"\n")
        cat("\n")
        if(x$NC >0){
            cat("\n")
            print(tmp,row.names=T)
        }
        cat("\n")
        prmatrix(tmp1)
        cat("\n")
        cat("----\nModel converged.\n")	
        cat("number of iterations: ", x$niter,"\n")
        cat("convergence criteria: parameters=", signif(x$cv[1],2), "\n")
        cat("                    : likelihood=", signif(x$cv[2],2), "\n") 
        cat("                    : second derivatives=", signif(x$cv[3],2), "\n")
    }else{
        switch(as.character(x$converged),
               "2"={ warning("Maximum number of iterations reached.",call.=FALSE)},
               "3"={ warning("Model did not converge.",call.=FALSE)})
    }

}
