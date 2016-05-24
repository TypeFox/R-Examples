#' Wrapper for Anderson-Darling tests
#' 
#' A set of Anderson-Darling tests (Anderson and Darling, 1952) are applied as
#' proposed by Aitchison (Aichison, 1986).
#' 
#' First, the data is transformed using the \sQuote{ilr}-transformation.  After
#' applying this transformation
#' 
#' - all (D-1)-dimensional marginal, univariate distributions are tested using
#' the univariate Anderson-Darling test for normality.
#' 
#' - all 0.5 (D-1)(D-2)-dimensional bivariate angle distributions are tested
#' using the Anderson-Darling angle test for normality.
#' 
#' - the (D-1)-dimensional radius distribution is tested using the
#' Anderson-Darling radius test for normality.
#' 
#' A print and a summary method are implemented. The latter one provides a similar output is proposed by (Pawlowsky-Glahn, et al. (2008). In addition
#' to that, p-values are provided.
#' 
#' @aliases adtestWrapper print.adtestWrapper summary.adtestWrapper
#' @param x compositional data of class data.frame or matrix
#' @param alpha significance level
#' @param R Number of Monte Carlo simulations in order to provide p-values.
#' @param robustEst logical
#' @param object an object of class adtestWrapper for the summary method
#' @param ... additional parameters for print and summary passed through
#' @return \item{res }{ a list including each test result } \item{check }{
#' information about the rejection of the null hypothesis} \item{alpha}{ the
#' underlying significance level } \item{info}{ further information which is
#' used by the print and summary method. } \item{est}{ \dQuote{standard} for
#' standard estimation and \dQuote{robust} for robust estimation }
#' @author Matthias Templ and Karel Hron
#' @seealso \code{\link{adtest}}, \code{\link{isomLR}}
#' @references Anderson, T.W. and Darling, D.A. (1952) \emph{Asymptotic theory
#' of certain goodness-of-fit criteria based on stochastic processes} Annals of
#' Mathematical Statistics, \bold{23} 193-212.
#' 
#' Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional Data}
#' Monographs on Statistics and Applied Probability. Chapman \& Hall Ltd.,
#' London (UK). 416p.
#' @keywords htest
#' @export
#' @examples
#' 
#' data(machineOperators)
#' a <- adtestWrapper(machineOperators, R=50) # choose higher value of R
#' a
#' summary(a)
#' 
adtestWrapper=function(x,alpha=0.05,R=1000, robustEst=FALSE){
  if(robustEst == TRUE ) robust <- "robust" else robust <- "standard"
  z=isomLR(x)
  n=ncol(z)
  if(ncol(z)==1){
    res<-info<-list()
    res[[1]]=adtest(z,R,locscatt=robust)
    info[[1]]=paste(1)
    check<- logical(1)
  }
  if(ncol(z)==2){
    res<-info<-list()
    res[[1]]=adtest(z[,1],R,locscatt=robust)
    res[[2]]=adtest(z[,2],R,locscatt=robust) 
    res[[3]]=adtest(z,R,locscatt=robust)
    info[[1]]=paste(1)
    info[[2]]=paste(2)
    info[[3]]=paste(3)
    check <- logical(3)
  }
  if(ncol(z)>2){
    res<-info<-list()
    for(i in 1:ncol(z)){
      res[[i]]=adtest(z[,i],R,locscatt=robust)
      info[[i]]=paste(i)
    }
    index=1
    for(i in 1:(ncol(z)-1)){
      for(j in (i+1):ncol(z)){     
        res[[n+index]]=adtest(z[,c(i,j)],R,locscatt=robust)
        info[[n+index]]=paste(i,j,collapse=":")
        index=index+1
      }
    }
    res[[n+index]]=adtest(z,R,locscatt=robust)
    info[[n+index]]=paste("all")
    check <- logical(n+index)  
  }
  
  for(i in 1:length(check)){  
    check[i] <- ifelse(res[[i]]$p.value > alpha, TRUE, FALSE)
  }
  output=list(res=res, check=check, alpha=alpha, info=info, est=robust)
  class(output)="adtestWrapper"
  invisible(output)
}

#' @rdname adtestWrapper
#' @export
print.adtestWrapper <- function(x, ...){
  if(all(x$check)){
    print(paste("The data follow the normal distribution on the simplex (alpha =",x$alpha,")",sep=""))
  } else { 
    print(paste("The data do not follow the normal distribution on the simplex (alpha =",x$alpha,")",sep=""))
    #print(x$check)
  }
}

#' @rdname adtestWrapper
#' @export
#' @method summary adtestWrapper
summary.adtestWrapper=function(object, ...){
  d=data.frame(ilrVars=unlist(object$info),
               testName=unlist(lapply(object$res,function(x) x$method)),
               testStat=unlist(lapply(object$res,function(x) x$statistic)),
               pvalue=unlist(lapply(object$res,function(x) x$p.value)),
               #alpha=object$alpha,
               check=object$check)
  #d=lapply(res,function(x) x$info)
  string <- paste("Anderson-Darling test results ( alpha =", object$alpha, "):")
  string2 <- paste("--> p-values and tests are obtained from", object$est, "estimates.")
  cat("\n  -----------------------------------------------")	
  cat("\n ", string)
  cat("\n  ----------------\n")
  print(d)
  cat("\n  -----------------------------------------------\n")
  cat("\n ", string2)
  cat("\n")
  invisible(d)
}
