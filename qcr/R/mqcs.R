#-----------------------------------------------------------------------------#
#                                                                             #
#            QUALITY CONTROL AND RELIABILITY IN R                             #
#                                                                             #
#  An R package for statistical in-line quality control.                      #
#                                                                             #
#  Written by: Miguel A. Flores Sánchez                                       #
#              Student Master of Statistical Techniques                       #
#              University of The Coruña, SPAIN                                #
#              mflores@outlook.com                                            #
#                                                                             #
#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
# Main function to create a 'mqcs' object
#-----------------------------------------------------------------------------#
##' Create an object of class 'mqcs' to perform statistical quality control.
##' This function is used to compute statistics required for to plot Multivariate Control Charts
##' 
##' @aliases mqcs summary.mqcs print.mqcs
##' 
##' @param x  Object mqcd (Multivariante Quality Control Data)
##' @param method Is the method employed to compute the covatiance matrix
##' in individual observation case. Two methods are used "sw" 
##' for compute according to (Sullivan,Woodall 1996a) and "hm" 
##' by (Holmes,Mergen 1993)
##' @param ... arguments passed to or from methods.
##' @export

mqcs <- function(x, method = "sw", ...)
  #.........................................................................  
  {
  
  p <- ncol(x) # quality characteristics
  m <- nrow(x) # number of samples or observations
  n <- dim(x)[3] # observations or sample size 
  
  x.jk <- matrix(0,m,p)  
  x.jk <- apply(x,1:2,mean)
  
  Xmv <- colMeans(x.jk)
  S <- covariance(x,method = method)
  
  result <- list (mqcd = x, mean = Xmv, S = S, mean.jk = x.jk) 
  
  oldClass(result) <- c("mqcs")
  
  return(result)
} # mqcs
#.........................................................................

##' @export
print.mqcs <- function(x, ...) str(x,1)
#.........................................................................
##' @export
summary.mqcs <- function(object, ...)
  #.........................................................................
{
  type <- object$type
  t2 <- object$statistics
  cat("\nSummary of group statistics:\n")
  print(summary(t2))
  cat("\nNumber of quality characteristics: ", ncol(object$mqcd))
  cat("\nNumber of samples or observations: ", nrow(object$mqcd))
  cat("\nNumber of observations or sample size: ", dim(object$mqcd)[3])
  
  center <- object$mean
  cat("\n\nMean Vector: \n", center)
  cat("\nCovariance Matrix:\n")
  S <-object$S
  print(S)
  
  limits <- object$limits
  if (!is.null(limits)) 
  { cat("\nControl limits:", "\n") 
    print(limits)
  }
  

  if (length(object$violations)== 0){
    cat("\nNumber beyond limits: 0", "\n") 
  } 
  else {cat("\nBeyond limits of control:", "\n")
        print(object$statistics[object$violations])
  }
  
  invisible()
  #.........................................................................
} # summary.mqcs



## ' Sample covariance
## ' 
## ' It allows to compute the sample covariance in presence of rational subgroups
## ' or for individuals according to (Sullivan,Woodall 1996) and (Holmes,Mergen
## ' 1993)
## ' 
## ' 
## ' @param x matrix or array of the quality characteristics.
## ' @param stat is the statistics
## ' @param method is the method used in individual observation case.
## ' @param \dots other parameters
## ' @note In individuals observation case (n = 1) use for default the
## ' (Sullivan,Woodall 1996) proposal
## ' @author Edgar Santos-Fernandez
## ' @references Holmes, D.S., Mergen, A.E.: Improving the performance of
## ' T-square control chart. Quality Engineering 5(4), 619-625 (1993)
## ' 
## ' Sullivan, J.H., Woodall, W.H.: A Comparison of Multivariate Quality Control
## ' Charts for Individual Observations. Journal of Quality Technology 28(4)
## ' (1996)
## ' @keywords ~kwd1 ~kwd2
## ' @examples
## ' 
## ' # individual case 
## ' data(dowel1)
## ' covariance(dowel1,method="sw")
## ' covariance(dowel1,method="hm")
 
covariance <- function(x, stat, method, ...){
    p <- ncol(x) # quality characteristics
    m <- nrow(x) # sample
    if (class(x) == "matrix" || class(x) == "data.frame") (x <- array(data.matrix(x),c(m,p,1)))
    n <- dim(x)[3] # observations or sample size
    
    s.jk <- matrix(0,m,p ^ 2) # matrix of the covariances of observations
    SS <- matrix(0,m,1) # matrix of /S/ statistic 
    
    if(n > 1){
      arrays <- expand.grid(1:p,1:p)
      
      for (i in 1 : m){
        for(j in 1 : p ^ 2){
          s.jk[i,j] <- cov(x[i,arrays[j,1],],x[i,arrays[j,2],])
        }
      } 
      
      S <- matrix(colMeans(s.jk),p,p)
      
      for (ii in 1 : m){
        SS[ii] <- det(matrix(s.jk[ii,],p,p))
      }
      
      if(missing(stat)) (return(S))
      else (return(SS))
      
    }    
    
    if(n == 1){
      if(missing(method))(method="sw")
      
      if(method == "sw"){
        B <- matrix(0,p,p)
        w <- sweep(x,2,(apply(x,2,mean))) #compute de value minus the mean
        for(i in 1:m){
          B <- B + w[i,,] %*% t(w[i,,])
        }
        S <- s1 <- B/(m - 1)
      }
      
      if(method == "hm"){
        V <- matrix(0,m-1,p)
        for(i in 1:m-1){
          V[i,] <- x[i+1,,] - x[i,,]
        }
        S <- s2 <- .5 * t(V) %*% V / (m - 1)
      }
      
      
      return(S)
    }
    
    
  }