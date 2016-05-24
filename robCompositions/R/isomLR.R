#' Isometric log-ratio transformation
#' 
#' An isometric log-ratio transformation (ilr1) and it's inverse transformation with a special choice of the balances
#' 
#' The isomLR transformation moves D-part compositional data from the simplex
#' into a (D-1)-dimensional real space isometrically. From our choice of the
#' balances (ilr1), all the relative information of the part \eqn{x_1} from the
#' remaining parts is separated. It is useful for estimating missing values in
#' \eqn{x_1} by regression of the remaining variables.
#' 
#' @aliases isomLR isomLRinv
#' @param x object of class data.frame or matrix. Positive values only for isomLR().
#' @param fast if TRUE, it is approx. 10 times faster but numerical problems in case of 
#' high-dimensional data may occur.
#' @return The transformed data.
#' @author Karel Hron, Matthias Templ
#' @references Egozcue J.J., V. Pawlowsky-Glahn, G. Mateu-Figueras and C.
#' Barcel'o-Vidal (2003) Isometric logratio transformations for compositional
#' data analysis. \emph{Mathematical Geology}, \bold{35}(3) 279-300. \
#' 
#' Hron, K. and Templ, M. and Filzmoser, P. (2010) Imputation of missing values
#' for compositional data using classical and robust methods
#' \emph{Computational Statistics and Data Analysis}, vol 54 (12), pages
#' 3095-3107.
#' @keywords math
#' @export
#' @examples
#' 
#' require(MASS)
#' Sigma <- matrix(c(5.05,4.95,4.95,5.05), ncol=2, byrow=TRUE)
#' z <- isomLRinv(mvrnorm(100, mu=c(0,2), Sigma=Sigma))
#' 
#' data(expenditures)
#' isomLR(expenditures)
#' 
#' x <- exp(mvrnorm(2000, mu=rep(1,10), diag(10)))
#' system.time(isomLR(x))
#' system.time(isomLR(x, fast=TRUE))
#' 
#' 
isomLR <- function(x, fast=FALSE){
	x.ilr=matrix(NA,nrow=nrow(x),ncol=ncol(x)-1)
	D=ncol(x)
        if(fast){
	  for (i in 1:ncol(x.ilr)){
	     x.ilr[,i]=sqrt((D-i)/(D-i+1))*log(((apply(as.matrix(x[,(i+1):D,drop=FALSE]),1,prod))^(1/(D-i)))/(x[,i]))
	  }
	} else {
  	  for (i in 1:ncol(x.ilr)){
#		x.ilr[,i]=sqrt((D-i)/(D-i+1))*log(((apply(as.matrix(x[,(i+1):D,drop=FALSE]),1,prod))^(1/(D-i)))/(x[,i]))
		x.ilr[,i]=sqrt((D-i)/(D-i+1))*log(apply(as.matrix(x[,(i+1):D]), 1, gm)/(x[,i]))	
	  } 
	}
	return(-x.ilr)
}
#' @rdname isomLR
#' @export
isomLRinv <- function(x){
  x <- -x
  y=matrix(0,nrow=nrow(x),ncol=ncol(x)+1)
  D=ncol(x)+1
  y[,1]=-sqrt((D-1)/D)*x[,1]
  for (i in 2:ncol(y)){
    for (j in 1:(i-1)){
      y[,i]=y[,i]+x[,j]/sqrt((D-j+1)*(D-j))
    }
  }
  for (i in 2:(ncol(y)-1)){
    y[,i]=y[,i]-sqrt((D-i)/(D-i+1))*x[,i]
  }
  yexp=exp(y)
  x.back=yexp/apply(yexp,1,sum) # * rowSums(derOriginaldaten)
  return(x.back)
  #return(yexp)
}
