#' Plot method for lb objects
#'
#' Produce a plot of an LB fit. The default is a complete coefficient path.
#' 
#' The default plot uses the fraction of L1 norm as the x. 
#' For multinomial case, the sum of absolute values of different class's 
#' coefficients are caculated to represent each variable.
#' The intercept term is not ploted
#' 
#' @param x lb object
#' @param xtype  The x-axis type. "t" or "norm". Default is "t".
#' @param omit.zeros When the number of variables  is much greater than
#' the number of observations, many coefficients will never be nonzero;
#' this logical (default \code{TRUE}) avoids plotting these zero coefficents
#' @param eps Definition of zero above, default is \code{1e-10}
#' @param \dots Additonal arguments for generic plot. Can be used to set xlims, change colors, line widths, etc
#' @return NULL
#' @author Feng Ruan, Jiechao Xiong and Yuan Yao
#' @keywords methods
#'

plot.lb <- function(x,xtype = c("t","norm"),omit.zeros=TRUE,eps= 1e-10,...){
  xtype <- match.arg(xtype)
  object <- x
	if (object$family=="multinomial")
	  coef <- sapply(1:dim(object$path)[3],function(x) colSums(abs(object$path[,,x])))
	else
	  coef <- object$path
	if(omit.zeros){
		c1 <- drop(abs(coef)%*%rep(1, ncol(coef)))
		nonzeros <- c1 > eps
		cnums <- seq(nonzeros)[nonzeros]
		coef <- coef[nonzeros,,drop=FALSE]
	}else {cnums <- seq(nrow(coef))}
	stepid <- seq(ncol(coef))
	if (xtype=="t")
	  s <- object$t
	else if(xtype=="norm"){
	  s <- apply(abs(coef),2,sum)
	  s <- s/max(s)
	}
	
	if (object$kappa==Inf)
	  matplot(s,t(coef),xlab="Solution-Path",ylab="Coefficients",ylim = range(coef),type="s",...)
	else
	  matplot(s,t(coef),xlab="Solution-Path",ylab="Coefficients",ylim = range(coef),type="l",...)
	#title(paste(object$type[1],object$type[2],sep="-"), line=3)
	abline(h=0, lty=2.5)
	axis(3, at=s, labels=paste(stepid), cex=.5)
	axis(4,at=coef[,ncol(coef)],labels=paste(cnums),cex=0.8)
}