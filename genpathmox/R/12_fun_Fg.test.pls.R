#' @title F-Global test. 
#' @details
#' Internal function. \code{Fg.test.pls} is called by \code{test.partition.pls}.
#' @param Y0 vector of the concatenate indipendent latent variables of null hypothesis global test.
#' @param X0 matrix of the concatenate predictor latent variables of null hypothesis global test.
#' @param Y1 vector of the concatenate indipendent latent variables of alternative hypothesis global test.
#' @param X1 matrix of the concatenate predictor latent variables of alternative hypothesis global test.
#' @param method string indicating the method: LM or LAD
#' @param \dots Further arguments passed on to \code{\link{Fg.test.pls}}. 
#' @return list containing the statistic and the p-value of global test
#' @keywords internal
#' @export

Fg.test.pls	<-	function(Y0,X0,Y1,X1,method,...)
{
	if (method == "lm")
	{
		reg0	=	lm(Y0~X0-1)                   	
		SSR0	=	sum(reg0$residuals^2)    
		df0	=	(nrow(X0) - ncol(X0))

		reg1	=	lm(Y1~X1-1)                   	
		SSR1	=	sum(reg1$residuals^2)    
		df1	=	(nrow(X1) - ncol(X1))                        
	}
	if (method == "lad")
	{
		reg0	=	rq(Y0~X0-1,method="fn")                   	
		SSR0	=	sum(abs(reg0$residuals))    
		df0	=	(nrow(X0) - ncol(X0))
			
		reg1	=	rq(Y1~X1-1,method="fn")                   	
		SSR1	=	sum(abs(reg1$residuals))    
		df1	=	(nrow(X1) - ncol(X1))                                              
	}
	Fg		=	((SSR0-SSR1)/(df0-df1))/(SSR1/df1)             
	pval.g	=	pf(Fg,(df0-df1),df1,lower.tail=FALSE)     

	list(Fg=Fg ,pvg=pval.g)
}
