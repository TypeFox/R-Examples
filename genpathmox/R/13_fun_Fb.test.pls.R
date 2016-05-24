#' @title Block test 
#' @details
#' Internal function. \code{Fb.test.pls} is called by \code{test.partition.pls}.
#' @param Y1 vector of the concatenate indipendent latent variables of alternative hypothesis block test.
#' @param X1 matrix of the concatenate predictor latent variables of alternative hypothesis block test.
#' @param info.block list contaning information about the endogenous equations of the pls model.
#' @param method string indicating the method: LM or LAD.
#' @param \dots Further arguments passed on to \code{\link{Fb.test.pls}}. 
#' @return list containing the statistic and the p-value of block test
#' @keywords internal
#' @export

Fb.test.pls	<-	function(Y1,X1,info.block,method,...)
{
	Fb			=	NULL
	pval.b		=	NULL
	k			=	ncol(X1)/2
   	
   	if (method == "lm")
	{
		reg1	=	lm(Y1~X1-1)                   	
		SSR1	=	sum(reg1$residuals^2)    
		df1	=	(nrow(X1) - ncol(X1))                                                
	}
   	if (method == "lad")
	{
		reg1	=	rq(Y1~X1-1,method ="fn")                   	
		SSR1	=	sum(abs(reg1$residuals))    
		df1	=	(nrow(X1) - ncol(X1))	                                                
	}
	
	lb	=	NULL
	for (i in 1:length(info.block)) lb[i] = length(info.block[[i]])				
	p.f = lb[1]
	for (i in 2:length(lb)){p.f[i]=p.f[i-1]+lb[i]}
	p.i = (p.f-lb)+1


	for (j in 1:length(lb))
	{
		A				=	X1    
		A[,p.i[j]:p.f[j]]	=	as.matrix(A[,p.i[j]:p.f[j]]+A[,(k+p.i[j]):(p.f[j]+k)])       
		X1.b				=	A[,-as.vector(seq(k+p.i[j],k+p.f[j],1))]
		
       	if (method == "lm")
		{
			df0.b	=	(nrow(X1.b)-ncol(X1.b)) 
			SSR0.b	=	sum(lm(Y1~X1.b-1)$residuals^2)
	    }
	    if (method == "lad")
		{
			df0.b	=	(nrow(X1.b)-ncol(X1.b)) 
			SSR0.b	=	sum(abs(rq(Y1~X1.b-1,method="fn")$residuals))
	    }        
		Fb[j]		=	((SSR0.b-SSR1)/(df0.b-df1))/(SSR1/df1)             
		pval.b[j]		=	pf(Fb[j],(df0.b-df1),df1,lower.tail=FALSE)    
	}
	
	names(Fb)	=	names(pval.b)	=	names(info.block)
	
	list(Fb=Fb,pvb=pval.b)
}
