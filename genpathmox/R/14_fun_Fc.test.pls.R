#' @title Coefficient test 
#' @details
#' Internal function. \code{Fc.test.pls} is called by \code{test.partition.pls}.
#' @param Y1 vector of the concatenate indipendent latent variables of alternative hypothesis coefficient test.
#' @param X1 matrix of the concatenate predictor latent variables of alternative hypothesis coefficient test.
#' @param path.name vector of path coefficients labels.
#' @param info.block list contaning information about the endogenous equations of the pls model.
#' @param method string indicating the method: LM or LAD.
#' @param \dots Further arguments passed on to \code{\link{Fc.test.pls}}. 
#' @return list containing the statistic and the p-value of coefficient test
#' @keywords internal
#' @export

Fc.test.pls	<-	function(Y1,X1,path.name,info.block,method,...)
{
	Fc			=	NULL
    pval.c		=	NULL
    new.Fc		=	list()
    new.pval.c	=	list()
    k			=	ncol(X1)/2
   	p.name=list()
   	
   	if (method == "lm")
	{
		reg1	=	lm(Y1~X1-1)                   	
		SSR1	=	sum(reg1$residuals^2)    
		df1	=	(nrow(X1) - ncol(X1))                                                
	}
   	if (method == "lad")
	{
		reg1	=	rq(Y1~X1-1)                   	
		SSR1	=	sum(abs(reg1$residuals))  
		df1	=	(nrow(X1) - ncol(X1))                                             
	}	
	
	for (j in 1:k)
	{

		A		= X1	      
		A[,j]	= as.matrix(A[,j]+A[,j+k])         
		X1.c		= A[,-(j+k)]
        
        if (method == "lm")
        {
        	df0.c	= (nrow(X1.c)-ncol(X1.c))                             
		SSR0.c	= sum(lm(Y1~X1.c-1)$residuals^2)
		}
	    if (method == "lad")
	    {         
        	df0.c	= (nrow(X1.c)-ncol(X1.c))                             
		SSR0.c	= sum(abs(rq(Y1~X1.c-1,method="fn")$residuals))
	    }     
		
		Fc[j]	 =	((SSR0.c-SSR1)/(df0.c-df1))/(SSR1/df1)             
		pval.c[j]	 =	pf(Fc[j],(df0.c-df1),df1,lower.tail=FALSE)     
	}
		
	
	for (j in 1:length(info.block))
	{
		p.name[[length(p.name)+1]]	=	paste(info.block[[j]],"->",names(info.block)[j])
	}

	names(Fc)=names(pval.c)=unlist(p.name)
	
	list(Fc=Fc,pvc=pval.c)
}
