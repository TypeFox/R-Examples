
#' @title Defining the candidates to the optimum partition for each of segmentation variables 
#' (regression model) 
#' @details
#' Internal function. \code{all.part.reg} is called by \code{partopt.reg}.
#' @param x matrix or data.frame with the data.
#' @param y matrix or data.frame of the segmentation variables.
#' @param method string indicating the method: LM or LAD
#' @param \dots Further arguments passed on to \code{\link{all.part.reg}}.
#' @return list containing information of the candidates to the optimum partition for each of segmentation variables  
#' @keywords internal
#' @export

all.part.reg <- function(x,y,method,...)
{
	part	=	partition(y)
	                                 
	p.bin	=	list()                                     
	modalidad=	list()                                    
	modtwo	=	list() 
	                                   
	length(p.bin)	=	length(modalidad)	=	length(modtwo)	=	ncol(y) 
	
	Ftest      	= NULL
	df0        	= NULL
	df1        	= NULL
	pvl        	= NULL 
	indg1      	= NULL
	indg2		= NULL
	
	for (j in 1:ncol(y))
	{                              
		if (is.null(part$split[[j]])) next
		
		jpart	=	splitopt.reg(x,splits=part$split[[j]],fact=y[,j],method) 
		
		if (is.null(jpart$pval)) next
		
		p.bin[[j]]			=	jpart$particion                   
		modalidad[[j]]		=	unlist(jpart$mod)
		pvl[j]				=	jpart$pval
		Ftest[j]			=	jpart$Fi
		indg1[j]			=	jpart$ind1
		indg2[j]			=	jpart$ind2
		modtwo[[j]]			=	jpart$modtwo
		df0[j]				=	jpart$df.num
		df1[j]				=	jpart$df.den
	}
	variable				=	names(y)
	
	list(p.bin=p.bin,variable=variable,modalidad=modalidad,Ftest=Ftest,pvl=pvl,indg1=indg1,indg2=indg2,df0=df0,df1=df1,modtwo=modtwo)
}
