#' @title Defining the optimum partition given a set of segmentation variables (regression model)
#' @details
#' Internal function. \code{partopt.reg} is called by \code{reg.pathmox}.
#' @param x matrix or data.frame with the data.
#' @param y matrix or data.frame of the segmentation variables.
#' @param method string indicating the method: LM or LAD
#' @param \dots Further arguments passed on to \code{\link{all.part.reg}}.
#' @return list containing information of the optimum partition given a set of segmentation variables
#' @keywords internal
#' @export

partopt.reg	<-	function(x,y,method,...)
{
	a.p	=	all.part.reg(x,y,method)
	
	if(any(!is.null(a.p$pvl))){
		
	fun.min	=	f.min(a.p$pvl)
	
		p.bin.opt		=	a.p$p.bin[[fun.min$p.min]]                  
		modalidad.opt	=	a.p$modalidad[[fun.min$p.min]] 
		variable.opt	=	a.p$variable[[fun.min$p.min]]  
		pvl.opt			=	a.p$pvl[fun.min$p.min]
		Ftest.opt		=	a.p$Ftest[fun.min$p.min]       
		indg1.opt		=	a.p$indg1[fun.min$p.min]      
		indg2.opt		=	a.p$indg2[fun.min$p.min]       
		modtwo.opt		=	a.p$modtwo[[fun.min$p.min]]		
		df0.opt			=	a.p$df0[fun.min$p.min]
		df1.opt			=	a.p$df1[fun.min$p.min]
		
		variable		=	a.p$variable[fun.min$all.v]   
		Ftest			=	a.p$Ftest[!is.na(a.p$Ftest)]
		pvl				=	a.p$pvl[!is.na(a.p$pvl)]
		df0				=	a.p$df0[!is.na(a.p$df0)]
		df1				=	a.p$df1[!is.na(a.p$df1)]
		indg1			=	a.p$indg1[!is.na(a.p$indg1)]
		indg2			=	a.p$indg2[!is.na(a.p$indg2)]
		
		i=1
		modalidad=NULL
		while(i<length(unlist(a.p$modalidad)))
		{
			modalidad=rbind(modalidad,unlist(a.p$modalidad)[c(i,i+1)])
			i=i+2
		}
		colnames(modalidad)	=	c("mod.g1","mod.g2")
		
		candidates	=	data.frame(variable,Ftest,pvl,df0,df1,indg1,indg2,modalidad)
		candidates	=	candidates[order(candidates[,2],decreasing=T),]

	list(candidates=candidates,
		p.bin.opt=p.bin.opt,
		variable.opt=variable.opt,
		modalidad.opt=modalidad.opt,
		Ftest.opt=Ftest.opt,
		pvl.opt=pvl.opt,
		indg1.opt=indg1.opt,
		indg2.opt=indg2.opt,
		df0.opt=df0.opt,
		df1.opt=df1.opt,
		modtwo.opt=modtwo.opt)
	}
	else
	{
		list(candidates=NULL,
		p.bin.opt=NULL,
		variable.opt=NULL,
		modalidad.opt=NULL,
		Ftest.opt=NULL,
		pvl.opt=NULL,
		indg1.opt=NULL,
		indg2.opt=NULL,
		df0.opt=NULL,
		df1.opt=NULL,
		modtwo.opt=NULL)	
	}
	
}
