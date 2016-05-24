#' @title Calculating the comparison tests for the regression model
#' @details
#' Internal function. \code{test.partition.reg} is called by \code{reg.pathmox}.
#' @param x matrix or data.frame with the data.
#' @param modtwo vector indicating the binary partition
#' @param signific string indicating a stop condition.
#' @param method string indicating the method: LM or LAD.
#' @param \dots Further arguments passed on to \code{\link{test.particion.reg}}.
#' @return list containing matrices needed for the comparison test
#' @keywords internal
#' @export

test.particion.reg<- function(x,modtwo,signif,method,...) 
{
	d.info	=	F.data.reg(x,modtwo,df)
	
	FG		=	Fg.test.reg(d.info$Y0,d.info$X0,d.info$Y1,d.info$X1,method)

	if(FG$pvg>signif)
	{
		list(Fg=FG$Fg ,pvg=FG$pvg,Fc=0,pvc=0)
	}
	else
	{
		FC	=	Fc.test.reg(d.info$Y1,d.info$X1,d.info$var.f,d.info$var.p,method)
		
		list(Fg=FG$Fg ,pvg=FG$pvg,Fc=FC$Fc,pvc=FC$pvc) 	
	}
}
