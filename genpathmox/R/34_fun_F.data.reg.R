
#' @title Defining the matrices needed for the comparison test (regression model)  
#' @details
#' Internal function. \code{F.data.reg} is called by \code{test.partition.reg}.
#' @param x matrix or data.frame with the data.
#' @param modtwo vector indicating the binary partition.
#' @param \dots Further arguments passed on to \code{\link{F.data.reg}}.
#' @return list containing matrices needed for the comparison test
#' @keywords internal
#' @export

F.data.reg	<-	function(x,modtwo,...)
{
	q  = ncol(x)                                          	
	n  = nrow(x)                                           	

	Y0	=	as.matrix(x[,1])
	X0	=	cbind(rep(1,nrow(as.matrix(x[,-1]))),as.matrix(x[,-1]))
	
	x1	=	subset(x,modtwo==1)                        
	x2	=	subset(x,modtwo==2)
	
	l1.resp	=	as.matrix(x1[, 1])
	l1.pred	=	cbind(rep(1,nrow(x1)),as.matrix(x1[,-1]))
	l2.resp	=	as.matrix(x2[, 1])
	l2.pred	=	cbind(rep(1,nrow(x2)),as.matrix(x2[,-1]))

	Y1		=	rbind(l1.resp,l2.resp)
	X1		=	blockdiag(l1.pred,l2.pred)

	var.f=NULL
	var.p=NULL
	
	for(i in 1: ncol(x)-1)
	{
		var.f[i]<- paste(dimnames(x)[[2]][i+1],sep="")
		var.p[i]<- paste(dimnames(x)[[2]][i+1],sep="")
	}
	var.f=c("Intercept",var.f)
	var.p=c("Intercept",var.p)



	list(Y0=Y0,X0=X0, Y1=Y1,X1=X1,var.f=var.f,var.p=var.p)

}

