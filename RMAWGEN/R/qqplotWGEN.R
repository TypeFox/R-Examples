
NULL

  
#' 
#' Makes a qqplot and Wilcoxon test between the two columns of \code{val}
#' 
#'  @param val a matrix with two columns containing the two samples to be compared
#'  @param xlab,ylab,main see \code{\link{plot.default}}
#'  @param xlim,ylim see \code{\link{plot.default}}
#'  @param diff logical variable, if \code{TRUE} the function is applied to \code{diff(val)} instead of \code{val}. See \code{\link{diff}}
#' 	@param quantile quantile value on which data samples in \code{val} are considered. Default is 0.
#'
#' @author  Emanuele Cordano, Emanuele Eccel
#' 
#' @export
#' 
#'        
#' @return  Wilcoxon test between the two columns of 'val'



qqplotWGEN <-
function (val,xlab="simulated",ylab="measured",main="title",ylim=c(min(val),max(val)),xlim=c(min(val),max(val)),diff=FALSE,quantile=0) { 
# =cbind(Tn_gen[,1]-SplineAdvTn[,1],Tn_mes[,1]-SplineAdvTn[,1])
	
	if(diff) val <- diff(val) 
	
	v1 <- val[,1]
	v2 <- val[,2]
	if (quantile>0) {
		
		v1 <- v1[v1>quantile(v1,probs=quantile,na.rm=TRUE)]
		v2 <- v2[v2>quantile(v2,probs=quantile,na.rm=TRUE)]
		
	} 
	
	qqplot(v1,v2,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,main=main)
	abline(0,1)
	
	out <- wilcox.test(v1,v2)
	
	return(out)
}

