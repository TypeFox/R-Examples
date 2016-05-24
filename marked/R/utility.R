#' Various utility functions
#' 
#' Several functions have been added to help visualize data including resight.matrix
#' which provides for each cohort the number of releases on the diagonal and the number 
#' resighted on each occasion in the upper-triangular matrix. naive.survival
#' provides a naive survival estimate at each time for each cohort.  The estimate for
#' time i is the number resighted at time i+1 or later divided by the number seen at time i or later.
#' These values should be interpreted cautiously because they are influenced by capture probability in
#' year i but it is useful to identify particularly high or low survival values. Functions Phi.mean
#' and p.mean compute average real parameter values Phi or p across time for a single age or across ages for a
#' single time.
#' @aliases resight.matrix naive.survival mcmc_mode
#' @export resight.matrix
#' @export naive.survival
#' @export mcmc_mode
#' @usage resight.matrix(x)
#' 
#'        naive.survival(x)
#'       
#'        mcmc_mode(x)
#' 
#' @param x processed data list - result from process.data in marked or real estimates from fitted model
#' @return matrix of values with cohort and year labels
#' @author Jeff Laake
#' @keywords utility
resight.matrix=function(x)
{
	zz=strsplit(x$data$ch,"")
	zz=t(sapply(zz,as.numeric))
	x$data$cohort=process.ch(x$data$ch)$first+x$begin-1
	resight.matrix=t(sapply(split(as.data.frame(zz),x$data$cohort),colSums))
	colnames(resight.matrix)=c(0,cumsum(x$time.intervals))+x$begin.time
	return(resight.matrix)
}
naive.survival=function(x)
{
	ll=process.ch(x$data$ch,all=TRUE)
	x$data$cohort=ll$first+x$begin-1
	surv.matrix=t(sapply(split(as.data.frame(ll$First-ll$Lplus),x$data$cohort),colSums))
	naive.S=matrix(NA,nrow=nrow(surv.matrix),ncol=ncol(surv.matrix))
	naive.S=surv.matrix[,2:(ncol(surv.matrix))]/surv.matrix[,1:(ncol(surv.matrix)-1)]
	naive.S[is.infinite(naive.S)]=NA
	naive.S[is.nan(naive.S)]=NA
	rownames(naive.S)=levels(x$data$cohort)
	colnames(naive.S)=x$begin.time+ c(0,cumsum(x$time.intervals[-length(x$time.intervals)]))
	return(naive.S)	
}
mcmc_mode <- function(x){
	dx <- density(x)
	dx$x[dx$y==max(dx$y)]
}



