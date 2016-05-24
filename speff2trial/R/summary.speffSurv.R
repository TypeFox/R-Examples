summary.speffSurv <- function(object,...){
	ci <- object$beta + tcrossprod(sqrt(object$varbeta)*qnorm(1-(1-object$conf.level)/2),c(-1,1))
	pval <- 2*pnorm(-abs(object$beta)/sqrt(object$varbeta))
	Tab <- cbind(object$beta, sqrt(object$varbeta), ci, pval)
	rownames(Tab) <- c("Prop Haz","Speff")
	colnames(Tab) <- c("Log HR","SE","LowerCI","UpperCI","p")
	out <- list(tab=Tab,method=object$method,formula=object$formula,fixed=object$fixed)
	class(out) <- "summary.speffSurv"
	out
}
