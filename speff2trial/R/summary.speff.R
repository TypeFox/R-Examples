summary.speff <- function(object,...){
	ci <- object$coef[,3] + tcrossprod(sqrt(object$varbeta)*qnorm(1-(1-object$conf.level)/2),c(-1,1))
	pval <- 2*pnorm(-abs(object$coef[,3])/sqrt(object$varbeta))
	Tab <- cbind(object$coef[,3], sqrt(object$varbeta), ci, pval)
	rownames(Tab) <- c("Naive","Speff")
	colnames(Tab) <- c(ifelse(object$endpoint=="quantitative","Mean Diff","Log OR"),"SE","LB","UB","p")
	out <- list(tab=Tab,method=object$method,rsq=object$rsq,formula=object$formula,predicted=object$predicted)
	class(out) <- "summary.speff"
	out
}
