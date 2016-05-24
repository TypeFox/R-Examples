summary.endogMNP <- function(object, CI=c(2.5, 97.5), discard=1, ...){

  p <- object$n.dim
  param <- object$param
n.cov <- ncol(param) - p*(p+1)/2
#	n.cov <- ncol(param) - object$n.dim 
  n.draws <- nrow(param)
	param.table <- cbind(apply(param[discard:n.draws,], 2, mean), apply(param[discard:n.draws,], 2, sd),
                       apply(param[discard:n.draws,], 2, quantile, min(CI)/100),
                       apply(param[discard:n.draws,], 2, quantile, max(CI)/100)) 
  colnames(param.table) <- c("mean", "std.dev.", paste(min(CI), "%", sep=""),
                             paste(max(CI), "%", sep=""))
  rownames(param.table) <- colnames(param)
  
  ans <- list(call = object$call, selBase = object$base1, outBase=object$base2,
              n.obs = object$n.obs,
              n.param = ncol(param), n.draws = n.draws,
              coef.table= if(n.cov > 1) param.table[1:n.cov,]
              else matrix(param.table[1,], nrow=1,
                          dimnames=list(rownames(param.table)[1], colnames(param.table))),
              cov.table=param.table[(n.cov+1):ncol(param),])  
  class(ans) <- "summary.endogMNP"
  return(ans)
}
