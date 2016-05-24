summary.bayesQR <- function(object, CI=c(2.5, 97.5),...){

#  p <- object$n.alt
  param <- object$param
  nReg <- object$betLen
  n.var <- ncol(param) - nReg
  n.draws <- nrow(param)
  param.table <- cbind(apply(param, 2, mean), apply(param, 2, sd),
                       apply(param, 2, quantile, min(CI)/100),
                       apply(param, 2, quantile, max(CI)/100)) 
  colnames(param.table) <- c("mean", "std.dev.", paste(min(CI), "%", sep=""),
                             paste(max(CI), "%", sep=""))
  rownames(param.table) <- colnames(param)
  
#  ans <- list(call = object$call,
#              n.obs = object$nObs,
#              n.param = ncol(param), n.draws = n.draws,
#              coef.table= param.table[1:nReg,], 
#              cov.table=data.frame(param.table[(nReg+1):ncol(param),]))  

  ans <- list(call = object$call,
              n.obs = object$nObs,
              n.param = ncol(param), n.draws = n.draws,
              coef.table= param.table)  

  class(ans) <- "summary.bayesQR"
  return(ans)
}
