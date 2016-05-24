#' Plot the nproc curve(s).
#' @export
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @param x fitted nproc object using \code{nproc}.
#' @param ... additional arguments.
#' @seealso \code{\link{npc}} and \code{\link{nproc}}
#' @examples
#' n = 1000
#' x = matrix(rnorm(n*2),n,2)
#' c = 1+3*x[,1]
#' y = rbinom(n,1,1/(1+exp(-c)))
#' fit = nproc(x, y, methods = 'svm')
#' plot(fit)
#' #fit = nproc(x, y, methods = c('svm','logistic'))
#' #plot(fit)
#' ##Compare Confidence Curves
#' #fit = nproc(x, y, methods = c('svm','logistic','lda'), conf = TRUE)
#' #plot(fit)
plot.nproc <- function(x, ...){

  roc.lo = x$roc.lo
  roc.up = x$roc.up
n.method = length(x$methods)
if(x$conf){
  plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab = "FPR", ylab = "TPR",
       main = paste("NP ROC: ", 1-x$delta, " Confidence Curve"))
} else {
  plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab = "FPR", ylab = "TPR",
       main = paste("NP ROC"))
}

for(i in 1:n.method){
  lines(x$alphalist, roc.lo[, 2*i],
       type = "s", col = i)
  if(x$conf){
  lines(x$alphalist, roc.up[, 2*i], col = i, type = 's')
  }
}
legend('bottomright',legend=x$methods, col=1:n.method, lty = rep(1,n.method))

}

