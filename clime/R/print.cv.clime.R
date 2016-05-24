print.cv.clime <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  print(cbind(lambda=signif(x$lambda,digits), mean=signif(x$loss.mean, digits), sd=signif(x$loss.sd, digits) ))
  cat("\n CV loss used:",x$loss,"\n CV optimal lambda=", signif(x$lambdaopt, digits), "\n")
}
