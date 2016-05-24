plot.cv.SGL = function(x,...){
  cvobj = x
  xlab = "log(Lambda)"
  ylab = "CV Negative Log Likelihood"
  cvUp <- cvobj$lldiff + cvobj$llSD
  cvDn <- cvobj$lldiff - cvobj$llSD
  plot.args = list(x=log(cvobj$lambdas), y = cvobj$lldiff, xlab = xlab, ylab = ylab, ylim = range(cvUp, cvDn))
  do.call("plot", plot.args)
  error.bars(log(cvobj$lambdas), cvUp, cvDn, width = 0.01, col ="darkgrey")
  points(log(cvobj$lambdas),cvobj$lldiff,pch=20,col="red")
}
