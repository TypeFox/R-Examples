# importFrom("graphics", "abline", "layout", "legend", "lines", "par", "plot", "plot.new")
# importFrom("stats", "IQR", "as.formula", "dist", "fft", "lm", "median", "rnorm", "runif", "sd", "ts")
#' @importFrom graphics abline layout legend lines par plot plot.new
#' @importFrom stats IQR as.formula dist fft lm median rnorm runif sd ts
NULL

vectorizePar = function(par, N, default=1:N){
  if (is.null(par)) par=default
  if (length(par) < N) par=rep(par,length.out=N)
  par
}
                        

propagateTakensAttr = function(x, takens){
  attr(x, "id") = attr(takens,"id")
  attr(x, "time.lag") = attr(takens, "time.lag")
  attr(x, "embedding.dim") = attr(takens, "embedding.dim")
  x
}
