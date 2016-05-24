summary.ellipsefit <-
function(object,boot=TRUE,N=1000,studentize=TRUE,center=FALSE,cbb=NULL,joint=FALSE,seed=NULL,...) {
if (!is.null(seed)) set.seed(seed)
thecall <- match.call()
if (boot==TRUE) {
  if (object$method=="harmonic2")
  {res <- harmonic2summary(object,N=N,studentize=studentize,cbb=cbb,joint=joint)}
else if (object$method=="lm")
{res <- lmsummary(object,N=N,studentize=studentize,center=center,cbb=cbb,joint=joint)}
else if (object$method=="geometric")
{res <- geometricsummary(object,N=N,studentize=studentize,center=center,cbb=cbb,joint=joint)}
else if (object$method=="direct")
{res <- directsummary(object,N=N,studentize=studentize,center=center,cbb=cbb,joint=joint)}
  else {res <- nlssummary(object,N=N,studentize=studentize,center=center,cbb=cbb,joint=joint,...)}
  res$boot <- TRUE
res$residuals <- sqrt((res$x-res$pred.x)^2+(res$y-res$pred.y)^2)
  res$period.time <- object$period.time
}
else res <- object
res$summarycall <- thecall
  class(res) <- "ellipsesummary"
  res
  }
