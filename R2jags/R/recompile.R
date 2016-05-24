recompile <- function(object, n.iter, refresh, progress.bar) UseMethod("recompile")


recompile.rjags <- function(object, n.iter = 100, refresh = n.iter/50, progress.bar = "text")
{
  object$model$recompile()
  adapt(object$model, n.iter = n.iter, by = refresh, progress.bar = progress.bar, end.adaptation=TRUE) 

}
  
recompile.rjags.parallel <- function(object, n.iter = 100, refresh = n.iter/50, progress.bar = "text")
{
  nchains <- length(object$model)
  for(i in 1:nchains){
    object$model[[i]]$recompile()
    adapt(object$model[[i]], n.iter = n.iter, by = refresh, progress.bar = progress.bar, end.adaptation=TRUE) 
  }
}
