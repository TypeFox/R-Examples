##' @export
model.frame.lvmfit <- function(formula, all=FALSE,...) {
  dots <- list(...)
  mydata <- formula$data$model.frame
  if (!is.data.frame(mydata) & !is.matrix(mydata))
    return(mydata)
  if (all) return(mydata)
##  xfix <- colnames(mydata)[(colnames(mydata)%in%parlabels(formula$model0,exo=TRUE))]
  xfix <- colnames(mydata)[(colnames(mydata)%in%parlabels(formula$model0))]
  return( mydata[,c(manifest(formula),xfix),drop=FALSE] )
}

##' @export
model.frame.multigroupfit <- function(formula,...) {
  dots <- list(...)
  mydata <- formula$model$data
  return(mydata)
}
