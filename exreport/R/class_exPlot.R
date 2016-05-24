is.exPlot <- function(x) {
  is(x, "exPlot")
}

#' @export
print.exPlot <- function (x, ...) {
  print(x$plot)
  print(x$tags$title)
}

#' @export
summary.exPlot <- function (object, ...) {
  print(object$tags$title)
  print(object$plot)
}

#Anonymous constructor
.exPlot <- function(ggplot2Obj, target=NULL, title,  alias = "plot", tags){
  
  if (!is.null(target))
    newTags <- .metaTags(title = title, target = target, alias = alias)
  else
    newTags <- .metaTags(title = title, alias = alias)
  
  tags <- .updateTags(tags, newTags)
  
  plo <- list(
    "plot"     = ggplot2Obj,
    "tags"     = tags
  )
  
  class(plo) <- c("exPlot", "reportable")
  
  plo
}