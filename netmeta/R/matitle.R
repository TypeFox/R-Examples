matitle <- function(x){
  tl <- options()$width-7
  ##
  if (!is.null(x$title))
    if (x$title!="")
      if (nchar(x$title) <= tl)
        cat("Title: ", x$title, "\n\n", sep="")
      else
        cat("Title: ", substring(x$title, 1, tl-4),
            " ...\n\n", sep="")
}
