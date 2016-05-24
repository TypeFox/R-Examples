summary.kohonen <- function(object, ...)
{
  cat(object$method, " map of size ",
      object$grid$xdim, "x", object$grid$ydim,
      " with a ", object$grid$topo,
      if (object$toroidal) "toroidal", " topology.", sep="")

  if (!is.null(object$data)) {
    switch(object$method,
           som = {
             cat("\nTraining data included; dimension is",
                 nrow(object$data), "by", ncol(object$data))
           },
           supersom = {
             cat("\nTraining data included of ",
                 nrow(object$data[[1]]), "objects")
             cat("\nThe number of layers is", length(object$data))
             if (length(object$data) > length(object$whatmap))
               cat(", of which", length(object$whatmap),
                   "have been used in training.")
           },
           {
             cat("\nTraining data included; dimension is",
                 nrow(object$data), "by", ncol(object$data))
             cat("\nDimension of Y:", nrow(object$Y), "by", ncol(object$Y))
             if (!is.null(object$predict.type)) {
               cat("\nPrediction type:",
                   ifelse(object$predict.type == "class",
                          "classification",
                          "regression"))
             }
           }
           )

    cat("\nMean distance to the closest unit in the map:",
        mean(object$distances))
  } else {
    cat("\nNo training data included in the object.")
  }

  cat("\n")

  invisible()
}

print.kohonen <- function(x, ...)
{
  cat(x$method, " map of size ", x$grid$xdim, "x", x$grid$ydim,
      " with a ", x$grid$topo, if (x$toroidal) " toroidal",
      " topology.", sep="")
  if (!is.null(x$data))
    cat("\nTraining data included.")
  cat("\n")
}

