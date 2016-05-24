summary.treeClust <- function(object, ...)
{
#
# Summary method for objects of class treeClust.
#
      cat("Call:\n")
      cat(deparse(object$call), "\n")
      if (any (names (object) == "mat"))
          cat("Structure is ", nrow(object$mat), "by", ncol(object$mat), "\n")
      if (object$final.algorithm == "None")
        cat ("No final clustering was performed\n")
      else        
        cat (paste ("Final clustering by", object$final.algorithm, "\n"))

}

