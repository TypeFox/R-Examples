print.treeClust <- function(x, ...)
{
#
# Print method for objects of class treeClust.
#
      cat("Call:\n")
      cat(deparse(x$call), "\n")
      if (any (names (x) == "mat"))
          cat("Structure is ", nrow(x$mat), "by", ncol(x$mat), "\n")
      if (x$final.algorithm == "None")
        cat ("No final clustering was performed\n")
      else        
        cat (paste ("Final clustering by", x$final.algorithm, "\n"))
      print(x$tbl)
}

