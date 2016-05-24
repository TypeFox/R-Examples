`print.draw.latent` <-
function(x, ...) {
   cat("Diagram has patterns with", x$which.npos, "items positive\n")
   cat("Original order was", x$new.order, "\n")
   if(!is.null(x$rootname))
      cat("Commands output to ", x$rootname, ".dt\n", sep = "")
   invisible(x)
} # end print method for draw.latent

