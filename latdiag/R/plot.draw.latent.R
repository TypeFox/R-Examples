`plot.draw.latent` <-
function(x, graphtype = "png", ...) {
   if(is.null(x$rootname)) {
      warning("Cannot plot if rootname NULL")
      res <- list(retval = 256, cmnd = NULL)
   } else {
      cmnd <- paste("dot -T", graphtype, " ", x$rootname, ".dt -o ",
         x$rootname, ".", graphtype, sep = "")
      retval <- system(cmnd)
      res <- list(retval = retval, cmnd = cmnd)
   }
   invisible(res)
} # end plot method for draw.latent

