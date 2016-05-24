## =============================================================================
## from the help page of rgl.setMouseCallbacks
## =============================================================================

 pan3d <- function(button) {
   start <- list()

   begin <- function(x, y) {
       start$userMatrix <<- par3d("userMatrix")
       start$viewport <<- par3d("viewport")
       start$scale <<- par3d("scale")
       start$projection <<- rgl.projection()
       start$pos <<- rgl.window2user( x/start$viewport[3], 1 - y/start$viewport[4], 0.5,
                                      projection=start$projection)
   }

   update <- function(x, y) {
        xlat <- (rgl.window2user( x/start$viewport[3], 1 - y/start$viewport[4], 0.5,
                                 projection = start$projection) - start$pos)*start$scale
        mouseMatrix <- translationMatrix(xlat[1], xlat[2], xlat[3])
        par3d(userMatrix = start$userMatrix %*% t(mouseMatrix) )
   }
   rgl.setMouseCallbacks(button, begin, update)
#   cat("Callbacks set on button", button, "of rgl device",rgl.cur(),"\n")
 }


