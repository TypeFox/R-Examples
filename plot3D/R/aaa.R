## ====================================================================
## A local environment for non user-visible data,
## ====================================================================

.plot3D <- new.env()

.plot3D$plist <- list()

.plot3D$refresh <- TRUE

refresh <- function(set = TRUE)     # internal function...
  .plot3D$refresh <- set
  
getplist <- function()
  .plot3D$plist

setplist <- function(plist)
  .plot3D$plist <- plist

initplist <- function(add) {
  if (add) 
    plist <- getplist()
  else
    plist <- NULL
# test for setting the correct plt parameters:
  if (!add & .plot3D$refresh) {
    plot.new()
    par(new = TRUE)
  }  
  plist
}    