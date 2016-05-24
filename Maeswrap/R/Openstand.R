#' @importFrom rgl lines3d
Openstand <- function(treesfile="trees.dat"){

    x0 <- readPAR(treesfile,"x0","plot",fail=FALSE)
    y0 <- readPAR(treesfile,"y0","plot",fail=FALSE)
    if(is.na(x0))x0 <- 0
    if(is.na(y0))y0 <- 0
    
    xmax <- readPAR(treesfile,"xmax","plot")
    ymax <- readPAR(treesfile,"ymax","plot")
    
  	M <- matrix(c(x0,y0,0,
  			      xmax,y0,0,
  				  xmax,ymax,0,
  				  x0,ymax,0,
  				  x0,y0,0), ncol=3, byrow=TRUE)
    lines3d(M, col="darkgrey")
    
}
