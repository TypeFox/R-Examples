cubecorners <- function(cdobj, fg="black",bg="white",size=0.5,circles=NULL, squares=NULL){
    ### function to draw custom symbols into the corners of the cube
    ### e.g. empty symbols
    ### that can be filled e.g. with averages or sample sizes
    ### and can be circles (default: all 8 corners are circles) or squares
    ### indicated by row numbers in cube data frame

    ### size refers to the side length of a square - circles have the same area as that square
    if (is.null(circles) & is.null(squares))
    stop("You need to specify at least one of circles or squares in function cubecorners.")
    data <- as.data.frame(cdobj$res3d$xyz.convert(cdobj$cub))
    if (!is.null(circles))
       symbols(data[circles,], circles=rep(size/sqrt(pi),length(circles)),
           bg=bg,fg=fg,add=TRUE, lwd=2, inches=FALSE)
    if (!is.null(squares))
       symbols(data[squares,], squares=rep(size,length(squares)),
           bg=bg,fg=fg,add=TRUE, lwd=2, inches=FALSE)
}