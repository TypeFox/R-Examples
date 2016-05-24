# draw axes in the data ellipse computed by ellipse3d
ellipse3d.axes <-
function (x, centre = c(0, 0, 0), scale = c(1, 1, 1), level = 0.95,
    t = sqrt(qchisq(level, 3)), which = 1:3, labels=TRUE, label.ends=c(2,4,6), ...) 
{
    stopifnot(is.matrix(x)) 
    stopifnot(dim(x)[1] ==  dim(x)[2])  # square matrix?
    
    cov <- x[which, which]
    eig <- eigen(cov)
    # coordinate axes, (-1, 1), in pairs, for X, Y, Z
    axes <- matrix(
      c(-1, 0, 0,   1, 0, 0,
        0, -1, 0,   0, 1, 0,
        0, 0, -1,   0, 0, 1),  6, 3, byrow=TRUE)
	rownames(axes)<- apply(expand.grid(c("min","max"),c("X","Y","Z"))[,2:1],1,paste,collapse="")

	# transform to PC axes
    axes <- axes %*% sqrt(diag(eig$values)) %*% t(eig$vectors)
    result <- rgl::scale3d(axes, t, t, t)
    if (!missing(scale)) {
        if (length(scale) != 3) scale <- rep(scale, length.out=3) 
        result <- rgl::scale3d(result, scale[1], scale[2], scale[3])
        }
    if (!missing(centre)) {
        if (length(centre) != 3) scale <- rep(centre, length.out=3) 
        result <- rgl::translate3d(result, centre[1], centre[2], centre[3])
        }
    rgl::segments3d(result, ...)
    if (!missing(labels)) {
    	if (is.logical(labels) & labels) labels <- paste("PC", 1:3, sep="")
    	if (length(labels)==1) labels <- paste(labels, 1:3, sep="")
    	rgl::texts3d(result[label.ends,], texts=labels, ...)
    }
    invisible(result)
}
