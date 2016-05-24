"krige" <-
function (s, point.obj, at, var.mod.obj, maxdist = NULL, extrap = FALSE, border=NULL) 
{
	    if (!inherits(point.obj, "point")) 
	    	    stop("point.obj must be of class, \"point\".\n")
	    # perform kriging on all given points:
	    if (!inherits(var.mod.obj, "variogram.model")) 
	    	    stop("var.mod.obj must be of class, \"variogram.model\".\n")
	    s$do <- c(rep(TRUE, length(s$x)))
	    # do nothing outside the convex hull?
	    # pull out the attribute vector...
	    if (!extrap) {
               if(is.null(border))
	    	    s$do <- in.chull(s$x,s$y,point.obj$x,point.obj$y) 
               else
                  if(is.null(border$x) | is.null(border$y) | 
                       length(border$x)!=length(border$y))
                     stop("border argument wrong!")
                  else
	    	    s$do <- in.polygon(s$x,s$y,border$x,border$y) 
	    }
	    at <- point.obj[[match(at, names(point.obj))]]
	    # if a maxdist hasn't been entered, then use all of the points...
	    if (is.null(maxdist)) 
	    	    krige.all(s, point.obj, at, var.mod.obj)
	    else krige.maxdist(s, point.obj, at, var.mod.obj, maxdist)
}
