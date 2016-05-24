## gpclib:  General Polygon Clipping library for R
## Copyright (C) 2003-2010 Roger D. Peng <rpeng@jhsph.edu>


## R functions for using GPC library and manipulating polygons

setClass("gpc.poly", representation(pts = "list"))
setClass("gpc.poly.nohole", "gpc.poly")

## setValidity("gpc.poly",
##             function(x) {
##                 pts <- x@pts
##                 lens <- sapply(pts, length)
##                 if(!identical(all(lens == 3), TRUE))
##                     "Not all contours have correct elements"
##                 ## correct names (x, y, hole)
##                 contour.names <- lapply(pts, names)
##                 correct.names <- c("x", "y", "hole")
##                 check <- lapply(contour.names, function(n) {
##                     n %in% correct.names
##                 })
##                 if(!identical(all(sapply(check, all)), TRUE))
##                     "Incorrect list names in contours"
##                 ## correct types
##                 ## equal lengths?
##             })

setMethod("show", "gpc.poly",
          function(object) {
              cat("GPC Polygon\n")
              cat("   Num. Contours: ", length(object@pts), "\n")
              if(length(object@pts) == 1)
                  cat("   Num. Vertices: ", length(object@pts[[1]]$x),"\n")
              bbox <- get.bbox(object)
              cat("   BBox (X): ", bbox$x[1], "-->", bbox$x[2], "\n")
              cat("   BBox (Y): ", bbox$y[1], "-->", bbox$y[2], "\n")
              invisible(object)
          })

setGeneric("get.bbox", function(x)
           standardGeneric("get.bbox"))

setMethod("get.bbox", signature(x = "gpc.poly"),
          function(x) {
              pts <- x@pts
              x <- unlist(lapply(pts, "[[", "x"))
              y <- unlist(lapply(pts, "[[", "y"))
              
              if(!is.null(x))
                  xlim <- range(x)
              else
                  xlim <- c(NA, NA)
              if(!is.null(y))
                  ylim <- range(y)
              else
                  ylim <- c(NA, NA)
              list(x = xlim, y = ylim)
          })

setGeneric("plot")

setMethod("plot", "gpc.poly",
          function(x, y, poly.args = list(), xlab = "X", ylab = "Y",
                   asp = 1, add = FALSE, ...) {
              if(!add) {
                  bbox <- get.bbox(x)
                  plot(0, 0, ylim = bbox$y, xlim = bbox$x, type="n",
                       xlab = xlab, ylab = ylab, asp = asp, ...)
              }
              invisible(lapply(x@pts, function(p) {
                  do.call("polygon", append(list(x = p), poly.args))
              }))
          })

setGeneric("intersect")
setGeneric("union")
setGeneric("setdiff")
if (!isGeneric("symdiff")) setGeneric("symdiff", function(x, y)
		standardGeneric("symdiff"))

setMethod("intersect", signature(x = "gpc.poly", y = "gpc.poly"),
          	function(x, y) {
		  		spx = as(x,"SpatialPolygons")
				spy = as(y,"SpatialPolygons")
				spres = gIntersection(spx, spy, drop_lower_td=TRUE)
		
                                if (is.null(spres))
                                    return(new("gpc.poly"))
				as(spres,"gpc.poly")
          	})


setMethod("setdiff", signature(x = "gpc.poly", y = "gpc.poly"),
          	function(x, y) {
		  		spx = as(x,"SpatialPolygons")
				spy = as(y,"SpatialPolygons")
				spres = gDifference(spx,spy, drop_lower_td=TRUE)
		
                                if (is.null(spres))
                                    return(new("gpc.poly"))
				as(spres,"gpc.poly")
          	})

setMethod("union", signature(x = "gpc.poly", y = "gpc.poly"),
          function(x, y) {
		  		spx = as(x,"SpatialPolygons")
				spy = as(y,"SpatialPolygons")
				spres = gUnion(spx,spy, drop_lower_td=TRUE)
		
                                if (is.null(spres))
                                    return(new("gpc.poly"))
				as(spres,"gpc.poly")
          })

setMethod("symdiff", signature(x = "gpc.poly", y = "gpc.poly"),
          function(x, y) {
		  		spx = as(x,"SpatialPolygons")
				spy = as(y,"SpatialPolygons")
				spres = gSymdifference(spx,spy, drop_lower_td=TRUE)

                                if (is.null(spres))
                                    return(new("gpc.poly"))
				as(spres,"gpc.poly")
          })

setGeneric("tristrip", function(x) 
           standardGeneric("tristrip"))
           
setMethod("tristrip", signature(x = "gpc.poly"),
		function(x) {
			stop("rgeos does not currently implement this feature of gpclib")
		
	      	#poly <- as(x, "numeric")
	      	#result <- .Call("Rgpc_polygon_to_tristrip", poly, PACKAGE = "gpclib")
	      	#result <- lapply(result, function(strip) matrix(strip, ncol=2, byrow=TRUE))
	  	})

setGeneric("triangulate", function(x) 
           standardGeneric("triangulate"))
 
setMethod("triangulate", signature(x = "gpc.poly"),
		function(x) {
			stop("rgeos does not currently implement this feature of gpclib")
	      	tristrip <- tristrip(x)
	      	triangles <- lapply(tristrip, 
	                     function(strip) {
	                          n <- nrow(strip)
	                          if (n > 3)
	                              result <- strip[c(1:3,4:2) + 2*rep(0:(n %/% 2 - 2), each=6), ]
	                          else
	                              result <- strip[0,]
	                          if (n %% 2 && n > 2) result <- rbind(result, strip[1:3 + n - 3,])
	                          return(result)
	                          })
	                          
	      	do.call(rbind, triangles)
	  	})

setMethod("[", "gpc.poly",
          function(x, i, j, ..., drop = FALSE) {
              new("gpc.poly", pts = x@pts[i])
          })

setAs("matrix", "gpc.poly",
      function(from, to) {
          if(ncol(from) > 2)
              stop("matrix must have 2 columns")
          p <- list(x = from[,1], y = from[,2], hole = FALSE)
          new("gpc.poly", pts = list(p))
      })

setAs("data.frame", "gpc.poly",
      function(from, to) {
          as(as.matrix(from), "gpc.poly")
      })

## Miscellaneous Utility Functions

setGeneric("append.poly", function(x, y)
           standardGeneric("append.poly"))

setMethod("append.poly",
          signature(x = "gpc.poly", y = "gpc.poly"),
          function(x, y) {
              newpts <- append(x@pts, y@pts)
              new("gpc.poly", pts = newpts)
          })

setGeneric("scale.poly", function(x, ...)
           standardGeneric("scale.poly"))

setMethod("scale.poly", signature(x = "gpc.poly"), 
          function(x, xscale, yscale = xscale, ...) {
              x@pts <- lapply(x@pts, function(p) {
                  p$x <- p$x / xscale
                  p$y <- p$y / yscale
                  p
              })
              x
          })

## Compute the area of each polygon in the polygon set contained in
## `object'.  

setGeneric("area.poly", function(object, ...)
           standardGeneric("area.poly"))

setMethod("area.poly", signature(object = "gpc.poly"),
          function(object, ...) {
              area <- function(x.mat) {
                  if(nrow(x.mat) < 3) 
                      return(0);   
                  x.segmat <- cbind(x.mat, rbind(x.mat[2:nrow(x.mat), ],
                                                 x.mat[1, ]));
                  abs(sum(x.segmat[,1] * x.segmat[,4] - x.segmat[,3]
                          * x.segmat[,2])) / 2
              }
              if(length(object@pts) == 0)
                  return(0)
              a <- sapply(object@pts, function(p) area(cbind(p$x, p$y)))
              holeflags <- sapply(object@pts, "[[", "hole")
              sum(a[!holeflags]) - sum(a[holeflags])
          })

## Added 2003/01/07
setGeneric("get.pts", function(object)
           standardGeneric("get.pts"))

setMethod("get.pts", signature(object = "gpc.poly"),
          function(object) {
              object@pts
          })

## These two functions are needed for the intersect/union/setdiff
## methods.  They basically serialize and unserialize the "gpc.poly"
## object.

## We need separate functions because polygons with holes and polygons
## without holes have different file formats

setAs("numeric", "gpc.poly.nohole", 
      function(from) {
          ## The shortest a vector can be is 8 numbers:  1. Num. Contours;
          ## 2. Num pts for first contour; and three vertices
          if(length(from) < 8)
              stop("numeric vector not long enough")
          expand.poly <- function(x) {
              ## `x' is just a long vector of numbers with a special format
              num.contours <- x[1]; x <- x[-1]
              polyfile <- x
              poly <- vector("list", length = num.contours)
              
              for(i in 1:num.contours) {
                  npts <- polyfile[1]; polyfile <- polyfile[-1]
                  m <- matrix(polyfile[1:(2*npts)], byrow = TRUE, ncol = 2)
                  poly[[i]] <- list(x = m[,1], y = m[,2], hole = FALSE)
                  polyfile <- polyfile[-(1:(2*npts))]
              }
              poly
          }
          new("gpc.poly.nohole", pts = expand.poly(from))
      })

setAs("numeric", "gpc.poly", 
      function(from) {
          ## The shortest a vector can be is 9 numbers:  1. Num. Contours;
          ## 2. Num pts for first contour; 3. hole flag; and three vertices
          if(length(from) < 9)
              stop("numeric vector not long enough")
          expand.poly <- function(x) {
              num.contours <- x[1]; x <- x[-1]
              polyfile <- x
              poly <- vector("list", length = num.contours)
              
              for(i in 1:num.contours) {
                  npts <- polyfile[1]
                  polyfile <- polyfile[-1]
                  hole <- as.logical(polyfile[1])
                  polyfile <- polyfile[-1]
                  
                  m <- matrix(polyfile[1:(2*npts)], byrow = TRUE, ncol = 2)
                  poly[[i]] <- list(x = m[,1], y = m[,2], hole = hole)
                  polyfile <- polyfile[-(1:(2*npts))]
              }
              poly
          }
          new("gpc.poly", pts = expand.poly(from))
      })

##

setAs("gpc.poly", "numeric",
      function(from) {
          flatten.poly <- function(poly) {
              num.contours <- length(poly@pts)
              flat <- lapply(poly@pts, function(p)
                         {
                             v <- as.vector(t(cbind(p$x, p$y)))
                             c(length(p$x), as.numeric(p$hole), v)
                         })
              c(num.contours, unlist(flat))
          }
          flatten.poly(from)
      })

setAs("gpc.poly", "matrix",
      function(from) {
          if(length(from@pts) > 1)
              stop("can only convert a single contour into a matrix")
          pts <- from@pts[[1]]
          m <- cbind(x = pts$x, y = pts$y)
      })

## 'from' is a list(x = ..., y = ...)

setAs("list", "gpc.poly",
      function(from) {
              if(!(all(c("x", "y") %in% names(from))))
                      stop("list should have names 'x' and 'y'")
              if(length(from$x) != length(from$y))
                      stop("'x' and 'y' elements should have the same length")
              as(cbind(from$x,from$y), "gpc.poly")
      })


## Read a polygon from a file

read.polyfile <- function(filename, nohole = TRUE) {
        polyfile <- scan(filename, quiet = TRUE)
        if(nohole) 
                as(polyfile, "gpc.poly.nohole")
        else
                as(polyfile, "gpc.poly")
}

## Write a "gpc.poly" object to a text file

write.polyfile <- function(poly, filename = "GPCpoly.txt") {    
    if(!is(poly, "gpc.poly"))
        stop("'poly' should be of class 'gpc.poly'")
    outfile <- file(filename, "w")
    on.exit(close(outfile))

    num.contours <- length(poly@pts)   
    cat(num.contours, "\n", file = outfile)

    for(i in 1:num.contours) {
        m <- as(poly[i], "matrix")
        cat(nrow(m), "\n", file = outfile, append = TRUE)

        if(!is(poly, "gpc.poly.nohole"))
            cat(as.numeric(poly@pts[[i]]$hole), "\n",
                file = outfile, append = TRUE)       
        write(t(m), file = outfile, ncolumns = 2, append = TRUE)
    }
}


