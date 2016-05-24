## added conditional recomputation of maximum
## NB: if there is nothing in chull,
## surface.OKrig will generate an error. Hence it is more logical
## to test the grid before running surface.OKrig (=>predict_surface)
## To that aim, I defined the safeSurface.OKrig() wrapper
predict_surface <- function (object, grid.list = NA,
                            extrap = ("extrapolateOutOfHull" %innc% blackbox.getOption("miscOptions")),
                            redundantpts = NULL, ## FR->FR this does not seem to be used anywhere in current code
                            nx = 80, ny = 80, xy = c(1, 2), order.variables = "xy", ...)
{ ## note that grid.list should not be NA and has info for *all* fittedNames
  if (is.na(grid.list)[1]) {
    if (is.null(object$x)) {
      stop.redef("Need an X matrix in the output object")
    }
    grid.list <- calcGridFromxy(object$x, nx = nx, ny = ny,
                                  xycols = xy)
  }
  xg <- make.surface.grid(grid.list)
  z <- predict(object, xg, ## predict.OKrig or Migraine's predict.OKriglist
               testhull=FALSE ##passed to purefn ## FALSE 25/05/2009
  )
  tmp <- as.surface(grid.list, z, order.variables = order.variables)
  if (!extrap) {
    if (!is.null(redundantpts)) {  ## then uses local redundantpts else global testspace
      ## it's just the set of parameter points in some adequate format
      message.redef("(!) local hull used in predictSurface. Consider redundant.mode argument there")
      constraints <- resetCHull(redundantpts, formats="constraints")[c("a", "b")] ## here local
    } else {
      constraints <- blackbox.getOption("hulls")$Kgtotal[c("a", "b")]
    }
    if (is.null(constraints)) {
      message.redef("(!) From predictSurface() : !extrap but is.null(constraints)")
    }
    inConvexHull <- apply(xg, 1, isPointInCHull, constraints=constraints) # a vector of T/F
    z[!inConvexHull] <- NA
    valsInConvex <- z[inConvexHull]
  } else { ## extrap
    valsInConvex <- z
  }
  if(length(valsInConvex)==0) {
    #  	if(extrap) {
    #			print("NB Plot range is completely outside convex envelope of data.")
    #			print("    but extrap=TRUE. -- (message from predict_surface)")
    #		} else {
    #			print("(!) Plot range completely outside convex envelope of data.")
    #			print("    Use the ExtrapolateSurfaces setting, or else")
    #			print("    surface.OKrig(..., extrap=TRUE, ...) [in R code]")
    #			print("    to plot extrapolated surface. -- (message from predict_surface)")
    #		}
  } else { # exists value(s) in convex hull
    if (sort(valsInConvex)[1] > blackbox.getOption("rosglobal")$value) {
      ## message operates in interactive mode and print's go to the R_out... file
      if (interactive()) {message.redef("predict_surface found a better maximum")}
      cat("predict_surface found a better maximum", "\n")
      ## as.data.frame(xg[which.max(z), , drop=FALSE]) to keep names as 'names', (vs 'colnames' in matrix )
      rosglobal <- findglobalMLE(initptinfK=as.data.frame(xg[which.max(z), , drop=FALSE])) ## < <- rosglobal, 2Ds2
      blackbox.options(rosglobal=rosglobal)

    }
  }
  out <- as.surface(grid.list, z, order.variables = order.variables) ## a may have only NA's
  out$zextrap <- tmp$z ## does not have only NA's
  out ## return value
}
