# Calculate the coordinate inside "pol" where a rectangle 
# of width "labw" and height "labh" centred on the coordinate
# has the largest possible distance from the polygon boundary.
# Uses a grid search with "gridpoints" coordinate pairs for the
# initial search, and refines the result using numerical 
# optimisation (the Nelder-Mead algorithm).
labpos.maxdist = function(pol, labw, labh, gridpoints) {

  # Create a rectangular polygon with a given size, centred at a given position
  makerect = function(x, y, wd, ht, rectID)
    Polygons(list(Polygon(cbind(c(x-wd/2, x+wd/2, x+wd/2, x-wd/2, x-wd/2), 
                                c(y-ht/2, y-ht/2, y+ht/2, y+ht/2, y-ht/2)))), rectID)

  # Calculate the distance from the given coordinates (2-column matrix)
  # to the polygon boundary. Returns 0 for any coordinates that results
  # in a rectangle not fully inside the polygon.
  getDists = function(co) {
    # Create a candidate label rectangle for each coordinate
    rects = SpatialPolygons(sapply(seq_len(nrow(co)),
                            function(i) makerect(co[i,1], co[i,2], labw, labh, i)),
                            proj4string=CRS(proj4string(pol)))

    # Check each rectangle to see if it"s inside the polygon
    inside = apply(gContains(pol, rects, byid = TRUE), 1, any)

    # Calculate and return the distance to the polygon boundary for
    # each rectangle (0 for any rectangle (partially) outside the polygon)
    di = numeric(nrow(co))
    if(any(inside))
      di[inside] = gDistance(pol.l, rects[inside], byid = TRUE)
    di
  }

  # Convert the polygon to lines, so that we can easily measure
  # the distance from the rectangles to the polygon boundary
  pol.l = as(pol, "SpatialLines")
  
  # Sample a regular grid of points inside the polygon
  co = coordinates(spsample(pol, n = gridpoints, type = "regular"))
  
  # Calculate the distances to the polygon boundary
  di = getDists(co)
  
  # Stop if no rectangles were inside the polygon
  if(all(di == 0))
    stop("Could not fit label inside polygon (with the given number of grid points)")
  
  # Use numerical optimisation to zero in on the optimal position,
  # with the coordinate found by grid search as initial values.
  # Unfortunately, the Nelder-Mead algorithm in "optim" has some
  # strange heuristics for the step size (size of the initial 
  # simplex) which is not possible to override, so, as a workaround, 
  # we do our own affine transformation inside a wrapper function.
  origin = co[which.max(di),]                                          # Initial value
  stepsize = min(apply(co, 2, function(x) diff(sort(unique(x))[1:2]))) # Distance between grid points
  
  getDistsScaled = function(co, origin, stepsize) {
    -getDists(matrix(co, ncol=2)*stepsize + origin)
  }
  optres = optim(c(0,0), getDistsScaled, method = "Nelder-Mead", 
               control = list(trace = FALSE, reltol = .001),
               origin = origin, stepsize = stepsize)
  optval = optres$par*stepsize + origin # "Optimal" label position
  
  # Only return the calculated position if the numerical optimisation succeeded
  if(optres$convergence != 0) {
    warning("Numerical optimisation did not converge. Returning result of grid search.")
    optval = origin
  }
  
  # Return the "optimal" label position
  optval
}


# Shrink the polygon (using a negative buffer) until
# the convex hull of the shrunken polygon is entirely
# contained inside the original polygon. Then calculate
# the label position as the centroid of the shrunken 
# polygon. (The centroid of the convex hull also gives
# similar results.)
labpos.buffer = function(pol, getLargestPolyPart) {
    init = 0                       # Initial amount to shrink
    estep = sqrt(gArea(pol)/pi)/10 # Additional amount to shrink for each step
    
    # Try repeatedly shrinking the polygon until we"re left
    # with a polygon whose convex hull fits inside the
    # original polygon.
    repeat {
      repeat {
        r = init + estep                 # Amount to shrink
        pol.b = gBuffer(pol, width = -r) # Shrink the polygon
        if (is.null(pol.b)) {
          estep = estep/2
        } else {
          if( gArea(pol.b) <= 0 )        # If the shrunken polygon is empty ...
            estep = estep/2 else break   # ... try again with a smaller value
        }
      }
      
      # If we"re left with more than one polygon, choose the largest one
      polb = getLargestPolyPart(pol.b)
      
      # Calculate the convex hull of the inner polygon.
      # If this is wholly contained in the original polygon,
      # break out of the loop and set the label point to
      # the centroid of the inner polygon.
      if( gContains(pol, gConvexHull(pol.b)) ) break else init=init+estep
    }
    coordinates(pol.b)
}


# Generate a random position inside the polygon
labpos.random = function(pol) {
  # Note that the "spsample" function sometimes fail to
  # find a random point inside the polygon (on the first try),
  # so we may need repeated sampling. Eventually a point is found.
  repeat {
    coord.rand = tryCatch(spsample(pol, n = 1, type = "random"), error = function(x) NULL)
    if(!is.null(coord.rand)) break
  }
  coordinates(coord.rand)
}


# Calculate optimal polygon labels.
#     pols: SpatialPolygons(DataFrame) object.
#   labels: The labels to apply to each polygon.
#           If NULL, the positition is calculated as if the
#           label is a square with one line long sides
#           (given the current value of "cex").
#   method: The method(s) used to calculate label positions.
#           Possible values: maxdist, buffer, centroid,
#                            random, pointonsurface
# polypart: Which parts each multipolygon to use ("all" or "largest").
#           Note that "largest" also removes any holes before
#           calculating the label position, so the labels are
#           no longer guaranteed not to overlap a hole.
#   doplot: Should the labels be plotted on the current graphics device?
#      ...: Additional arguments sent to "text"
polygonsLabel = function(pols, labels = NULL, method = c("maxdist",
                         "buffer", "centroid", "random", "inpolygon")[1],
                         gridpoints = 60, polypart = c("all", "largest")[1],
                         cex = 1, doPlot = TRUE, ...) {

  # Return a SpatialPolygons object containing the
  # biggest polygon part of a multipolygon
  getLargestPolyPart = function(multipol) {
    pols = multipol@polygons[[1]]@Polygons
    
    if(length(pols) > 1) {
      areas = lapply(pols, function(x) x@area) # List of areas of the polygon parts
      pols = pols[which.max(areas)]            # The largest polygon part
    
    # Create new SpatialPolygons object only containing the largest polygon part
    multipol = SpatialPolygons(list(Polygons(pols, ID="1")),
                             proj4string=CRS(as.character(proj4string(multipol))))
    }
    
    multipol
  }
  
  # If no labels are given, assume empty labels
  if(is.null(labels)) labels = ""
  
  # If ncessary, convert the labels into character labels
  # (needed for factor labels)
  if( !is.character(labels))
    labels=as.character(labels)
  
  # Recycle all arguments so they contain the same
  # number of elements as there are polygons
  n = length(pols)
  labels = rep(labels, length = n)
  method = rep(method, length = n)
  polypart = rep(polypart, length = n)
  
  # For each polygon, calculate the optimal label position
  ret = matrix(ncol=2, nrow = n)
  for (i in seq_along(labels)) {
      labw = strwidth(labels[i], cex = cex)  # Label width
      labh = strheight(labels[i], cex = cex) # Label height
      if(labw == 0) labw = labh # For the empty string, use a one-line high square
      pol = pols[i,]            # The polygon of interest
      
      # Optionally remove smaller polygon parts
      if(polypart[i] == "largest") pol = getLargestPolyPart(pol)
          
      # Calculate the label optimal label position
      ret[i,] = switch(method[i], "maxdist" = labpos.maxdist(pol, labw, labh, gridpoints),
                                  "buffer" = labpos.buffer(pol, getLargestPolyPart=getLargestPolyPart),
                                  "centroid" = coordinates(pol),
                                  "random" = labpos.random(pol),
                                  "inpolygon" = coordinates(gPointOnSurface(pol)),
                                  stop(paste("Unknown method:", method[i])))
  }
  
  # Optionally plot the labels
  if(doPlot)
    text(ret, labels, cex=cex, ...)
  
  ret
}

# TODO:
#
#   Automatic reprojection for unprojected maps
#     (simple quirectangular, to match what "plot" does).
#
#   Option to auto-rotate labels (default < +/-30 degrees)
#     to make them fit (better).
#
#   Better error handling (option to stop or to use
#   centroids when no valid label positions are found).
#
#   A "removeholes" option?
