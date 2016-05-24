# FField.R ####################################################################
# FField Package
# Author: Grigori Kapoustin, 2013
# License: GPL-3
# Force field simulation for mutual repulsion by set of points.
# Very useful for placing text labels on graphs, such as scatterplots.

# 'wt' should be regarded as defined globally 
# when the check tool is applied to this package. 
if (getRversion() >= '2.15.1') {
  globalVariables("wt")
}

FFieldPtRep <- function(coords,
                        rep.fact = 20,
                        rep.dist.lmt = 10,
                        attr.fact = 0.2,
                        adj.max = 0.1,
                        adj.lmt = 0.5,
                        iter.max = 10000) {
  # Performs force field simulation for mutual repulsion by set of points.
  # Points experience repulsion from one another and attraction to
  # their original positions.
  # Repulsion is inversely proportional to the
  # square of the distance.
  # Attraction is directly proportional to the distance.
  # Very useful for placing text labels on graphs, such as scatterplots.
  # Depending on the nature of the plot, parameters may need to be masaged
  # for the simulation to converge.
  # Assumes 1x1 coordinate aspect ration and re-scaling of inputs
  # may be needed.
  #
  # Args:
  #   coords: matrix or data.frame consisting of two columns 
  #     (x and y coordinates).
  #   rep.fact: repulsion force factor.
  #   rep.dist.lmt: repulsion distance limit.
  #   attr.fact: attraction force factor.
  #   adj.max: maximum position adjustment at each iteration.
  #   adj.lmt: position adjustment limit at which the simulation stops.
  #   iter.max: the maximum number of iterations beyond which simulation
  #     will end and a warning will be reported.  
  #
  # Returns:
  #   coordinates of the points at completion of the simulation
  
  if (length(dim(coords)) != 2) {
    stop("FFieldPtRep: dim(coords) must be 2\n")    
  }
  if (ncol(coords) < 2) {
    stop("FFieldPtRep: ncol(coords) must be >= 2\n")    
  }
  if (nrow(coords) < 2) {
    stop("FFieldPtRep: nrow(coords) must be >= 2\n")    
  }
  
  coords <- as.data.frame(coords)
  colnames(coords)[(1:2)] <- c("x", "y")  
  coords.orig <- coords
  
  FVCalc <- function(vects.x, 
                     vects.y, 
                     f.fact, 
                     f.type = "invsq") {
    # Force vector calculation common code    
    
    d.sq <- (vects.x ^ 2 + vects.y ^ 2)
    d <- sqrt(d.sq)
    
    # Normalize the vectors
    vects.x <- vects.x / d
    vects.y <- vects.y / d
    
    # Get the force vector matrices
    if (f.type == "invsq") {
      d.sq[d >= rep.dist.lmt] <- Inf
      vect.f.x <- vects.x / d.sq * f.fact
      vect.f.y <- vects.y / d.sq * f.fact  
    } else if (f.type == "lin") {
      vect.f.x <- vects.x * d * f.fact
      vect.f.y <- vects.y * d * f.fact    
    } else {
      stop("FFieldPtRep: Unexpected f.type\n")
    }
    
    # Remove NaNs that occur when d == 0 
    # (occuring when calculating the repulsion of point
    # and itself and the attraction of point at the origin and the origin).
    vect.f.x[is.na(vect.f.x)] <- 0
    vect.f.y[is.na(vect.f.y)] <- 0
    
    # Combine the force vectors acting upon each point
    f.vect <- cbind(colSums(vect.f.x), colSums(vect.f.y))
    return (f.vect)
  }
  
  iter <- 0
  
  repeat {
  
    # Calculate repulsion forces.
    # Direction is from other points to a given point.
    vects.x <- apply(coords, 1, function(c) (c[1] - coords$x))
    vects.y <- apply(coords, 1, function(c) (c[2] - coords$y))    
    f.rep.v <- FVCalc(vects.x = vects.x, 
                      vects.y = vects.y, 
                      f.fact = rep.fact, 
                      f.type = "invsq")
    
    # Calculate attraction forces.
    # Direction is from each point to its original position.
    vects.orig <- coords.orig - coords
    f.attr.v <- FVCalc(vects.x = t(as.matrix(vects.orig$x)), 
                       vects.y = t(as.matrix(vects.orig$y)), 
                       f.fact = attr.fact, 
                       f.type = "lin")
    
    # Combine the forces.
    f.v <- f.rep.v + f.attr.v
    if (all(abs(f.v) <= adj.lmt)) {
      break()
    }
    
    # Adjust the coordinates    
    mv.vect <- apply(f.v, 
                     c(1, 2), 
                     function(x) sign(x) * min(abs(x), adj.max))    
    coords <- coords + mv.vect    
    
    if ((iter <- iter + 1) > iter.max) {
      warning("FFieldPtRep: Maximum iterations exceeded ",
              "without convergence.\n")
      break()
    }
  }
  
  return(coords)
}

FFieldPtRepDemo <- function() {
  # Demonstrates the utility of FFieldPtRep for placing labels
  # in a scatterplot.
  #
  # Args:
  #   none
  #
  # Returns:
  #   none  
  
  library(ggplot2)
  library(gridExtra)
  
  # Sample plot with crowded labels
  p1 <- 
    (ggplot(mtcars, aes(x = wt, 
                        y = mpg, 
                        label = rownames(mtcars)))  
     + geom_point()
     + geom_text()
     + ggtitle("Before"))
  
  # Normalize coordinates to maintain constant aspect ratio
  x.fact <- 100 / max(mtcars$wt)
  y.fact <- 100 / max(mtcars$mpg)
  
  # Repel points
  coords <-
    FFieldPtRep(coords = cbind(mtcars$wt * x.fact, 
                               mtcars$mpg * y.fact),
                rep.fact = 40)
  
  # Convert back to plot coordinates
  x.t <- coords$x / x.fact
  y.t <- coords$y / y.fact
  
  # Sample plot with repelled labels
  p2 <- 
    (ggplot(mtcars, aes(x = wt, 
                        y = mpg, 
                        label = rownames(mtcars)))  
     + geom_point()
     + geom_text(x = x.t,
                 y = y.t)
     + geom_segment(data = mtcars,
                    xend = x.t,
                    yend = y.t)
     + ggtitle("After"))
  
  grid.arrange(p1, p2)
}
