# New generic functions
# nNode = function(obj) UseMethod("nNode")
# nChild = function(obj) UseMethod("nChild")
# totalLength = function(obj) UseMethod("totalLength")
# coords = function(obj) UseMethod("coords")
# yrange = function(obj) UseMethod("yrange")
# xrange = function(obj) UseMethod("xrange")
# yrange = function(obj) UseMethod("zrange")
# insertionPosition = function(obj) UseMethod("insertionPosition")
# meanInsertionAngle = function(obj) UseMethod("meanInsertionAngle")
#meanInterbranch = function(obj, allroot=T) UseMethod("meanInterbranch")
# nLatRoot = function(obj) UseMethod("nLatRoot")
# nPrimRoot = function(obj) UseMethod("nPrimRoot")
# nRoot = function(obj) UseMethod("nRoot")  

############################################################

#' Plant object, containing a root system, composed of roots
#' @param roots    the root object contained in the plant. Can be null and incremented afterward
#' @keywords rsml
#' @author Guillaume Lobet - guillaume.lobet(at)ulg.ac.be
#' @return the plant
#' @export
#' @examples
#' pl <- plant()
plant = 
  function(roots = NULL)
  {
    pls = list(roots = roots)
    class(pls) = "plant"
    pls
  }

############################################################

#' Add a root to an existing plant. Returns the plant with the added root
#' @param pl  the plant to add the root to
#' @param ro  the root object to add to the plant.
#' @keywords rsml
#' @author Guillaume Lobet - guillaume.lobet(at)ulg.ac.be
#' @return the new plant, with the added root
#' @export
#' @examples
#' data(lupin)
#' r <- root()
#' lupin <- addRootToPlant(lupin, r)
addRootToPlant = 
  function(pl, ro)
  {
    if(class(pl) != "plant" || class(ro) != "root" )
      stop("need objects of class plant and root")
    pl$roots[[length(pl$roots) + 1]] <- ro
    pl
  }

############################################################

#' Compute the mean interbranch distance of all the primary roots in the image
#' @keywords rsml
#' @param obj of class plant
#' @param allroot if true, compute the interbanch distance on the whole root
#' @author Guillaume Lobet - guillaume.lobet(at)ulg.ac.be
#' @return the mean interbranch distance of the root system
#' @export
#' @examples
#' data(lupin)
#' meanInterbranchPlant(lupin)
meanInterbranchPlant = 
  function(obj, allroot=F) 
  {
    dens <- 0
    for(i in 1:nPrimRoot(obj)){
      dens <- dens + meanInterbranch(obj$root[[i]], allroot)
    }
    dens / nPrimRoot(obj)
  }

############################################################

#' Compute the mean insertion angle of all the laterals in the plant
#' @keywords rsml
#' @param obj of class plant
#' @author Guillaume Lobet - guillaume.lobet(at)ulg.ac.be
#' @return the mean insertion angle of the root system
#' @export
#' @examples
#' data(lupin)
#' meanInsertionAnglePlant(lupin)
meanInsertionAnglePlant = 
  function(obj) 
  {
    ang <- 0
    for(i in 1:nPrimRoot(obj)){
      ang <- ang + meanInsertionAngle(obj$root[[i]])
    }
    ang / nPrimRoot(obj)
  }


############################################################
#' Compute the total number of roots in the plant
#' @param obj of class plant
#' @keywords rsml
#' @author Guillaume Lobet - guillaume.lobet(at)ulg.ac.be
#' @return the number of root in the plant
#' @export
#' @examples
#' data(lupin)
#' nRoot(lupin)
nRoot = 
  function(obj) 
  {
    l <- length(obj$roots)
    for(i in 1:length(obj$roots)){
      l <- l + length(obj$root[[i]]$children)
      for(j in 1:length(obj$root[[i]]$children)){
        l <- l + length(obj$root[[i]]$children[[j]]$children)
        for(k in 1:length(obj$root[[i]]$children[[j]]$children)){
          l <- l + length(obj$root[[i]]$children[[j]]$children[[k]]$children)
        }
      }
    }
    l
  }


############################################################
#' Compute the number of primary roots in the plant
#' @param obj of class plant
#' @keywords rsml
#' @author Guillaume Lobet - guillaume.lobet(at)ulg.ac.be
#' @return the number of primary root in the plant
#' @export
#' @examples
#' data(lupin)
#' nPrimRoot(lupin)
nPrimRoot = 
  function(obj) 
  {
    length(obj$roots)
  }

############################################################

#' Compute the number of lateral roots in the plant
#' @param obj of class plant
#' @keywords rsml
#' @return the number of lateral root in the plant
#' @export
#' @examples
#' data(lupin)
#' nLatRoot(lupin)
nLatRoot = 
  function(obj) 
  {
    nRoot(obj) - nPrimRoot(obj)
  }


############################################################

#' Compute the length of the primary root based on the coordinates of its nodes
#' @param obj of class plant
#' @keywords rsml
#' @return the total length of the primary roots
#' @export
#' @examples
#' data(lupin)
#' primLength(lupin)
primLength = 
  function(obj)
  {
    l = 0
    for(j in 1:nPrimRoot(obj)){
      r <- obj$roots[[j]]
      l <- l + length(r)
    }
    l
  }

############################################################

#' Compute the length of the lateral root based on the coordinates of its nodes
#' @param obj of class plant
#' @keywords rsml
#' @return the total length of the lateral roots
#' @export
#' @examples
#' data(lupin)
#' latLength(lupin)
latLength = 
  function(obj)
  {
    length(obj) - primLength(obj)
  }
############################################################

#' Compute the length of the root based on the coordinates of its nodes
#' @param x object of class plant
#' @keywords rsml
#' @return the total length of the plant roots
#' @export
#' @examples
#' data(lupin)
#' length(lupin)
length.plant = 
  function(x)
  {
    l = 0
    for(i in 1:nPrimRoot(x)){
      l <- l + length(x$roots[[i]])
      for(j in 1:length(x$roots[[i]]$children)){
        l <- l + length(x$roots[[i]]$children[[j]])
        for(k in 1:length(x$roots[[i]]$children[[j]]$children)){
          l <- l + length(x$roots[[i]]$children[[j]]$children[[k]])
        }
      }
    }
    l
  }

############################################################

#' Compute the yrange of the plant
#' @param obj of class plant
#' @keywords rsml
#' @return c(y1,y2) where y1 and y2 are the y limits of the plant
#' @export
#' @examples
#' data(lupin)
#' yrangePlant(lupin)
yrangePlant = 
  function(obj)
  {
    ymin = 1e9;
    ymax = -1e9;
    for(i in 1:nPrimRoot(obj)){
      r <-  obj$roots[[i]]
      ymin <- min(c(ymin, -yrange(r)))
      ymax <- max(c(ymax, -yrange(r)))
      if(nChild(r) > 0){
        for(j in 1:nChild(r)){
          rr <-  r$children[[j]]
          ymin <- min(c(ymin, -yrange(rr)))
          ymax <- max(c(ymax, -yrange(rr)))
        } 
      }
    }   
    c(ymin, ymax)
  }

############################################################
#' Compute the xrange of the plant
#' @param obj of class plant
#' @keywords rsml
#' @return c(x1,x2) where x1 and x2 are the x limits of the plant
#' @export
#' @examples
#' data(lupin)
#' xrangePlant(lupin)
xrangePlant = 
  function(obj)
  {
    xmin = 1e9;
    xmax = -1e9;
    for(i in 1:nPrimRoot(obj)){
      r <-  obj$roots[[i]]
      xmin <- min(c(xmin, xrange(r)))
      xmax <- max(c(xmax, xrange(r)))
      if(nChild(r) > 0){
        for(j in 1:nChild(r)){
          rr <-  r$children[[j]]
          xmin <- min(c(xmin, xrange(rr)))
          xmax <- max(c(xmax, xrange(rr)))
        } 
      }
    }   
    c(xmin, xmax)
  }

############################################################
#' Compute the zrange of the plant
#' @param obj of class plant
#' @keywords rsml
#' @return c(z1,z2) where z1 and z2 are the z limits of the plant
#' @export
#' @examples
#' data(lupin)
#' zrangePlant(lupin)
zrangePlant = 
  function(obj)
  {
    zmin = 1e9;
    zmax = -1e9;
    for(i in 1:nPrimRoot(obj)){
      r <-  obj$roots[[i]]
      zmin <- min(c(zmin, zrange(r)))
      zmax <- max(c(zmax, zrange(r)))
      if(nChild(r) > 0){
        for(j in 1:nChild(r)){
          rr <-  r$children[[j]]
          zmin <- min(c(zmin, zrange(rr)))
          zmax <- max(c(zmax, zrange(rr)))
        } 
      }
    }   
    c(zmin, zmax)
  }
############################################################

#' Plot the root system
#' @param x object of class plant
#' @param threed make a 3D plot for the plant
#' @param ... plot options
#' @keywords rsml
#' @import rgl
#' @return null
#' @export
#' @examples
#' # Plot 2D plant
#' data(lupin)
#' plot(lupin, threed=FALSE)
#' 
#' # Plot 3D plant
#' data(anagallis)
#' plot(anagallis, threed=TRUE)
plot.plant = 
  function(x, threed = F, ...)
  {
    obj <- x
    if(threed){
      plot3d(1, 1, 1, type="n", xlim=xrangePlant(obj), ylim=yrangePlant(obj) , zlim=zrangePlant(obj), ylab="", xlab="")
      for(j in 1:nPrimRoot(obj)){
        r <- obj$roots[[j]]
        plot3d(coords(r)$x, coords(r)$y, coords(r)$z, lwd=2, col="red", type="l", add=T)
        if(nChild(r) > 0){
          for(i in 1:nChild(r)){
            rr <- r$children[[i]]
            plot3d(coords(rr)$x, coords(rr)$y, coords(rr)$z, lwd=2, col="green", type='l', add=T)
            if(nChild(rr) > 0){
              for(k in 1:nChild(rr)){
                rrr <- rr$children[[k]]
                plot3d(coords(rrr)$x, coords(rrr)$y, coords(rrr)$z, lwd=2, col="yellow", type='l', add=T)
              }
            }
          }
        }
      }
    } else {
      plot(1, 1, type="n", xlim=xrangePlant(obj), ylim=yrangePlant(obj), ylab="", xlab="")
      for(j in 1:nPrimRoot(obj)){
        r <- obj$roots[[j]]
        lines(coords(r)$x, -coords(r)$y, lwd=2, col="red")
        if(nChild(r) > 0){
          for(i in 1:nChild(r)){
            rr <- r$children[[i]]
            lines(coords(rr)$x, -coords(rr)$y, lwd=2, col="green")
            if(nChild(rr) > 0){
              for(k in 1:nChild(rr)){
                rrr <- rr$children[[k]]
                lines(coords(rrr)$x, -coords(rrr)$y, lwd=2, col="yellow")
              }
            }
          }
        }
      }
    } 
  }

############################################################
#' Summary of the plant
#' @param object object of class node
#' @param ... summary options
#' @author Guillaume Lobet - guillaume.lobet(at)ulg.ac.be
#' @export
#' @examples
#' data(lupin)
#' sum.lup <- summary(lupin)
#' sum.lup$total.length$value # Get total length
summary.plant = 
  function(object, ...)
  {
    obj <- object
    list(
      total.length = list(value = round(length(obj), 2), unit = "cm"),
      prim.length = list(value = round(primLength(obj), 2), unit = "cm"),      
      lat.length = list(value = round(latLength(obj), 2), unit = "cm"),
      n.root = list(value = nRoot(obj), unit = "-"),
      n.prim = list(value = nPrimRoot(obj), unit = "-"),
      n.lat = list(value = nLatRoot(obj), unit = "-"),      
      mean.insertion.angle = list(value = round(meanInsertionAnglePlant(obj), 2), unit = "degree"),
      mean.interbranch = list(value = round(meanInterbranchPlant(obj), 2), unit = "root/cm")
      )

  }

###########################################################
#' Print the plant
#' @param x  object of class node
#' @param ... print options
#' @author Guillaume Lobet - guillaume.lobet(at)ulg.ac.be
#' @export
#' @examples
#' data(lupin)
#' print(lupin)
print.plant = 
  function(x, ...)
  {
    obj <- x
    variable <- c("Total length", 
                  "Primary length", 
                  "Lateral length", 
                  "Number of roots", 
                  "Number of primary roots", 
                  "Number of lateral roots", 
                  "Mean insertion angle",
                  "Mean interbranch distance")
    value <- c(round(length(obj), 2),
               round(primLength(obj), 2),
               round(latLength(obj), 2),
               nRoot(obj),
               nPrimRoot(obj),
               nLatRoot(obj),
               round(meanInsertionAnglePlant(obj), 2),
               round(meanInterbranchPlant(obj), 2))
    units <- c("cm","cm","cm", "-", "-", "-", "degree", "root/cm")
    data.frame(variable, value, units)  
  }

