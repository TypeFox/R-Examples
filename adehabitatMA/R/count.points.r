"count.points" <- function(xy, w)
  {
      ## Verifications
      if (is(w, "SpatialGrid"))
          fullgrid(w) = FALSE
      if (!inherits(w, "SpatialPixels"))
          stop("w should inherit the class SpatialPixels")
      if (!inherits(xy, "SpatialPoints"))
          stop("xy should inherit the class SpatialPoints")
      pfsx <- proj4string(w)
      pfsxy <- proj4string(xy)
      if (!identical(pfsx, pfsxy))
          stop("different proj4string in w and xy")

      gr <- gridparameters(w)
      if (nrow(gr) > 2)
          stop("w should be defined in two dimensions")
      if ((gr[1, 2] - gr[2, 2])> get(".adeoptions", envir=.adehabitatMAEnv)$epsilon)
          stop("the cellsize should be the same in x and y directions")
      if (ncol(coordinates(xy))>2)
          stop("xy should be defined in two dimensions")


      ## Prepares a vector containing the boundaries of the pixels
      meth <- "one"
      if (inherits(xy, "SpatialPointsDataFrame")) {
          if (ncol(xy)==1) {
              meth="sev"
          } else {
              meth="one"
              warning("several columns in the SpatialPointsDataFrame, no id considered")
          }
      }
      if (meth=="one") {
          ov <- over(xy, geometry(w))
          oo <- table(ov)
          repo <- rep(0, length(w[[1]]))
          repo[as.numeric(names(oo))] <- oo
          repo <- data.frame(x=repo)
          coordinates(repo) <- coordinates(w)
          gridded(repo) <- TRUE
          if (!is.na(pfsx))
              proj4string(repo) <- CRS(pfsx)
          return(repo)
      } else {
          id <- factor(xy[[1]])
          xy2 <- as.data.frame(coordinates(xy))
          lixy <- split(xy2, id)
          cp <- lapply(lixy, function(x) {
              count.points(SpatialPoints(x), w)
          })
          cp <- do.call("data.frame", lapply(cp, function(x) x[[1]]))
          coordinates(cp) <- coordinates(w)
          gridded(cp) <- TRUE
          if (!is.na(pfsx))
              proj4string(cp) <- CRS(pfsx)
          return(cp)
      }
  }

