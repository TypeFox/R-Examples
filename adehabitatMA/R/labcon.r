"labcon" <- function(x)
  {
      ## Verifications
      if (is(x, "SpatialGrid"))
          fullgrid(x) = FALSE
      if (!inherits(x, "SpatialPixelsDataFrame"))
          stop("should be an object of class SpatialPixelsDataFrame")
      pfs <- proj4string(x)
      gr <- gridparameters(x)
      if (nrow(gr) > 2)
          stop("x should be defined in two dimensions")
      if ((gr[1, 2] - gr[2, 2])> get(".adeoptions", envir=.adehabitatMAEnv)$epsilon)
          stop("the cellsize should be the same in x and y directions")
      maa<-as.image.SpatialGridDataFrame(x)
      y <- maa$z

      ## rajfond adds empty lines and columns on the borders of a map
      rajfond <- function(x) {
          nr <- nrow(x)
          nc <- ncol(x)
          f <- rep(0, nr)
          x <- cbind(f, x, f)
          f <- rep(0, nc + 2)
          x <- rbind(f, x, f)
      }

      ## The map is transformed so that it takes either
      ## the value 0 (NA) or 1 (mapped value)
      y[!is.na(y)] <- 1
      y[is.na(y)] <- 0
      y <- rajfond(y)

      ## sequential labelling of connex components
      ## with the C function "seqeticorr"
      toto <- .C("seqeticorr", as.double(t(y)), as.integer(nrow(y)),
                 as.integer(ncol(y)), PACKAGE="adehabitatMA")

      ## output
      etiquete <- matrix(toto[[1]], nrow = nrow(y), byrow = TRUE)
      ## and we delete the empty lines and columns added
      etiquete <- etiquete[-c(1, nrow(etiquete)), -c(1, ncol(etiquete))]
      etiquete[etiquete==0]<-NA
      maa$z <- etiquete
      grid <- image2Grid(maa)
      grid <- as(grid, "SpatialPixelsDataFrame")
      gridded(grid) <- TRUE

      if (!is.na(pfs))
          proj4string(grid) <- CRS(pfs)

      return(grid)
  }

