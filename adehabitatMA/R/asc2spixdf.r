"asc2spixdf" <- function(a)
  {
      ## Verifications
      if (!inherits(a, "asc"))
          stop("a should be of class \"asc\"")

      ## creates the data frame of coordinates
      xyc <- .getXYcoords(a)
      xc <- rep(xyc$x, times=length(xyc$y))
      yc <- rep(xyc$y, each=length(xyc$x))
      xyc<-data.frame(x=xc,y=yc)

      ## keep only the mapped areas for the variable
      cons <- (1:length(c(a)))[!is.na(c(a))]
      var <- c(a)[cons]
      xyc <- xyc[cons,]
      names(xyc) <- c("x","y")

      ## created the spatial data frame
      df1 <- data.frame(xyc, var)
      coordinates(df1) <- c("x","y")
      gridded(df1) <- TRUE

      ## Output
      return(df1)
  }

