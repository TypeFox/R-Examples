"kasc2spixdf" <- function(ka)
{
    ## Verifications
    if (!inherits(ka, "kasc"))
      stop("ka should be of class \"kasc\"")

    ## gets the coordinates of the kasc
    xyc <- .getXYcoords(ka)
    xc <- rep(xyc$x, times=length(xyc$y))
    yc <- rep(xyc$y, each=length(xyc$x))
    xyc<-data.frame(x=xc,y=yc)

    ## remove the NA
    ka <- .managNAkasc(ka)
    cons <- (1:nrow(ka))[!is.na(ka[,1])]
    df <- ka[cons,]

    ## The spatial data frame
    class(df) <- "data.frame"
    df <- as.data.frame(lapply(df, function(x) as.numeric(x)))
    xyc <- xyc[cons,]
    names(xyc) <- c("x","y")
    df1 <- data.frame(xyc, df)
    coordinates(df1) <- c("x","y")
    gridded(df1) <- TRUE
    return(df1)
  }

