"traj2spdf" <- function(tr)
{
    ## Verifications
    if (!inherits(tr, "traj"))
      stop("tr should be of class \"traj\"")

    ## Conversion to data frame
    class(tr) <- "data.frame"
    xy <- tr[,c("x","y")]
    tr$y <- tr$x <- NULL

    ## Conversion to SpatialPointsDataFrame
    res <- SpatialPointsDataFrame(xy, tr)

    ## Output
    return(res)
  }

