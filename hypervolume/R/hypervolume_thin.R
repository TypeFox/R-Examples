hypervolume_thin <- function(hv, factor=NULL, npoints=NULL)
{
  rp <- hv@RandomUniformPointsThresholded
  nrp <- nrow(rp)  
  
  if (!is.null(factor))
  {
    if (factor <= 0 | factor >=1)
    {
      stop("Thinning factor must be in (0,1).")
    }
  }
  else if (!is.null(npoints))
  {
    if (npoints < 1)
    {
      stop("Number of points must be greater than zero.")
    }

    if (!is.null(factor))
    {
      stop("Cannot specify both factor and npoints.")
    }
    else
    {
      # make sure we don't take more points than exist in the dataset
      npoints <- min(nrp, npoints)
      
      # recast in terms of a factor
      factor <- npoints / nrp
    }
  }
  else
  {
    stop("Must specify either factor or npoints.")
  }
  
  hv_out <- hv
  
  hv_out@RandomUniformPointsThresholded <- rp[sample(1:nrow(rp),nrow(rp)*factor),]
  
  hv_out@PointDensity <- hv@PointDensity * factor
  
  hv_out@RepsPerPoint <- hv@RepsPerPoint * factor
  
  return(hv_out)
}