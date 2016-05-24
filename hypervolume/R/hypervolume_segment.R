hypervolume_segment <- function(hv,distancefactor=hv@Dimensionality, npmax=NULL)
{
  if(class(hv) != "Hypervolume")
  {
    stop("Input must be a Hypervolume class object")
  }
  
  # thin the input
  if (!is.null(npmax))
  {
    hv <- hypervolume_thin(hv, npoints=npmax)
  }
  
  hvrp <- hv@RandomUniformPointsThresholded
  
  characteristicdistance <- (1/hv@PointDensity)^(1/hv@Dimensionality)
  
  hc <- fastcluster::hclust.vector(hvrp, method='single')
  
  membership <- cutree(hc, h = characteristicdistance * distancefactor)
  
  ngroups <- max(membership)
  
  hvs <- vector("list",ngroups)
  
  for (i in 1:ngroups)
  {
    hv_temp <- hv
    hv_temp@RandomUniformPointsThresholded <- hvrp[membership==i,,drop=FALSE]
    hv_temp@ProbabilityDensityAtRandomUniformPoints <- hv@ProbabilityDensityAtRandomUniformPoints[membership==i]
    hv_temp@Volume <- hv@Volume * nrow(hv_temp@RandomUniformPointsThresholded) / nrow(hv@RandomUniformPointsThresholded)
    
    hvs[[i]] <- hv_temp
  }
  
  hvs_segmented <- do.call("hypervolume_join",hvs)
  
  return(hvs_segmented)
}