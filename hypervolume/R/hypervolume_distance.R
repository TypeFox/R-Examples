hypervolume_distance <- function(hv1, hv2, type="centroid", npmax=1000, check_memory=TRUE)
{
  hv1p <- hv1@RandomUniformPointsThresholded
  hv2p <- hv2@RandomUniformPointsThresholded
  
  if (type=="centroid")
  {
    hv1p_center <- colMeans(hv1p, na.rm=T)
    hv2p_center <- colMeans(hv2p, na.rm=T)
    
    centroid_distance <- sqrt(sum((hv1p_center - hv2p_center)^2))
    
    return(centroid_distance)
  }
  else if (type=="minimum")
  {    
    hv1p_ss <- hv1p[ sample(1:nrow(hv1p), min(npmax, nrow(hv1p)))  ,]
    hv2p_ss <- hv2p[ sample(1:nrow(hv2p), min(npmax, nrow(hv2p)))  ,]
    
    message(sprintf('Calculation will require %d pairwise distance calculations.',nrow(hv1p_ss)*nrow(hv2p_ss)))

    if (check_memory==TRUE)
    {
      message('Re-run with check_memory=FALSE to continue.')
      stop()
    }

    crossdistances <- fastPdist2(hv1p_ss, hv2p_ss)
    
    minimum_distance <- min(as.numeric(as.matrix(crossdistances)),na.rm=T)
    
    return(minimum_distance)
  }
  else
  {
    stop('Argument \'type\' takes unrecognized value.')
  }
}