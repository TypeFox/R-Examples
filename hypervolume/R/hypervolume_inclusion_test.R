hypervolume_inclusion_test <- function(hv, points, reduction_factor=1, verbose=T, distance_factor=1.0)
{  
  np = nrow(hv@RandomUniformPointsThresholded)
  dimhv = ncol(hv@RandomUniformPointsThresholded)
  dimp = ncol(points)
  
  if (dimp != dimhv)
  {
    stop('Dimensionality of hypervolume and points is not the same.')
  }
  
  if (reduction_factor <= 0 | reduction_factor > 1)
  {
    stop('Reduction factor is not in (0,1].')
  }

  # now pick a uniformly random subset of these points
  # assuming that the set of points is already uniformly random
  numpointstokeep_hv = floor(np * reduction_factor)        
  if (reduction_factor < 1)
  {
    hv_points_ss = hv@RandomUniformPointsThresholded[sample(1:np,size=numpointstokeep_hv),]
  }
  else
  {
    hv_points_ss = hv@RandomUniformPointsThresholded
  }
  
  if (verbose==TRUE)
  {
    cat(sprintf('Retaining %d points for %d inclusion tests.\n', numpointstokeep_hv, nrow(points)))
  }
  
  # determine the reduced hypervolume's point density
  point_density = nrow(hv_points_ss) / hv@Volume
  
  # calculate characteristic distances
  cutoff_dist = point_density^(-1/dimhv) * distance_factor
  
  # figure out which points are 'close enough' to other points
  # (within a n-ball of the critical distance)
  points_in_hv_all = evalfspherical(hv_points_ss, cutoff_dist, points, verbose=verbose)
  
  # flag those points and return them
  points_in = points_in_hv_all > 0
  
  # clean up memory
  gc()
  
  return(points_in)
}