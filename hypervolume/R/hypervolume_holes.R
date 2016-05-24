hypervolume_holes <- function(hv_obs, hv_exp, set_npoints_max=NULL, set_check_memory=TRUE)
{		
  # initialize result
  finalresult <- NULL
  
  if (is.null(set_npoints_max))
  {
    set_npoints_max = floor(100*10^sqrt(hv_obs@Dimensionality))
    cat(sprintf('Choosing set_npoints_max=%.0f (choose a larger value for more accuracy.)\n',set_npoints_max))    
  }
  
  # make sure we are using a hypervolume with data (i.e. not the output of set operations)
  if(all(is.nan(hv_obs@Data)))
  {
    stop('Hypervolume must be associated with datapoints')
  }
  
  if (hv_obs@Dimensionality != hv_exp@Dimensionality)
  {
    stop('Observed and expected hypervolumes must have same dimensionality.')
  }
  
  # FIND THE DIFFERENCE between the convex hull shape and the real hypervolume
  cat("Beginning set operations (resampling to minimum density)...")
  hvs_overlap <- hypervolume_set(hv_obs, hv_exp, check_memory=set_check_memory, npoints_max=set_npoints_max)
  if (set_check_memory)
  {
    stop('Set set_check_memory=F to continue.\n')
  }
  cat("Finished set operations.\n")
  
  # find the distance between all points in the difference
  randompoints <- hvs_overlap@HVList$Unique_2@RandomUniformPointsThresholded
  if (is.null(randompoints) || nrow(randompoints) == 0)
  {		
    cat('No holes found.\n');
  }
  else
  {  	
    cat(sprintf("Retaining %.0f random points in set difference.\n", nrow(randompoints)))

    criticaldistance <- hvs_overlap@HVList$Unique_2@PointDensity ^(-1/hvs_overlap@HVList$Unique_2@Dimensionality)    
    
    # find points with minimum neighbor distance less than threshold
    distances <- as.matrix(dist(randompoints, method="euclidean"))
    diag(distances) <- NA
    isin <- (apply(distances, 1, min, na.rm=T) < criticaldistance)
    
    cat(sprintf("Removing %d stray points...\n", length(which(isin==0))))
    
    randompoints_trimmed <- randompoints[isin,]
    
    thishv <- hvs_overlap@HVList$Unique_2 # copy base information
    thishv@RandomUniformPointsThresholded <- randompoints_trimmed
    thishv@Volume <- hvs_overlap@HVList$Unique_2@Volume * nrow(randompoints_trimmed) / nrow(randompoints)
    thishv@Name <- sprintf("Hole in %s relative to %s", hv_obs@Name, hv_exp@Name)
    
    finalresult <- thishv

    # return the final hypervolumelist
    cat('Returning all holes.\n');
    
    if (is.null(finalresult))
    {
      cat('No holes found. Function will return NULL.\n')
    }
        
  }
  
  return(finalresult)
}