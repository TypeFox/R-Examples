hypervolume <- function(data, repsperpoint=NULL, bandwidth, quantile=0.0, name=NULL, verbose=T, warnings=T)
{
  data <- as.data.frame(data)
  
  dim = ncol(data)
  np = nrow(data)
  
  if (is.null(repsperpoint))
  {
    repsperpoint = floor(500*ncol(data))
    cat(sprintf('Choosing repsperpoint=%.0f (use a larger value for more accuracy.)\n',repsperpoint))
  }
  
  if (nrow(data) == 0)
  {
    stop('Hypervolume cannot be computed with empty input data.')
  }
  
  if (any(is.na(data)))
  {
    stop('Hypervolume cannot be computed with missing values. Reduce dimensionality or remove observations with missing values.')
  }
  
  if (length(bandwidth) == 1)
  {
    bandwidth <- rep(bandwidth, ncol(data))
  }  
    
  if (length(bandwidth) != ncol(data))
  {
    stop('Input bandwidth vector not same length as dimensionality of dataset.')
  }
  
  if (any(bandwidth==0))
  {
    stop('Bandwidth must be non-zero.')
  }
  
  ### CHECK FOR WARNINGS
  if (warnings==TRUE)
  {
    corMatrix = cor(data)
    corMatrix[lower.tri(corMatrix,diag=TRUE)] <- NA
    correlationThreshold = 0.8
    corInd = which(abs(corMatrix) > correlationThreshold,arr.ind=T)
    finalstring = ""
    if (length(corInd) >= 1)
    {
      for (i in 1:nrow(corInd))
      {
        finalstring <- c(finalstring, sprintf('\n\tDimensions %s and %s are highly correlated (r=%.2f)', 
                                              names(data)[corInd[i,1]], names(data)[corInd[i,2]], corMatrix[corInd[i,1], corInd[i,2]]))
      }
      finalstring <- c(finalstring, "\nConsider removing some axes.")
      warning(finalstring)
    }
    
    if (log(np) <= dim)
    {
      warning(sprintf('Log number of observations (%.2f) is less than or equal to the number of dimensions (%d).\nResults may be inaccurate and sensitive to bandwidth choice.\nConsider reducing the dimensionality of the analysis.',
                      log(np), dim))
    }
    
    if (nrow(data) > 2)
    {
      sds = apply(data, 2, sd)
      sdsds <- sd(sds)
      if (!any(is.na(sdsds)))
      {
			if (sd(sds) > 5)
   			{
   	     		finalstring = 'Some dimensions have much more higher standard deviations than others:\n'
   	     		for (i in 1:length(sds))
       	 		{
        			finalstring <- c(finalstring, sprintf('\t%s %.2f\n',names(sds)[i], sds[i]))
        		}
        		warning(c(finalstring, 'Consider rescaling axes before analysis.'))
		 	}
		}
    }
    else
    {
      warning('Data contain two or fewer observations - cannot perform scaling checks. ')
    }
	}
  
  ### ERROR CHECKING FINISHED - BEGIN ALGORITHM

  
  # figure out the hypervolume of one kernel and the random point density within it
  hyperbox_volume = prod(2*bandwidth) # hyperbox is in both + and - dimensions
  point_density = repsperpoint / hyperbox_volume
 
  # generate random point cloud around each data point
  random_points = 2*(matrix(runif(repsperpoint*np*dim,min=0,max=1),nrow=repsperpoint*np, ncol=dim) - 0.5) * repmat(t(as.matrix(bandwidth)), repsperpoint*np, 1)
  offset = repmat(as.matrix(data), repsperpoint, 1)  
  data_points = random_points + offset 

  # determine the probability density at each random point
  if (verbose == TRUE)
  {
    cat('Evaluating probability density...\n')
  }
  point_counts = evalfrectangular(data, bandwidth, data_points,verbose=verbose)
  if (verbose == TRUE)
  {
    cat('Finished evaluating probability density.\n')
  }
  
  # infer the total volume based on the random point counts and the quantile
  # note that it may not be possible to achieve the exact quantile specified
  vc = volume_calculation(point_counts, point_density, quantile,verbose)
  
  disjunctfactor <- vc$final_volume / (nrow(as.matrix(data)) * prod(2*bandwidth))
  
  
  # see if the points are disjunct (if the total volume equals the summed volume of the hyperbox around each point)
  if (disjunctfactor > 0.9)
  {
    if (warnings==TRUE)
    {
      warning(sprintf('Disjunct factor: %.2f\nRatio of inferred volume to summed volume of hyperbox kernels around each data point is close to unity. Points are likely disjunct. Consider increasing the bandwidth.',disjunctfactor))
  
    }
  }
  # threshold the random points to include only those above the chosen quantile threshold
  point_counts_final = cbind(data_points, point_counts)[point_counts >= vc$index,]
  

  # downweight by the number of times the point intersect and standardize the point density
  weights = 1 / point_counts_final[,ncol(point_counts_final)]
  ow <- getOption('warn')
  options(warn=-1)
  weightedsample = sample(x=1:nrow(point_counts_final),size=floor(vc$final_volume * point_density),replace=T,prob=weights)
  options(warn=ow)
  # keep only unique points
  weightedsample = unique(weightedsample)

  # prepare object for output
  points_uniform_final = as.data.frame(point_counts_final[weightedsample,1:(ncol(point_counts_final)-1)])
  names(points_uniform_final) = names(data)
  density_uniform_final = point_counts_final[weightedsample,ncol(point_counts_final)]
  point_density_final = nrow(points_uniform_final) / vc$final_volume

  # clean up memory
  gc()
  

  result = new("Hypervolume", Name=ifelse(is.null(name), "untitled", toString(name)))
  result@Data = as.matrix(data)
  result@Dimensionality = dim
  result@Volume = vc$final_volume
  result@PointDensity = point_density_final
  result@Bandwidth = bandwidth
  result@RepsPerPoint = repsperpoint
  result@DisjunctFactor = disjunctfactor
  result@QuantileThresholdDesired = quantile
  result@QuantileThresholdObtained = vc$quantile_obtained
  result@RandomUniformPointsThresholded = as.matrix(points_uniform_final);  
  result@ProbabilityDensityAtRandomUniformPoints = density_uniform_final
  
  if (verbose==TRUE)
  {
    cat(sprintf('Quantile requested: %.2f   obtained: %.2f\n', result@QuantileThresholdDesired, result@QuantileThresholdObtained))
  }
  
  return(result)  

}
