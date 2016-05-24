hypervolume_importance <- function(hv)
{
  data <- hv@Data
  bandwidth <- hv@Bandwidth
  quantile <- hv@QuantileThresholdDesired
  reps <- hv@RepsPerPoint
  
  hv_others <- rep(NA, ncol(data))
  
  for (i in 1:length(hv_others) )
  {
    vars <- setdiff(1:ncol(data),i)

    hv_this <- hypervolume(data[,vars], repsperpoint=reps, bandwidth=bandwidth[vars],quantile=quantile,warnings=F,verbose=F)
    hv_others[i] <- hv_this@Volume
  }
  
  scaledvalues <- hv@Volume / hv_others
  
  names(scaledvalues) <- colnames(data)
  
  return(scaledvalues)
}
