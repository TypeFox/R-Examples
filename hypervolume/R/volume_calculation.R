volume_calculation <- function(point_counts_raw, point_density, quantile, verbose=T)
{
  # count up the number 
  if (verbose==TRUE) {cat('Beginning volume calculation... ')}
  tabulated_counts <- cbind((1:max(point_counts_raw)), (tabulate(point_counts_raw)))
  if (verbose==TRUE) {cat('done. \n')}
  
  weighted_count = tabulated_counts[,2] / tabulated_counts[,1]

  cumulatesum = cumsum(weighted_count)/sum(weighted_count)

  qindex = head(which(cumulatesum >= quantile),n=1)
  qval = c(0,cumulatesum)[qindex]
  
  final_volume = sum(weighted_count[qindex:length(weighted_count)]) / point_density
  
  return(list(final_volume = final_volume, index_out = qindex, quantile_obtained = qval))
}
