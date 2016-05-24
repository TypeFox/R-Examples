# Distribute rainfall (pe) input for this time step and across HRU groupings 
# ==============================================================================
# Inputs
# ==============================================================================
# rain: rainfall or pe in time step
# ==============================================================================
# returns
# updated input flows at time step
# adding recharge to river hsus 
# ==============================================================================

allocate.PE <- function(pe, groups, it)
{
	# currently as for rainfall 
	# partition rain equally across groups 
	return(allocate.rain(pe, groups, it))  
}

is.datetime <- function(obj)
{
	return(inherits(obj, "POSIXt"))
	
}

# return a vector with the correct rainfall allocated from multiple gauges to appropriate area,
# identified by gauge.id 
# this needs optimisation: the rainfall vectors could be allocated at initialisation, one for each
# HRU (which might lead to them being duplicated...) and the index automatically 
allocate.rain <- function(rain, groups, it)
{
	# can index by time instead of integer (performance implications?)
	if(is.datetime(it)){it <- which(index(rain)==it)}
# 	if(all(colnames(rain)==groups$id))
# 	{
# 		# already split the rainfall input up by response unit
# 		rain.dist <- rain[it,]   
# 	}
# 	else
# 	{	
# 	
# 	  if(it>nrow(rain))
# 	  {
# 	    warning("Rainfall record truncated at ", tail(index(rain), 1))
# 	    it <- nrow(rain)
# 	  }
# 	
	#ncol <- ncol(rain)
	# gauge id, limited to nunber of data columns	
	rain.dist <- rain[it,pmin(groups$gauge.id, ncol(rain))]   
	
  return(as.vector(rain.dist))
}

