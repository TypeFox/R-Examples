simmr_load =
function(mixtures, source_names, source_means, source_sds, correction_means=NULL, correction_sds=NULL, concentration_means=NULL, group=NULL) {

# Function to load in data for simmr and check whether it's appropriate for running through simmr_mcmc

# Go through each object and check that it matches the requirements

# Mixtures must be a matrix - the number of rows is the number of observations and the number of columns is the number of tracers
if(!is.matrix(mixtures)) stop("mixtures object must be a matrix")
n_obs = nrow(mixtures)
n_tracers = ncol(mixtures)

# Add column names if they're not there
if(is.null(colnames(mixtures))) colnames(mixtures) = paste0('tracer',1:n_tracers)

# source_names must be a character vector - the length of it is the number of sources
if(!is.vector(source_names)) stop("source_names must be a vector")
source_names = as.character(source_names)
n_sources = length(source_names)

# source_means and source_sds must both be matrices where the number of rows is n_sources (in the same order as source_names) and the number of columns is n_tracers
if(length(dim(source_means))!=2) stop("source_means must have two dimensions")
if(length(dim(source_sds))!=2) stop("source_sds must have two dimensions")
if(nrow(source_means)!=n_sources) stop('Number of rows in source_means does not match length(source_names)')
if(ncol(source_means)!=n_tracers) stop('Number of columns in source_means does not match ncol(mixtures)')
if(nrow(source_sds)!=n_sources) stop('Number of rows in source_sds does not match length(source_names)')
if(ncol(source_sds)!=n_tracers) stop('Number of columns in source_sds does not match ncol(mixtures)')

# Check that either neither or both of corrections_means and sds are given
if(is.null(correction_means) & !is.null(correction_sds)) stop("Both correction_means and correction_sds must be supplied")
if(!is.null(correction_means) & is.null(correction_sds)) stop("Both correction_means and correction_sds must be supplied")

# correction_means and correction_sds must be matrix and of same dimension as n_sources
if(!is.null(correction_means)) {
  if(length(dim(correction_means))!=2) stop("correction_means must have two dimensions")
  if(length(dim(correction_sds))!=2) stop("correction_sds must have two dimensions")
  if(nrow(correction_means)!=n_sources) stop('Number of rows in correction_means does not match length(source_names)')
  if(ncol(correction_means)!=n_tracers) stop('Number of columns in correction_means does not match ncol(mixtures)')
  if(nrow(correction_sds)!=n_sources) stop('Number of rows in correction_sds does not match length(source_names)')
  if(ncol(correction_sds)!=n_tracers) stop('Number of columns in correction_sds does not match ncol(mixtures)')
} else {
  correction_means = matrix(0,ncol=n_tracers,nrow=n_sources)
  correction_sds = matrix(0,ncol=n_tracers,nrow=n_sources)
}

# concentration_means must be a matrix where all elements are less than 1
if(!is.null(concentration_means)) {
  if(length(dim(concentration_means))!=2) stop("concentration_means must have two dimensions")
  if(nrow(concentration_means)!=n_sources) stop('Number of rows in concentration_means does not match length(source_names)')
  if(ncol(concentration_means)!=n_tracers) stop('Number of columns in concentration_means does not match ncol(mixtures)')
  if(any(concentration_means>1)) stop("Invalid values in concentration_means; must all be less than 1")
} else {
  concentration_means = matrix(1,ncol=n_tracers,nrow=n_sources)
}

# Check the groups are the right length and structure if given
if(!is.null(group)) {
  if(!is.integer(group)) stop("group variable needs to be of integer type. Perhaps use as.integer?")
  if(min(group)!=1) stop("Group number must start at 1")
  if(!all(diff(sort(unique(group)))==1)) stop("Group numbers must proceed sequentially from 1, 2, 3, etc")
  if(length(group)!=n_obs) stop("Number of group values not equal to number of observations")
} else {
  group = rep(1,n_obs)
}
n_groups = length(unique(group))

# Prepare output and give class
out = list(mixtures=mixtures, source_names=source_names, source_means=source_means, source_sds=source_sds, correction_means=correction_means, correction_sds=correction_sds, concentration_means=concentration_means, group=group, n_obs=n_obs, n_tracers=n_tracers, n_sources=n_sources, n_groups=n_groups)

# Look through to see whether there are any missing values in anything
if(any(unlist(lapply(out,'is.na')))) stop("Missing values provided for some values. Check your inputs")

class(out) = 'simmr_input'

return(out)

}
