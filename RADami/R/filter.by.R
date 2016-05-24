filter.by <-
function(dat, taxa, threshold = 'all') {
  ## returns just loci for which the requested taxa are present at some threshold
  ## default to returning 'all'
  if(!class(dat) %in% c('summary.pyRAD.loci', 'pyRAD.loci')) stop("This function only works with summary.pyRAD.loci datatypes")
  if(class(dat) == 'pyRAD.loci') dat.mat <- dat$radSummary$inds.mat[taxa, ]
    else dat.mat <- dat$inds.mat[taxa, ] # this is the case if you pass in a summary.pyRAD.loci object
  if(threshold == 'all') threshold <- length(taxa)
  return(names(which(apply(dat.mat, 2, sum) >= threshold)))
  }
