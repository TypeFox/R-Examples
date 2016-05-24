match.lnL.to.trees <- function(directory = 'getwd()', 
							   lnLprefix = 'RAxML_info.', lnLsuffix = '.lnL', 
							   treeIndexFile = 'tree.index.lines.txt', locus.names = NULL, ...) {
## updated 2014-01-23 to prune out files where no trees written 
  treeIndex <- read.delim(paste(directory, '/', treeIndexFile, sep = ''), as.is = T, header = F, row.names = 1)
  full.lnL <- try(get.raxml.treeLikelihoods(paste(directory, '/RAxML_info.fullMatrixOut.lnL', sep = ''))) # added 2014-02-13
  if(is.null(locus.names)) locus.names <- row.names(treeIndex)
  logfile = file(format(Sys.time(), "match.lnL.%Y-%m-%d.log.txt"), open = "a")
  lnL.list <- lapply(paste(directory, '/', lnLprefix, locus.names, lnLsuffix, sep = ''), get.raxml.treeLikelihoods)
  close(logfile)
  names(lnL.list) <- locus.names
  raxml.worked <- names(which(lnL.list != 'FAIL')) ## added 2014-01-23
  lnL.list <- lnL.list[raxml.worked] ## added 2014-01-23
  treeIndex <- treeIndex[raxml.worked, ] ## added 2014-01-23
  locus.names <- raxml.worked ## added 2014-01-23
  out.mat <- matrix(NA, nrow = length(locus.names), ncol = dim(treeIndex)[2], dimnames = list(locus.names, NULL))
  for(i in locus.names) {
	error.out = try(names(lnL.list[[i]]) <- try(unique(as.character(treeIndex[i,]))))
	if(class(error.out) == 'try-error') message(paste('assigning names to locus', i, 'failed'))
    else out.mat[i, ] <- lnL.list[[i]][as.character(treeIndex[i,])]
	}
  out.mat <- out.mat[-c(which(apply(out.mat, 1, function(x) any(is.na(x))))), ]
  if(class(full.lnL) != 'try-error') attr(out.mat, 'full.lnL') <- full.lnL
  class(out.mat) <- 'partitionedRAD'
  return(out.mat)
  }
