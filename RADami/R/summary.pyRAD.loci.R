summary.pyRAD.loci <-
function(object, ...) {
# Arguments:
#  object = a pyRAD.loci object
#  var.only = if T, only includes variable loci; as written, the function assumes a "*" if there is
#  REVISED 2014-02-04, so var.only option is excluded; doesn't make much sense in a summary method 
  message("\nDoing a pyRAD summary")
  reportInterval <- 2000 # this is just for screen reporting... only matters with really long files
  ## currently, locus.names includes a null (""), b/c the break lines have no locus name
  locus.names <- as.character(unique(object$locus.index)) # this slows things down by a factor of 2 or 3, but it seems to prevent a subscript-out-of-bounds error
  locus.names <- locus.names[locus.names != ""]
  ## REWRITE TO LOOK FOR BREAKLINES THAT HAVE * OR - IN THEM
  variable.loci <- locus.names[!is.na(object$seqs[object$breaks])] # note that this works with the pyRAD output we are currently using... should be checked
  # if(var.only) locus.names <- variable.loci
  num.loci <- length(locus.names)
  tip.names <- as.character(unique(object$tips[-c(object$breaks, object$cons)]))
  message("Splitting tips by locus name...")
  tips.per.locus <- split(object$tips, object$locus.index)[locus.names]
  seqs.per.locus <- split(object$seqs, object$locus.index)[locus.names]
  num.inds.per.locus <- sapply(tips.per.locus, length)
  inds.mat <- matrix(NA, nrow = length(tip.names), ncol = num.loci, dimnames = list(tip.names, locus.names))
  message("Making tips matrix...")
  start.time <- Sys.time()
  ## is there some way to vectorize the following:
  for(i in seq(num.loci)) {
	temp <- try(tip.names %in% tips.per.locus[[locus.names[i]]])
	if(class(temp) != "try-error") {
	  inds.mat[ , locus.names[i]] <- temp
	  names(seqs.per.locus[[i]]) <- tips.per.locus[[i]]
	  }
	else(message(paste("Error occurred with locus", locus.names[i])))
    if(i / reportInterval - i %/% reportInterval == 0) {
  	   message(paste('...', i, 'of', num.loci, 
 	   '-- Estimated time remaining =', ((Sys.time() - start.time) / i) * (num.loci - i), attr(Sys.time() - start.time, 'units')
  	   ))
	   }
	 }

  out <- list(num.loci = num.loci, tips.per.locus = tips.per.locus, break.vectors = object$break.vectors, seqs.per.locus = seqs.per.locus, num.inds.per.locus = num.inds.per.locus, variable.loci = variable.loci, inds.mat = inds.mat, locus.lengths = lengths.report(object, 0))
  class(out) <- 'summary.pyRAD.loci'
  out
  }
