# The NanoStringNorm package is copyright (c) 2012 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

other.normalization.quantile <- function(x, anno, OtherNorm = 'none', verbose = TRUE, genes.to.fit = NA, genes.to.predict = NA) { 

	genes.to.fit <- if(is.na(genes.to.fit)) "all" else genes.to.fit;
	genes.to.predict <- if(is.na(genes.to.predict)) "all" else genes.to.predict;
	
	# parse the genes.to.fit options
	if (any(genes.to.fit == 'all')) {
		genes.to.fit <- unique(anno$Code.Class);
		}
	if (any(genes.to.fit == 'controls')) {
		genes.to.fit <- c('Housekeeping','Negative','Positive','Control');
		}

	# which genes need to be fit
	genes.to.fit = anno$Name %in% genes.to.fit | grepl(paste(genes.to.fit, collapse = "|"), anno$Code.Class, ignore.case = TRUE) ;
	
	# set values not to be fit to NA i.e. ignored
	x.fit <- x;
	x.fit[!genes.to.fit,] <- NA;

	# zeros are considered missing so ignored for empirical distribution generation and application
	x.fit[x.fit == 0] <- NA;

	# list of Code.Classes or list of genes to apply the method to.  
	if (verbose == TRUE) {
		if (any(!genes.to.fit)) {
			cat("OtherNorm.quantile: The following genes will not be processed:\n");
			print(anno[!genes.to.fit,c("Code.Class","Name")]);
			}
		}

	# sort each sample independently into a new dataset
	x.sort <- apply(
		X = x.fit, 
		MARGIN = 2,
		FUN = sort,
		na.last = TRUE
		);

	# take median across samples for each rank and use this as an empirical distribution.
	empirical.distribution <- apply(
		X = x.sort, 
		MARGIN = 1,
		FUN = median,
		na.rm = TRUE
		);

	# get the ranks in the unsorted data
	x.rank <- apply(
		X = x.fit, 
		MARGIN = 2,
		FUN = rank,
		na.last = TRUE,
		ties = 'first'
		);

	# set ranks to NA if NA in the original data
	x.rank[is.na(x.fit)] <- NA;

	# convert the ranks into probabilities to use quantile function
	get.prob.from.rank <- function(xval) ( xval - 1 ) / ( length(na.omit(xval)) - 1 );

	x.rank.prob <- apply(
		X = x.rank,
		MARGIN = 2,
		get.prob.from.rank
		);

	# use the empirical distribution to peg the ranks too
	get.quantile <- function(y, dist) { quantile(x = dist, probs = y, na.rm = TRUE) }

	# only run the quantile function on the endogenous but add back the controls
	x.fit <- apply(
		X = x.rank.prob,
		MARGIN = 2, 
		FUN = get.quantile, 
		dist = empirical.distribution
		);

	# NA's back to 0
	x.fit[is.na(x.fit)] <- 0;

	# add back the original counts to genes that should be ignored 
	x.fit <- data.frame(x.fit);
	
	if (any(!genes.to.fit)) {
		x.fit[!genes.to.fit, ] <- x[!genes.to.fit,];
		}
	
	rownames(x.fit) <- anno$Name;
	return(x.fit);
	}
