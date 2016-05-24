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

other.normalization.vsn <- function(x, anno, OtherNorm = 'none', verbose = TRUE, genes.to.fit = NA, genes.to.predict = NA, ...) {

	# set the defaults to use all genes
	genes.to.fit = if(is.na(genes.to.fit)) "all" else genes.to.fit;
	genes.to.predict = if(is.na(genes.to.predict)) "all" else genes.to.predict;

	# do not proceed if vsn cannot be loaded
	if (!requireNamespace("vsn")) {
		stop("OtherNorm.vsn:  VSN is not available.");
		}

	# parse the genes.to.fit options
	if (any(genes.to.fit == 'all')) {
		genes.to.fit <- unique(anno$Code.Class);
		}
	else if (any(genes.to.fit == 'controls')) {
		genes.to.fit <- c('Housekeeping','Negative','Positive','Control');
		}	

	# parse the genes.to.predict options
	if (any(genes.to.predict == 'all')) {
		genes.to.predict <- unique(anno$Code.Class);
		}
	else if (any(genes.to.predict == 'controls')) {
		genes.to.predict <- c('Housekeeping','Negative','Positive','Control');
		}

	# which genes need to use
	genes.to.fit = anno$Name %in% genes.to.fit | grepl(paste(genes.to.fit, collapse = "|"), anno$Code.Class, ignore.case = TRUE) ;
	genes.to.predict = anno$Name %in% genes.to.predict | grepl(paste(genes.to.predict, collapse = "|"), anno$Code.Class, ignore.case = TRUE) ;

	# set values not to be fit to NA i.e. ignored
	x.fit <- x;
	x.fit[!genes.to.fit,] <- NA;

	x.predict <- x;
	x.predict[!genes.to.predict,] <- NA;

	# list of Code.Classes or list of genes to apply the method to.  
	if (verbose == TRUE) {
		#cat("OtherNorm.vsn: The following genes will used to fit the model:\n");
		#print(anno[genes.to.fit,c("Code.Class","Name")]);
		}

	#### START METHOD #########################################################################
		
	# vsn normalization
	if ( all(genes.to.fit == genes.to.predict) ) {
		# using common data
		# x.fit <- justvsn(x.fit);
		x.fit <- vsn::vsn2(x.fit, verbose = verbose, ...);
		x.predict <- predict(x.fit, newdata = x.predict, useDataInFit = TRUE);
		}
	else {
		# fit model on controls and predict on endogenous
		x.fit <- vsn::vsn2(x.fit, minDataPointsPerStratum = 16, verbose = verbose, ...);
		x.predict <- predict(x.fit, newdata = x.predict, useDataInFit = FALSE);
		}
	
	# extract the expression matrx
	#x.predict.exprs <- exprs(x.predict);

	### END METHOD ############################################################################	

	# add back the original counts to genes that should be ignored 
	x.predict <- data.frame(x.predict);
	
	if (any(!genes.to.predict)) {
		x.predict[!genes.to.predict, ] <- x[!genes.to.predict,];
		}

	rownames(x.predict) <- anno$Name;

	return(x.predict);
	}
