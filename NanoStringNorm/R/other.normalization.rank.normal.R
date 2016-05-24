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

other.normalization.rank.normal <- function(x, anno, OtherNorm = 'none', verbose = TRUE, genes.to.fit = NA) { 

	genes.to.fit <- if (all(is.na(genes.to.fit))) { "all"; } else { genes.to.fit; }

	# parse the genes.to.fit options
	if (any(genes.to.fit == 'all')) {
		genes.to.fit <- unique(anno$Code.Class);
		}

	if (any(genes.to.fit == 'controls')) {
		genes.to.fit <- c('Housekeeping','Negative','Positive','Control');
		}

	# set the genes not to be normalized to NA
	genes.to.fit <- anno$Name %in% genes.to.fit | grepl(paste(genes.to.fit, collapse = "|"), anno$Code.Class, ignore.case = TRUE);

	# list of Code.Classes or list of genes to apply the method to.  
	if (verbose == TRUE) {
		cat("OtherNorm.rank.normal: Genes to Fit:\n");
		print(anno[genes.to.fit, 'Name']);
		}

	#set values not to be fit to NA i.e. ignored
	x.fit <- x;
	x.fit[!genes.to.fit,] <- NA;

	# list of Code.Classes or list of genes to apply the method to.  
	if (verbose == TRUE) {
		if (any(!genes.to.fit)) {
			cat("OtherNorm.rank.normal: The following genes will not be processed:\n");
			print(anno[!genes.to.fit,c("Code.Class","Name")]);
			}
		}

	# rank (forced normal) normalization
	x.fit <- apply(
		X = x.fit,
		MARGIN = 2,
		FUN = function(y) qnorm (p = (rank(y, na.last = 'keep') -.5)/length(y))
		);

	# add back the original counts to genes that should be ignored 
	x.fit <- data.frame(x.fit);

	if (any(!genes.to.fit)) {
		x.fit[!genes.to.fit,] <- x[!genes.to.fit,];
		}

	rownames(x.fit) <- anno$Name;

	return(x.fit);
	}
