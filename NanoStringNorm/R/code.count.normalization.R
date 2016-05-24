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

code.count.normalization <- function(x, anno, CodeCount = 'none', logged = FALSE, verbose = TRUE) {

	# check if missing
	if (is.na(CodeCount)) {
		stop('CodeCount: CodeCount normalization method cannot be missing.  Try setting to *none*');
		}

	# Code Count Normalization using sum: take sum of positive controls per sample
	if ( CodeCount == 'sum' ) { 
		pos.sample <- apply(
			X = x[anno$Code.Class == 'Positive',,drop=FALSE], 
			MARGIN = 2, 
			FUN = sum
			);
		}
	
	# Code Count Normalization using geometric mean
	else if ( CodeCount == 'geo.mean' ) { 
		pos.sample <- apply(
			X = x[anno$Code.Class == 'Positive',,drop=FALSE],
			MARGIN = 2, 
			FUN = get.geo.mean,
			logged = logged
			);
		}
		
	else {
		stop('CodeCount: Unimplemented CodeCount normalization method');
		}
	
	# calculate the normalization factor as the ratio of mean vs sample values
	pos.norm.factor <- mean(pos.sample) / pos.sample;

	# warning if norm factor is outside expected range
	if ( any(pos.norm.factor > 3) | any(pos.norm.factor < 0.3) ) {
	
		pos.norm.factor.flagged.samples <- as.data.frame(pos.norm.factor[(pos.norm.factor > 3) | (pos.norm.factor < 0.3)]);
		colnames(pos.norm.factor.flagged.samples) <- 'pos.norm.factor';

		if (verbose) {
			cat('CodeCount: The following samples have positive normalization factors outside the \n\t recommended range of (0.3 to 3).  Consider removing them.\n\n');
			print(signif(pos.norm.factor.flagged.samples,3));
			cat('\n');
			}
		}

	# multiply normalization factor to data values.  output needs to be transposed to get correct orientation
	x <- t(
		apply(
			X = x, 
			MARGIN = 1, 
			FUN = '*', 
			pos.norm.factor
			)
		);

	return(
		list(
			x = x,
			pos.norm.factor = pos.norm.factor,
			pos.sample = pos.sample
			)
		);
	}
