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

output.formatting <- function(x, anno, OtherNorm = 'none', round.values = FALSE, is.log = FALSE, take.log = FALSE, verbose = TRUE) {

	if ( (!is.logical(round.values) | !is.logical(take.log)) | is.na(round.values) | is.na(take.log) | !is.logical(is.log) | is.na(is.log) ) {
		stop('OUTPUT.FORMATTING: Round and take.log and is.log need to be TRUE or FALSE\n\n');
		}

	if ( is.log == TRUE & take.log == TRUE) {
		warning("OUTPUT.FORMATTING:  Are you sure you want to take the log of logged data?  Both is.log and take.log are TRUE.");
		}

	# shouldn't be log-transforming z-scores
	if (OtherNorm %in% c('rank.normal','zscore') & verbose == TRUE) {
		cat('OUTPUT.FORMATTING: log/rounding of zscores/rank.normal is intentionally not implemented.\n\n');
		}
	
	# round the data to discrete values.
	if (round.values == TRUE) {
		x <- round(x, digits = 0);
		}

	# log2 transformation
	if (take.log == TRUE & !OtherNorm %in% c('rank.normal','zscore')) {

		# handle negative values
		if (any(na.omit(x) < 1)) {

			if (verbose) {
				cat('log: Setting values less than 1 to 1 in order to calculate the log in positive space.\n\n');
				}

			# set all values less than one on the raw scale to 0 on the log scale i.e. undetected
			x[x < 1] <- 1;
			}

		# log-transform (base2 by default)
		x <- signif(log2(x), 4);
		}

	return(x);
	}
