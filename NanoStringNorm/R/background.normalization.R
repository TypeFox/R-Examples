 #The NanoStringNorm package is copyright (c) 2012 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

background.normalization <- function(x, anno, Background = 'none', verbose = TRUE) {

	# check if missing
	if (is.na(Background)) {
		stop('Background: Background normalization method cannot be missing.  Try setting to *none*');		
		}

	# Background Correction
	if (Background != 'none') {

		# take the mean of negative controls
		if (Background == 'mean') {
			background.level <- apply(
				X = x[anno$Code.Class == 'Negative',,drop=FALSE], 
				MARGIN = 2, 
				FUN = mean, 
				na.rm = TRUE
				);
			}

		# take the mean of negative controls + 2 sd
		else if (Background == 'mean.2sd') {
			mean.plus.2sd <- function(y) mean(y, na.rm = TRUE) + 2 * sd(y, na.rm = TRUE);
			background.level <- apply(
				X = x[anno$Code.Class == 'Negative',,drop=FALSE], 
				MARGIN = 2, 
				FUN = mean.plus.2sd
				);
			}

		# take the maximum of negative controls
		else if (Background == 'max'){
			background.level <- apply(
				X = x[anno$Code.Class == 'Negative',,drop=FALSE], 
				MARGIN = 2, 
				FUN = max, 
				na.rm = TRUE
				);
			}

		# give an error if an unknown method is used
		else {
			stop('Background: Unimplemented Background method');
			}

		# flag any samples that are outliers
		background.level.sd.from.mean <- data.frame(background.zscore = (background.level - mean(background.level)) / sd(background.level));
		#colnames(background.level.sd.from.mean) <- 'background.zscore';

		if (verbose & any(abs(background.level.sd.from.mean) > 3)) {
			cat('Background: The following samples have an estimated background greater than \n\t 3 standard deviations from the mean.\n\n');
			print(signif(subset(background.level.sd.from.mean, abs(background.level.sd.from.mean) > 3),3));
			cat('\n');
			}

		# subtract backgound
		x <- t(apply(
			X = x, 
			MARGIN = 1, 
			FUN = '-', 
			background.level
			));

		# set negative values to zero
		x[x < 0] <- 0;

		# calculate count of missing values (DW missing could be moved to own function)
		genes.proportion.missing <- apply(
			X = x[grep('Endogenous', anno$Code.Class),,drop=FALSE],
			MARGIN = 1,
			FUN = function(y) { sum(y <= 0, na.rm = TRUE) / length(y) }
			);

		samples.proportion.missing <- apply(
			X = x[grep('Endogenous', anno$Code.Class),,drop=FALSE],
			MARGIN = 2,
			FUN = function(y) { sum(y <= 0, na.rm = TRUE) / length(y) }
			);

		samples.proportion.missing <- data.frame(row.names = colnames(x), proportion.missing = samples.proportion.missing);
		#colnames(samples.proportion.missing) <- 'proportion.missing';

		if (verbose == TRUE) {
			cat(paste('Background: After correction' , sum(samples.proportion.missing <= .90), 'samples and', sum(genes.proportion.missing <= .90), '\n\tEndogenous genes have less than 90% missing. \n\n'));
			}

		if (verbose & any(samples.proportion.missing > 0.9) ) {
			print(signif(subset(samples.proportion.missing, samples.proportion.missing > 0.9,),3));
			cat('\n');
			}

		}

	return(
		list(
			x = x,
			background.level = background.level
			)
		);
	}
