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

predict.concentration <- function(x, anno, is.log, take.log, verbose = TRUE) {

	# concentration in fM
	if ( take.log == TRUE | is.log == TRUE) { 
		positive.controls.concentration <- c(.125, .5, 2, 8, 32, 128);
		}
	else {
		positive.controls.concentration <- log2(c(.125, .5, 2, 8, 32, 128));
		}

	# generate average model based on mean of postive controls
	endogenous.gene.mean <- apply(
		X = x[grepl('Endogenous', anno$Code.Class),],
		MARGIN = 1,
		FUN = mean,
		na.rm = TRUE
		);

	positive.control.mean <- apply(
		X = x[anno$Code.Class == 'Positive', ],
		MARGIN = 1,
		FUN = mean,
		na.rm = TRUE
		);

	fm1.intercept <- lm(positive.control.mean ~ positive.controls.concentration)$coefficients[1];
	fm1.slope     <- lm(positive.control.mean ~ positive.controls.concentration)$coefficients[2];

	mean.endogenous.gene.concentration <- (endogenous.gene.mean - fm1.intercept) / fm1.slope;

	# generate individual specific model of postive controls
	all.endogenous.gene.concentration <- NULL;
	for ( i in 1:ncol(x) ) {
		sample <- x[,i];
		fm2.intercept                     <- lm(sample[anno$Code.Class == 'Positive'] ~ positive.controls.concentration)$coefficients[1];
		fm2.slope                         <- lm(sample[anno$Code.Class == 'Positive'] ~ positive.controls.concentration)$coefficients[2];
		endogenous.gene.concentration     <- (sample[grepl('Endogenous', anno$Code.Class)] - fm2.intercept) / fm2.slope;
		all.endogenous.gene.concentration <- cbind(all.endogenous.gene.concentration, endogenous.gene.concentration);
		}

	colnames(all.endogenous.gene.concentration) <- colnames(x);

	gene.concentration <- data.frame(
		row.names = anno[grepl('Endogenous', anno$Code.Class),'Name'], 
		mean.concentration = round(mean.endogenous.gene.concentration, 3), 
		round(all.endogenous.gene.concentration, 3),
		stringsAsFactors = FALSE
		)

	return(gene.concentration);
	}
