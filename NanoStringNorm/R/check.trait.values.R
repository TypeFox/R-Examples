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

check.trait.values <- function(x, anno = NA, traits = NA) {

	# check if NA
	if (all(is.na(traits))) {
		return(FALSE);
		}

	# attempt to convert traits into a matrix
	traits <- as.matrix(traits[,!colnames(traits) %in% 'pair.ids',drop = FALSE]);

	# if a vector add a dummy colname
	if ( is.null(colnames(traits)) ) {
		colnames(traits) <- paste('trait', 1:ncol(traits), sep = '');
		}

	# initialize variables
	check.na <- check.values <- check.samples <- NA;

	# check if input is a single value specifically if NA
	if ( all(dim(traits) == 1) ) {
		if (is.na(traits)) {
			traits.ok <- FALSE;
			}
		else {
			stop("Trait: Unrecognized trait input.");
			}
		}
	else {
		# check if the traits contain only 1 and 2
		if ( is.numeric(traits) & all(traits %in% c(NA,1,2)) ) {
			traits.ok <- TRUE;
			}
		else {
			stop("Trait: Only numeric variables with values NA, 1 or 2 are accepted.  The effect is terms of the second level i.e. disease.");
			}
		# check if traits have the right ids in the right order
		if ( nrow(traits) == ncol(x) ) {
			traits.ok <- TRUE;
			}
		else {
			stop("Trait: The number of traits is different form the number of samples.");
			}
		# check that that order of the sample and tratis matches	
		if (all(rownames(traits) == colnames(x)) ) {
			traits.ok <- TRUE;
			}
		else {
			stop("Trait: Values must include the same samples as the NanoString input data.  Confirm that traits have rownames in the same order as the columns in the expression data.")
			}

		# count the number of unique values per trait
		n.values <- apply(X = traits, MARGIN = 2, FUN = function(y) {length(unique(na.omit(y)))});

		if ( any(n.values == 1) ) {
			stop("Trait: Some of your traits only have one value.");
			}
		else {
			traits.ok <- TRUE;
			}

		}

	return(traits.ok)
	}
