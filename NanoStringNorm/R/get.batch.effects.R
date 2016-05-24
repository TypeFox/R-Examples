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

get.batch.effects <- function(x, anno = NA, take.log, traits, sample.summary.stats, guess.cartridge = TRUE) {

	# convert the traits to matrix
	traits <- as.matrix(traits);
	if ( is.null(colnames(traits)) ) colnames(traits) <- paste('col', 1:ncol(traits), sep = '');

	traits <- traits[,!grepl('pair.ids',colnames(traits)),drop=FALSE];

	if (guess.cartridge) {
		# a best guess of how many cartridges are there.
		cartridge.max <- floor(ncol(x)/12) + 1;

		# create a design matrix that has dummy traits for each cartridge vs the rest.  do not create if there is already a supplied cartridge variable
		if ( all(!grepl('cartridge', colnames(traits),ignore.case=TRUE)) & cartridge.max > 1 ) {

			# create dummy matrix to hold the values
			cartridge.id <- rep(1:cartridge.max, each = 12)[1:ncol(x)];
			cartridge.design.matrix <- matrix(
				data = 1, 
				nrow = ncol(x), 
				ncol = cartridge.max,
				dimnames = list(colnames(x), paste('cartridge',1:cartridge.max,sep='')) 
				);

			# assign values to each cartridge variable
			cartridge.remove <- NULL;
			for (i in 1:cartridge.max) {
				if ( sum(cartridge.id == i) < 6 ) {
					cartridge.remove <- c(cartridge.remove, i);
					}
				cartridge.design.matrix[cartridge.id == i, i] <- 2;
				}

			# remove cartridge variables if less that 3 counts
			if ( !is.null(cartridge.remove) ) {
				cartridge.design.matrix <- cartridge.design.matrix[,-cartridge.remove];
				}

			# add the cartridge variable to the traits matrix 
			if ( all(dim(traits) == 1) ) {
				traits <- cartridge.design.matrix;
				}
			else {
				traits <- cbind(traits, cartridge.design.matrix);
				}
			}
		}

	# if no traits are supplied just exit
	if ( all(dim(traits) == 1) ) {
		return(data.frame(NA));
		}

	# use the probe summary values not the normalization factors which are ratio's
	sample.summary.stats <- sample.summary.stats[,!colnames(sample.summary.stats) %in% c('pos.norm.factor','sampleContent.norm.factor')];

	# deterine summary stats with no values
	sample.summary.stats.remove <- NULL;
	for (i in 1:ncol(sample.summary.stats) ) {
		if ( all(is.na(sample.summary.stats[,i])) ) {
			# add the column names to a character vector
			sample.summary.stats.remove <- c(sample.summary.stats.remove, i);
			}
		}

	# remove the columns that have no values
	if ( !is.null(sample.summary.stats.remove) ) {
		sample.summary.stats <- sample.summary.stats[,-sample.summary.stats.remove];
		}

	# initialize values
	batch.effects <- NULL;

	# loop over each trait and get the ttest pvalue and foldchange
	for ( i in 1:ncol(traits) ) {
		trait.value <- as.numeric(traits[,i]);
		trait.name  <- colnames(traits)[i];

		# get the mean and ttest of each summary stat grouped by trait
		batch.effect.mean.grp <- apply(
			X = sample.summary.stats, 
			MARGIN = 2,
			FUN = function(xval,trait.value) c(mean(xval[trait.value == 1 & !is.na(trait.value)]),mean(xval[trait.value == 2 & !is.na(trait.value)])),
			trait.value
			);

		# do a ttest comparing the distribution of the sample statistics for each of the trait groups
		 t.test.pvalue <- function(xval, trait, ...) {
		 	x <- xval[trait == 1];
			y <- xval[trait == 2];
			obj <- try(t.test(x, y, ...), silent=TRUE);
			if (is(obj, "try-error")) return(NA) else return(obj$p.value);
		 	}

		batch.effect.ttest <- apply(
			X = sample.summary.stats, 
			MARGIN = 2,
			FUN = t.test.pvalue,
			trait = trait.value
			);

		# combine the results
		batch.effect  <- cbind(trait.name, colnames(sample.summary.stats), signif(t(batch.effect.mean.grp),3), signif(batch.effect.ttest,3));
		batch.effects <- rbind(batch.effects, batch.effect);

		}

	# format the output
	rownames(batch.effects)  <- NULL;
	batch.effects <- data.frame(
		trait.name = batch.effects[,1],
		sample.statistics = batch.effects[,2],
		mean.grp1 = as.numeric(batch.effects[,3]),
		mean.grp2 = as.numeric(batch.effects[,4]),
		p.ttest = as.numeric(batch.effects[,5])
		);
	names(batch.effects) <- c("trait.name", "sample.statistics", "mean.grp1", "mean.grp2", "p.ttest");
	batch.effects$sample.statistics <- gsub('Sample.', '', batch.effects$sample.statistics, fixed = TRUE);
	batch.effects[batch.effects$sample.statistics == 'pos.controls', 'sample.statistics'] <- 'PositiveControls';
	batch.effects[batch.effects$sample.statistics == 'background.level', 'sample.statistics'] <- 'NegativeControls';
	batch.effects[batch.effects$sample.statistics == 'rna.content', 'sample.statistics'] <- 'RNA Content';

	return(batch.effects);
	} 


