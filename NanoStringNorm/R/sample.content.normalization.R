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

sample.content.normalization <- function(x, anno, SampleContent = 'none', logged = FALSE, verbose = TRUE) {

	# check if missing
	if (is.na(SampleContent)) {
		stop('SampleContent: SampleContent normalization method cannot be missing.  Try setting to *none*');
		}

	# Sample Content Normalization
	if (SampleContent != 'none') {

		# get list of endogenous probes excluding viral genes from normalization factors due to potential associatio with phenotype
		endogenous.genes <- grepl('Endogenous', anno$Code.Class) & !grepl('bkv-|ebv-|hcmv-|hiv1-|hsv1-|hsv2-|kshv-|mcv-',anno$Name);

		# Check for Endogenous: take sum of housekeeping genes.
		if (SampleContent == 'housekeeping.sum') {
			rna.content <- apply(
				X = x[anno$Code.Class %in% c('Control', 'Housekeeping', 'housekeeping'),,drop=FALSE ], 
				MARGIN = 2, 
				FUN = sum
				);
			}
	
		else if (SampleContent == 'housekeeping.geo.mean') {
			
			# take the geometric mean of the housekeeping genes.
			rna.content <- apply(
				X = x[anno$Code.Class %in% c('Control', 'Housekeeping', 'housekeeping'),,drop=FALSE ], 
				MARGIN = 2, 
				FUN = get.geo.mean,
				logged = logged
				);
			}

		# take mean of all counts
		else if (SampleContent == 'total.sum') {
			rna.content <- apply(
				X = x[endogenous.genes,,drop=FALSE],
				MARGIN = 2, 
				FUN = sum
				);
			}

		# take geometric mean of any endogenous or HK genes that seem stable
		else if (SampleContent == 'low.cv.geo.mean') {

		# get a list of expressed stable genes
			gene.summary.stats  <- as.data.frame(get.gene.summary.stats(x, anno));
			low.cv.genes <- (
				grepl('endogenous|housekeeping|control', anno$Code.Class, ignore.case=TRUE) & 
#				gene.summary.stats$Mean > 10 & 
				gene.summary.stats$Mean > quantile(gene.summary.stats$Mean, 0.5, na.rm = TRUE) &
				gene.summary.stats$CV < quantile(gene.summary.stats$CV, 0.5, na.rm = TRUE)
				);

			low.cv.genes.rank.lt10 <- rank(gene.summary.stats[low.cv.genes,'CV'], ties.method = 'random') <= 10;

			if (verbose) {
				cat('SampleContent: The following genes have been chosen to adjust for RNA Content.\n\n');
				print(data.frame("HK.Candidates" = rownames(x[low.cv.genes,][low.cv.genes.rank.lt10,])));
				cat('\n');
				}

			rna.content <- apply(
				X = x[low.cv.genes,][low.cv.genes.rank.lt10,,drop=FALSE],
				MARGIN = 2, 
				FUN = get.geo.mean,
				logged = logged
				);
			}

		else if (SampleContent == 'top.mean') {

			# sum each RNA 
			sum.rna <- apply(
				X = x[endogenous.genes,,drop=FALSE],
				MARGIN = 1,
				FUN = sum
				);

			# get the ranks of the RNA sums
			rank.rna <- (length(sum.rna) + 1) - rank(sum.rna, ties.method = 'first');

			rna.content <- apply(
				X = x[endogenous.genes, ][rank.rna <= 75,,drop=FALSE],
				MARGIN = 2,
				FUN = mean,
				na.rm = TRUE
				);
			}

		else if (SampleContent == 'top.geo.mean') {

			# sum each RNA
			sum.rna <- apply(
				X = x[endogenous.genes,,drop=FALSE],
				MARGIN = 1,
				FUN = sum
				);

			# get the ranks of the RNA sums
			rank.rna <- (length(sum.rna) + 1) - rank(sum.rna, ties.method = 'first');

			rna.content <- apply(
				X = x[endogenous.genes,,drop=FALSE][rank.rna <= 75,,drop=FALSE],
				MARGIN = 2,
				FUN = get.geo.mean,
				logged = logged
				);
			}

		else {
			stop('SampleContent: Unimplemented SampleContent method');
			}

		# handle cases where rna.content is estimated to be zero
		if (any(rna.content < 1)) {

			if (verbose) {
				# print warning 
				cat('SampleContent: The following samples have estimates of zero RNA content.  Consider removing them.\n\n');
				print(names(rna.content[rna.content < 1]));
				cat('\n');
				}

			# fudge rna.content to the expected minimum of 1.
			rna.content[rna.content < 1] <- 1; 
			}

		# calc normalization factor
		sampleContent.norm.factor <- mean(rna.content) / rna.content;

		# flag any samples that are outliers
		rna.content.sd.from.mean <- data.frame(rna.zscore = (rna.content - mean(rna.content)) / sd(rna.content));

		if (verbose & any(abs(rna.content.sd.from.mean) > 3)) {
			cat('SampleContent: The following samples have sample/rna content greater than \n\t3 standard deviations from the mean.\n\n');
			print(signif(subset(rna.content.sd.from.mean, abs(rna.content.sd.from.mean) > 3),3));
			cat('\n');
			}

		# adjust the data based on the normalization factor
		x <- t(apply(
			X = x,
			MARGIN = 1,
			FUN = '*',
			sampleContent.norm.factor
			));
		}
	return(
		list(
			x = x,
			sampleContent.norm.factor = sampleContent.norm.factor,
			rna.content = rna.content
			)
		);
	}
