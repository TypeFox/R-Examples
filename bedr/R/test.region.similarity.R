# The bedr package is copyright (c) 2014 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

test.region.similarity <- function (x, y, n = 1e3, stratify.by.chr = FALSE, species = "human", build = "hg19", mask.gaps = FALSE, mask.repeats = FALSE, check.zero.based = TRUE, check.chr = TRUE, check.valid = TRUE, verbose = TRUE) {

#	is.valid.region(x, check.zero.based = check.zero.based, check.chr = check.chr, verbose = verbose);
#	is.valid.region(y, check.zero.based = check.zero.based, check.chr = check.chr, verbose = verbose);

	x <- convert2bed(x, check.zero.based = check.zero.based, check.chr = check.chr, check.valid = check.valid, verbose = verbose);
	y <- convert2bed(y, check.zero.based = check.zero.based, check.chr = check.chr, check.valid = check.valid, verbose = verbose);

	jaccard.orig     <- jaccard(x,y, check.chr = FALSE, check.valid = FALSE, check.sort = FALSE, check.merge = FALSE)[1,3]; # larger more similar (0-1)
	reldist.orig     <- reldist(x,y, check.chr = FALSE, check.valid = FALSE, check.sort = FALSE, check.merge = FALSE);# smaller more similar (0-0.5). should be uniform
	
	if (nrow(reldist.orig)!=0) {
		reldist.orig.median           <- tail(reldist.orig[cumsum(reldist.orig$count/sum(reldist.orig$count))  <= 0.5,"reldist"], 1);
		}
	else {
		reldist.orig.median <- 0.25;
		}

	if (0 %in% reldist.orig$reldist) {
		reldist.orig.fraction.zero    <- reldist.orig[reldist.orig$reldist == 0,"fraction"]; # the fraction at reldist 0.  if uniform all 50 fractions = 0.02
		}
	else {
		reldist.orig.fraction.zero <- 0;
		}
		
	if (0.5 %in% reldist.orig$reldist) {
		reldist.orig.fraction.fifty    <- reldist.orig[reldist.orig$reldist == 50,"fraction"]; # the fraction at reldist 0.  if uniform all 50 fractions = 0.02
		}
	else {
		reldist.orig.fraction.fifty <- 0;
		}
	
	catv("PERMUTATION TEST\n")

	#catv(paste0(" * iterator: 0"))
	#backspace <- "\b";

	batch.size <- 20;
	n.batches <- ceiling(n/batch.size);

	#if (batch.size > n) batch.size <- n;

	for (i in 1:n.batches) {

		if (1 == n.batches) {
			batch.size <- n;
			}
		else if  (i == n.batches) {
			batch.size <- batch.size %% n;
			}

	#	catv(paste0(backspace, i));
	#	backspace <- rep("\b", nchar(i));

		# send batches of permutations

		perm.results <- mclapply(list(1:n.batches), iterate.perm, x.region = x, y.region = y, batch.size = batch.size, jaccard.orig = jaccard.orig, reldist.orig.median = reldist.orig.median, reldist.orig.fraction.zero = reldist.orig.fraction.zero, reldist.orig.fraction.fifty = reldist.orig.fraction.fifty, stratify.by.chr = stratify.by.chr, species = species, build = build, mask.gaps = mask.gaps, mask.repeats = mask.repeats);
  
#		if (i == 100 ) {
#			if (jaccard.perm.gt > 20 && reldist.perm.median.lt > 20) break;
#			}
		}
	
	perm.results <- colSums(perm.results[[1]]);

	jaccard.p                 <- perm.results["jaccard.perm.gt"] / n;
	reldist.median.p          <- perm.results["reldist.perm.median.lt"] / n;
	reldist.fraction.zero.p   <- perm.results["reldist.perm.fraction.zero.gt"] / n;
	reldist.fraction.fifty.p  <- perm.results["reldist.perm.fraction.fifty.gt"] / n;

	tests <- c("jaccard.stat", "reldist.median","reldist.fraction.zero","reldist.fraction.50")
	effect.sizes <- c(jaccard.orig, reldist.orig.median, reldist.orig.fraction.zero, reldist.orig.fraction.fifty) 
	empirical.pvalues <- c(jaccard.p, reldist.median.p, reldist.fraction.zero.p, reldist.fraction.fifty.p);
	perm.results <- data.frame(test = tests, effect = effect.sizes, p = empirical.pvalues, stringsAsFactors = FALSE);

	return(list(results = perm.results, n = n, relative.distance = reldist.orig))
	}

# permute region and compare with original
# function used for parrallelization
iterate.perm <- function(x, x.region, y.region, batch.size, jaccard.orig, reldist.orig.median, reldist.orig.fraction.zero, reldist.orig.fraction.fifty, stratify.by.chr, species, build, mask.gaps, mask.repeats) {

	jaccard.perm.gt <- integer(batch.size);
	reldist.perm.median.lt <- integer(batch.size);  
	reldist.perm.fraction.zero.gt <- integer(batch.size); 
	reldist.perm.fraction.fifty.gt <- integer(batch.size);
	
	for (i in 1:batch.size) {
		x.perm       <- permute.region(x.region, stratify.by.chr = stratify.by.chr, species = species, build = build, mask.gaps = mask.gaps, mask.repeats = mask.repeats, is.checked = TRUE, sort.output = TRUE);
		jaccard.perm <- jaccard(x.perm, y.region, check.chr = FALSE, check.valid = FALSE, check.sort = FALSE, check.merge = FALSE, verbose = FALSE)[1,3];
		reldist.perm <- reldist(x.perm, y.region, check.chr = FALSE, check.valid = FALSE, check.sort = FALSE, check.merge = FALSE, verbose = FALSE);

		# if reldist fails
		if (nrow(reldist.perm)!=0) {
			reldist.perm.median           <- tail(reldist.perm[cumsum(reldist.perm$count/sum(reldist.perm$count))  <= 0.5,"reldist"], 1);
			}
		else {
			reldist.perm.median <- 0.25;
			}

		if (0 %in% reldist.perm$reldist) {
			reldist.perm.fraction.zero    <- reldist.perm[reldist.perm$reldist == 0,"fraction"]; # the fraction at reldist 0.  if uniform all 50 fractions = 0.02
			}
		else {
			reldist.perm.fraction.zero <- 0;
			}
			
		if (0.5 %in% reldist.perm$reldist) {
			reldist.perm.fraction.fifty    <- reldist.perm[reldist.perm$reldist == 50,"fraction"]; # the fraction at reldist 0.  if uniform all 50 fractions = 0.02
			}
		else {
			reldist.perm.fraction.fifty <- 0;
			}
		
		if (jaccard.perm >= jaccard.orig) jaccard.perm.gt <- jaccard.perm.gt + 1;
		if (reldist.perm.median <= reldist.orig.median)  reldist.perm.median.lt <- reldist.perm.median.lt + 1;
		if (reldist.perm.fraction.zero >= reldist.orig.fraction.zero)  reldist.perm.fraction.zero.gt <- reldist.perm.fraction.zero.gt + 1;
		if (reldist.perm.fraction.fifty >= reldist.orig.fraction.fifty)  reldist.perm.fraction.fifty.gt <- reldist.perm.fraction.fifty.gt + 1;
		}

	return(data.frame(batch.size, jaccard.perm.gt, reldist.perm.median.lt, reldist.perm.fraction.zero.gt, reldist.perm.fraction.fifty.gt))
	}
