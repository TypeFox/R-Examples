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

convert2bed <- function(x, set.type = TRUE, check.zero.based = TRUE, check.chr = TRUE, check.valid = TRUE, check.sort = TRUE, check.merge = TRUE, verbose = TRUE) {
	# vector index (0), bed (1), index in first column (2), rownmames are index (3), unrecognized(4)

	catv("CONVERT TO BED\n");

	input.type <- determine.input(x, verbose = verbose);

	# extract the region and convert to bed
	if (input.type == 0) {
		x <- index2bed(x, set.type = set.type);
		}
	else if (input.type == 1) {
		colnames(x)[1:3] <- c("chr","start","end");
		}
	else if (input.type == 2) {
		bed <- index2bed(x[,1], set.type = set.type);
		x <- data.frame(bed, x[,-1,drop=FALSE], stringsAsFactors = FALSE);
		}
	else if (input.type == 3) {
		bed <- index2bed(rownames(x), set.type = set.type);
		x <- data.frame(bed, x, stringsAsFactors = FALSE);
		}
	else {
		catv("  ERROR: Not sure what the input format is!\n");
		stop();
		}

	# check if the region is valid
	if (check.valid) {
		attr(x, "input.type") <- 1; # temporarily assign type to bed b/c it was just converted
		is.valid <- is.valid.region(x, throw.error = TRUE, check.zero.based = check.zero.based, check.chr = check.chr, verbose = verbose);
		}

	# check if sorted
	if (check.sort) {
		is.sorted <- is.sorted.region(x, method = "lexicographical", check.valid = FALSE, check.zero.based = check.zero.based, check.chr = check.chr, check.merge = FALSE, verbose = verbose);
		
		#catv(" * Checking sort order... ");
		#if (!is.sorted) {
		#	catv("FAIL\n");
		#	catv(paste0("   The input for object is not *lexographically* ordered!\n   This can cause unexpected results for some set operations.\n   try: x <- bedr.sort.region(x)\n"));
		#	}
		#else {
		#	catv("PASS\n");
		#	}
		}
	
	# check if merged
	if (check.merge) {
		is.merged <- is.merged.region(x, check.valid = FALSE, check.zero.based = check.zero.based, check.chr = check.chr, check.sort = FALSE, verbose = verbose);
		#if (!is.merged) {
		#	catv(paste0("   The input for object has overlapping features!\n   This can cause unexpected results for some set operations.\n   i.e. x <- bedr.merge.region(x)\n", sep = ""));
		#	}
		}

	if (set.type) {
		input.type <- attr(x, "input.type");
		}

	attr(x, "input.type") <- input.type;
	attr(x, "is.index")   <- input.type %in% c(0,2,3);
	attr(x, "is.vector")  <- input.type == 0;

	return(x);
	}

