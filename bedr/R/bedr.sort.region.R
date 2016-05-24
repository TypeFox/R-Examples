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

bedr.sort.region <- function(x, method = "lexicographical", engine = "R", chr.to.num = c("X" = 23, "Y" = 24, "M" = 25),  check.zero.based = TRUE, check.chr = TRUE, check.valid = TRUE, check.merge = TRUE, verbose = TRUE) {
# engine bedtools, bedops, unix, R

catv("SORTING\n");

	if (!is.vector(x)) {
		header  <- colnames(x);
		}
	else {
		header <- "index";
		}
	
	# run validation first
	if (check.valid) {
		is.valid <- is.valid.region(x, check.zero.based = check.zero.based, check.chr = check.chr, throw.error = TRUE, verbose = verbose);
		}

	x       <- convert2bed(x, set.type = FALSE, check.zero.based = FALSE, check.chr = FALSE, check.valid = FALSE, check.sort = FALSE, check.merge = FALSE, verbose = FALSE);
	is.unix <- .Platform$OS.type == "unix";
	if (is.unix) {
		# handle Solaris as its non GNU, and hence no --version option
		if (grepl('SunOS', Sys.info()['sysname'])) {
			sort.version <- 1;
			}
		else{
			sort.version <- as.numeric(
				strsplit(system("sort --version", intern = TRUE), split = " ")[[1]][4]
				);
			}
		}

	# check engine
	if (!engine %in% c("unix", "bedtools", "bedops", "R")) {
		catv("    ERROR: Sort engine must be one of 'unix, bedtools, bedops, R'\n");
		stop();
		}

	# check if unix is the os
	if (!is.unix && engine == "unix") {
		catv("    This is not unix. Defaulting to R for natural sorting.\n");
		engine <- "R";
		}

	# check if the version of gnu sort supports natural/version sorting
	if (method == "natural" && is.unix && engine == "unix" && sort.version < 7.0 ) {
		catv("    Version of gnu sort is too low for natural sorting.  Defaulting to R.\n");
		engine <- "R";
		}

	# check if R/unix were chosen for natural sorting 
	if (method == "natural" && !engine %in% c("R","unix")) {
		if (is.unix &&  sort.version >= 7.0) {
			catv("    Bedtools/bedops does not support natural sorting.  Defaulting to unix.\n");
			engine <- "unix";
			}
		else {
			catv("    Bedtools/bedops does not support natural sorting.  Defaulting to R.\n");
			engine <- "R";
			}
		}

	if (method == "natural") {
		if (engine == "unix") {
			tmp.file   <- create.tmp.bed.file(x, "sort");
			command    <- paste("sort -k1,1V -k2,2n", tmp.file);
			output     <- system(command, intern = TRUE);
			output     <- strsplit2matrix(output, split = "\t");
			output[,2] <- as.numeric(output[,2]);
			output[,3] <- as.numeric(output[,3]);
			}
		else if (engine == "R"){
			catv("    Natural sorting is done in R which could be memory intensive for large files\n");
			chr.integer <- gsub("^chr", "", x[,1], ignore.case = TRUE);
			if (!is.null(chr.to.num)) {
				for (chr.char in names(chr.to.num)) {
					chr.integer[chr.integer == chr.char] <- chr.to.num[chr.char];
					}
				}
			chr.integer <- as.numeric(chr.integer);
			x[,2] <- as.numeric(x[,2]);
			x[,3] <- as.numeric(x[,3]);
			output <- x[order(chr.integer, x[,2]),];
			}
		}
	else if (method %in% c("lexicographical","lex")) {

		if (engine == "unix") {
			tmp.file   <- create.tmp.bed.file(x, "sort");
			command    <- paste("sort -k1,1 -k2,2n", tmp.file);
			output     <- system(command, intern = TRUE);
			output     <- strsplit2matrix(output, split = "\t");
			output[,2] <- as.numeric(output[,2]);
			output[,3] <- as.numeric(output[,3]);
			}
		else if (engine == "bedtools") {
			output <- bedr(method = "sort", engine = "bedtools", input = list(i=x), check.zero.based = FALSE, check.chr = FALSE, check.valid = FALSE, check.sort = FALSE, check.merge = FALSE, verbose = verbose);
			output[,2] <- as.numeric(output[,2]);
			output[,3] <- as.numeric(output[,3]);
			}
		else if (engine == "bedops") {
			output <- bedr(method = "sort", engine = "bedops", input = list(i=x), check.zero.based = FALSE, check.chr = FALSE, check.valid = FALSE, check.sort = FALSE, check.merge = FALSE, verbose = verbose);
			output[,2] <- as.numeric(output[,2]);
			output[,3] <- as.numeric(output[,3]);
			}
		else {
			output <- x[order(x[,1],x[,2]),];
			}
		}
	else {
		catv("   ERROR: Sort method can be natural or lexicographical\n")
		stop();
		}

	# recreate original index
	if (attr(x,"is.index")) {
		
		old.scipen <- getOption("scipen")
		options(scipen = 999);
		index <- paste(output[,1],":",output[,2],"-",output[,3], sep = "");
		options(scipen = old.scipen);

		# delete the chr,start,stop columns and return a vector that was the input
		if (attr(x,"is.vector")) {
			output <- as.vector(index); 
			}
		else {
			# set rownames to index if possible
			if (nrow(output) == length(unique(index))) rownames(output) <- index;
			output <- output[,!colnames(output) %in% c("chr","start","end"),drop = FALSE];
			output <- data.frame(index, output, stringsAsFactors = FALSE);
			}
		}

	# check merging after sorting b/c bedtools doesn't work if not sorted!
	if (check.merge && any(is.merged.region(output, check.zero.based = FALSE, check.chr = FALSE, check.valid = FALSE, check.sort = FALSE))) {
		catv(" * Overlapping regions can cause unexpected results.\n")
		}

	if (!is.vector(output)) colnames(output) <- header;

	return(output);
	}
