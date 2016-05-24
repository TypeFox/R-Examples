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

is.valid.region <- function(x, check.zero.based = TRUE, check.chr = TRUE, throw.error = FALSE , verbose = TRUE) {
	# vector index (0), bed (1), index in first column (2), rownmames are index (3), unrecognized(4)

	catv("VALIDATE REGIONS\n");
	is.error <- FALSE;
	# check if input type has already been determined
	input.type <- attr(x, "input.type"); 

	# determine input format
	if (is.null(input.type)) {
		input.type <- determine.input(x, verbose = verbose);
		attr(x, "input.type") <- input.type;
		}

	if (input.type == 0) { # index vector
		is.index <- TRUE;
		}
	else if (input.type == 1) { # bed
		# check the data types.  this can screw up sorting in R
		is.correct.datatype <- all(sapply(x[,1:3], mode) == c("character", "numeric", "numeric"));
		if (!is.correct.datatype) {
			catv(" * Check data types... FAIL\n   Make sure the the chr column is character and the start and end positions are numbers\n");
			if (is.factor(x[,1])) catv("   your chr column is a factor!\n")
			is.error <- TRUE;
			}

		old.scipen <- getOption("scipen");
		options(scipen = 999);
		x <- paste0(x[,1],":",x[,2],"-",x[,3]);
		options(scipen = old.scipen);

		}
	else if (input.type == 2){ # index tabular
		x <- x[,1];
		}
	else if (input.type == 3) { # index tabular as rownames
		x <- rownames(x);
		}
	else {
		catv("ERROR: Not sure what the input format is!\n");
		stop();
		}

	is.string <- is.character(x);
	if (!is.string) {
		catv(" * Check if index is a string... FAIL\n   Is the input a factor?\n");
		is.error <- TRUE;
		}
	else {
		catv(" * Check if index is a string... PASS\n");
		}

	if (check.chr) {
		pattern <- "(^chr[0-9XYMTxymt]{1,2}|^chrGL.*):\\d*-\\d*$";
		}
	else {
		pattern <- "^.*:\\d*-\\d*$";
		}
	
	is.valid.pattern <- grepl(pattern, x);

	if (any(!is.valid.pattern)) {
		catv(" * Check index pattern... FAIL\n   Use check.chr = FALSE if no 'chr' prefix\n");
		if(verbose) print(head(x[!is.valid.pattern]));
		is.error <- TRUE;
		}
	else {
		catv(" * Check index pattern... PASS\n");
		}

	# split for more checks
	x <- index2bed(x);
	is.missing <- is.na(x[,1]) | is.na(x[,2]) | is.na(x[,3]);
	if (any(is.missing)) {
		catv("  * Check for missing values... FAIL\n");
		if(verbose) if (sum(is.missing)>5) {print(head(x[is.missing,]))} else {print(x[is.missing,])};
		is.error <- TRUE;
		}
	else {
		catv(" * Check for missing values... PASS\n");
		}

	is.valid.position <- x$start <= x$end;
	if (any(is.missing) )  is.valid.position[is.missing] <- TRUE;

	if (any(!is.valid.position)) {
		catv(" * Check if start is larger than end position... FAIL.\n");
		if(verbose) if (sum(!is.valid.position)>5) {print(head(x[!is.valid.position,]))} else {print(x[!is.valid.position,])};
		is.error <- TRUE;
		}
	else {
		catv(" * Check for larger start position... PASS.\n");
		}

	if (check.zero.based) {
		is.zero.based <- x$start != x$end;
		is.zero.based[is.na(is.zero.based)] <- TRUE;
		if (any(!is.zero.based)) {
			catv(" * Check if zero based... FAIL\n   There are identical start and stop positions.\n   This usually indicates 1 based point variants.\n");
			if(verbose) if (sum(!is.zero.based)>5) {print(head(x[!is.zero.based,]))}else {print(x[!is.zero.based,])};
			}
		else {
			catv(" * Check if zero based... PASS\n");
			}
		}
	else {
		is.zero.based <- TRUE;
		}

	if(throw.error && is.error) {
		catv("\nERROR: Sorry, the program is stopping until problem features are fixed.  \n    You can try using is.valid.region and fix.region functions to debug\n");
		stop();
		}

	is.valid <- !(!is.string | is.missing | !is.zero.based | !is.valid.position | !is.valid.pattern);

	return(is.valid);
	
	}

fix.region <- function(x) {
	x <- lapply(x, function(y) y[is.valid.region(y)])
	}
