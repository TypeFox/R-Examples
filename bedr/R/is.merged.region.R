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

is.merged.region <- function(x, check.zero.based = TRUE, check.chr = TRUE, check.valid = TRUE, check.sort = TRUE, verbose = FALSE) {

	x.merge   <- bedr.merge.region(x, check.zero.based = check.zero.based, check.chr = check.chr, check.valid = check.valid, check.sort = check.sort, verbose = FALSE);
	if (is.vector(x)) {
		is.merged <- length(x) == length(x.merge);
		}
	else {
		if (is.null(nrow(x)) || is.null(nrow(x.merge))) {
			is.merged <- FALSE;
			}
		else {
			is.merged <- nrow(x) == nrow(x.merge);
			}
		}
	
	catv(" * Checking for overlapping 'contiguous' regions... ")
	if (is.merged) {
		catv("PASS\n");
		}
	else {
		catv("FAIL\n")
		catv(paste0("   The input for object has overlapping features!\n   This can cause unexpected results for some set operations.\n   i.e. x <- bedr.merge.region(x)\n", sep = ""));
		}

	return(is.merged);
	}
