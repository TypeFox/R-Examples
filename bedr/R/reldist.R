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

reldist <- function(x, y, detail = FALSE, check.zero.based = TRUE, check.chr = TRUE, check.valid = TRUE, check.sort = TRUE, check.merge = TRUE, verbose = TRUE){

	catv("RELATIVE DISTANCE\n");
	if (detail) {
		params <- "-detail";
		}
	else {
		params <- "";
		}

	x <- bedr(engine = "bedtools", input = list(a = x, b=y), method = "reldist", params = params, check.zero.based = check.zero.based, check.chr = check.chr, check.valid = check.valid, check.sort = check.sort, check.merge = check.merge, verbose = verbose);

	if (detail) {
		colnames(x) <- c("chr","start","stop","reldist");
		x$start <- as.integer(x$start);
		x$end <- as.integer(x$end);
		x$reldist <- as.integer(x$reldist);
		}
	else {
		x[,1] <- as.numeric(x[,1]);
		x[,2] <- as.numeric(x[,2]);
		x[,3] <- as.numeric(x[,3]);
		x[,4] <- as.numeric(x[,4]);
		}
#	colnames(x) <- c("reldist", "count", "total", "fraction");;

	return(x);
	}
