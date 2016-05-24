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

cluster.region <- function(x, distance = 0, check.zero.based = TRUE, check.chr = TRUE, check.valid = TRUE, check.sort = TRUE, verbose = TRUE) {

catv("CLUSTERING\n")

	if (is.vector(x) || ncol(x)<4) {list.names=FALSE}

	if(is.vector(x, mode = "character")) {
		n.rec.before <- length(x)
		}
	else {
		n.rec.before <- nrow(x)
		}

	x <- bedr(engine = "bedtools", input = list(i = x), method = "cluster", params = paste("-d", distance), check.zero.based = check.zero.based, check.chr = check.chr, check.valid = check.valid, check.merge = FALSE, check.sort = check.sort, verbose = verbose);

	if(is.vector(x, mode = "character")) {
		n.rec.after <- length(unique(x[,2]));
		colnames(x)[2] <- "regionIndex";
		}
	else {
		n.rec.after <- length(unique(x[,ncol(x)]))
		colnames(x)[ncol(x)] <- "regionIndex";
		}

	catv(paste0("  * Collapsing ", n.rec.before, " --> ", n.rec.after, " regions... NOTE\n"))

	return(x);
	}
