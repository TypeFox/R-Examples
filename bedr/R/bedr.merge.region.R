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

bedr.merge.region <- function(x, distance = 0, list.names = TRUE, number = FALSE, stratify.by = NULL, check.zero.based = TRUE, check.chr = TRUE, check.valid = TRUE, check.sort = TRUE, verbose = TRUE) {

catv("MERGING\n")

	if (is.vector(x) || ncol(x)<4) {list.names=FALSE}

	if (list.names && number) {
		list.names.param <- "-c 4,1 -o collapse,count"
		}
	else if (list.names) {
		list.names.param <- "-c 4 -o collapse"
		}
	else if (number) {
		list.names.param <- "-c 1 -o count"
		}
	else {
		list.names.param <- ""
		}

	# deprecated
	# list.names.param <- ifelse(list.names, "-nms", "");
	# number     <- ifelse(number, "-n", "");

	if(is.vector(x, mode = "character")) {
		n.rec.before <- length(x)
		}
	else {
		n.rec.before <- nrow(x)
		}
	
	# run validation first 
	if (check.valid) {
		is.valid <- is.valid.region(x, check.zero.based = check.zero.based, check.chr = check.chr, throw.error = TRUE, verbose = verbose);
		}

	# need to check sort first to avoid error in bedtools merge
	if (check.sort) {
		if (!is.sorted.region(x, check.zero.based = FALSE, check.chr = FALSE, check.valid = FALSE, check.merge = FALSE, verbose = FALSE)) {
			catv(" * Bedtools requires sorted input for merging!\n");
			catv("   Your data is being automatically sorted\n");
			x <- bedr.sort.region(x, check.zero.based = FALSE, check.chr = FALSE, check.valid = FALSE, check.merge = FALSE, verbose = FALSE);
			}
		}

	# should the collapsing be done within groups i.e. genes
	if (is.null(stratify.by)) {
		x <- bedr(engine = "bedtools", input = list(i = x), method = "merge", params = paste("-d", distance, list.names.param), check.zero.based = FALSE, check.chr = FALSE, check.valid = FALSE, check.merge = FALSE, check.sort = FALSE, verbose = verbose);
		}
	else {
		if (!stratify.by %in% colnames(x)) {
			catv("ERROR: the statified column does not exist \n");
			return(1);
			}
		x <- by(x, x[,stratify.by], function(x) {bedr(engine = "bedtools", input = list(i = x), method = "merge", params = paste("-d", distance, list.names.param), check.zero.based = FALSE, check.chr = FALSE, check.valid = FALSE, check.merge = FALSE, check.sort = FALSE, verbose = verbose)});
		x <- do.call(rbind, x);
		}

	if(is.vector(x, mode = "character")) {
		n.rec.after <- length(x)
		}
	else {
		n.rec.after <- nrow(x)
		}

	catv(paste0(" * Collapsing ", n.rec.before, " --> ", n.rec.after, " regions... NOTE\n"))

	if (length(colnames(x))>=3) {
		colnames(x)[1:3] <- c("chr","start","end");
		}

	# replace repeating merged names
	if (list.names) {
		x[,4] <- unlist(lapply(lapply(strsplit(x[,4], split=","), "unique"), paste, collapse=",")); # yuck
		colnames(x)[4] <- "names";
		}

	return(x);
	}
