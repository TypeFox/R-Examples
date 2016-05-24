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

bedr.join.multiple.region <- function(x = list(), fraction.overlap = 1/1e9, empty = FALSE, missing.values = ".", cluster = FALSE, species = "human", build = "hg19", check.zero.based = TRUE, check.chr = TRUE, check.valid = TRUE, check.sort = TRUE, check.merge = TRUE, verbose = TRUE) {

	catv("JOINING\n")

	if (is.null(names(x))) names(x) <- letters[1:length(x)];
	dataset.names.orig    <- names(x);
	dataset.names         <- paste("-names ", paste(dataset.names.orig, collapse = " "));

	fraction.overlap <- ifelse(fraction.overlap == 1/1e9, "", paste("-f ", fraction.overlap));
	missing.values   <- paste("-filler", missing.values, sep = " ");
	empty            <- ifelse(empty, "-empty", "");
	cluster          <- ifelse(!cluster, "", "-cluster");
	genome.file      <- paste0("genomes/", species, ".", build, ".genome");
	genome.file      <- system.file(genome.file, package = "bedr");
	genome.file      <- paste0("-g ", genome.file);
	
	if (!check.sort || !check.merge) {
		catv("  * bedr.join.multiple.region assumes sorted regions.\n    Also, overlapping regions cause unexpected results.\n")
		}

	# only the first region should be names i
	names(x)[1] <- "i";
	for (i in 2:length(x)) {
		names(x)[i] <- "";
		}

	xyz <- bedr(input = x, method = "multiinter", params = paste(fraction.overlap, dataset.names, genome.file, missing.values, cluster), check.zero.based = check.zero.based, check.chr = check.chr, check.valid = check.valid, check.sort = check.sort, check.merge = check.merge, verbose = verbose);

	colnames(xyz)[(ncol(xyz)-(1+length(x))):ncol(xyz)] <- c("n.overlaps","names", dataset.names.orig);

	return(xyz);

	}
