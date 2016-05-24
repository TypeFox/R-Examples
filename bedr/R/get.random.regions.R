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

get.random.regions <- function(n = 10, chr = NULL, species = "human", build = "hg19", size.mean = 10, size.sd = 0.25, mask.gaps = FALSE, mask.repeats = FALSE, sort.output = TRUE, verbose = TRUE) {

	region.size  <- ceiling(rlnorm(n, meanlog = size.mean, sdlog = size.sd));

	if (!is.null(chr)) {
		stratify.by.chr = TRUE;
		x <- data.frame(chr = chr, start = 1, end = region.size, stringsAsFactors = FALSE); # don't worry about the length b/c correct at end
		}
	else {
		x <- data.frame(chr = "chr1",start =  1, end = region.size, stringsAsFActors = FALSE);
		stratify.by.chr = FALSE;
		}

	x <- permute.region(x, stratify.by.chr = stratify.by.chr, species = species, build = build, mask.gaps = mask.gaps, mask.repeats = mask.repeats, sort.output = sort.output, is.checked = TRUE, verbose = verbose);

	if (length(x) > length(n)) {
		x <- head(x, n);
		}

#	x <- paste0(x$chr,":",x$start,"-",x$end);

	return(x);
	}
