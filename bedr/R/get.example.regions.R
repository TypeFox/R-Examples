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

get.example.regions <- function() {
	# unordered and overlapping
	a <- c("chr1:10-100", "chr1:101-200","chr1:200-210", "chr1:211-212", "chr2:10-50","chr10:50-100", "chr2:40-60","chr20:1-5");

	# nothing special but should overlap with (a)  
	b <- c("chr1:1-10", "chr1:111-250","chr1:2000-2010", "chr2:1-5","chr10:100-150", "chr2:40-60","chr20:6-10");

	# incorrect format start>end, strange chr, longer than expected chr
	d <-  c("chr1:10-100", "chr1:101:200", "chr10:100-101", "1:200-210", "chr2:50-40","chr10:50-50", "chr2:o-60","chr20:1-5x");

	x <- list(a = a, b = b, d = d)

	return(x);
	}

