# The VennDiagram package is copyright (c) 2012 Ontario Institute for Cancer Research (OICR)
# This package and its accompanying libraries is free software; you can redistribute it and/or modify it under the terms of the GPL
# (either version 1, or at your option, any later version) or the Artistic License 2.0.  Refer to LICENSE for the full license text.
# OICR makes no representations whatsoever as to the SOFTWARE contained herein.  It is experimental in nature and is provided WITHOUT
# WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR IMPLIED. OICR MAKES NO REPRESENTATION
# OR WARRANTY THAT THE USE OF THIS SOFTWARE WILL NOT INFRINGE ANY PATENT OR OTHER PROPRIETARY RIGHT.
# By downloading this SOFTWARE, your Institution hereby indemnifies OICR against any loss, claim, damage or liability, of whatsoever kind or
# nature, which may arise from your Institution's respective use, handling or storage of the SOFTWARE.
# If publications result from research using this SOFTWARE, we ask that the Ontario Institute for Cancer Research be acknowledged and/or
# credit be given to OICR scientists, as scientifically appropriate.

### FUNCTION TO FIND DISTANCE BETWEEN CIRCLES #####################################################
find.dist <- function(area1, area2, cross.area, inverted = FALSE) {

	if (inverted) {
		r2 <- sqrt(area1 / pi);
		r1 <- sqrt(area2 / pi);
		}
	else {
		r1 <- sqrt(area1 / pi);
		r2 <- sqrt(area2 / pi);
		}

	# set up a sequence of distances corresponding to full intersection to 0 intersection with set resolution (step)
	d <- r1 + r2;
	resolution <- 0.001;
	d.list <- seq(r1 - r2 + resolution, d, resolution);
	int.list <- sapply(d.list, find.intersect, r1, r2);
	match.list <- (int.list >= cross.area);
	index.true <- length(match.list[match.list]);
	index.false <- index.true + 1;

	if (0 == index.true) {
		return(d.list[index.false]);
		}
	else {
		return(mean(d.list[index.true], d.list[index.false]));
		}
	}
