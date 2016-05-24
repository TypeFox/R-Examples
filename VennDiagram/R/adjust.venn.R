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

### FUNCTION TO FIT VENN DIAGRAM TO SIZE ##########################################################
adjust.venn <- function(gList1, margin = 0.01, ...) {
	# adjusts the positions of the ellipses to account for changes made to the plot so far

	# vectors containing all the data sequentially
	x.vect <- vector();
	y.vect <- vector();

	# list containing each vector in x.vect and y.vect
	x.list <- list();
	y.list <- list();

	for (i in 1:length(gList1)) {
		x.vect <- c(x.vect, as.vector(gList1[i][[1]]$x, mode = 'numeric'));
		y.vect <- c(y.vect, as.vector(gList1[i][[1]]$y, mode = 'numeric'));
		x.list[[i]] <- as.vector(gList1[i][[1]]$x, mode = 'numeric');
		y.list[[i]] <- as.vector(gList1[i][[1]]$y, mode = 'numeric');
		}

	# get dimensions of ellipses
	max.x <- max(x.vect) + margin;
	min.x <- min(x.vect) - margin;
	max.y <- max(y.vect) + margin;
	min.y <- min(y.vect) - margin;
	x.centre <- (max.x + min.x) / 2;
	y.centre <- (max.y + min.y) / 2;

	size <- 0.99;

	# wider than tall
	if (max.x - min.x >= max.y - min.y) {
		for (i in 1:length(x.list)) {
			x.list[[i]] <- unit((x.list[[i]] - x.centre) * (size / (max.x - min.x)) + 0.5, 'npc');
			y.list[[i]] <- unit((y.list[[i]] - y.centre) * (size / (max.x - min.x)) + 0.5, 'npc');
			}
		}
	else {
		for (i in 1:length(x.list)) {
			x.list[[i]] <- unit((x.list[[i]] - x.centre) * (size / (max.y - min.y)) + 0.5, 'npc');
			y.list[[i]] <- unit((y.list[[i]] - y.centre) * (size / (max.y - min.y)) + 0.5, 'npc');
			}
		}

	for (i in 1:length(gList1)) {
		gList1[i][[1]]$x <- x.list[[i]];
		gList1[i][[1]]$y <- y.list[[i]];
		}

	return(gList1);
	}
