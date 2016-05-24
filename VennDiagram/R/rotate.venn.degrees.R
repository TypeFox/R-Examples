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

### FUNCTION TO ROTATE DIAGRAM BY DEGREES #########################################################
rotate.venn.degrees <- function(gList1, angle = 90, x.centre = 0.5, y.centre = 0.5) {

	x.vect <- vector();
	y.vect <- vector();
	x.list <- list();
	y.list <- list();
	x2.list <- list();
	y2.list <- list();

	for (i in 1:length(gList1)) {
		x.vect <- c(x.vect, as.vector(gList1[i][[1]]$x, mode = 'numeric'));
		x.list[[i]] <- as.vector(gList1[i][[1]]$x, mode = 'numeric');
		}
	for (i in 1:length(gList1)) {
		y.vect <- c(y.vect, as.vector(gList1[i][[1]]$y, mode = 'numeric'));
		y.list[[i]] <- as.vector(gList1[i][[1]]$y, mode = 'numeric');
		}
	for (i in 1:length(x.list)) {
		x2.list[[i]] <- (x.list[[i]] - x.centre) * cos(angle * pi / 180) - (y.list[[i]] - y.centre) * sin(angle * pi / 180) + x.centre;
		y2.list[[i]] <- (x.list[[i]] - x.centre) * sin(angle * pi / 180) + (y.list[[i]] - y.centre) * cos(angle * pi / 180) + y.centre;
		}
	for (i in 1:length(gList1)) {
		gList1[i][[1]]$y <- unit(y2.list[[i]], 'npc');
		gList1[i][[1]]$x <- unit(x2.list[[i]], 'npc');
		}

	return(gList1);
	}

