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

# flip a Venn diagram
flip.venn <- function(gList1, axis = 'v') {

	x.vect <- vector();
	y.vect <- vector();
	x.list <- list();
	y.list <- list();

	if ('v' == axis) {
		for (i in 1:length(gList1)) {
			x.vect <- c(x.vect, as.vector(gList1[i][[1]]$x, mode = 'numeric'));
			x.list[[i]] <- as.vector(gList1[i][[1]]$x, mode = 'numeric');
			}
		for (i in 1:length(x.list)) { x.list[[i]] <- unit(1 - x.list[[i]], 'npc') }
		for (i in 1:length(gList1)) { gList1[i][[1]]$x <- x.list[[i]] }

		return(gList1);
		}

	else if ('h' == axis) {
		for (i in 1:length(gList1)) {
			y.vect <- c(y.vect, as.vector(gList1[i][[1]]$y, mode = 'numeric'));
			y.list[[i]] <- as.vector(gList1[i][[1]]$y, mode = 'numeric');
			}
		for (i in 1:length(y.list)) { y.list[[i]] <- unit(1 - y.list[[i]], 'npc') }
		for (i in 1:length(gList1)) { gList1[i][[1]]$y <- y.list[[i]] }

		return(gList1);
		}
	else {
		flog.error('Unknown axis type',name="VennDiagramLogger")
stop('Unknown axis type');
		}
	}
