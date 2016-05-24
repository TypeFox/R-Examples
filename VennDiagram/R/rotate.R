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

# rotate & reverse triple Venn diagrams
rotate <- function(
	area.vector,
	category.vector,
	lwd,
	lty,
	col,
	label.col,
	cex,
	fontface,
	fontfamily,
	cat.col,
	cat.cex,
	cat.fontface,
	cat.fontfamily,
	alpha,
	rotation,
	reverse,
	fill = NULL
	) {

	rot.1.f.7 <- 1:7;
	rot.1.f.3 <- 1:3;
	rot.1.r.7 <- c(3,2,1,6,5,4,7);
	rot.1.r.3 <- c(2,1,3);
	rot.2.f.7 <- c(3,6,7,2,5,4,1);
	rot.2.f.3 <- c(2,3,1);
	rot.2.r.7 <- c(7,6,3,4,5,2,1);
	rot.2.r.3 <- c(3,2,1);
	rot.3.f.7 <- c(7,4,1,6,5,2,3);
	rot.3.f.3 <- c(3,1,2);
	rot.3.r.7 <- c(1,4,7,2,5,6,3);
	rot.3.r.3 <- c(1,3,2);

	if (reverse) {
		if (1 == rotation) {
			order.7 <- rot.1.r.7;
			order.3 <- rot.1.r.3;
			}
		else if (2 == rotation) {
			order.7 <- rot.2.r.7;
			order.3 <- rot.2.r.3;
			}
		else if (3 == rotation) {
			order.7 <- rot.3.r.7;
			order.3 <- rot.3.r.3;
			}
		}
	else {
		if (1 == rotation) {
			order.7 <- rot.1.f.7;
			order.3 <- rot.1.f.3;
			}
		else if (2 == rotation) {
			order.7 <- rot.2.f.7;
			order.3 <- rot.2.f.3;
			}
		else if (3 == rotation) {
			order.7 <- rot.3.f.7;
			order.3 <- rot.3.f.3;
			}
		}

	return(
		list(
			a = area.vector[order.7],
			c = category.vector[order.3],
			lwd = lwd[order.3],
			lty = lty[order.3],
			col = col[order.3],
			label.col = label.col[order.7],
			cex = cex[order.7],
			fontface = fontface[order.7],
			fontfamily = fontfamily[order.7],
			cat.col = cat.col[order.3],
			cat.cex = cat.cex[order.3],
			cat.fontface = cat.fontface[order.3],
			cat.fontfamily = cat.fontfamily[order.3],
			fill = fill[order.3],
			alpha = alpha[order.3]
			)
		);

	}
