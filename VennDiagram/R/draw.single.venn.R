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

### FUNCTION TO DRAW VENN DIAGRAM WITH A SINGLE SET ###############################################
draw.single.venn <- function(
	area,
	category = "",
	lwd = 2,
	lty = "solid",
	col = "black",
	fill = NULL,
	alpha = 0.5,
	label.col = "black",
	cex = 1,
	fontface = "plain",
	fontfamily = "serif",
	cat.pos = 0,
	cat.dist = 0.025,
	cat.cex = 1,
	cat.col = "black",
	cat.fontface = "plain",
	cat.fontfamily = "serif",
	cat.just = list(c(0.5, 0.5)),
	cat.default.pos = "outer",
	cat.prompts = FALSE,
	rotation.degree = 0,
	rotation.centre = c(0.5, 0.5),
	ind = TRUE,
	...
	) {

	# check parameter lengths
	if (length(category) != 1) { flog.error("Unexpected parameter length for 'category'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'category'"); }
	if (length(lwd) != 1) { flog.error("Unexpected parameter length for 'lwd'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'lwd'"); }
	if (length(lty) != 1) { flog.error("Unexpected parameter length for 'lty'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'lty'"); }
	if (length(col) != 1) { flog.error("Unexpected parameter length for 'col'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'col'"); }
	if (length(label.col) != 1) { flog.error("Unexpected parameter length for 'label.col'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'label.col'"); }
	if (length(cex) != 1) { flog.error("Unexpected parameter length for 'cex'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'cex'"); }
	if (length(fontface) != 1) { flog.error("Unexpected parameter length for 'fontface'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'fontface'"); }
	if (length(fontfamily) != 1) { flog.error("Unexpected parameter length for 'fontfamily'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'fontfamily'"); }
	if (length(fill) != 1 & length(fill) != 0) { flog.error("Unexpected parameter length for 'fill'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'fill'"); }
	if (length(alpha) != 1 & length(alpha) != 0) { flog.error("Unexpected parameter length for 'alpha'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'alpha'"); }
	if (length(cat.pos) != 1) { flog.error("Unexpected parameter length for 'cat.pos'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'cat.pos'"); }
	if (length(cat.dist) != 1) { flog.error("Unexpected parameter length for 'cat.dist'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'cat.dist'"); }
	if (length(cat.col) != 1) { flog.error("Unexpected parameter length for 'cat.col'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'cat.col'"); }
	if (length(cat.cex) != 1) { flog.error("Unexpected parameter length for 'cat.cex'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'cat.cex'"); }
	if (length(cat.fontface) != 1) { flog.error("Unexpected parameter length for 'cat.fontface'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'cat.fontface'"); }
	if (length(cat.fontfamily) != 1) { flog.error("Unexpected parameter length for 'cat.fontfamily'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'cat.fontfamily'"); }
	if (!(class(cat.just) == "list" & length(cat.just) == 1 & length(cat.just[[1]]) == 2)) { flog.error("Unexpected parameter format for 'cat.just'",name="VennDiagramLogger")
stop("Unexpected parameter format for 'cat.just'"); }

	cat.pos <- cat.pos + rotation.degree;

	# check category label defaults
	if (cat.default.pos != 'outer' & cat.default.pos != 'text' & category != '' & cat.prompts) {
		flog.info("No default location recognized.  Automatically changing to 'outer'",name="VennDiagramLogger");
		cat.default.pos <- 'outer';
		}
	if (cat.default.pos == "outer" & category != "" & cat.prompts) {
		flog.info("Placing category labels at default outer locations.  Use 'cat.pos' and 'cat.dist' to modify location.",name="VennDiagramLogger");
		flog.info(paste("Current 'cat.pos':", cat.pos, 'degrees'),name="VennDiagramLogger");
		flog.info(paste("Current 'cat.dist':", cat.dist),name="VennDiagramLogger");
		}
	if (cat.default.pos == "text" & category != "" & cat.prompts) {
		flog.info("Placing category labels at default text locations.  Use 'cat.pos' and 'cat.dist' to modify location.",name="VennDiagramLogger");
		flog.info(paste("Current 'cat.pos':", cat.pos, 'degrees'),name="VennDiagramLogger");
		flog.info(paste("Current 'cat.dist':", cat.dist),name="VennDiagramLogger");
		}

	max.circle.size = 0.2;

	# obtain radius corresponding to the circle with given area and convert it to Grid dimensions
	r1 <- sqrt(area / pi);
	shrink.factor <- max.circle.size / r1;
	r1 <- r1 * shrink.factor;

	# initialize gList to hold all Grobs generated
	grob.list <- gList();

	# plot Venn diagram
	tmp <- VennDiagram::ellipse(
		x = 0.5,
		y = 0.5,
		a = r1,
		b = r1,
		gp = gpar(
			lty = 0,
			fill = fill,
			alpha = alpha
			)
		);
	grob.list <- gList(grob.list, tmp);
	tmp <- VennDiagram::ellipse(
		x = 0.5,
		y = 0.5,
		a = r1,
		b = r1,
		gp = gpar(
			lwd = lwd,
			lty = lty,
			col = col,
			fill = 'transparent'
			)
		);
	grob.list <- gList(grob.list, tmp);
	tmp <- textGrob(
		label = area,
		x = 0.5,
		y = 0.5,
		gp = gpar(
			col = label.col,
			cex = cex,
			fontface = fontface,
			fontfamily = fontfamily
			)
		);
	grob.list <- gList(grob.list, tmp);

	if (cat.default.pos == 'outer') { cat.pos.1 <- find.cat.pos(0.5, 0.5, cat.pos, cat.dist, r1) }
	if (cat.default.pos == 'text') { cat.pos.1 <- find.cat.pos(0.5, 0.5, cat.pos, cat.dist) }
	tmp <- textGrob(
		label = category,
		x = cat.pos.1$x,
		y = cat.pos.1$y,
		just = cat.just[[1]],
		gp = gpar(
			col = cat.col,
			cex = cat.cex,
			fontface = cat.fontface,
			fontfamily = cat.fontfamily
			)
		);
	grob.list <- gList(grob.list, tmp);
	grob.list <- VennDiagram::adjust.venn(VennDiagram::rotate.venn.degrees(grob.list, rotation.degree, rotation.centre[1], rotation.centre[2]), ...);
	if (ind) { grid.draw(grob.list); }
	return(grob.list);
	}
