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

### FUNCTION TO ADD TEXT TO VENN DIAGRAM ##########################################################
add.title <- function(
	gList,
	x,
	pos = c(0.5, 1.05),
	cex = 1,
	fontface = 'plain',
	fontfamily = 'serif',
	col = 'black',
	just = c(0.5, 1),
	...
	) {

	tmp <- textGrob(
		label = x,
		x = pos[1],
		y = pos[2],
		just = just,
		gp = gpar(
			col = col,
			cex = cex,
			fontface = fontface,
			fontfamily = fontfamily
			)
		);

	grob.list <- gList(gList, tmp);
	return(VennDiagram::adjust.venn(grob.list, ...))
	}
