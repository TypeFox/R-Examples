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

### FUNCTION TO DRAW SPECIAL CASES ################################################################
draw.sp.case <- function(
	area.list,
	enabled.areas,
	area.x,
	area.y,
	attach.label.to,
	x.centres,
	y.centres,
	a.list,
	b.list,
	straight.reverse,
	reverse = FALSE,
	category,
	cat.default.pos = "outer",
	lwd = rep(2, 3),
	lty = rep("solid", 3),
	col = rep("black", 3),
	label.col = rep("black", 7),
	cex = rep(1, 7),
	fontface = rep("plain", 7),
	fontfamily = rep("serif", 7),
	cat.pos = c(-40, 40, 0),
	cat.dist = c(0.05, 0.05, 0.025),
	cat.col = rep("black", 3),
	cat.cex = rep(1, 3),
	cat.fontface = rep("plain", 3),
	cat.fontfamily = rep("serif", 3),
	cat.just = list(c(0.5, 1), c(0.5, 1), c(0.5, 0)),
	cat.prompts = FALSE,
	fill = NULL,
	alpha = rep(0.5, 3),
	print.mode = "raw",
    sigdigs=3,
	...
	) {

	grob.list <- gList();

	# create the three ellipses
	for (i in 1:3) {
		grob.list <- gList(
			grob.list,
			VennDiagram::ellipse(
				x = x.centres[i],
				y = y.centres[i],
				a = a.list[i],
				b = b.list[i],
				gp = gpar(
					lty = 0,
					fill = fill[i],
					alpha = alpha[i]
					)
				)
			);
		}

	# create the three ellipse-borders
	for (i in 1:3) {
		grob.list <- gList(
			grob.list,
			VennDiagram::ellipse(
				x = x.centres[i],
				y = y.centres[i],
				a = a.list[i],
				b = b.list[i],
				gp = gpar(
					lwd = lwd[i],
					lty = lty[i],
					col = col[i],
					fill = 'transparent'
					)
				)
			);
		}

	# add the text labels
	# make it percents if it is enabled
	# else give the count number
	processedLabels <- rep("",length(area.list));
    if(print.mode[1] == "percent"){
			processedLabels <- paste(signif(area.list/sum(area.list)*100,digits=sigdigs),"%",sep="");
			if(isTRUE(print.mode[2] == "raw"))
			{
				processedLabels <- paste(processedLabels,"\n(",area.list,")",sep="");
			}
		}
	if(print.mode[1] == "raw"){
			processedLabels <- area.list;
			if(isTRUE(print.mode[2] == "percent"))
			{
				processedLabels <- paste(processedLabels,"\n(",paste(signif(area.list/sum(area.list)*100,digits=sigdigs),"%)",sep=""),sep="");
			}
		}
	
	for (i in 1:7) {
			if (i %in% enabled.areas) {
				grob.list <- gList(
					grob.list,
					textGrob(
						label = processedLabels[i],
						x = area.x[i],
						y = area.y[i],
						just = c('centre', 'centre'),
						gp = gpar(
							col = label.col[i],
							cex = cex[i],
							fontface = fontface[i],
							fontfamily = fontfamily[i]
							)
						)
					);
				}
			}
		

	# create category labels
	for (i in 1:3) {

		# try to auto-assign category position labels
		if ('outer' == cat.default.pos) {
			this.cat.pos <- VennDiagram::find.cat.pos(
				x = x.centres[i],
				y = y.centres[i],
				pos = cat.pos[i],
				dist = cat.dist[i],
				r = a.list[i]
				);
			}
		else if ('text' == cat.default.pos) {
			this.cat.pos <- VennDiagram::find.cat.pos(
				x = area.x[attach.label.to[i]],
				y = area.y[attach.label.to[i]],
				pos = cat.pos[i],
				dist = cat.dist[i]
				);
			}
		else {
			flog.error('Invalid cat.default.pos setting',name="VennDiagramLogger")
stop('Invalid cat.default.pos setting');
			}

		# create the label
		grob.list <- gList(
			grob.list,
			textGrob(
				label = category[i],
				x = this.cat.pos$x,
				y = this.cat.pos$y,
				just = cat.just[[i]],
				gp = gpar(
					col = cat.col[i],
					cex = cat.cex[i],
					fontface = cat.fontface[i],
					fontfamily = cat.fontfamily[i]
					)
				)
			);
		}

	# fit Venn diagram to size
	grob.list <- VennDiagram::adjust.venn(grob.list, ...);

	if (straight.reverse) {
		if (reverse) {
			return(VennDiagram::flip.venn(grob.list, axis = 'v'));
			}
		}
	return(grob.list);
	}
