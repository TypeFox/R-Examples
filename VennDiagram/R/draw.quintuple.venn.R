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

### FUNCTION TO DRAW VENN DIAGRAM WITH FIVE SETS #################################################
draw.quintuple.venn <- function(
	area1,
	area2,
	area3,
	area4,
	area5,
	n12,
	n13,
	n14,
	n15,
	n23,
	n24,
	n25,
	n34,
	n35,
	n45,
	n123,
	n124,
	n125,
	n134,
	n135,
	n145,
	n234,
	n235,
	n245,
	n345,
	n1234,
	n1235,
	n1245,
	n1345,
	n2345,
	n12345,
	category = rep("", 5),
	lwd = rep(2, 5),
	lty = rep("solid", 5),
	col = rep("black", 5),
	fill = NULL,
	alpha = rep(0.5, 5),
	label.col = rep("black", 31),
	cex = rep(1, 31),
	fontface = rep("plain", 31),
	fontfamily = rep("serif", 31),
	cat.pos = c(0, 287.5, 215, 145, 70),
	cat.dist = rep(0.2, 5),
	cat.col = rep("black", 5),
	cat.cex = rep(1, 5),
	cat.fontface = rep("plain", 5),
	cat.fontfamily = rep("serif", 5),
	cat.just = rep(list(c(0.5, 0.5)), 5),
	rotation.degree = 0,
	rotation.centre = c(0.5, 0.5),
	ind = TRUE,
	cex.prop=NULL,
	print.mode = "raw",
    sigdigs=3,
    direct.area=FALSE,
    area.vector=0,
	...
	) {

	#area1 > area2 > area3 > area4 > area5
	# check parameter lengths
	if (length(category) == 1) { cat <- rep(category, 5); }
	else if (length(category) != 5) { flog.error("Unexpected parameter length for 'category'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'category'"); }

	if (length(lwd) == 1) { lwd <- rep(lwd, 5); }
	else if (length(lwd) != 5) { flog.error("Unexpected parameter length for 'lwd'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'lwd'"); }

	if (length(lty) == 1) { lty <- rep(lty, 5); }
	else if (length(lty) != 5) { flog.error("Unexpected parameter length for 'lty'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'lty'"); }

	if (length(col) == 1) { col <- rep(col, 5); }
	else if (length(col) != 5) { flog.error("Unexpected parameter length for 'col'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'col'"); }

	if (length(label.col) == 1) { label.col <- rep(label.col, 31); }
	else if (length(label.col) != 31) { flog.error("Unexpected parameter length for 'label.col'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'label.col'"); }

	if (length(cex) == 1) { cex <- rep(cex, 31); }
	else if (length(cex) != 31) { flog.error("Unexpected parameter length for 'cex'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'cex'"); }

	if (length(fontface) == 1) { fontface <- rep(fontface, 31); }
	else if (length(fontface) != 31) { flog.error("Unexpected parameter length for 'fontface'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'fontface'"); }

	if (length(fontfamily) == 1) { fontfamily <- rep(fontfamily, 31); }
	else if (length(fontfamily) != 31) { flog.error("Unexpected parameter length for 'fontfamily'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'fontfamily'"); }

	if (length(fill) == 1) { fill <- rep(fill, 5); }
	else if (length(fill) != 5 & length(fill) != 0) { flog.error("Unexpected parameter length for 'fill'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'fill'"); }

	if (length(alpha) == 1) { alpha <- rep(alpha, 5); }
	else if (length(alpha) != 5 & length(alpha) != 0) { flog.error("Unexpected parameter length for 'alpha'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'alpha'"); }

	if (length(cat.pos) == 1) { cat.pos <- rep(cat.pos, 5); }
	else if (length(cat.pos) != 5) { flog.error("Unexpected parameter length for 'cat.pos'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'cat.pos'"); }

	if (length(cat.dist) == 1) { cat.dist <- rep(cat.dist, 5); }
	else if (length(cat.dist) != 5) { flog.error("Unexpected parameter length for 'cat.dist'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'cat.dist'"); }

	if (length(cat.col) == 1) { cat.col <- rep(cat.col, 5); }
	else if (length(cat.col) != 5) { flog.error("Unexpected parameter length for 'cat.col'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'cat.col'"); }

	if (length(cat.cex) == 1) { cat.cex <- rep(cat.cex, 5); }
	else if (length(cat.cex) != 5) { flog.error("Unexpected parameter length for 'cat.cex'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'cat.cex'"); }

	if (length(cat.fontface) == 1) { cat.fontface <- rep(cat.fontface, 5); }
	else if (length(cat.fontface) != 5) { flog.error("Unexpected parameter length for 'cat.fontface'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'cat.fontface'"); }

	if (length(cat.fontfamily) == 1) { cat.fontfamily <- rep(cat.fontfamily, 5); }
	else if (length(cat.fontfamily) != 5) { flog.error("Unexpected parameter length for 'cat.fontfamily'",name="VennDiagramLogger")
stop("Unexpected parameter length for 'cat.fontfamily'"); }

	if (!(class(cat.just) == "list" & length(cat.just) == 5 & length(cat.just[[1]]) == 2 & length(cat.just[[2]]) == 2 & length(cat.just[[3]]) == 2 & length(cat.just[[4]]) == 2 & length(cat.just[[5]]) == 2)) { flog.error("Unexpected parameter format for 'cat.just'",name="VennDiagramLogger")
stop("Unexpected parameter format for 'cat.just'"); }
	cat.pos <- cat.pos + rotation.degree;
	
	if(direct.area){
		areas <- area.vector;
		#create the variables and assign their values from the area vector
		for(i in 1:31)
		{
			assign(paste("a",i,sep=""),area.vector[i]);
		}
	}
	else{
		# generate partial areas from given arguments
		a31 <- n12345;
		a30 <- n1234 - a31;
		a29 <- n1235 - a31;
		a28 <- n1245 - a31;
		a27 <- n1345 - a31;
		a26 <- n2345 - a31;
		a25 <- n245 - a26 - a28 - a31;
		a24 <- n234 - a26 - a30 - a31;
		a23 <- n134 - a27 - a30 - a31;
		a22 <- n123 - a29 - a30 - a31;
		a21 <- n235 - a26 - a29 - a31;
		a20 <- n125 - a28 - a29 - a31;
		a19 <- n124 - a28 - a30 - a31;
		a18 <- n145 - a27 - a28 - a31;
		a17 <- n135 - a27 - a29 - a31;
		a16 <- n345 - a26 - a27 - a31;
		a15 <- n45 - a18 - a25 - a16 - a28 - a27 - a26 - a31;
		a14 <- n24 - a19 - a24 - a25 - a30 - a28 - a26 - a31;
		a13 <- n34 - a16 - a23 - a24 - a26 - a27 - a30 - a31;
		a12 <- n13 - a17 - a22 - a23 - a27 - a29 - a30 - a31;
		a11 <- n23 - a21 - a22 - a24 - a26 - a29 - a30 - a31;
		a10 <- n25 - a20 - a21 - a25 - a26 - a28 - a29 - a31;
		a9 <- n12 - a19 - a20 - a22 - a28 - a29 - a30 - a31;
		a8 <- n14 - a18 - a19 - a23 - a27 - a28 - a30 - a31;
		a7 <- n15 - a17 - a18 - a20 - a27 - a28 - a29 - a31;
		a6 <- n35 - a16 - a17 - a21 - a26 - a27 - a29 - a31;
		a5 <- area5 - a6 - a7 - a15 - a16 - a17 - a18 - a25 - a26 - a27 - a28 - a31 - a20 - a29 - a21 - a10;
		a4 <- area4 - a13 - a14 - a15 - a16 - a23 - a24 - a25 - a26 - a27 - a28 - a31 - a18 - a19 - a8 - a30;
		a3 <- area3 - a21 - a11 - a12 - a13 - a29 - a22 - a23 - a24 - a30 - a31 - a26 - a27 - a16 - a6 - a17;
		a2 <- area2 - a9 - a10 - a19 - a20 - a21 - a11 - a28 - a29 - a31 - a22 - a30 - a26 - a25 - a24 - a14;
		a1 <- area1 - a7 - a8 - a18 - a17 - a19 - a9 - a27 - a28 - a31 - a20 - a30 - a29 - a22 - a23 - a12;

		# check plausibility and 0 partial areas
		areas <- c(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, a24, a25, a26, a27, a28, a29, a30, a31);
	}
	areas.error <- c(
		"a1 <- area1 - a7 - a8 - a18 - a17 - a19 - a9 - a27 - a28 - a31 - a20 - a30 - a29 - a22 - a23 - a12",
		"a2 <- area2 - a9 - a10 - a19 - a20 - a21 - a11 - a28 - a29 - a31 - a22 - a30 - a26 - a25 - a24 - a14",
		"a3 <- area3 - a21 - a11 - a12 - a13 - a29 - a22 - a23 - a24 - a30 - a31 - a26 - a27 - a16 - a6 - a17",
		"a4 <- area4 - a13 - a14 - a15 - a16 - a23 - a24 - a25 - a26 - a27 - a28 - a31 - a18 - a19 - a8 - a30",
		"a5 <- area5 - a6 - a7 - a15 - a16 - a17 - a18 - a25 - a26 - a27 - a28 - a31 - a20 - a29 - a21 - a10",
		"a6 <- n35 - a16 - a17 - a21 - a26 - a27 - a29 - a31",
		"a7 <- n15 - a17 - a18 - a20 - a27 - a28 - a29 - a31",
		"a8 <- n14 - a18 - a19 - a23 - a27 - a28 - a30 - a31",
		"a9 <- n12 - a19 - a20 - a22 - a28 - a29 - a30 - a31",
		"a10 <- n25 - a20 - a21 - a25 - a26 - a28 - a29 - a31",
		"a11 <- n23 - a21 - a22 - a24 - a26 - a29 - a30 - a31",
		"a12 <- n13 - a17 - a22 - a23 - a27 - a29 - a30 - a31",
		"a13 <- n34 - a16 - a23 - a24 - a26 - a27 - a30 - a31",
		"a14 <- n24 - a19 - a24 - a25 - a30 - a28 - a26 - a31",
		"a15 <- n45 - a18 - a25 - a16 - a28 - a27 - a26 - a31",
		"a16 <- n345 - a26 - a27 - a31",
		"a17 <- n135 - a27 - a29 - a31",
		"a18 <- n145 - a27 - a28 - a31",
		"a19 <- n124 - a28 - a30 - a31",
		"a20 <- n125 - a28 - a29 - a31",
		"a21 <- n235 - a26 - a29 - a31",
		"a22 <- n123 - a29 - a30 - a31",
		"a23 <- n134 - a27 - a30 - a31",
		"a24 <- n234 - a26 - a30 - a31",
		"a25 <- n245 - a26 - a28 - a31",
		"a26 <- n2345 - a31",
		"a27 <- n1345 - a31",
		"a28 <- n1245 - a31",
		"a29 <- n1235 - a31",
		"a30 <- n1234 - a31",
		"a31 <- n12345"
		);
	for (i in 1:length(areas)) {
		if (areas[i] < 0) {
			flog.error(paste("Impossible:", areas.error[i], "produces negative area"),name="VennDiagramLogger")
stop(paste("Impossible:", areas.error[i], "produces negative area"));
			}
		}

        ## rescaling area labels to be proportional to area
        if(length(cex.prop) > 0){

            if(length(cex.prop) != 1) flog.error("Value passed to cex.prop is not length 1",name="VennDiagramLogger")
stop("Value passed to cex.prop is not length 1")

            ## figure out what function to use
            func = cex.prop
            if(class(cex.prop) != "function"){
                if(cex.prop == "lin"){
                    func = function(x) x
                }
                else if(cex.prop == "log10"){
                    func = log10
                }
                else flog.error(paste0("Unknown value passed to cex.prop: ", cex.prop),name="VennDiagramLogger")
stop(paste0("Unknown value passed to cex.prop: ", cex.prop))
            }

            ## rescale areas
            maxArea = max(areas)            
            for(i in 1:length(areas)){                
                cex[i] = cex[i] * func(areas[i]) / func(maxArea)
                if(cex[i] <= 0) stop(paste0("Error in rescaling of area labels: the label of area ",
                          i, " is less than or equal to zero"))
            }
        }


	# initialize gList to hold all Grobs generated
	grob.list <- gList();

	# plot the ellipses of the Venn diagram
	dist <- 0.13;
	a <- 0.24;
	b <- 0.46;
	init.angle <- -20;

	ellipse.positions <- matrix(
		nrow = 5,
		ncol = 3
		);
	colnames(ellipse.positions) <- c('x', 'y', 'rotation');

	ellipse.positions[1,] <- c(
		0.5 + dist * sin(init.angle * pi / 180),
		0.5 + dist * cos(init.angle * pi / 180),
		0
		);

	ellipse.positions[2,] <- c(
		0.5 - dist * cos((288 + init.angle - 270) * pi / 180),
		0.5 + dist * sin((288 + init.angle - 270) * pi / 180),
		-110
		);

	ellipse.positions[3,] <- c(
		0.5 - dist * sin((216 + init.angle - 180) * pi / 180),
		0.5 - dist * cos((216 + init.angle - 180) * pi / 180),
		145
		);

	ellipse.positions[4,] <- c(
		0.5 + dist * sin((180 - 144 - init.angle) * pi / 180),
		0.5 - dist * cos((180 - 144 - init.angle) * pi / 180),
		35
		);

	ellipse.positions[5,] <- c(
		0.5 + dist * cos((init.angle + 72 - 90) * pi / 180),
		0.5 - dist * sin((init.angle + 72 - 90) * pi / 180),
		-72.5
		);

	# create the ellipses
	for (i in 1:5) {
		grob.list <- gList(
			grob.list,
			VennDiagram::ellipse(
				x = ellipse.positions[i,'x'],
				y = ellipse.positions[i,'y'],
				a = a,
				b = b,
				rotation = ellipse.positions[i,'rotation'],
				gp = gpar(
					lty = 0,
					fill = fill[i],
					alpha = alpha[i]
					)
				)
			);
		}

	# create the ellipse borders
	for (i in 1:5) {
		grob.list <- gList(
			grob.list,
			VennDiagram::ellipse(
				x = ellipse.positions[i,'x'],
				y = ellipse.positions[i,'y'],
				a = a,
				b = b,
				rotation = ellipse.positions[i,'rotation'],
				gp = gpar(
					lwd = lwd[i],
					lty = lty[i],
					col = col[i],
					fill = 'transparent'
					)
				)
			);
		}

	# add area labels
	label.matrix <- matrix(
		nrow = 31,
		ncol = 3
		);
	colnames(label.matrix) <- c('label', 'x', 'y');

	label.matrix[ 1,] <- c(a1,  0.4555, 0.9322);
	label.matrix[ 2,] <- c(a2,  0.0800, 0.6000);
	label.matrix[ 3,] <- c(a3,  0.3000, 0.1000);
	label.matrix[ 4,] <- c(a4,  0.7900, 0.1700);
	label.matrix[ 5,] <- c(a5,  0.9000, 0.6800);
	label.matrix[ 6,] <- c(a6,  0.7400, 0.6950);
	label.matrix[ 7,] <- c(a7,  0.6300, 0.8050);
	label.matrix[ 8,] <- c(a8,  0.4000, 0.7950);
	label.matrix[ 9,] <- c(a9,  0.2550, 0.7150);
	label.matrix[10,] <- c(a10, 0.1930, 0.4800);
	label.matrix[11,] <- c(a11, 0.2250, 0.3330);
	label.matrix[12,] <- c(a12, 0.4200, 0.2050);
	label.matrix[13,] <- c(a13, 0.5720, 0.1800);
	label.matrix[14,] <- c(a14, 0.7530, 0.3200);
	label.matrix[15,] <- c(a15, 0.8230, 0.4700);
	label.matrix[16,] <- c(a16, 0.7470, 0.5820);
	label.matrix[17,] <- c(a17, 0.6620, 0.7500);
	label.matrix[18,] <- c(a18, 0.4880, 0.7610);
	label.matrix[19,] <- c(a19, 0.3230, 0.7370);
	label.matrix[20,] <- c(a20, 0.2530, 0.5730);
	label.matrix[21,] <- c(a21, 0.2250, 0.3950);
	label.matrix[22,] <- c(a22, 0.3550, 0.2900);
	label.matrix[23,] <- c(a23, 0.5150, 0.2050);
	label.matrix[24,] <- c(a24, 0.6550, 0.2900);
	label.matrix[25,] <- c(a25, 0.7830, 0.4200);
	label.matrix[26,] <- c(a26, 0.7200, 0.4450);
	label.matrix[27,] <- c(a27, 0.6050, 0.7010);
	label.matrix[28,] <- c(a28, 0.3420, 0.6680);
	label.matrix[29,] <- c(a29, 0.2940, 0.4100);
	label.matrix[30,] <- c(a30, 0.5220, 0.2730);
	label.matrix[31,] <- c(a31, 0.5000, 0.5000);

	processedLabels <- rep("",length(label.matrix[,'label']));
    if(print.mode[1] == "percent"){
			processedLabels <- paste(signif(label.matrix[,'label']/sum(label.matrix[,'label'])*100,digits=sigdigs),"%",sep="");
			if(isTRUE(print.mode[2] == "raw"))
			{
				processedLabels <- paste(processedLabels,"\n(",label.matrix[,'label'],")",sep="");
			}
		}
	if(print.mode[1] == "raw"){
			processedLabels <- label.matrix[,'label'];
			if(isTRUE(print.mode[2] == "percent"))
			{
				processedLabels <- paste(processedLabels,"\n(",paste(signif(label.matrix[,'label']/sum(label.matrix[,'label'])*100,digits=sigdigs),"%)",sep=""),sep="");
			}
		}

	for (i in 1:nrow(label.matrix)) {
		tmp <- textGrob(
			label = processedLabels[i],
			x = label.matrix[i,'x'],
			y = label.matrix[i,'y'],
			gp = gpar(
				col = label.col[i],
				cex = cex[i],
				fontface = fontface[i],
				fontfamily = fontfamily[i]
				)
			);
		grob.list <- gList(grob.list, tmp);
		}


	# find the location and plot all the category names
	cat.pos.x <- c(0.4555, 0.08, 0.3, 0.79, 0.90);
	cat.pos.y <- c(0.9322, 0.60, 0.1, 0.17, 0.68);

	for (i in 1:5) {

		# find the label position
		this.cat.pos <- find.cat.pos(
			x = cat.pos.x[i],
			y = cat.pos.y[i],
			pos = cat.pos[i],
			dist = cat.dist[i]
			);

		# and plot it
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

	#if (exists("margin")) {margin <- margin + 0.05} else {margin <- 0.05}

	# adjust grob.list to fit and return grob.list
	grob.list <- VennDiagram::adjust.venn(VennDiagram::rotate.venn.degrees(grob.list, rotation.degree, rotation.centre[1], rotation.centre[2]), ...);
	if (ind) { grid.draw(grob.list); }
	return(grob.list);
	}
