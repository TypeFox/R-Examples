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

bedr.plot.region <- function(input, filename = NULL, type = "venn", feature = "interval", fraction.overlap = 1e-9, group = NULL, params = list(), verbose = TRUE) {
	# Venn, Box or BarPlot, heatmap, circos 
	# feature is region/bp for intersect
	# https://groups.google.com/forum/#!msg/bedtools-discuss/93n2WUP0MME/HvCN7kyjroEJ

	input.overlap <- list();
	main <- "";
	sub  <- "";
	n <- length(input);
	
	# plot venn diagrams for region intersections 
	if (type == "venn") {

		# if only one dataset compare with a collapsed/merged
		if (n == 1) {
			input[[1]]   <- bedr.sort.region(input[[1]]);
			input[[2]]   <- bedr.merge.region(input[[1]], number = TRUE);
			last.var     <- names(input[[2]])[ncol(input[[2]])]
			input[[2]]   <- input[[2]][input[[2]][,last.var] >  1,];
			names(input) <- c("Original","Merged");
			main <- "Merged Intervals"
			sub <- paste0("n.before=", nrow(as.data.frame(input[[1]])), "  n.after=",nrow(as.data.frame(input[[2]])));
			n <- 2;
			}

		# map input to genes before doing overlap
		if (feature == "gene") { # gene overlap
			for (i in 1:n) {
				#if (!exists(refgene)) {
					gene.region <- query.ucsc("refGene");
				#	}
				input.overlap[[names(input)[i]]] <- gene.region[in.region(gene.region, input[[i]]), "symbol"];
				}
			}
		else if (feature == "reference") { # map input to reference before doing overlap.  reference being the first set of regions.
			if (n <= 2) {catv("Plot.region: you need at least three regions to plot by reference\n"); stop()}
			for (i in 2:n) {
				input.overlap[[names(input)[i]]] <- bed2index(input[[1]][in.region(input[[1]], input[[i]])]);
				}
			}
		else if (feature == "interval") { # interval overlap
			input.overlap <- bedr.join.multiple.region(input);
			}
		else if (feature == "cluster") { # use clustering algorithm to identify overlaps
			input.overlap <- bedr.join.multiple.region(input, cluster = TRUE);
			}
		else if (feature == "bp") {
			input.overlap <- bedr.join.multiple.region(input);
			size          <- size.region(input.overlap);
			}
		else {
			catv("Plot.regin: unrecognized feature type.\n")
			stop();
			}

		# add names of input to output
		names(input.overlap)[(ncol(input.overlap)-n+1):ncol(input.overlap)] <- names(input);

		# summarize counts for each of the intersections 
		if (feature == "bp") {
			# size by intersection
			overlap.summary <- tapply(size, input.overlap[,ncol(input.overlap)-n], sum);
			}
		else {
			# convert df into list for venn plotting
			overlap.summary <- table(input.overlap[ncol(input.overlap)-n]);
			# don't use the list of labels input
			#input.list <- df2list(input.overlap[,names(input.overlap) %in% names(input)]);
			}
		
		# convert table of counts to venn format	
		venn.input <- table2venn(x = overlap.summary, var.names = names(input)); 
	
		# setup plotting defaults
		venn.2   <- alist(
			fill = c("cornflowerblue", "yellow"),
			cat.col = c("cornflowerblue", "yellow"),
			cat.dist = c(0.03, 0.03),
			cat.pos = c(-20, 14),
			cat.cex = 3,
			cex = 4,
			euler.d = TRUE,
			scaled = TRUE
			);

		venn.3   <- alist(
			fill = c("red", "blue", "green"),
			label.col = c("darkred", "white", "darkblue", "white", "white", "white", "darkgreen"),
			cat.col = c("darkred", "darkblue", "darkgreen"),
			cat.dist = c(0.06, 0.06, 0.03),
			cat.pos = 0,
			euler.d = FALSE,
			scaled = FALSE
			);

		venn.4   <- alist(
			fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
			label.col = c("orange", "white", "darkorchid4", "white", "white", "white", "white", "white", "darkblue", "white", "white", "white", "white", "darkgreen", "white"),
			cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4")
			);

		venn.5   <- alist(
			fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
			cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
			cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),
			cat.cex = 1.5,
			margin = 0.05
			);
		
		if (n == 1) {
			venn.n <- venn.2;
			venn.type <- "pairwise";
			}	
		else if (n == 2) {
			venn.n <- venn.2;
			venn.type <- "pairwise";
			}
		else if (n == 3) {
			venn.n <- venn.3;
			venn.type <- "triple";
			}
		else if (n == 4) {
			venn.n <- venn.4;
			venn.type <- "quad"
			}
		else if (n == 5) {
			venn.n <- venn.5;
			venn.type <- "quintuple"
			}
		else {
			catv("Plot: Max 5 beds for a Venn\n")
			}
		
		# decide what venn function
		venn.function <- paste("draw", venn.type, "venn", sep = ".");
		
		# set up custom plotting parameters
		venn.default <- formals(venn.function);
		
		venn.new   <- alist(
			filename = NULL,
			category = names(input),
			col = "black",
			#lty = "dashed",
			lwd = 2,
			alpha = 0.50,
			cex = 2.5,
			fontfamily = "serif",
			fontface = "bold",
			cat.fontface = "bold",
			cat.cex = 2.5,
			cat.fontfamily = "serif",
			main = main,
			main.cex = 3,
			sub = sub,
			sub.cex = 1.2
			);

		# update with any user specified venn plotting customizations
		# venn.new <- modifyList(venn.new, ...)

		venn.new['filename'] <- list(filename);

		# create the paramater list.  modifyList2 allows assignment of NULL without deletion
		
		# start with generic defaults
		venn.new  <-  modifyList2(venn.default, venn.new);
		
		# then add defaults specific to size of the input
		venn.new  <-  modifyList2(venn.new, venn.n);
		venn.new  <-  modifyList2(venn.new, venn.input);
		
		# include additonal plotting parameters specified by the calling function
		venn.new  <-  modifyList2(venn.new, params);
		
		# dump the ellipses
		venn.new[["..."]] <- NULL;
		#formals(venn.diagram) <- venn.new;
	
		# try and do.call
		grid.newpage();
		venn <- do.call(venn.function, venn.new);

		#grid.draw(venn);
		} 
	else if ("barplot") {
		
		}
	else if ("boxplot") {
		
		}
	else if ("circos") {
		}
		
	return(invisible(input.overlap))
	}

