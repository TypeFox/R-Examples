linkage.blocks.network <-
function(data.obj, geno.obj = NULL, collapse.linked.markers = TRUE, threshold.power = 1, plot.blocks = FALSE){
	
	net.data <- data.obj$var.to.var.p.val
	pheno.net.data <- data.obj$max.var.to.pheno.influence
	
	if(length(net.data) == 0){
		stop("calc.p() must be run to calculate variant-to-variant influences.")
		}
	
	full.geno <- get.geno.with.covar(data.obj, geno.obj)
	
	is.char <- as.logical(is.na(suppressWarnings(as.numeric(colnames(full.geno)[1]))))
			
	
	#find all the chromosomes that were used in the pairwise scan and sort them
	all.marker.chr <- get.marker.chr(data.obj, sort(unique(as.numeric(c(net.data[,1], net.data[,2])))))
	all.marker.names <- get.marker.name(data.obj, sort(unique(as.numeric(c(net.data[,1], net.data[,2])))))
	all.marker.num <- get.marker.num(data.obj, all.marker.names)
	
	u_chr <- sort(unique(as.numeric(all.marker.chr)))
	if(u_chr[1] == 0){u_chr <- c(u_chr[-1], 0)}
	
	#========================================================================================
	# internal functions
	#========================================================================================
	my.palette <- colorRampPalette(c("lightblue2", "green4"),space = "rgb")
	
	#if we are not collapsing the markers into blocks, 
	#or a chromosome only has one marker, or we are
	#adding the covariate chromosome, just add all 
	#individual markers to the list of blocks.
	add.ind.markers <- function(link.blocks, ch, chr.markers){
		chr.block.num <- 1
		if(is.null(link.blocks[[1]])){
			num.blocks = 1
			}else{
			num.blocks = length(link.blocks) + 1	
			}
		for(i in 1:length(chr.markers)){
			link.blocks[[num.blocks]] <- chr.markers[i]
			if(ch == 0){
				names(link.blocks)[num.blocks] <- chr.markers[i]
				}else{
				names(link.blocks)[num.blocks] <- paste("Chr", ch, "_", chr.block.num, sep = "")
				}
			num.blocks <- num.blocks + 1
			chr.block.num <- chr.block.num + 1
			}
		return(link.blocks)
		}
	

	#get the recombination data for the markers on a given chromosome
	get.chr.cor <- function(chr.markers){		
		just.markers <- unique(chr.markers)
		marker.locale <- match(just.markers, colnames(full.geno))
		chr.geno <- full.geno[,marker.locale]
		chr.cor <- cor(chr.geno, use = "pairwise.complete.obs")
		# chr.cor[lower.tri(chr.cor, diag = FALSE)] <- NA
		return(chr.cor)
		}
		

	
	add.chr.blocks <- function(link.blocks, new.blocks){
		chr.block.num <- 1
		if(is.null(link.blocks[[1]])){
			num.blocks = 1
			}else{
			num.blocks = length(link.blocks) + 1	
			}		
		for(i in 1:length(new.blocks)){
			link.blocks[[num.blocks]] <- new.blocks[[i]]
			names(link.blocks)[num.blocks] <- names(new.blocks)[i]
			num.blocks <- num.blocks + 1
			}
		return(link.blocks)
		}
	
	
	#========================================================================================
	# end internal functions
	#========================================================================================


	if(plot.blocks){pdf(paste("Recomb.Images.Genotype.Net.Thresh.", threshold.power, ".pdf", sep = ""), width = 10, height = 5)}
	#go through each chromosome separately and find the linkage blocks on each chromosome
	link.blocks <- vector(mode = "list", length = 1)
	num.blocks <- 1
	for(ch in u_chr){
		chr.blocks = 1
		if(is.char){
		chr.markers <- all.marker.names[which(as.numeric(all.marker.chr) == as.numeric(ch))]
		}else{
		chr.markers <- all.marker.num[which(as.numeric(all.marker.chr) == as.numeric(ch))]	
		}
		
		if(!collapse.linked.markers || length(chr.markers) == 1 || ch == 0){
			if(ch == 0){
			link.blocks <- add.ind.markers(link.blocks, ch, get.marker.name(data.obj, chr.markers))	
			}else{
			link.blocks <- add.ind.markers(link.blocks, ch, get.marker.num(data.obj, chr.markers))	
			}
			num.blocks <- num.blocks + 1
			}else{
			all.cor <- get.chr.cor(chr.markers)
			diag(all.cor) <- 0
			thresh.mat <- abs(all.cor^threshold.power)
			net <- graph.adjacency(thresh.mat, mode = "undirected", weighted = TRUE)
			comm <- fastgreedy.community(net)$membership
			
			# comm.num <- 1
			adj.comm <- consec.pairs(comm)
			cm.changes <- which(!apply(adj.comm, 1, function(x) x[1] == x[2])) #find everywhere the community number changes
			if(length(cm.changes) == 0){ #if there are no changes, put the whole chromosome into the block
				link.blocks[[num.blocks]] <- get.marker.num(data.obj, chr.markers)
				names(link.blocks)[num.blocks] <- paste("Chr", ch, "_", chr.blocks, sep = "")
				num.blocks <- num.blocks + 1
				}else{ #otherwise, step through the communities and add each one as a block
					chr.marker.locale <- which(all.marker.chr == ch)
					chr.marker.names <- all.marker.names[chr.marker.locale]

					for(cm in 1:(length(cm.changes)+1)){
						if(chr.blocks == 1){
							marker.names <- V(net)$name[1:cm.changes[cm]]
							}
						if(cm > length(cm.changes)){
							marker.names <- V(net)$name[(cm.changes[(cm-1)]+1):length(comm)]
							}
							
						if(cm <= length(cm.changes) && chr.blocks > 1){
							marker.names <- V(net)$name[(cm.changes[cm-1]+1):cm.changes[cm]]
							}
							
						link.blocks[[num.blocks]] <- get.marker.num(data.obj, marker.names)
						names(link.blocks)[num.blocks] <- paste("Chr", ch, "_", chr.blocks, sep = "")
						num.blocks <- num.blocks + 1
						chr.blocks <- chr.blocks + 1
						} #end looping through communities
					} #end adding blocks based on communities
				} #end case for if we are not dealing with the covariate chromosome or not collapasing markers


			if(plot.blocks && ch != 0){
										
											
						layout(matrix(c(1,2), nrow = 1))
						image(1:dim(all.cor)[1], 1:dim(all.cor)[2], all.cor, main = paste("Marker Correlation Chr", ch), xlim = c(0,(dim(all.cor)[1]+1)), ylim = c(0,(dim(all.cor)[1]+1)), col = my.palette(50), axes = FALSE, ylab = "", xlab = "")
						
						block.names <- sapply(strsplit(names(link.blocks), "_"), function(x) x[1])
						chr.names <- sapply(strsplit(block.names, "Chr"), function(x) x[2])
						final.block.locale <- which(chr.names == ch)
						start.block = 0.5
						#outline each block
						for(b in 1:length(final.block.locale)){
							end.block <- start.block + length(link.blocks[[final.block.locale[b]]])
							segments(x0 = start.block, y0 = start.block, x1 = start.block, y1 = end.block, lwd = 3)
							segments(x0 = start.block, y0 = start.block, x1 = end.block, y1 = start.block, lwd = 3)
							segments(x0 = end.block, y0 = start.block, x1 = end.block, y1 = end.block, lwd = 3)
							segments(x0 = start.block, y0 = end.block, x1 = end.block, y1 = end.block, lwd = 3)						
							start.block <- end.block
							} #end outlining blocks

						image(1:dim(thresh.mat)[1], 1:dim(thresh.mat)[2], thresh.mat, main = "Thresholded Correlations", xlim = c(0,(dim(thresh.mat)[1]+1)), ylim = c(0,(dim(thresh.mat)[1]+1)), col = my.palette(50), axes = FALSE, ylab = "", xlab = "")						
						block.names <- sapply(strsplit(names(link.blocks), "_"), function(x) x[1])
						chr.names <- sapply(strsplit(block.names, "Chr"), function(x) x[2])
						final.block.locale <- which(chr.names == ch)
						start.block = 0.5
						#outline each block
						for(b in 1:length(final.block.locale)){
							end.block <- start.block + length(link.blocks[[final.block.locale[b]]])
							segments(x0 = start.block, y0 = start.block, x1 = start.block, y1 = end.block, lwd = 3)
							segments(x0 = start.block, y0 = start.block, x1 = end.block, y1 = start.block, lwd = 3)
							segments(x0 = end.block, y0 = start.block, x1 = end.block, y1 = end.block, lwd = 3)
							segments(x0 = start.block, y0 = end.block, x1 = end.block, y1 = end.block, lwd = 3)						
							start.block <- end.block
							} #end outlining blocks		
							
			} #end plotting
	

			} #end looping through chromosomes	
	
	if(plot.blocks){dev.off()}
	
	if(collapse.linked.markers){
		data.obj$linkage.blocks.collapsed <- link.blocks
		}else{
		data.obj$linkage.blocks.full <- link.blocks
		}

	return(data.obj)
	
}
