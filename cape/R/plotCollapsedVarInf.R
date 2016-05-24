plotCollapsedVarInf <-
function(data.obj, expand.labels = FALSE, all.markers = FALSE, scale.effects = c("none", "log10", "sqrt")){
	
	adj.mat <- data.obj$collapsed.net
	if(is.null(adj.mat)){
		stop("This function operates on the collapsed network. get.network() with collapsed.linked.markers = TRUE must be run first.")
		}
	
	if(length(grep("n", scale.effects)) > 0){
		scale.effects <- "none"
		}
	if(length(scale.effects) == 1){
		if(scale.effects != "log10" & scale.effects != "sqrt" & scale.effects != "none"){
			stop("scale.effects must be 'log10', 'sqrt' or 'none.'")
			}
		}

	
	blocks <- data.obj$linkage.blocks.collapsed
	

	#replace marker numbers with names
	get.marker.names <- function(block){
		marker.locale <- which(colnames(data.obj$geno) %in% as.character(block))
		marker.names <- data.obj$marker.names[marker.locale]
		return(marker.names)
		}


	if(!all.markers){	

		adj.mat[which(adj.mat == 0)] <- NA
		
		if(expand.labels){
			named.blocks <- lapply(blocks, get.marker.names)			
			marker.names <- sapply(named.blocks, function(x) paste(x, collapse = ", "))
			pheno.names <- names(data.obj$max.var.to.pheno.influence)
			rownames(adj.mat) <- marker.names
			colnames(adj.mat) <- c(marker.names, pheno.names)
			}

		main <- "Condensed Variant Influences"
	
		if(scale.effects == "log10"){
			neg.locale <- which(adj.mat < 0)
			scaled.effects <- log10(abs(adj.mat))
			scaled.effects[neg.locale] <- scaled.effects[neg.locale]*-1	
			adj.mat <- scaled.effects
			main <- "log10 of Condensed Variant Influences"
			}
		if(scale.effects == "sqrt"){
			neg.locale <- which(adj.mat < 0)
			scaled.effects <- sqrt(abs(adj.mat))
			scaled.effects[neg.locale] <- scaled.effects[neg.locale]*-1	
			adj.mat <- scaled.effects
			main <- "Square Root of Condensed Variant Influences"
			}


		myImagePlot(adj.mat, min.x = (max(abs(adj.mat), na.rm = TRUE)*-1), max.x = max(abs(adj.mat), na.rm = TRUE), main = "Condensed Variant Influences", xlab = "Target", ylab = "Source")

		}else{ #if we are including all markers, we need to expand the adjacency matrix
			all.markers <- colnames(data.obj$geno)
			pheno.names <- names(data.obj$max.var.to.pheno.influence)
			expanded.adj.mat <- matrix(0, ncol = length(all.markers), nrow = length(all.markers))
			colnames(expanded.adj.mat) <- rownames(expanded.adj.mat) <- all.markers
			expanded.pheno.mat <- matrix(0, nrow = length(all.markers), ncol = length(pheno.names))
			colnames(expanded.pheno.mat) <- pheno.names; rownames(expanded.pheno.mat) <- all.markers
			expanded.adj.mat <- cbind(expanded.adj.mat, expanded.pheno.mat)
					
			#go through the blocks. Replace markers in blocks with the block number
			extra.markers.in.blocks <- sapply(blocks, length)
			big.block.locale <- which(extra.markers.in.blocks > 1)
			if(length(big.block.locale) > 0){
				for(i in 1:length(big.block.locale)){
					block.markers <- blocks[[big.block.locale[i]]]
					marker.locale <- which(colnames(expanded.adj.mat) %in% block.markers)
					col.to.remove <- marker.locale[-1]
					expanded.adj.mat <- expanded.adj.mat[,-col.to.remove]
					expanded.adj.mat <- expanded.adj.mat[-col.to.remove,]
					colnames(expanded.adj.mat)[marker.locale[1]] <- names(big.block.locale)[i]
					rownames(expanded.adj.mat)[marker.locale[1]] <- names(big.block.locale)[i]
					}
				}
			
							
			small.block.locale <- which(extra.markers.in.blocks == 1)		
			small.block.markers <- sapply(blocks, function(x) x[[1]])[small.block.locale]
			exp.block.locale <- which(rownames(expanded.adj.mat) %in% small.block.markers)
			rownames(expanded.adj.mat)[exp.block.locale] <- colnames(expanded.adj.mat)[exp.block.locale] <- names(small.block.markers)
			
			get.marker.name <- function(marker.number){
				marker.locale <- which(colnames(data.obj$geno) == marker.number)
				return(data.obj$marker.names[marker.locale])
				}
				
			#rename the rest of the markers with their names
			marker.rows <- rownames(expanded.adj.mat)
			block.locale <- suppressWarnings(which(is.na(as.numeric(marker.rows))))
			not.block <- suppressWarnings(which(!is.na(as.numeric(marker.rows))))
			not.block.locale <- marker.rows[not.block]
			row.marker.names <- apply(matrix(not.block, nrow = 1), 2, get.marker.name)
			marker.rows[not.block] <- row.marker.names

			rownames(expanded.adj.mat) <- marker.rows
			colnames(expanded.adj.mat)[1:length(marker.rows)] <- marker.rows

			
			sig.locale <- which(adj.mat != 0, arr.ind = TRUE)
			if(length(sig.locale) > 0){
				for(i in 1:length(sig.locale[,1])){
					location <- sig.locale[i,]
					adj.rowname <- rownames(adj.mat)[location[1]]
					adj.colname <- colnames(adj.mat)[location[2]]
					expanded.adj.mat[as.character(adj.rowname), as.character(adj.colname)] <- adj.mat[as.numeric(location[1]), as.numeric(location[2])]
					}
				}
						
		if(expand.labels){
			block.locale <- which(colnames(expanded.adj.mat) %in% names(blocks))
			if(length(block.locale) > 0){
				for(i in 1:length(block.locale)){
					block.name <- colnames(expanded.adj.mat)[block.locale[i]]
					block.markers <- blocks[[block.name]]
					block.marker.locale <- which(colnames(data.obj$geno) %in% block.markers)
					block.marker.name <- paste(data.obj$marker.names[block.marker.locale], collapse = ", ")
					colnames(expanded.adj.mat)[block.locale[i]] <- block.marker.name
					rownames(expanded.adj.mat)[block.locale[i]] <- block.marker.name
					}
				}	
			}
			
		expanded.adj.mat[which(expanded.adj.mat == 0)] <- NA
		
			main <- "Condensed Variant Influences"
	
	if(scale.effects == "log10"){
		neg.locale <- which(expanded.adj.mat < 0)
		scaled.effects <- log10(abs(expanded.adj.mat))
		scaled.effects[neg.locale] <- scaled.effects[neg.locale]*-1	
		expanded.adj.mat <- scaled.effects
		main <- "log10 of Condensed Variant Influences"
		}
	if(scale.effects == "sqrt"){
		neg.locale <- which(expanded.adj.mat < 0)
		scaled.effects <- sqrt(abs(expanded.adj.mat))
		scaled.effects[neg.locale] <- scaled.effects[neg.locale]*-1	
		expanded.adj.mat <- scaled.effects
		main <- "Square Root of Condensed Variant Influences"
		}

		myImagePlot(expanded.adj.mat, min.x = (max(abs(expanded.adj.mat), na.rm = TRUE)*-1), max.x = max(abs(expanded.adj.mat), na.rm = TRUE), main = "Condensed Variant Influences", xlab = "Target", ylab = "Source")

			
		}
	
	
	}
