get.network <-
function(data.obj, p.or.q = 0.05, standardized = TRUE, min.std.effect = 0, collapse.linked.markers = TRUE, verbose = FALSE, plot.linkage.blocks = FALSE, threshold.power = 1){
	
	# linkage.method = c("genotype.comm", "genotype.step", "effects", "prominence")
	# r.thresh = 0.5
	# effect.drop = 1

	linkage.method = "genotype.comm"
	
	
	#get the linkage blocks based on the significant markers
	# data.obj <- linkage.blocks(data.obj, collapse.linked.markers = collapse.linked.markers, threshold.power = threshold.power, verbose = verbose, plot.results = plot.linkage.blocks)
	# method.check <- grep("comm", linkage.method)
	# if(length(method.check) > 0){linkage.method <- "genotype.comm"}
	# if(linkage.method == "effects"){
		# data.obj <- linkage.blocks.cormat(data.obj, p.or.q = p.or.q, effect.drop = effect.drop, collapse.linked.markers = collapse.linked.markers, threshold.power = threshold.power, verbose = verbose, plot.results = plot.linkage.blocks)
		# }
	# if(linkage.method == "genotype.step"){
		# data.obj <- linkage.blocks.stepwise(data.obj, collapse.linked.markers = collapse.linked.markers, r.thresh = r.thresh, plot.blocks = plot.linkage.blocks)
		# }
	# if(linkage.method == "genotype.comm"){
		data.obj <- linkage.blocks.network(data.obj, collapse.linked.markers = collapse.linked.markers, threshold.power = threshold.power, plot.blocks = plot.linkage.blocks)
		# }
	# if(linkage.method == "prominence"){
		# data.obj <- linkage.blocks.prominence(data.obj, p.or.q = p.or.q, collapse.linked.markers = collapse.linked.markers, threshold.power = threshold.power, verbose = verbose, plot.results = plot.linkage.blocks)
		# }
	
	
	if(collapse.linked.markers){
		blocks <- data.obj$linkage.blocks.collapsed
		}else{
		blocks <- data.obj$linkage.blocks.full
		}
		
	#make sure all block entries are numeric
	blocks <- lapply(blocks, function(x) get.marker.num(data.obj, x))
	
	#build a network based on the block structure
	all.net.data <- data.obj$var.to.var.p.val
	
	
	var.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(all.net.data))))
	#filter to just the significant interactions
	net.data <- all.net.data[which(as.numeric(all.net.data[,var.sig.col]) <= p.or.q),,drop = FALSE]
	#and filter by effect size if appropriate
	net.data <- net.data[which(as.numeric(net.data[,5]) >= min.std.effect),,drop = FALSE]
	pheno.tables <- data.obj$max.var.to.pheno.influence
	phenotypes <- names(pheno.tables)	


	adj.mat <- matrix(0, ncol = length(blocks), nrow = length(blocks))
	colnames(adj.mat) <- rownames(adj.mat) <- names(blocks)
	
	block.markers <- NULL
	for(i in 1:length(blocks)){
		block.markers <- rbind(block.markers, cbind(rep(names(blocks)[i], length(blocks[[i]])), blocks[[i]]))	
		}
	
	#get the block pairs for the interactions with significant effects
	get.sig.block.pair <- function(marker.pair){
		#find the block that marker 1 is in
		block1 <- which(block.markers[,2] == marker.pair[1])
		#find the block that marker 2 is in
		block2 <- which(block.markers[,2] == marker.pair[2])
		if(length(block1) > 0 && length(block2) > 0){
			return(c(block.markers[block1,1], block.markers[block2,1]))
			}else{
				return(c(0,0))
				}
		}

	sig.block.pairs <- unique(t(apply(matrix(net.data[,1:2], ncol = 2), 1, get.sig.block.pair)))

	if(length(sig.block.pairs) > 0){
	
		sig.block.pairs <- sig.block.pairs[which(sig.block.pairs[,1] != 0),,drop=FALSE]
	
		#for each pair of blocks
		get.adj.weight <- function(block.pair){
			#get all the markers in the two blocks
			all.markers1 <- blocks[[block.pair[1]]]
			all.markers2 <- blocks[[block.pair[2]]]
			
			#find the maximum weight between markers in the blocks
			#in both directions
			block1.source.locale <- which(net.data[,"Source"] %in% all.markers1)
			block2.target.locale <- which(net.data[,"Target"] %in% all.markers2)
			block1.to.block2 <- intersect(block1.source.locale, block2.target.locale)
			
			if(length(block1.to.block2) > 0){
				if(standardized){
					all.effects <- as.numeric(net.data[block1.to.block2,"Effect"])/as.numeric(net.data[block1.to.block2,"SE"])
					}else{
					all.effects <- as.numeric(net.data[block1.to.block2,"Effect"])	
					}
				adj.mat[block.pair[1], block.pair[2]] <- all.effects[which(abs(all.effects) == max(abs(all.effects)))][1]
				}
			
	
			block2.source.locale <- which(net.data[,"Source"] %in% all.markers2)
			block1.target.locale <- which(net.data[,"Target"] %in% all.markers1)
			block2.to.block1 <- intersect(block2.source.locale, block1.target.locale)
			
			if(length(block2.to.block1) > 0){
				if(standardized){
				all.effects <- as.numeric(net.data[block2.to.block1,"Effect"])/as.numeric(net.data[block2.to.block1,"SE"])
				}else{
					all.effects <- as.numeric(net.data[block2.to.block1,"Effect"])
					}
				adj.mat[block.pair[2], block.pair[1]] <- all.effects[which(abs(all.effects) == max(abs(all.effects)))][1]
				}
				
			return(adj.mat)
			}
		
	
	
		for(i in 1:length(sig.block.pairs[,1])){
			adj.mat <- get.adj.weight(sig.block.pairs[i,])
			}
		}	
	
	#Now add the phenotypic effects continuing to use the maximum significant effect from each block
	pheno.mat <- matrix(0, nrow = length(blocks), ncol = length(phenotypes))
	colnames(pheno.mat) <- phenotypes
	rownames(pheno.mat) <- names(blocks)
	
	pheno.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(pheno.tables[[1]]))))
	get.block.inf <- function(block){
		all.markers <- blocks[[block]]
		for(i in 1:length(pheno.tables)){
			sig.inf <- pheno.tables[[i]][which(pheno.tables[[i]][,pheno.sig.col] <= p.or.q),,drop = FALSE]
			# sig.inf[,1] <- get.marker.name(data.obj, as.numeric(sig.inf[,1]))
			if(length(sig.inf) > 0){
				block.locale <- which(sig.inf[,"marker"] %in% all.markers)
				# if(length(block.locale) == length(all.markers)){
				if(length(block.locale) >= 1){
					if(standardized){
					all.effects <- as.numeric(sig.inf[block.locale,"coef"])/as.numeric(sig.inf[block.locale,"se"])
					}else{
						all.effects <- as.numeric(sig.inf[block.locale,"coef"])
						}
					pheno.mat[block,i] <- all.effects[which(abs(all.effects) == max(abs(all.effects)))][1]
					}
				}
			}
		return(pheno.mat)	
		}

	for(i in 1:length(blocks)){
		pheno.mat <- get.block.inf(block = names(blocks)[i])
		}
	
	final.mat <- cbind(adj.mat, pheno.mat)

	if(collapse.linked.markers){
		data.obj$collapsed.net <- final.mat
		}else{
		data.obj$full.net <- final.mat	
		}
	
	return(data.obj)
	
	}
