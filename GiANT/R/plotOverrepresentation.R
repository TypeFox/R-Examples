plotOverrepresentation <- function(object,
		signLevel = object$signLevel,
		subset = NULL,
		aggregate = FALSE,
		ask = FALSE,
		...){

	if(class(object) != "gsaResult"){
		stop("'object' mut be of class gsaResult.")
	}
	if(object$analysis$name != "overrepresentation"){
		stop("'object' must contain results of an overrepresentation analysis.")
	}
	
	requireNamespace("limma")

	oldAsk <- par("ask")
	nms <- names(object$res.all)

	if(is.null(subset)){
		sel <- object$adjustedPValues < signLevel
	}else{
		sel <- seq_along(object$adjustedPValues) %in% subset
	}

	getTab <- function(nSets){
		t(sapply(0:(2^nSets -1),function(x){as.integer(intToBits(x))}))[,nSets:1]
	}

	nSel <- sum(sel)

	if(nSel == 0){
		stop("No significant gene sets found for the given significance level.")
	}

	if(aggregate){
		if(nSel > 4){
			stop("Cannot plot a venn diagram for more than 4 selected sets.")
		}

		ask <-  FALSE

		genes <- unique(c(unlist(sapply(which(sel), function(i){object$res.all[[i]]$geneSet})),
			object$res.all[[1]]$geneSetValues$coreSet))

		mat <- do.call(rbind,
			c(list(coreSet = genes %in% object$res.all[[1]]$geneSetValues$coreSet),
			lapply(object$res.all, function(x){
					genes %in% x$geneSet
				})))

		counts <- getTab(nSel+1)
		colnames(counts) <- c("coreSet", nms[sel])

		aktCounts <- cbind(counts, apply(counts,1,function(x){
				if(sum(x) == 0){
					return(object$res.all[[1]]$geneSetValues$nAllGenes)
				}else if(sum(x) == 2 && x[1] == 1 && object$adjustedPValues[colnames(counts)[which(x==1)[2]]] < signLevel){
					return(paste(sum(colSums(mat[colnames(counts)[as.logical(x)],,drop=FALSE]) == sum(x)), " *", sep =""))
				}
				return(sum(colSums(mat[colnames(counts)[as.logical(x)],,drop=FALSE]) == sum(x)))
			}))

		colnames(aktCounts) <- c("coreSet", nms[sel], "Counts")
		class(aktCounts) <- "VennCounts"

		#vennDiagramm in limma uses wrong color order if 5 sets are given
		if(nSel == 4){
			colors <- c(rep("black", nSel), "red")
		}else{
			colors <- c("red", rep("black", nSel))
		}
		#print(aktCounts)

		limma::vennDiagram(aktCounts, circle.col = colors, cex = c(1,0.75,0.7), ...)
	}else{
		ask <- ask

		counts <- getTab(2)

		lapply(1:sum(sel), function(j){

				i <- which(sel)[j]
				aktCounts <- cbind(counts, c(
						object$res.all[[i]]$geneSetValues$nAllGenes,
						length(object$res.all[[i]]$geneSet),
						length(object$res.all[[i]]$geneSetValues$coreSet),
						length(object$res.all[[i]]$geneSetValues$intersectGeneSetCoreSet)
					))
				if(object$adjustedPValues[i] < signLevel){
					aktCounts[nrow(aktCounts),ncol(aktCounts)] <- paste(aktCounts[nrow(aktCounts),ncol(aktCounts)], " *", sep ="")
				}

				colnames(aktCounts) <- c("coreSet", nms[i], "Counts")
				class(aktCounts) <- "VennCounts"

				par(ask=ask && j > 1 && dev.interactive())

				colors <- c("blue", rep("black", nSel))

				#print(aktCounts)

				limma::vennDiagram(aktCounts, circle.col = colors, cex = 1, ...)
		})
	}

	par(ask=oldAsk)

	return(invisible(NULL))
}
