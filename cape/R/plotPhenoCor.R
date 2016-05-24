plotPhenoCor <-
function(data.obj, pheno.which = NULL, color.by = NULL, group.labels = NULL, text.cex = 1, pheno.labels = NULL){
			
	all.pheno <- data.obj$pheno

	if(is.null(pheno.which)){
		pheno.names <- colnames(all.pheno)
		}else{
			if(is.numeric(pheno.which)[1]){
				pheno.names <- colnames(all.pheno)[pheno.which]
				}else{
				pheno.names <- pheno.which	
				}
			}
	
	if(is.null(pheno.labels)){
		pheno.labels <- pheno.names
		}
		
	pairs.matrix <- pair.matrix(pheno.names)

	all.cols <- brewer.pal(8, "Accent")


	if(!is.null(color.by)){
		covar.info <- get.covar(data.obj)
		cols <- rep(NA, dim(all.pheno)[1])
		group.col <- which(covar.info$covar.names %in% color.by)
			
		if(length(group.col) == 0){
			stop(paste("I couldn't find the", color.by, "column. Please check the case and the spelling."))
			}
			
		group.mem <- covar.info$covar.table[,group.col]
		groups <- sort(unique(group.mem))
		num.groups <- length(groups)
		if(num.groups > 8){
			stop("There cannot be more than 8 groups")
			}
		for(i in 1:num.groups){
			cols[which(group.mem == groups[i])] <- all.cols[i]
			}
		}else{
		cols <- rep("black", dim(all.pheno)[1])
		groups <- NULL
		}

	if(is.null(group.labels)){
		group.labels <- groups
		}


			
	plot.cor <- function(x,y){
		points(x,y, col = cols)
		if(!is.null(color.by)){
			legend("topleft", legend = group.labels, fill = all.cols[1:length(group.labels)], cex = 0.7)
			}
		}

	
	write.cor <- function(x,y){
		if(length(which(!is.na(x+y))) > 0){
	
			if(!is.null(groups)){
				group.cor <- try(apply(matrix(groups, ncol = 1), 1, function(a) cor(x[which(group.mem == a)], y[which(group.mem == a)], use = "complete")), silent = TRUE)
				}else{group.cor = NULL}
				pheno.cor <- try(cor(x, y, use = "complete"), silent = TRUE)
				}

			cors <- c(group.cor, pheno.cor)
			group.labels <- c(group.labels, "R")
			x.shrinkage = 0.4
			y.shrinkage = 0.3
			x.range <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
			y.range <- max(y, na.rm = TRUE) - min(y, na.rm = TRUE)
			x.shift <- x.shrinkage*x.range
			y.shift <- y.shrinkage*y.range
			y.placement <- segment.region(min(y, na.rm = TRUE)+y.shift, max(y, na.rm = TRUE)-y.shift, max(c(length(cors),2)), alignment = "center")[1:length(cors)]
			x.placement <- segment.region(min(x, na.rm = TRUE)+x.shift, max(x, na.rm = TRUE)-x.shift,3, alignment = "center")
			text(rep(x.placement[1], length(cors)), y.placement, group.labels, cex = text.cex, adj = 1)
			text(rep(x.placement[2], length(cors)), y.placement, rep("=", length(cors)), cex = text.cex, adj = 0.5)
			text(rep(x.placement[3], length(cors)), y.placement, signif(cors, 3), cex = text.cex, adj = 0)	
			}
	
	panel.hist <- function(x, ...){
    		usr <- par("usr"); on.exit(par(usr))
    		par(usr = c(usr[1:2], 0, 1.5) )
	    	h <- hist(x, plot = FALSE)
    		breaks <- h$breaks; nB <- length(breaks)
	    	y <- h$counts; y <- y/max(y)
    		rect(breaks[-nB], 0, breaks[-1], y, ...)
		}


	pairs(all.pheno, lower.panel = plot.cor, upper.panel = write.cor, diag.panel = panel.hist, labels = pheno.labels)
	# legend("topleft", legend = group.labels, fill = all.cols[1:length(group.labels)])
		

}
