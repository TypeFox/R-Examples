plotPheno <-
function(data.obj, pheno.which = NULL, color.by = NULL, group.labels = NULL, pheno.labels = NULL){
	
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

		layout.mat <- get.layout.mat(dim(data.obj$pheno)[2])
		layout(layout.mat)
		for(p in 1:length(pheno.labels)){
			plot(data.obj$pheno[,p], col = cols, pch = 16, xlab = "Individual", ylab = pheno.labels[p], main = pheno.labels[p])
			if(!is.null(groups)){
				legend("topleft", pch = 16, col = all.cols[1:length(groups)], legend = group.labels)			
				}
			}
		
	}
