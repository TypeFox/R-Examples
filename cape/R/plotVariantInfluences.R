plotVariantInfluences <-
function(data.obj, p.or.q = 0.05, min.std.effect = 0, plot.all.vals = FALSE, standardize = TRUE, pos.col = "brown", neg.col = "blue", not.tested.col = "lightgray", light.dark = "light", show.marker.labels = FALSE, show.chr = TRUE, label.chr = TRUE, scale.effects = c("log10", "sqrt", "none"), pheno.width = 11, covar.width = 11, covar.labels = NULL, phenotype.labels = NULL){
	
	if(plot.all.vals){
		light.dark = "light"
		}
	
	if(length(grep("n", scale.effects)) > 0){
		scale.effects <- "none"
		}
	if(length(scale.effects) == 1){
		if(scale.effects != "log10" & scale.effects != "sqrt" & scale.effects != "none"){
			stop("scale.effects must be 'log10', 'sqrt' or 'none.'")
			}
		}
	
	if(covar.width < 1 || !is.numeric(covar.width)){
		stop("covar.width must be a whole positive number.")
		}

	if(pheno.width < 1 || !is.numeric(pheno.width)){
		stop("pheno.width must be a whole positive number.")
		}

	var.influences <- data.obj$var.to.var.p.val
		
	pheno.inf <- data.obj$max.var.to.pheno.influence
	if(is.null(phenotype.labels)){
		pheno.names <- names(data.obj$max.var.to.pheno.influence)
		}else{
		pheno.names <- phenotype.labels
		if(length(pheno.names) != length(names(data.obj$max.var.to.pheno.influence))){
			stop("I am detecting the wrong number of phenotype labels for the phenotypes present.")
			}
		}
	num.pheno <- length(pheno.names)
	
	if(not.tested.col == TRUE){
		not.tested.col = "lightgray"
		}
	
	if(is.null(var.influences)){
		stop("calc.p() must be run to calculate variant-to-variant influences.")
		}

	if(is.null(pheno.inf)){
		stop("direct.influence() must be run to calculate variant-to-trait influences.")
		}

	#This function expands the given rows or columns of 
	#a matrix by a given amount
	expand.matrix <- function(mat, row.col.num, row.or.col, expansion.factor){
		if(row.or.col == "row"){
			row.labels <- 1:nrow(mat)
			for(i in 1:length(row.col.num)){
				mat.before <- mat[which(row.labels < row.col.num[i]),,drop=FALSE]
				mat.after <- mat[which(row.labels > row.col.num[i]),,drop=FALSE]
				row.to.expand <- mat[which(row.labels == row.col.num[i]),]
				mat.to.add <- matrix(row.to.expand, nrow = expansion.factor, ncol = ncol(mat), byrow = TRUE)
				row.labels <- c(row.labels[which(row.labels < row.col.num[i])], rep(row.col.num[i], expansion.factor), row.labels[which(row.labels > row.col.num[i])])
				mat <- rbind(mat.before, mat.to.add, mat.after)
				}
			return(mat)
			}


		if(row.or.col == "col"){
			col.labels <- 1:ncol(mat)
			for(i in 1:length(row.col.num)){
				mat.before <- mat[,which(col.labels < row.col.num[i]),drop=FALSE]
				mat.after <- mat[,which(col.labels > row.col.num[i]),drop=FALSE]
				col.to.expand <- mat[,which(col.labels == row.col.num[i])]
				mat.to.add <- matrix(col.to.expand, ncol = expansion.factor, nrow = nrow(mat), byrow = FALSE)
				col.labels <- c(col.labels[which(col.labels < row.col.num[i])], rep(row.col.num[i], expansion.factor), col.labels[which(col.labels > row.col.num[i])])
				mat <- cbind(mat.before, mat.to.add, mat.after)
				}
			return(mat)
			}
		}

	#This function replaces a row or column in a 
	#matrix with another matrix. If adding rows,
	#the two matrices must match in the number 
	#of columns
	replace.row.col <- function(orig.matrix, replace.matrix, row.col.num, row.or.col){
		if(row.or.col == "row"){
			orig.before <- orig.matrix[1:(min(row.col.num)-1),]
			if(max(row.col.num) < nrow(orig.matrix)){
				orig.after <- orig.matrix[(max(row.col.num)+1):nrow(orig.matrix), ]
				}else{
				orig.after <- NULL
				}
			new.mat <- rbind(orig.before, replace.matrix, orig.after)
			return(new.mat)
			}
		if(row.or.col == "col"){
			orig.before <- orig.matrix[,1:(min(row.col.num)-1)]
			if(max(row.col.num) < nrow(orig.matrix)){
				orig.after <- orig.matrix[,(max(row.col.num)+1):nrow(orig.matrix)]
				}else{
				orig.after <- NULL
				}
			new.mat <- cbind(orig.before, replace.matrix, orig.after)
			return(new.mat)
			}
	}
	
	unique.markers <- sort(as.numeric(unique(c(as.vector(var.influences[,"Source"]), as.vector(var.influences[,"Target"]), pheno.inf[[1]][,"marker"]))))
		
		
	marker.chr <- get.marker.chr(data.obj, unique.markers)
	covar.locale <- which(marker.chr == 0)
	num.covar <- length(covar.locale)
	if(num.covar > 0){
		expanded.chr <- c(marker.chr[-covar.locale], rep(0, (num.covar*covar.width)-num.covar))
		expanded.covar <- unlist(lapply(unique.markers[covar.locale], function(x) rep(x, covar.width)))
		expanded.markers <- c(unique.markers[-covar.locale], expanded.covar)
		}else{
		expanded.chr <- marker.chr
		expanded.markers <- unique.markers
		}

	covar.info <- get.covar(data.obj)
	covar.names <- covar.info$covar.names[which(covar.info$covar.type == "p")]
	
	#get coordinates of the chromosome boundaries
	if(show.chr){
		just.marker.chr <- marker.chr[which(marker.chr != 0)]
		if(num.covar > 0){
			for(i in 1:length(covar.names)){
				just.marker.chr <- c(just.marker.chr, rep(covar.names[i], covar.width))
				}
			}
		u_chr <- unique(just.marker.chr[which(!is.na(just.marker.chr))])
		chr.boundaries <- apply(matrix(u_chr, ncol = 1), 1, function(x) max(which(just.marker.chr == x))) + 0.5
		chr.boundaries <- c(0, chr.boundaries)
		if(label.chr){
			chr.names <- unique(just.marker.chr)
			}else{
			chr.names <- NULL
			}
		}else{
		chr.boundaries <- NULL
		chr.names <- NULL
		}
	
			
	var.influence.mat <- matrix(NA, nrow = length(unique.markers), ncol = length(unique.markers))
	var.pval.mat <- matrix(NA, nrow = length(unique.markers), ncol = length(unique.markers))
	colnames(var.influence.mat) <- rownames(var.influence.mat) <- colnames(var.pval.mat) <- rownames(var.pval.mat) <- unique.markers

	pheno.influence.mat <- matrix(NA, nrow = length(unique.markers), ncol = num.pheno)
	pheno.pval.mat <- matrix(NA, nrow = length(unique.markers), ncol = num.pheno)
	colnames(pheno.influence.mat) <- colnames(pheno.pval.mat) <- pheno.names
	rownames(pheno.influence.mat) <- rownames(pheno.pval.mat) <- unique.markers
	
	
	#fill the variant-to-variant matrix with test statistics with sources in rows and targets in columns
	var.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(var.influences))))
	for(i in 1:length(var.influences[,1])){
		if(standardize){
			var.influence.mat[as.character(var.influences[i,"Source"]), as.character(var.influences[i,"Target"])] <- as.numeric(as.vector(var.influences[i,"Effect"]))/as.numeric(as.vector(var.influences[i,"SE"]))
			}else{
			var.influence.mat[as.character(var.influences[i,"Source"]), as.character(var.influences[i,"Target"])] <- as.numeric(as.vector(var.influences[i,"Effect"]))		
				}
		var.pval.mat[as.character(var.influences[i,"Source"]), as.character(var.influences[i,"Target"])] <- as.numeric(as.vector(var.influences[i,var.sig.col]))
		}
	var.pval.mat[which(!is.finite(var.pval.mat))] <- 1
	
	#expand the covariate rows and columns to the specified width
	if(length(expanded.chr) > ncol(var.influence.mat)){
		covar.locale <- which(marker.chr == 0)
		
		new.var.inf <- expand.matrix(var.influence.mat, covar.locale, "row", covar.width)
		new.var.inf <- expand.matrix(new.var.inf, covar.locale, "col", covar.width)
		
		new.var.pval <- expand.matrix(var.pval.mat, covar.locale, "row", covar.width)
		new.var.pval <- expand.matrix(new.var.pval, covar.locale, "col", covar.width)		
		
		var.influence.mat <- new.var.inf
		var.pval.mat <- new.var.pval
		}
		

	
	#fill the variant-to-phenotype matrix with test statistics 
	#(still with sources in rows and targets in columns)
	#use phenotypes or eigentraits based on user input
	pheno.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(pheno.inf[[1]]))))
	for(i in 1:length(unique.markers)){
			for(j in 1:length(pheno.names)){
				marker.locale <- which(pheno.inf[[j]][,"marker"] == unique.markers[i])
				if(length(marker.locale) > 0){
					if(standardize){	
						pheno.influence.mat[as.character(unique.markers[i]), pheno.names[j]] <- pheno.inf[[j]][marker.locale, "t.stat"]
						}else{
						pheno.influence.mat[as.character(unique.markers[i]), pheno.names[j]] <- pheno.inf[[j]][marker.locale, "coef"]
						}
				pheno.pval.mat[as.character(unique.markers[i]), pheno.names[j]] <- pheno.inf[[j]][marker.locale, pheno.sig.col]
				}else{
					pheno.influence.mat[as.character(unique.markers[i]), pheno.names[j]] <- NA
					pheno.pval.mat[as.character(unique.markers[i]), pheno.names[j]] <- NA
					}
			}
		}
	
	#expand the phenotype influence matrix to give it more visual weight in the plot
	expanded.pheno.mat <- expand.matrix(pheno.influence.mat, 1:ncol(pheno.influence.mat), "col", pheno.width)
	expanded.pheno.pval.mat <- expand.matrix(pheno.pval.mat, 1:ncol(pheno.pval.mat), "col", pheno.width)
	
	#also expand the regions where the covariates are if there are covariates
	if(length(covar.locale) > 0){
		expanded.pheno.mat <- expand.matrix(expanded.pheno.mat, covar.locale, "row", covar.width)
		expanded.pheno.pval.mat <- expand.matrix(expanded.pheno.pval.mat, covar.locale, "row", covar.width)
		}
	
	expanded.pheno.names <- rep("", dim(pheno.influence.mat)[2]*pheno.width)
	
	start.col = 1
	for(i in 1:dim(pheno.influence.mat)[2]){
		expanded.pheno.names[median(start.col:(start.col+pheno.width-1))] <- colnames(pheno.pval.mat)[i]
		start.col = start.col + pheno.width
		}

	
	full.inf.mat <- cbind(var.influence.mat, expanded.pheno.mat)
	full.pval.mat <- cbind(var.pval.mat, expanded.pheno.pval.mat) 

	full.inf.mat <- apply(full.inf.mat, 2, as.numeric)
	full.pval.mat <- apply(full.pval.mat, 2, as.numeric)
	
	marker.names  <- get.marker.name(data.obj, expanded.markers)
	rownames(full.inf.mat) <- marker.names
	colnames(full.inf.mat) <- c(marker.names, expanded.pheno.names)
	
	if(num.covar > 0){
		for(i in 1:length(covar.names)){
			covar.locale <- which(rownames(full.inf.mat) == covar.names[i])
			new.covar.names <- rep("", length(covar.locale))
			new.covar.names[median(1:length(new.covar.names))] <- covar.names[i]
			rownames(full.inf.mat)[covar.locale] <- new.covar.names
			colnames(full.inf.mat)[covar.locale] <- new.covar.names
			}
		}	
	
	
	#get the coordinates for all pairs not tested
	# not.tested.locale <- which(is.na(rotate.mat(full.inf.mat)), arr.ind = TRUE)
	not.tested.locale <- which(is.na(rotate.mat(full.inf.mat)), arr.ind = TRUE)
	
	if(not.tested.col == FALSE || is.na(not.tested.col)){
		not.tested.locale <- NULL
		}
	
	#take out any values that aren't significant according
	#to the user cutoff, and do not have a high enough
	#effect size
	#if we are not plotting all value
	#use an extra color matrix to highlight significant interactions
	#if we are plotting all markers
	extra.col.mat <- NULL
	if(!plot.all.vals){
		full.inf.mat[which(full.pval.mat > p.or.q)] <- NA
		full.inf.mat[which(abs(full.inf.mat) < min.std.effect)] <- NA
		}else{
		extra.col.mat <- matrix(NA, nrow = nrow(full.pval.mat), ncol = ncol(full.pval.mat))

		min.effect <- which(abs(full.inf.mat) > min.std.effect)
		neg.effect <- intersect(min.effect, which(full.inf.mat < 0))
		pos.effect <- intersect(min.effect, which(full.inf.mat > 0))

		sig.neg <- intersect(neg.effect, which(full.pval.mat < p.or.q))
		sig.pos <- intersect(pos.effect, which(full.pval.mat < p.or.q))

		extra.col.mat[sig.neg] <- get.col(neg.col, "dark")[3]
		extra.col.mat[sig.pos] <- get.col(pos.col, "dark")[3]
		
		}
	main <- "Variant Influences"
	
	if(scale.effects == "log10"){
		neg.locale <- which(full.inf.mat < 0)
		pos.locale <- which(full.inf.mat > 0)
		scaled.effects <- log10(abs(full.inf.mat))
		scaled.effects[neg.locale] <- scaled.effects[neg.locale]*-1	
		scaled.effects[pos.locale] <- abs(scaled.effects[pos.locale])
		full.inf.mat <- scaled.effects
		main <- "log10 Variant Influences"
		}
	if(scale.effects == "sqrt"){
		neg.locale <- which(full.inf.mat < 0)
		scaled.effects <- sqrt(abs(full.inf.mat))
		scaled.effects[neg.locale] <- scaled.effects[neg.locale]*-1	
		full.inf.mat <- scaled.effects
		main <- "Square Root of Variant Influences"
		}
	
		
	if(length(which(na.omit(as.vector(full.inf.mat)) != 0)) == 0){
			plot.new()
			plot.window(xlim = c(0,1), ylim = c(0,1))
			text(0.5, 0.5, "No Significant Interactions")
		}else{
		myImagePlot(x = full.inf.mat, min.x = min(full.inf.mat, na.rm = TRUE), max.x = max(full.inf.mat, na.rm = TRUE), main = main, xlab = "Target", ylab = "Source", mark.coords = not.tested.locale, mark.col = not.tested.col, show.labels = show.marker.labels, chromosome.coordinates = chr.boundaries, chr.names = chr.names, show.pheno.labels = TRUE, extra.col.mat = extra.col.mat, pos.col = pos.col, neg.col = neg.col, col.pal = light.dark)


		#add phenotype names
		if(!is.null(not.tested.locale)){
			legend("topright", legend = "not testable", col = not.tested.col, pch = 16)
			}
		}
	
	invisible(full.inf.mat)

	
	}
