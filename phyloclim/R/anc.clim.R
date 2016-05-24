anc.clim <- function(target, posterior = NULL, pno, n = 100, method = "GLS"){
	
	if (is.null(posterior))
		posterior <- list(target, target)
	
	# sample 'pno' 'n' times
	# ----------------------
	ntax <- dim(pno)[2] - 1
	x <- matrix(nrow = n, ncol = ntax)
	colnames(x) <- colnames(pno)[-1]
	for (i in 2:dim(pno)[2])
		x[, i - 1] <- sample(pno[, 1], size = n, replace = TRUE, 			prob = pno[, i])
	
	# core function: fits the model for all samples x to 
	# one tree sample from the posterior distribution
	# --------------------------------------------------
	core <- function(x, phy){
		
		# order x according to tiplabels in phy:
		# -----------------------
		x <- x[, match(phy$tip.label, colnames(x))]
		
		# ancestral character estimation:
		# -----------------------
		if (method == "GLS")
		out <- apply(x, 1, FUN = ace, 						phy = phy, type = "continuous", 				method = "GLS", CI = TRUE, model = "BM", 			kappa = 1, 	corStruct = corBrownian(1,  phy))
			
		if (method == "ML")	
		out <- apply(x, 1, FUN = ace, phy = phy,			type = "continuous", method = "ML", 			CI = FALSE, model = "BM", kappa = 1)
			
		# turn list elements into matrix:
		# number of samples x number of internal nodes
		# -------------------------------------------
		id <- which(c("GLS", "ML") %in% method)
		out <- lapply(out, function(x){x <- x[[id]]})
		out <- matrix(unlist(out), nrow = length(out), 			ncol = phy$Nnode, byrow = TRUE)
		out.mean <- apply(out, 2, mean)
	}
	out <- lapply(posterior, core, x = x)
	out <- matrix(unlist(out), nrow = length(out), 		ncol = posterior[[1]]$Nnode, byrow = TRUE)
	
	# calculate mean for nodes of target tree
	# ---------------------------------------
	for (i in seq(along = posterior)){

		# index vector for reordering
		id <- AC.node.trans(target, posterior[[i]], index = TRUE)
		nomatch <- which(!seq(along = id) %in% id)
		id[is.na(id)] <- nomatch
		
		# reorder ith row
		out[i, ] <- out[i, id]
		out[i, nomatch] <- NA
	}
	out.mean <- apply(out, 2, mean, na.rm = TRUE)
	
	# calculate weighted means for tips of target tree
	# -------------------------------------------------
	tips <- apply(x, 2, mean)
	tips <- tips[match(target$tip.label, names(tips))]
	
	# calculate 80% central density
	# -------------------------------------------------
	cd80 <- apply(x, 2, quantile, probs = c(0.1, 0.9))
	cd80 <- cd80[, match(target$tip.label, colnames(cd80))]
	
	out <- c(tips, out.mean)
	# check if reconstruction failed for some nodes:
	if (any(is.na(out))){
		nas <- which(is.na(out))
		nas <- paste(nas, collapse = ", ")
		stop("Reconstruction failed for nodes: ", nas)
	}
		
	out <- list(tree = target, means = out, 
		central.density = cd80)
}