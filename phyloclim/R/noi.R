noi <- function(tr, group, regex = NULL, monophyletic = FALSE){
	
	# test to filter out trees with
	# non-consecutive tip labels
	# --------------------------
	canonical <- seq(along = tr$tip.label)
	canonical <- as.numeric(canonical)
    given <- tr$edge[, 2][tr$edge[, 2] %in% canonical]
    given <- as.numeric(given)
    if (!identical(canonical, given)) 
        stop("tips are not numbered consecutively.",
	    " Type '?fixTips' for help.")
	
	# matrix of pairwise MRCA
	# -----------------------
	x <- mrca(tr)
	
	# core function
	# -------------
	foo <- function(group, tr, regex){
		if (is.null(regex)){
			y <- which(tr$tip.label %in% group)
			group <- tr$tip.label[y]	
			mintax <- group[y == min(y)]
			maxtax <- group[y == max(y)]
		}												                   else {
		    regex<- paste(group, collapse = "|")
		    y <- grep(regex, tr$tip.label)
			mintax <- tr$tip.label[min(y)]
			maxtax <- tr$tip.label[max(y)]
		}
		x <- x[rownames(x) == mintax, colnames(x) == maxtax]
		if (monophyletic){
			test <- tr$tip.label[descendants(tr, x)]
			if (!all(test %in% group))
				x <- NA
		}
		x
	}
	
	# turn 'group' to list
	if (!is.list(group)) group <- list(group)
	
	# Convert tip numbers to tip labels
	# ---------------------------------
	if (mode(group[[1]]) == "numeric" & is.null(regex))
		group <- lapply(group, function(x, phy) 
	        phy$tip.label[x], phy = tr)
	    
	# check tip labels
	# ----------------
	chk <- unlist(group)[!unlist(group) %in% tr$tip.label]
	if (length(chk) > 0 & is.null(regex))						
	    stop(paste(chk, collapse = "\n"), 
	        "\nare not valid tip labels.")
	
	
	nodes <- sapply(group, foo, tr = tr, regex = regex)
	nodes
}
