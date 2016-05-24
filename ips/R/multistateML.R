multistateML <- 								
function(phy, traits, model = "ARD", anc.states = TRUE, 
         path = "/Applications/BayesTraits", dir = NULL){
  
  ## checks
  ## ------
  if (!inherits(phy, "phylo")){
    stop("'phy' is not of class 'phylo'")
  }
  if (!inherits(traits, "data.frame")){
    stop("'traits' is not of class 'data.frame'")
  }
	
	# delete node.label, etc ....
	canonical <- c("edge", "Nnode", "tip.label", "edge.length")
	ind <- names(phy) %in% canonical
	if (!all(ind))
		phy[which(!ind)] <- NULL
		
	# integers as symbols for character states
	# ----------------------------------------
	c2i <- function(x, states){
		which(states == x)
	}
	states <- unique(traits[, 2])
	nbstates <- length(states)
	max.nb.rates <- nbstates^2 - nbstates
	traits[, 2] <- sapply(traits[, 2], c2i, states = states)
		
	# construct Addnode command
	# -------------------------
  if (anc.states){
    nodes <- length(phy$tip.label) + 2:phy$Nnode
    addnode <- function(node, phy){
      phy$tip.label[descendants(phy, node)]  	
    }
    anc.states <- lapply(nodes, addnode, phy = phy)
    names(anc.states) <- paste("node", nodes, sep = "")
  }
	if (is.list(anc.states)){
		if (class(phy) == "multiPhylo")
			g <- phy[[1]]								
    else 
      g <- phy
		check.added.node <- function(x, g){
			if(!all(x %in% g$tip.label))
				stop("incompatible taxon sets")
			paste(x, collapse = " ")
		}
		addednodes <- sapply(anc.states, check.added.node, g = g)
		addednodes <- paste("AddMRCA", names(anc.states), addednodes)
	}
	
	# set model
	# ------------
	m <- matrix(1, ncol = nbstates, nrow = nbstates)
	diag(m) <- 0
	m[lower.tri(m)] <- 1:(max.nb.rates/2)
	m <- t(m)
	m[lower.tri(m)] <- 1:(max.nb.rates/2) + (max.nb.rates/2)
	
	rts <- max.nb.rates
	
	if (model == "ER"){
		m <- matrix(1, ncol = nbstates, nrow = nbstates)
		diag(m) <- 0
	}
	if (model == "FB"){
		m[lower.tri(m)] <- 1
		m <- t(m)
		m[lower.tri(m)] <- 2
	}
	if (model == "ROW"){
		m <- matrix(sort(rep(1:nbstates, nbstates)), ncol = nbstates, nrow = nbstates)
		diag(m) <- 0
	}
	if (model == "SYM"){
		m[lower.tri(m)] <- 1:(max.nb.rates/2)
		m <- t(m)
		m[lower.tri(m)] <- 1:(max.nb.rates/2)
	}
	
	# create Restrict command
	# -----------------------
	stat <- sort(unique(as.vector(m)))
	stat <- stat[stat != 0]
	rts <- length(stat)

	cR <- function(x){
		x <- paste(x, collapse = "")
		paste("q", x, sep = "")
	}

	restrict <- extract <- vector(length = length(stat))
	for (i in seq(along = stat)){
		ind <- which(m == i, arr.ind = TRUE)
		r <- apply(ind, 1, cR)
		extract[i] <- tail(r, 1)
		restrict[i] <- paste(apply(ind, 1, cR), collapse = " ")
	}
	restrict <- paste("Restrict", restrict)
  
	# assemble file names
	# -----------------------
	fns <- paste("ML", paste(nbstates, "states", sep = ""), 
		paste(rts, "rates", sep = ""), sep = "_")
	fns <- paste(fns, "log", sep = ".")
  fns <- c("inputtree", "inputdata", "commands", fns)
  if (!is.null(dir)){
    dir.create(dir)
    fns <- paste(dir, fns, sep = "/")
  }
  
  # assamble commands file
  # -------------------
  cmds <- c(1, 1)
  if (is.list(anc.states))
    cmds <- c(cmds, addednodes)
  if (model != "ARD")
    cmds <- c(cmds, restrict)
  cmds <- c(cmds, paste("lf", fns[4]), "run")
  
  ## write input files and command file
  ## ----------------------------------
  write.nexus(phy, file = fns[1], translate = TRUE)
  write.table(traits, file = fns[2], quote = FALSE,   	
              row.names = FALSE, col.names = FALSE)
  write(cmds, file = fns[3])
	
	## execute BayesTraits:
	## --------------------
  ARGS <- paste(fns[1], fns[2], "<", fns[3])
  CALL <- paste(path, "/BayesTraits ", ARGS, sep = "")
  if (.Platform$OS.type == "unix") {
    system2(CALL)  # unix-alikes
  }
  else {
    system(CALL)   # windows
  }
	
	# read BayesTrait output
	# -----------------------
	x <- scan(fns[4], what = "c", quiet = TRUE, sep = "\n")
	sk <- grep("Tree No", x) - 1
	x <- read.table(fns[4], skip = sk, sep = "\t", header = TRUE)
	
	# write table
	# -----------
	xx <- x
	names(xx) <- gsub(" ", "", names(xx))
	names(xx)[1] <- "state"
	fn <- gsub("[.]log", ".table.log", fns[4])
	write.table(x, fn, quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
		
	# probabilities of ancestral states
	# ---------------------------------
	state.freq <- grep("[.]P[.]", colnames(x))
	anclik <- as.matrix(x[, state.freq])	
	if (class(phy) == "multiPhylo")	{
	  anclik <- apply(anclik, 2, mean)
	}									
	nds <- colnames(anclik)
	nds <- unique(gsub("[.]P[.][[:digit:]][.]", "", nds))
	anclik <- matrix(anclik, ncol = nbstates, byrow = TRUE)
	rownames(anclik) <- nds
	if (nrow(anclik) > 1) {
    anclik <- anclik[c(1, nrow(anclik):2), ]
	}
  colnames(anclik) <- states
	anclik <- anclik[, match(sort(states), states)]
		
	# calculate mean rates
	# --------------------
	rates <- x[ , grep("^q[[:digit:]]{2}", names(x))]
	if (class(phy) == "multiPhylo") rates <- mean(rates)
	if (model != "ARD") rates <- rates[names(rates) %in% extract]
  ## This hack is necessary because class 'ace' requires standard
  ## errors for each rate parameter: creation of of a NA-dummy-vector
  rates <- unlist(rates)
  se <- rep(NA, length(rates))
	
  ## assemble output object of class 'ace'
  ## see ?ace in ape
  ## ---------------
	obj <- list(x[, 2], rates, se, m, anclik)
	names(obj) <- c("loglik", "rates", "se", "index.matrix", "lik.anc")
	class(obj) <- "ace"
  return(obj)
}
