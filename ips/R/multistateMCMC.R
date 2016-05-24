multistateMCMC <- 
  function(phy, traits, model = "ARD", anc.states = TRUE, rd = 2.0, 
            rjhp = NULL, fixNodes = NULL, it = 100000, bi = 10000, 
            sa = 1000, path = "/Applications/BayesTraits", dir = NULL){
	
	# delete node.label, etc ....
	canonical <- c("edge", "Nnode", "tip.label", "edge.length")
	ind <- names(phy) %in% canonical
	if (!all(ind))
		phy[which(!ind)] <- NULL
		
		
	# integers as symbols for character states
	# ----------------------------------------
	foo <- function(x, states){
		which(states == x)
	}
	states <- unique(traits[, 2])
	nbstates <- length(states)
	max.nb.rates <- nbstates^2 - nbstates
	traits[, 2] <- sapply(traits[, 2], foo, states = states)
	
	# construct Addnode command
	# -------------------------
	if (is.list(anc.states)){
		g <- phy[[1]]
		addnode <- function(x, g){
			if(!all(x %in% g$tip.label))
				stop("incompatible taxon sets")
			paste(x, collapse = " ")		
		}
		addednodes <- sapply(anc.states, addnode, g = g)
		addednodes <- paste("AddMRCA", names(anc.states), addednodes)
	}
	
	# construct Fossil command
	# -------------------------
	if (is.list(fixNodes)){
		g <- phy[[1]]
		fixnode <- function(x, g, s){
			if(!all(x[-1] %in% g$tip.label))
				stop("incompatible taxon sets")
			x[1] <- which(s %in% x[1]) - 1
			paste(x, collapse = " ")
		}
		fixednodes <- sapply(fixNodes, fixnode, g = g, s = states)
		fixednodes <- paste("Fossil", names(fixNodes), fixednodes)
	}		
	
	#traits <- traits[traits[,1] %in% phy$tip.label,]
	
	m <- matrix(1, ncol = nbstates, nrow = nbstates)
	diag(m) <- 0
	m[lower.tri(m)] <- 1:(max.nb.rates/2)
	m <- t(m)
	m[lower.tri(m)] <- 1:(max.nb.rates/2) + (max.nb.rates/2)
	
	# set model
	if (model == "SYM"){
		restrict <- vector()
		extract <- vector()
		for(i in 1:(nbstates-1)){
			for(j in (i + 1):nbstates){
				restrict <- c(restrict, paste("Restrict q", i, j, " q", j, i, sep = ""))
				extract <- c(extract, paste("q", i, j, sep = ""))			}
		}
		m[lower.tri(m)] <- 1:(max.nb.rates/2)
		m <- t(m)
		m[lower.tri(m)] <- 1:(max.nb.rates/2)
	}
	if (model == "ER"){
		restrict <- vector()
		for(i in 1:(nbstates-1)){
			for(j in (i + 1):nbstates){
				restrict <- c(restrict, paste("q", i, j, 					" q", j, i, sep = ""))					}
		}
		restrict <- paste(restrict, collapse = " ")
		restrict <- paste("Restrict", restrict)
		extract <- "q12"
		m <- matrix(1, ncol = nbstates, nrow = nbstates)
		diag(m) <- 0
	}

	# assemble file names
	# -----------------------
	rts <- max.nb.rates
	if (model == "SYM") rts <- rts / 2
	if (model == "ER") rts <- 1
	if (!is.null(rjhp)) {
    rts <- "rjhp"
	}
	else {
    rts <- paste(rts, "rates", sep = "")
	}
	fns <- paste("MCMC", paste(nbstates, "states", sep = ""), 
	                 rts, sep = "_")
	fns <- paste(fns, "log", sep = ".")
	fns <- c("inputtree", "inputdata", "commands", fns)
	if (!is.null(dir)){
	  dir.create(dir)
	  fns <- paste(dir, fns, sep = "/")
	}
	
	# write commands file
	# -------------------
	cmds <- c(1, 2)
	if (is.list(anc.states))
		cmds <- c(cmds, addednodes)
	if (is.list(fixNodes))
		cmds <- c(cmds, fixednodes)
	if (model != "ARD")
		cmds <- c(cmds, restrict)
	cmds <- c(cmds, paste("rd", rd))
	if (!is.null(rjhp))
		cmds <- c(cmds, paste("rjhp", rjhp))
	cmds <- c(cmds, paste("it", it), paste("bi", bi), paste("sa", sa), 
          paste("lf", fns[4]), "run")
  
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
	  system2(CALL)  
	}
	if (.Platform$OS.type == "unix") {
	  system(CALL)  
	}
		
	# read BayesTrait output
	# -----------------------
	if (is.null(rjhp)){
		x <- scan(fns[4], what = "c", quiet = TRUE, sep = "\n")
		sk <- grep("Tree No", x)
		header <- scan(fns[4], skip = sk-1, nlines = 1, 			what = "c", sep = "\t", quiet = TRUE)
		x <- read.table(fns[4], skip = sk, sep = "\t")
		colnames(x) <- gsub(" |[(]|[)]", ".", header)
	}												else {
		x <- scan(fns[4], what = "c", quiet = TRUE, sep = "\n")
		log <- grep("\t'", x)
		log <- c(min(log) - 1, log)
		x <- x[log]
		x <- gsub("\t'", "\t", x)
		write(x, "temp.txt")
		x <- read.table("temp.txt", sep = "\t", header = TRUE)
	}
		
	xx <- x
	names(xx) <- gsub(" ", "", names(xx))
	names(xx)[1] <- "state"
	fn <- gsub("[.]log", ".table.log", fns[4])
	write.table(x, fn, quote = FALSE, col.names = 		TRUE, row.names = FALSE, sep = "\t")
	
	# probabilities of ancestral states
	# ---------------------------------
	state.freq <- grep("[.]P[.]", colnames(x))
	anclik <- as.matrix(x[, state.freq])
	
	mu <- apply(anclik, 2, mean)
	
	nds <- names(mu)
	nds <- unique(gsub("[.]P[.][[:digit:]][.]", "", nds))
	mu <- matrix(mu, ncol = nbstates, byrow = TRUE)
	rownames(mu) <- nds
	if (dim(mu)[1] > 1) mu <- mu[c(1, dim(mu)[1]:2), ]
	colnames(mu) <- states
	mu <- mu[, match(sort(states), states)]
	
	# calculate mean rates
	# --------------------
	rates <- x[ , grep("^q[[:digit:]]{2}", names(x))]
	rates <- mean(rates)
	if (model != "ARD") rates <- rates[names(rates) %in% extract]
	
	# harmonic mean
	hm <- tail(x$Harmonic.Mean, 1)
	
	# create output object
	# --------------------
	obj <- list(mean(x$Lh), hm, rates, m, mu)
	names(obj) <- c("loglik", "harmonic.mean", "rates", "index.matrix", "lik.anc")
	class(obj) <- "ace"
	obj	
}

