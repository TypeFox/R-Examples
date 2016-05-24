age.range.correlation <- function(phy, overlap, tri = "upper", 
                                  n = 10000){
	
	# check input
	# -----------
	if ( !inherits(phy, "phylo") ) 
	  stop("object 'phy' is not of class 'phylo'")
	if ( !is.ultrametric(phy) )
    stop("object 'phy' must be ultrametric")
    		
	# ages:
	# -----
	age <- branching.times(phy)	
	
	# make matrix symmetrical
	# -----------------------
	ovlap <- overlap
	if ( tri == "upper" )
		ovlap[lower.tri(ovlap)] <- t(ovlap)[lower.tri(ovlap)]
	if ( tri == "lower" )
		ovlap[upper.tri(ovlap)] <- t(ovlap)[upper.tri(ovlap)]
	
	# match matrix to tree
	# --------------------
	id <- match(phy$tip.label, rownames(ovlap))
	ovlap <- ovlap[id, id]
	
	# calculate 'nested mean overlap'
	# -------------------------------
	overlap <- sapply(names(age), nested.mean.overlap, phy = phy, 		
                    olap = ovlap)

	x <- cbind(age, overlap)	
	x.lm <- lm(overlap ~ age)
	
	# randomization:
	# --------------
	randomization <- function(phy, o, n, age){
		id <- sample(seq(along = o[, 1]))
		rownames(o) <- colnames(o) <- colnames(o)[id]
		o <- sapply(names(age), nested.mean.overlap, phy = phy, 
                olap = o)
		o <- cbind(age, o)
		o <- lm(o[, 2] ~ age)
		o$coefficients
	}
	random.x <- lapply(1:n, randomization, o = ovlap, 
                     phy = phy, age = age)
	random.x <- do.call(rbind, random.x) # do conversion explicitly!
	
  ## fraction of intercepts and slopes from the randomization,
  ## which are greater than the observed values
  ## ------------------------------------------
	f.intercept <- random.x[, "(Intercept)" ] > x.lm$coefficients["(Intercept)"]
  f.intercept <- length(which((f.intercept))) / n
  
  f.slope <- random.x[, "age" ] > x.lm$coefficients["age"]
	f.slope <- length(which(f.slope)) / n
	
  ## 2-sided p-values
  ## ----------------
	f <- c(f.intercept, f.slope)
	p <- sapply(f, function(x) 2 * min(x, 1 - x))
	sig <- cbind(f, p)
	rownames(sig) <- c("(intercept)", "age")
	
	list(age.range.correlation = x, 
		linear.regression = x.lm, 
		sig = sig,
		MonteCarlo.replicates = random.x	
	)
}