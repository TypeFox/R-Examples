pwc <- function(levels,abbrev.names){
	m <- ncol(comb <- combn(length(levels),2))
	M <- matrix(,nrow=2,ncol=m)
	M[1,] <- levels[comb[1,]]
	M[2,] <- levels[comb[2,]]
	if (abbrev.names){
		if (is.logical(abbrev.names)){
			levels <- abbreviate(levels)
		}else  levels <- abbreviate(levels,minlength=abbrev.names)
		colnames(M) <- paste(levels[comb[1,]],levels[comb[2,]],sep="-")
	}else colnames(M) <- paste(M[1,],M[2,],sep="-")
	return(M)
}

flm <- function(levels,abbrev.names){
	M <- matrix(levels,1)
	if (abbrev.names){
		if (is.logical(abbrev.names)){
			levels <- abbreviate(levels)
		}else  levels <- abbreviate(levels,minlength=abbrev.names)
	}
	colnames(M) <- levels
	return(M)
}

residcontr <- function(levels,abbrev.names){
	n <- length(levels)
	M <- diag(n)-matrix(1/n,n,n)
	if (abbrev.names){
		if (is.logical(abbrev.names)){
			levels <- abbreviate(levels)
		}else  levels <- abbreviate(levels,minlength=abbrev.names)
	}
	colnames(M) <- paste(levels,"(resid.)")
	return(M)
}

testInteractions <- function(model, pairwise=NULL, fixed=NULL, residual=NULL, across=NULL, custom=NULL, slope=NULL, adjustment=NULL, label.factors=FALSE, abbrev.levels=FALSE, ...){
	# Levels in existing factors
	xlevels <-  getXLevels(model)
	dots <- list(...)
	# Include within-subjects factors, if they are defined
	if ("idata" %in% names(dots)) xlevels <- c(xlevels, lapply(dots$idata, levels))
	if (missing(pairwise) && missing(fixed) && missing(residual) && missing (across) && missing(custom)) pairwise=names(xlevels)
	# Check the factors specified in the arguments
	if (!all(c(pairwise,fixed,across,names(custom)) %in% names(xlevels))){
		warning("Some factors with specified contrasts are not in the model and will be ignored.")
		pairwise <- pairwise[pairwise %in% names(xlevels)]
		fixed <- fixed[fixed %in% names(xlevels)]
		residual <- residual[residual %in% names(xlevels)]
		across <- across[across %in% names(xlevels)]
		custom <- custom[names(custom) %in% names(xlevels)]
	}
	# Set of tested contrasts
	pairwise.contrasts <- lapply(xlevels[pairwise],pwc,abbrev.levels)
	fixed.levels <- lapply(xlevels[fixed],flm,abbrev.levels)
	residual.effects <- lapply(xlevels[residual],residcontr,abbrev.levels)
	custom <- lapply(custom,as.matrix)
	tested.contrasts <- c(pairwise.contrasts,fixed.levels,residual.effects,custom)
	# Check that there are no repeated names in the list of matrices
	repeatedfactors <- (table(names(tested.contrasts)) > 1L)
	if (any(repeatedfactors)) stop("Redundant definition of factor contrasts")
	# All possible combinations of columns of different matrices
	combinations <- as.matrix(expand.grid(lapply(tested.contrasts,function(x) seq(1,ncol(x)))))
	combination.labels <- character()
	# Define the rownames for the different tests, from column names of contrast matrices
	for (fname in names(tested.contrasts)){
		cont.combinations <- combinations[,fname]
		cont.names <- colnames(tested.contrasts[[fname]])
		# Append factor name if required
		if (is.null(cont.names)){
			column.labels <- paste(fname,cont.combinations,sep="")
		}else{
			if (label.factors) cont.names <- paste(fname,cont.names,sep="")
			column.labels <- cont.names[cont.combinations]
		}
		combination.labels <- paste(combination.labels, format(column.labels,justify="right"), sep=" : ")
	}
	combination.labels <- substr(combination.labels,4,nchar(combination.labels[1]))
	nc <- nrow(combinations)
	# The result will be an anova table built by testFactors
	anova.table <- NULL
	# Argument test.formula for testFactors
	if (is.null(across) && is.null(slope)){
		terms.formula <- ~1
		terms.message <- ""
	}else{
		terms.formula <- paste(c(across,slope),collapse=":")
		terms.formula <- as.formula(paste("~0 +",terms.formula))
		terms.message <- if (is.null(slope)) "" else paste("Adjusted slope for",paste(slope,collapse=":"),"\n")
	}
	for (i in seq(1,max(1,nc))){
		# Run testFactors to obtain adjustes terms and the anova table
		if (length(tested.contrasts) > 0){
			levels <- mapply(function(lab,n) tested.contrasts[[lab]][,n], names(tested.contrasts), combinations[i,], SIMPLIFY=FALSE)
		}else{
			levels <- NULL
			combination.labels <- "Mean"
		}
		test <- testFactors(model,levels=levels,terms.formula=terms.formula,...)
		# Get adjusted values (stripped into a vector, if they are a matrix),
		# and attach them to the anova table
		adjusted.values <- summary(test)$adjusted.values[[1]]
		if (nrow(test.table <- summary(test)$anova.table) == 0) stop(test$terms[[1]]$test)
		anova.table <- rbind(anova.table,c(as.numeric(adjusted.values),unlist(test.table[1,])))
	}
	# Adjust p-values
	nc <- ncol(anova.table)
	has.P <- grepl("^(P|Pr)\\(", dimnames(anova.table)[[2]][nc])
	if (has.P) anova.table[,nc] <- p.adjust(anova.table[,nc],adjustment)
	anova.table <- as.data.frame(anova.table)
	# Name for the adjusted values column (take last test as pattern)
	if (length(adjusted.values) > 1L){
		adjval.names <- attr(adjusted.values,"dimnames")
		adjval.names[sapply(adjval.names,is.null)] <- NULL
		adjval.names <- expand.grid(adjval.names)
		adjval.names <- if (length(adjval.names)==2L) paste(adjval.names[1],adjval.names[2],sep=":") else as.character(unlist(adjval.names))
		names(anova.table)[1:length(adjusted.values)] <- adjval.names
	}else names(anova.table)[1] <- "Value"
	# Row names
	rownames(anova.table) <- combination.labels
	# Add row for residuals, if existed in the anova table
	if (nrow(test.table) > 1){
		anova.table <- rbind(anova.table,Residuals=c(matrix(,ncol=length(adjusted.values)),unlist(test.table[2,])))
	}
	# Attach combination labels and make anova table class
	heading <- attr(summary(test)$anova.table,"heading")
	# Redefine message about calculated term, if it is the link function in a glm
	if ("glm" %in% class(model) && attr(test,"means")=="link" && terms.message=="") terms.message <- "Adjusted link function\n"
	adjustment.message <- paste("\nP-value adjustment method:", if (is.null(adjustment)) p.adjust.methods[1] else adjustment)
	heading <- paste(terms.message,heading,adjustment.message,sep="")
	structure(anova.table, heading=heading, class=c("anova", "data.frame"))
}

