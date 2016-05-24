`dlmapdet` <-
function(input, algorithm, s.chr, chrSet, prevLoc=NULL, ...)
{
  dfMerged <- input$dfMerged
  n.chr <- length(input$map)
  nphe <- input$nphe
  type <- attr(input, "type") 
  results <- list()
  results$converge <- TRUE
  formula <- list()
  wald <- rep(0, length(input$map[[s.chr]]))

  if (algorithm=="asreml") {
    dfMrk <- input$dfMrk
    envModel <- input$envModel
    formula <- envModel
    groups <- list()

    ### calculate the minimum number of markers per chromosome
    nmrkchr <- vector(length=n.chr)
    for (i in 1:n.chr) nmrkchr[i] <- length(grep(paste("C", i, "M", sep=""), colnames(dfMrk)))
## Set up new dfMerged based on grouped random effects ##
###########################################
  if (min(nmrkchr) > nrow(dfMrk)) {
   dfm1 <- dfMerged[, c(1:nphe, match(prevLoc, names(input$dfMerged)))]

   # Set up random effects for markers on each chromosome	
   mat <- list()
   index <- list()
   ncolm <- vector(length=n.chr)
   for (kk in 1:n.chr) 
   {
    	index[[kk]] <- setdiff(grep(paste("C", kk, "M", sep=""), colnames(input$dfMrk)[2:ncol(input$dfMrk)]), match(prevLoc, colnames(input$dfMrk)[2:ncol(input$dfMrk)]))+1
	mat[[kk]] <- as.matrix(input$dfMrk[, index[[kk]]])
	mat[[kk]] <- mroot(mat[[kk]] %*% t(mat[[kk]]))
	ncolm[kk] <- ncol(mat[[kk]])
   }
   dfm2 <- do.call("cbind", mat)
   dfm2 <- as.data.frame(dfm2)
   ## need to set up idname to match
   dfm2 <- cbind(input$dfMrk[,1], dfm2)
   names(dfm2)[1] <- names(input$dfMrk)[1]
   dfMerged2 <- merge(dfm1, dfm2, by=names(dfm2)[1], all.x=TRUE, sort=FALSE)

   cumind <- c(0, cumsum(ncolm))
   for (kk in 1:n.chr)
   formula$group[[paste("g_", kk, "chr", sep="")]] <- ncol(dfm1) + (cumind[kk]+1):cumind[kk+1]
  } else {
  for (kk in 1:n.chr)
   formula$group[[paste("g_", kk, "chr", sep="")]] <- setdiff(grep(paste("C", kk, "M", sep=""),colnames(dfMerged)[(nphe+1):ncol(dfMerged)]), match(prevLoc, colnames(dfMerged)[(nphe+1):ncol(dfMerged)])) + nphe
  
   dfMerged2 <- dfMerged
  }
############################################
  groups <- formula$group

  ### in fact, this may need to be from a separate data frame 
  # Loop over positions on the selected chromosome
  mrkloop <- unlist(lapply(strsplit(names(dfMerged)[setdiff(grep(paste("C", s.chr, "M", sep=""), colnames(dfMerged)[(nphe+1):ncol(dfMerged)]), match(prevLoc, colnames(dfMerged)[(nphe+1):ncol(dfMerged)]))+nphe], "M"), function(x) return(x[2])))

  if (type=="f2")
	mrkloop <- unique(substr(mrkloop, 1, nchar(mrkloop)-1))
  mrkloop <- as.numeric(mrkloop)

  for (jj in 1:length(mrkloop))
  {
  	### add the selected marker onto dfMerged2
 	if (type=="f2")
	mrknam <- paste("C", s.chr, "M", mrkloop[jj], C("D", "A"), sep="")
	else 
	mrknam <- paste("C", s.chr, "M", mrkloop[jj], sep="")

	dfMerged3 <- dfMerged2
	if (length(intersect(mrknam, colnames(dfMerged2)))==0) {
	dfMerged3 <- cbind(dfMerged[, match(mrknam, colnames(dfMerged))], dfMerged2) 
	for (kk in 1:n.chr) formula$group[[paste("g_", kk, "chr", sep="")]] <- groups[[paste("g_", kk, "chr", sep="")]]+1} else {
	m <- match(mrknam, colnames(dfMerged2))
	dfMerged3 <- dfMerged3[,c(m, setdiff(1:ncol(dfMerged3),m))] 
	for (kk in 1:(s.chr-1)) formula$group[[paste("g_", kk, "chr", sep="")]] <- groups[[paste("g_", kk, "chr", sep="")]]+1 }	
	names(dfMerged3)[1] <- mrknam
 
	chrnam <- paste("idv(grp(g_", setdiff(chrSet,s.chr), "chr))", sep="")
	formula$random <- paste("~", paste(chrnam, collapse="+"))

	if (!is.null(envModel$random))
	  formula$random <- paste(formula$random, "+", as.character(envModel$random[2]), sep="")
	formula$random <- as.formula(formula$random)

	formula$fixed <- paste(as.character(envModel$fixed)[2], "~", as.character(envModel$fixed[3]), sep="")

	if (length(prevLoc) >0)
	formula$fixed <- paste(formula$fixed, "+", paste(prevLoc, collapse="+"), sep="")
	if (type=="f2")
	formula$fixed <- paste(formula$fixed, "+", paste(paste("C", s.chr, "M", mrkloop[jj], c("D", "A"), sep=""), collapse="+"), sep="") else
	formula$fixed <- paste(formula$fixed, "+", paste("C", s.chr, "M", mrkloop[jj], sep=""), sep="")

	formula$fixed <- as.formula(formula$fixed)

	formula$control <- envModel$control
	formula$eqorder <- 3
	formula$data <- "dfMerged3"
   	formula$Cfixed <- TRUE
  	formula <- c(formula, ...)
   	formula <- formula[!duplicated(formula)]
	formula <- formula[!sapply(formula, is.null)]
	if (length(chrSet)>1)  
	model <- do.call("asreml", formula)

	if (length(chrSet)==1)
	{
	formula1 <- formula
	formula1$random <- envModel$random
   	formula1$group <- envModel$group
	formula1 <- formula1[!sapply(formula1, is.null)]
	model <- do.call("asreml", formula1)
	}

	if (model$converge==FALSE) 	results$converge <- FALSE

	names <- paste("C", s.chr, "M", mrkloop[jj], c("", "D", "A"), sep="")
	if (model$coefficients$fixed[names(model$coefficients$fixed) %in% names]!=0) 	
	if (type=="f2")
	wald[mrkloop[jj]] <- waldtest.asreml(model, list(list(which(model$coefficients$fixed[names(model$coefficients$fixed) %in% names]!=0), "zero")))$zres$zwald
	else	{
		cf <- summary(model, all=TRUE)$coef.fixed
		wald[mrkloop[jj]] <- (cf[which(rownames(cf) %in% names),3])^2
		}
  }

  }  ## end of algorithm==asreml

  if (algorithm=="lme") {
	  fixed <- input$envModel$fixed
  chrRE <- vector()

  # Set up random effects for markers on each chromosome	
  for (kk in 1:n.chr)
  chrRE[kk] <- paste("pdIdent(~", paste(setdiff(names(dfMerged)[grep(paste("C", kk, "M", sep=""), names(dfMerged))], prevLoc), collapse="+"), "-1)", sep="")

  # Loop over positions on the selected chromosome
  mrkloop <- unlist(lapply(strsplit(names(dfMerged)[setdiff(grep(paste("C", s.chr, "M", sep=""), names(dfMerged)[(nphe+1):ncol(dfMerged)]), match(prevLoc, colnames(dfMerged)[(nphe+1):ncol(dfMerged)]))+nphe], "M"), function(x) return(x[2])))

  if (type=="f2")
	mrkloop <- unique(substr(mrkloop, 1, nchar(mrkloop)-1))
  mrkloop <- as.numeric(mrkloop)

  for (jj in 1:length(mrkloop))
  {
    # Only include random effects for chromosomes not being scanned
    if (length(chrSet)>2)	
	formula$random <- paste("pdBlocked(list(", paste(chrRE[setdiff(chrSet, s.chr)], collapse=","), "))", sep="")

    if (length(chrSet)==2)
	formula$random <- chrRE[setdiff(chrSet,s.chr)]

    formula$fixed <- paste(as.character(fixed)[2], "~", as.character(fixed)[3], sep="")

	# If there are markers already mapped on other chromosomes, add in as fixed effects 
    if (length(prevLoc) >0)
	formula$fixed <- paste(formula$fixed, "+",paste(prevLoc, collapse="+"), sep="")

    	effectnames <- paste("C", s.chr, "M", mrkloop[jj], sep="")
    	if (type=="f2") effectnames <- paste(effectnames, c("A", "D"), sep="")
    	
	formula$fixed <- paste(formula$fixed, "+", paste(effectnames, collapse="+"), sep="")

	formula$fixgrp <- paste(formula$fixed, "| grp1", sep="")
	formula$fixed <- as.formula(formula$fixed)
	formula$fixgrp <- as.formula(formula$fixgrp)

	gd <- groupedData(formula$fixgrp, data=dfMerged)

	# Fit model - different forms depending on relevant terms
	if (length(chrSet)>1)
	model <- lme(fixed=formula$fixed, random=eval(parse(text=formula$random)), data=gd, control=lmeControl(maxIter=input$maxit), na.action=na.omit)

	if (length(chrSet)==1)
	model <- lme(fixed=formula$fixed, random=~1|grp1, data=gd, control=lmeControl(maxIter=input$maxit), na.action=na.omit)

	# Get output - Wald statistic for position fixed effect
	wald[mrkloop[jj]] <- anova(model, Terms=effectnames)[1,3]
  }
} ## end of algorithm==lme  
 
  results$wald <- wald
  return(results)
}

