`dltest` <-
function(input, algorithm, chrSet, prevLoc=NULL, ...)
{
  dfMerged <- input$dfMerged
  n.chrSet <- length(chrSet)
  null.ll <- vector(length=n.chrSet)
  null.forma <- list()
  random.forma <- list() 
  results <- list()
  results$converge <- TRUE
  formula <- list() 

  if (algorithm=="asreml") {
    dfMrk <- input$dfMrk
    envModel <- input$envModel
    n.perm <- input$nperm
    n.chr <- length(input$map)
    nphe <- input$nphe
    idname <- input$idname

    formula <- envModel 
    npop <- ngen(input)

    # Create permutation matrices  
    perm.test <- matrix(nrow=n.perm+1, ncol=n.chrSet)
    maxp <- vector(length=n.perm+1)
    perm.mat <- matrix(nrow=npop, ncol=n.perm+1)
    perm.mat[,1] <- c(1:npop)

    if (n.perm>0)
    for (kk in 2:(n.perm+1))
	perm.mat[,kk] <- sample(npop)

    # Fit full model, with all chromosomes having separate VCs
    formula$fixed <- paste(as.character(envModel$fixed)[2], "~",as.character(envModel$fixed[3]), sep="")
    # Include fixed effects for all markers which have already been mapped
    if (length(prevLoc)>0)
    formula$fixed <- paste(formula$fixed, "+", paste(prevLoc, collapse="+"), sep="")
    formula$fixed <- as.formula(formula$fixed) 
 
    nmrkchr <- vector(length=n.chr)
    for (i in 1:n.chr) nmrkchr[i] <- length(grep(paste("C", i, "M", sep=""), colnames(dfMrk)))
############################################
# Random effects for all markers on a chromosome excluding those which
## Set up new dfMerged based on grouped random effects ##
###########################################
  if (min(nmrkchr) > nrow(dfMrk)) {
   dfm1 <- dfMerged[,c(1:nphe, match(prevLoc, names(dfMerged)))]
   index <- list()
   mat <- list()
   ncolm <- vector(length=n.chr)
   for (ii in 1:n.chr) {
    index[[ii]] <- 1+setdiff(grep(paste("C", ii, "M", sep=""), names(dfMrk)[2:ncol(dfMrk)]), match(prevLoc, names(dfMrk)[2:ncol(dfMrk)]))
    mat[[ii]] <- as.matrix(dfMrk[, index[[ii]]])
    mat[[ii]] <- mroot(mat[[ii]] %*% t(mat[[ii]]))
    ncolm[ii] <- ncol(mat[[ii]])
   }
   dfm2 <- do.call("cbind", mat)
   dfm2 <- as.data.frame(dfm2)
   dfm2 <- cbind(dfMrk[,1], dfm2)
   names(dfm2)[1] <- colnames(dfMrk)[1]
   dfMerged2 <- merge(dfm1, dfm2, by=names(dfm2)[1], all.x=TRUE, sort=FALSE)

   colnames(dfMerged2)[(ncol(dfm1)+1):ncol(dfMerged2)] <- paste("var", 1:(ncol(dfMerged2)-ncol(dfm1)), sep="")

   cumind <- c(0, cumsum(ncolm)) 

   for (ii in 1:n.chr)
	formula$group[[paste("g_", ii, "chr", sep="")]] <- ncol(dfm1) + (cumind[ii]+1):cumind[ii+1]
  } else {
  for (ii in 1:n.chr)
    formula$group[[paste("g_", ii, "chr", sep="")]] <- nphe+setdiff(grep(paste("C", ii, "M", sep=""), names(dfMerged)[(nphe+1):ncol(dfMerged)]), match(prevLoc, names(dfMerged)[(nphe+1):ncol(dfMerged)]))
    dfMerged2 <- dfMerged
  } 

############################################

  # Random effects for each chromosome in selected subset
  # Markers are modelled as independent and same variance within chromosomes
  chrnam <- paste("idv(grp(g_", chrSet, "chr))", sep="")
  formula$random <- paste("~", paste(chrnam, collapse="+"))

  if (length(envModel$random)>0)
	formula$random <- paste(formula$random, "+", as.character(envModel$random[2]), sep="")

  formula$random <- as.formula(formula$random)

  formula$dump.model <- TRUE
  formula$data <- "dfMerged2"
  formula$control <- envModel$control
  formula$eqorder <- 3
  formula <- c(formula, ...)
  formula <- formula[!duplicated(formula)]
  formula <- formula[!sapply(formula, is.null)]
 
  full <- do.call("asreml", formula)

  if (n.chrSet==1)
  {
     form.null <- formula
     form.null$random <- envModel$random
     form.null <- form.null[!sapply(form.null, is.null)]
     null.forma[[1]] <- do.call("asreml", form.null)
     null.forma[[1]]$control$eqorder <- 3
  }
 
  if (n.chrSet>1)
  for (cc in 1:n.chrSet)
  {
 	# fit model leaving out each chromosome to test VC
	chrnam <- paste("idv(grp(g_", setdiff(chrSet, chrSet[cc]), "chr))", sep="")
	rndf <- paste("~", paste(chrnam, collapse="+"), sep="")
	if (!is.null(envModel$random))
	rndf <- paste(rndf, "+", as.character(envModel$random[2]), sep="")

	rndf <- as.formula(rndf)
  	form.null <- formula	
	form.null$random <- rndf
	null.forma[[cc]] <- do.call("asreml", form.null)
	null.forma[[cc]]$control$eqorder <- 3
  }

  # Vector of observed test statistics from LRTs
  if (n.perm==0)
  {
	run <- asreml(model=full) 
	full.ll <- run$loglik
	if (run$converge==FALSE) results$converge <- FALSE

	for (cc in 1:n.chrSet)
	{
	  run <- asreml(model=null.forma[[cc]])
	  null.ll[cc] <- run$loglik
	  if (run$converge==FALSE) results$converge <- FALSE
 	}
  	perm.test[1,] <- 2*(full.ll-null.ll)
	results$obs <- perm.test[1,]

	results$raw.pval <- sapply(perm.test[1,], pvfx)
	if (input$multtest=="bon")
	results$adj.pval <- sapply(results$raw.pval, function(x) return(min(x*n.chrSet,1))) else {
	pval <- as.matrix(rbind(c(1:length(results$raw.pval)), results$raw.pval))
	pval <- as.matrix(pval[,order(pval[2,])])
   	pval[2,] <- sapply(pval[2,]*(n.chrSet:1), function(x) return(min(x,1)))
	results$adj.pval <- pval[2, order(pval[1,])] }
	results$thresh <- qchibar(input$alpha/n.chrSet)
  }

  if (n.perm>0)
  {
    namesrnd <- setdiff(names(dfMrk)[2:ncol(dfMrk)], prevLoc)

    for (ii in 1:(n.perm+1))
    {
	if (ncol(dfMrk) > 20*nrow(dfMrk)) {
	 ### now instead of merging the final matrix on, reconstruct
	 names(dfmp)[1] <- names(dfMrk)[1]
	 df3 <- merge(dfm1, dfmp, by=names(dfMrk)[1], all.x=TRUE, sort=FALSE)
	 dfmp <- cbind(dfMrk[,1], dfm2[perm.mat[,ii],])
    	} else {
         df2 <- cbind(dfMrk[,1],dfm2[perm.mat[,ii],which(names(dfMrk)%in%namesrnd)])
         names(df2)[1] <- names(dfMrk)[1]
         df4 <- dfMerged[, match(c(idname, setdiff(names(dfMerged), names(df2))), names(dfMerged))]
         df3 <- merge(df4, df2, by=idname, all.x=TRUE, sort=FALSE)
	}

	df3 <- df3[,match(names(dfMerged2), names(df3))]
	full <- update(full, data=df3)	

        # replace data in model for random marker effects
#        index <- match(namesrnd, names(full$data))
#        index <- index[!is.na(index)]
#        index2 <- match(names(full$data)[index], names(df3))
#        full$data[,index] <- df3[, index2]

	# run full model
	run <- asreml(model=full)
	full.ll <- run$loglik
        if (run$converge==FALSE)
	  results$converge <- FALSE

	for (cc in 1:n.chrSet)
	{
	   # replace data for random marker effects for each of the null models
#          index <- match(namesrnd, names(null.forma[[cc]]$data))
#          index <- index[!is.na(index)]
#          index2 <- match(names(null.forma[[cc]]$data)[index], names(df3))
	   null.forma[[cc]] <- update(null.forma[[cc]], data=df3)
#          null.forma[[cc]]$data[, index] <- df3[, index2]
           run <- asreml(model=null.forma[[cc]])
           if (run$converge==FALSE) results$converge <- FALSE
                null.ll[cc] <- run$loglik
	}

	perm.test[ii,] <- 2*(full.ll-null.ll)
    }	# end of loops over permutations

  	# For each permutation, store the maximum (over chromosomes) LRT
	maxp <- apply(perm.test, 1, max)

	# Permutation threshold is the (1-alpha) percentile of max values
  	results$thresh <- sort(maxp[2:(n.perm+1)])[floor((1-input$alpha)*n.perm)]
	results$raw.pval <- sapply(perm.test[1,], pvfx)
  	results$adj.pval <- sapply(perm.test[1,], function(x) sum(x<=maxp[2:(n.perm+1)])/n.perm)

	results$obs <- perm.test[1,]
	results$perm.ts <- perm.test

  } # end of check whether n.perm>0
  } # end of algorithm==asreml
  
  if (algorithm=="lme") {
    fixed <- input$envModel$fixed
    f.mrk <- vector()
    chrRE <- vector()

    LRTStats <- vector() 

  # Construct vector of already mapped markers (f.mrk) 
  formula$fixed <- paste(as.character(fixed)[2], as.character(fixed)[1], as.character(fixed)[3], sep="")

  # Include fixed effects for all markers which have already been mapped
  if (length(prevLoc)>0)
  formula$fixed <- paste(formula$fixed, "+", paste(prevLoc, collapse="+"))

  formula$fixgrp <- paste(formula$fixed, "|grp1", sep="")
  formula$fixed <- as.formula(formula$fixed)
  formula$fixgrp <- as.formula(formula$fixgrp)

  gd <- groupedData(formula$fixgrp, data=dfMerged)

  # Random effects for all markers on a chromosome excluding those which
  # enter the model as fixed effects
  for (ii in 1:length(input$map))
   	chrRE[ii] <- paste("pdIdent(~", paste(setdiff(names(dfMerged)[grep(paste("C", ii, "M", sep=""), names(dfMerged))], prevLoc), collapse="+"), "-1)", sep="")

  formula$random <- paste("pdBlocked(list(", paste(chrRE[chrSet], collapse=","), "))", sep="")

  if (length(chrSet)==1)
  formula$random <- chrRE[chrSet]

  full <- lme(fixed=formula$fixed, random=eval(parse(text=formula$random)), data=gd, control=lmeControl(maxIter=input$maxit), na.action=na.omit)

  full.ll <- full$logLik

  # If there is only one chromosome in the subset, compare a full model to the model
  # with no random effects
  if (n.chrSet==1)
    null.forma[[1]] <- lme(fixed=formula$fixed, data=gd, control=lmeControl(maxIter=input$maxit), na.action=na.omit)
 
  # Otherwise, compare the full model to leave-one-VC-out models, removing 
  # each chromosome effect one at a time
  if (n.chrSet>1)
  for (cc in 1:n.chrSet)
  {
    random.forma[[cc]] <- paste("pdBlocked(list(", paste(chrRE[setdiff(chrSet, chrSet[cc])], collapse=","), "))", sep="")

    if (n.chrSet==2)
    random.forma[[cc]] <- chrRE[setdiff(chrSet, chrSet[cc])]
  
    # Fit the null model, where we omit the specified chromosome random effect
    null.forma[[cc]] <- lme(fixed=formula$fixed, random=eval(parse(text=random.forma[[cc]])), data=gd, control=lmeControl(maxIter=input$maxit), na.action=na.omit)

    null.ll[cc] <- null.forma[[cc]]$logLik
  }
 
  LRTStats <- 2*(full.ll-null.ll)
  results$obs <- LRTStats
  results$raw.pval <- sapply(LRTStats, pvfx)

  # depends on multtest value
  if (input$multtest=="bon")
  results$adj.pval <- sapply(results$raw.pval*n.chrSet, function(x) return(min(x, 1))) else {
	pval <- as.matrix(rbind(1:length(results$raw.pval), results$raw.pval))
	pval <- as.matrix(pval[,order(pval[2,])])
   	pval[2,] <- sapply(pval[2,]*(n.chrSet:1), function(x) return(min(x,1)))
	results$adj.pval <- pval[2,order(pval[1,])] }
  results$thresh <- qchibar(input$alpha/n.chrSet)

  }

  return(results)
}

