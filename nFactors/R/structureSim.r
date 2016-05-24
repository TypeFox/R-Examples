structureSim <-
function(fload, reppar=30, repsim=100, N, quantile=0.95, model="components",
         adequacy=FALSE, details=TRUE, r2limen=0.75, all=FALSE) {
 simulation   <- sim.structure(fx=fload, n=N, raw=TRUE)
 if (adequacy == TRUE) print(factanal(covmat=simulation$model, factors=dim(fload)[2])) # Verification of the adequacy of the model
 eigenvalues  <- eigenComputes(simulation$r, cor=TRUE, model=model)
 variables    <- length(eigenvalues) # Compute the number of variables
 aparallel    <- parallel(var=dim(fload)[1],subject=N,rep=reppar,cent=quantile,model=model)$eigen$qevpea  # The percentile
 components   <- matrix(NA, ncol=15,nrow=repsim)
 analysis     <- NA
 values       <- matrix(NA, ncol=length(eigenvalues),nrow=repsim)
 for (i in 1:repsim) {
  simulation             <- sim.structure(fx=fload, n=N, raw=TRUE)
  aparallel              <- parallel(var=dim(fload)[1],subject=N,rep=reppar,cent=quantile,model=model)$eigen$qevpea
  eigenvalues            <- eigenComputes(simulation$r, cor=TRUE, model=model)
  values[i,]             <- eigenvalues
  results                <- nScree(x=eigenvalues,aparallel = aparallel, cor=TRUE, model=model)
  components[i,(1:4)]    <- t(results$Components)
  ### PERMUTATIONS
  if (eigenFrom(data.frame(simulation$observed)) == "data")  {
   permutation <- eigenBootParallel(x=data.frame(simulation$observed), quantile=quantile, model=model)$quantile
   }
  results                <- nScree(x=eigenvalues,aparallel = permutation, cor=TRUE, model=model)
  components[i, 5]       <- results$Components$nparallel
  ### ...
  components[i, 6]       <- nCng(x=eigenvalues, model=model)$nFactors
  components[i, (7:9)]   <- nMreg(x=eigenvalues, model=model)$nFactors
  components[i, (10:11)] <- nSeScree(x=eigenvalues, model=model, r2limen=r2limen)$nFactors

  if (model == "components") {
   components[i, (12:14)] <- nBartlett(x=eigenvalues, N=N, alpha=1-quantile, cor=TRUE, correction=TRUE)$nFactors
   if (all == TRUE) {
    cat(paste("-- repsim = ", i, "/",repsim,"\n", sep=""))
    components[i, (15)] <- nBentler(x=eigenvalues, N=N, alpha=1-quantile, cor=TRUE)$nFactors
    }
   }
  # analysis       <- rbind(analysis, results$Analysis)
  #components[2,] <- t(results$Components);components
  }

 names                <- colnames(results$Components)
 names                <- c("oc", "af", "par", "mean.eig", "per")
 components           <- data.frame(components)
 colnames(components) <- c(names,"cng","b","t.b","p.b","sescree","R2","Bartlett","Anderson","Lawley","Bentler")
 if (details == TRUE) analysis <- list(components=components, eigenvalues=values)
 if (repsim > 1)      components <- moreStats(components, quantile=quantile) else components <- NA

 res <- list(details=analysis, nFactors=components)
 class(res) <- 'structureSim'
 return(res)
 }
## LIGNE 21 MODIFIEE: ETAIT quantile=0.95
## LIGNE 42 MODIFIEE: EAIT c("oc", "af", "par", "per", "mean.eig")