# This set of function extracts the significant test restults, using different methods.

# pvalues: The vector of p-values
# method: The multiple testing adjustment method. Options are
#	plain: take simply the alpha as cut-off
#	bonferroni: Take alpha/length(p) as cut-off
#	simes: The simes extension of the bonferroni approach
#	bh: The benjamini-hochberg approach to control the FDR
#	csD: This one extracts all the significant tests for which the distance between expected and observed is largest
#	csR: This one extracts all the significant tests for which the ratio between expected and observed is largest
# alpha: Parameter for the multiple testing

# Version: 21-05-2013:
#	* Added the option for 'csD'
#	* Added the option for 'csR'
# Version: 30-06-2013:
#	* Start with the Westfall & Young approach, DF 
# Version: 03-07-2013:
#	* Finnished the Westfall & Young (maxT) approach, DF
#	* Adjusted the input for the remaining methods, according to the changes in the object-class gmw
# Version: 19-09-2013:
#       * Fixed a problem with the isVector flag

getSigTests <- function(pvalues, alpha=0.05, method="plain"){

# Input checks
  method <- match.arg(method,c("plain","bonferroni","simes","bh","csD","csR","maxT"))

# Set constants
  ifelse((is.numeric(pvalues)) & (!is.matrix(pvalues)),isVector <- TRUE, isVector <- FALSE   ) 

# Internal function, needed in case of matrix input
  getSigTests.internal <- function(pvalues, alpha, method){

  # Create slots for the return
    multAlpha <- NULL
    sigPvalues <- NULL    
    sigTests <- NULL     
	
  # No multiple adjustment, jut extraact all p-values that are smaller then a certain thresholds
    if(method=="plain")
    {
      # Extract the plain p-values that are smaller than a certain threshold alpha
	sigTests <- which((pvalues<=alpha)==TRUE)
      # Get the significant p-values of these positions:
	sigPvalues <- pvalues[sigTests]
      
  # Now the Bonferroni correction:
    } else if(method=="bonferroni") {
      # adjust the Alpha
	bonAlpha <- alpha/length(pvalues)
      # Extract the plain p-values that are smaller than the adjusted threshold bonAlpha
	sigTests <- which((pvalues<=bonAlpha)==TRUE)
      # Get the significant p-values of these positions:
	sigPvalues <- pvalues[sigTests]    
	multAlpha <- bonAlpha
      
  # The Simes extension of the Bonferroni correction:
    } else if (method=="simes"){
      # Sort the p-values
	pvalues.sorted <- sort(pvalues)
      # Adjust the Alpha
	bonAlpha <- alpha/length(pvalues)
      # Adjust the Bonferroni alpha to the adjusted values of the Simes approach:
	simAlpha <- bonAlpha * 1:length(pvalues)
      # get the i_max value
	imax <- which((pvalues.sorted <= simAlpha)==TRUE)
	if(sum(imax)>0)
	{
	    imax <- max(imax)
	  # Extract the plain p-values that are smaller than the adjusted threshold simAlpha
	    sigTests <- which((pvalues<=imax*bonAlpha)==TRUE)
	  # Get the significant p-values of these positions:
	    sigPvalues <- pvalues[sigTests]    
	    multAlpha <- imax*bonAlpha
	}
	
    } else if(method=="bh"){
      # Sort the p-values
	pvalues.sorted <- sort(pvalues)
      # Calculated the largest i, as in the formular to control the FDR, depending on alpha:
	kmax <- which((pvalues.sorted <= alpha*(1:length(pvalues))/length(pvalues))==TRUE)
	if(sum(kmax)>0)
	{
	    kmax <- max(kmax)
	  # Get the significant p-values of these positions:
	    sigPvalues <- pvalues.sorted[1:kmax]    
	    multAlpha <- pvalues.sorted[kmax]
	  # Get the original positions of the significant tests:
	    sigTests <- which(is.element(pvalues,sigPvalues))     
	} 
    } else if(method=="csD"){
      ## WARNING, I THINK WE SHOULD TAKE HERE FIRST THE ALPHA, TEST WITH THE CS CRITERIA IF THERE IS A DIFFERENCE IN THE P-VALUES
      ## AND THEN CHECK ONLY THOSE POSSIBLE ALPHA.GRID CUTOFF WHICH HAVEN BEEN SMALLER THAN THE POSITION WERE THE CS TEST REJECTS!
    
      # Set possible cut-off values
	alphaGrid <- seq(0,alpha,0.001)
	csDMat <- matrix(NA,ncol=length(alphaGrid),nrow=2)
      # check for each alpha the suitable value
	for(i in 1:length(alphaGrid)){
	  csDMat[1,i] <- sum(pvalues<=alphaGrid[i])
	  csDMat[2,i] <- alphaGrid[i]*length(pvalues)
	}
      # Now check the optimal cut-off:
	csD <- csDMat[1,]-csDMat[2,]
      # Get the optimal alpha:
	multAlpha <- alphaGrid[which(csD==max(csD))]
      # Get the significant tests and p-values:
	sigTests <- which(pvalues<=multAlpha)
	sigPvalues <- pvalues[pvalues<=multAlpha]
	
    } else if(method=="csR"){
      ## WARNING, I THINK WE SHOULD TAKE HERE FIRST THE ALPHA, TEST WITH THE CS CRITERIA IF THERE IS A DIFFERENCE IN THE P-VALUES
      ## AND THEN CHECK ONLY THOSE POSSIBLE ALPHA.GRID CUTOFF WHICH HAVEN BEEN SMALLER THAN THE POSITION WERE THE CS TEST REJECTS!
    
      # Set possible cut-off values
	alphaGrid <- seq(0.001,alpha,0.001)
	csRMat <- matrix(NA,ncol=length(alphaGrid),nrow=2)
      # check for each alpha the suitable value
	for(i in 1:length(alphaGrid)){
	  csRMat[1,i] <- sum(pvalues<=alphaGrid[i])
	  csRMat[2,i] <- alphaGrid[i]*length(pvalues)
	}
      # Now check the optimal cut-off:
	csR <- csRMat[1,]/csRMat[2,]
      # Get the optimal alpha:
	multAlpha <- alphaGrid[which(csR==max(csR))]
      # Get the significant tests and p-values:
	sigTests <- which(pvalues<=min(multAlpha))
	sigPvalues <- pvalues[pvalues<=multAlpha]
    } 
    
    result <- list(sigTests=sigTests, sigPvalues=sigPvalues, pvalues=pvalues, method=method, alpha=alpha, multAlpha=multAlpha)
    result
  }

  getSigTests.WY <- function(obs, pm, alpha, alternative, pvalues){
    Ntests <- length(obs)
    # In case the alternative is smaller we have to use some kind of minT approach, otherwise the regular maxT
    if(alternative=="smaller"){
      obs.sorted <- sort(obs,decreasing=FALSE)
      pm.sorted <- pm[,order(obs,decreasing=FALSE)]
      multAlpha <- c()
      for(i in 1:(Ntests-1)){
	temp <- apply(pm.sorted[,(i:Ntests)],1,min)
	multAlpha[i] <- sum(temp<obs.sorted[i])/nrow(pm)
      }
      multAlpha[Ntests] <- sum(pm.sorted[,Ntests] < obs.sorted[Ntests])/nrow(pm)
      multAlpha <- cummax(multAlpha)[rank(obs)]
    } else {
      obs.sorted <- sort(obs,decreasing=TRUE)
      pm.sorted <- pm[,order(obs,decreasing=TRUE)]
      multAlpha <- c()
      for(i in 1:(Ntests-1)){
	temp <- apply(pm.sorted[,(i:Ntests)],1,max)
	multAlpha[i] <- sum(temp>obs.sorted[i])/nrow(pm)
      }
      multAlpha[Ntests] <- sum(pm.sorted[,Ntests] > obs.sorted[Ntests])/nrow(pm)
      multAlpha <- cummax(multAlpha)[rank(-obs)]
    }
    sigTests <- which(multAlpha<=alpha)
    sigPvalues <- pvalues[multAlpha<=alpha]
    names(multAlpha) <- colnames(pvalues) 
    result <- list(sigTests=sigTests, sigPvalues=sigPvalues, pvalues=pvalues, method=method, alpha=alpha, multAlpha=multAlpha)
    result
  }

# In case that we give a matrix of p-values, run this function for each row and create a capsulated list
  if(method=="maxT"){
    if(class(pvalues)!="gmw") stop ("The pvalues object has to be of class gmw, please put here the outpt of the function gmw.")
    if(attr(pvalues,"keepPM")!=TRUE) stop ("For the Westfall & Young approach the permutation matrix has to be available, please choose the option 'keepPM=TRUE' in the gmw call.")
    alternative <- attr(pvalues,"alternative")

    if(is.matrix(pvalues$p.values) && nrow(pvalues$p.values)==1) isVector <- TRUE # pvalues <- as.vector(pvalues$p.values)
    
    if(!isVector){
      result <- list()
      for(i in 1:nrow(pvalues$p.values))
      {
	result[[i]] <- getSigTests.WY(pvalues$obsValue[[i]],pvalues$nullDist[[i]], alpha=alpha, alternative=alternative, pvalues=t(as.matrix(pvalues$p.values[i,])))
      }
      result$inputN <- nrow(pvalues$p.values)
    } else {
      result <- getSigTests.WY(pvalues$obsValue[[1]],pvalues$nullDist[[1]], alpha=alpha, alternative=alternative, pvalues=pvalues$p.values)
      result$inputN <- 1
    }
  } else {
    if(class(pvalues)!="gmw"){
      temp <- pvalues
      pvalues <- list()
      pvalues$p.values <- temp
    }
    if(is.matrix(pvalues$p.values) && nrow(pvalues$p.values)==1) isVector <- TRUE # pvalues <- as.vector(pvalues$p.values) 
    if(!isVector){
      result <- list()
      for(i in 1:nrow(pvalues$p.values))
      {
	 result[[i]] <- getSigTests.internal(pvalues$p.values[i,], alpha=alpha,method)
      }
      result$inputN <- nrow(pvalues$p.values)
    } else {
      result <- getSigTests.internal(pvalues$p.values,alpha,method)
      result$inputN <- 1
    }
  }
  class(result) <- "re"
  result
}
