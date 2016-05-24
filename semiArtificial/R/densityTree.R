
# creates random tree
# all features are treated equally, there is no class/predicted value
densityTree<-function(dataset, minNodeWeight=2, noSelectedAttr=2, densitySplitMethod, densityData){
  
  
  noInst=nrow(dataset)
  noAttr = ncol(dataset)

  # stopping criterions: leaf size 
  if (noInst <= minNodeWeight) {
    return(createDensityLeaf(noInst))
  }
  
  uniqueVal <- list()
  
  if (densitySplitMethod == "maxVariance") {
    
    # select attribute with maximal variance among randomly selected noSelectedAttr
    attrs <- c()
    available <- sample(1:noAttr, size=noAttr, replace=FALSE)
    while (length(attrs) < noSelectedAttr && length(available)>0) {
      
      sel <- available[1]
      
      un <- unique(dataset[,sel])
      un <- un[!is.na(un)]
      if (length(un) > 1) {
        attrs <- c(attrs, sel)
      }
      else if (length(un) == 1)
        uniqueVal[[1+length(uniqueVal)]] <- list(attr=sel,val=un[1])
      
      available <- available[-1]
    }    
    #select the attribute with maximal variance from selected
    maxVar <- 0
    selAttr <- 0
    for (i in seq(along=attrs)){
      vals <- dataset[,attrs[i]]
      if (is.factor(vals)) 
        vals <- as.integer(vals)
      # normalize to [0,1] and compute variance
      minA <- min(vals)
      maxA <- max(vals)
      intLen <- maxA -minA
      if (intLen > 0) {
        vals <-(vals - minA)/intLen
        avar <- var(vals, na.rm=TRUE)
        if (avar > maxVar) {
          maxVar <- avar
          selAttr <- attrs[i]
        }
      }
    }
    if (selAttr != 0)
      node <-  makeDensityNodeSplit(dataset[,selAttr], minNodeWeight, densitySplitMethod)   
  }
  else if (densitySplitMethod %in% c("balancedSplit","randomSplit")) {
    # as node split randomly select an attribute which does not have all values equal
    attrOrder <- sample(1:noAttr, size=noAttr, replace=FALSE)
    selAttr <- 0
    sel <- 1
    while (selAttr == 0 && sel <= noAttr) {
      un <- unique(dataset[,attrOrder[sel]])
      un <- un[!is.na(un)]
      if (length(un) > 1) {  
        selAttr <- attrOrder[sel]      
        node <-  makeDensityNodeSplit(dataset[,selAttr], minNodeWeight, densitySplitMethod)
        if (node$nodeType == "leaf") {
          selAttr <- 0
          sel <- sel + 1
        }       
      }
      else {
        if (length(un)==1)
          uniqueVal[[1+length(uniqueVal)]] <- list(attr=attrOrder[sel],val=un[1])
        sel <- sel + 1;  
      }
    }
  } 
  if (selAttr == 0)
    return(createDensityLeaf(noInst)) ;  
  
  if (node$nodeType=="leaf") {
    return(node)
  }
  else { 
    node$splitAttr <- selAttr
  
    if (densityData == "topDown" || densityData == "bottomUp") {
      node$unique <- uniqueVal
      node$pLeft <- length(node$leftInstances) / nrow(dataset)
    }
    
     # recursive splitting
  
    node$leftTree <- densityTree(dataset[node$leftInstances,], minNodeWeight,  noSelectedAttr, densitySplitMethod, densityData)
    node$leftInstances <- NULL
    node$rightTree <- densityTree(dataset[node$rightInstances,], minNodeWeight, noSelectedAttr, densitySplitMethod, densityData)
    node$rightInstances <- NULL
    
    return(node)
  }
}


createDensityLeaf <- function(noInst) {
  list(nodeType = "leaf", noInst = noInst)
}


# create a node of a density tree and create a split
makeDensityNodeSplit <- function(attrData,  minNodeWeight, densitySplitMethod) {
  node <- list(nodeType = "interior")
  
  # splitData  contains the split attribute, 
  splitData <- attrData # to keep original       
  
  ctable <- table(splitData)
  
  if (length(ctable)<=1) {
    # a single non NA value, attribute not suitable for splitting
    return(createDensityLeaf(length(splitData)))
  }
  allSum <- sum(ctable)
  if (is.factor(splitData)) {
    # discrete attribute
    node$splitType = "discrete"
    attrLevels = levels(splitData)
    
    if (!is.ordered(splitData))
      ctable<-sample(ctable) # radom shuffle of values
    
    if (densitySplitMethod %in% c("balancedSplit","maxVariance")) {
      # the split shall select the number of elements closest to the half
      leftSum <- 0
      iLeft <- 0
      halfSum  <- allSum/2 ;
      while (leftSum < halfSum ) {
        iLeft <- iLeft +1
        leftSum <- leftSum + ctable[iLeft];
      }
      l1Sum = leftSum - ctable[iLeft]
      if (leftSum - halfSum > halfSum - l1Sum ) {
        leftSum = l1Sum ;
        iLeft <- iLeft - 1
      }
    } 
    else if (densitySplitMethod == "randomSplit") {
      if (length(ctable) == 2)
        iLeft <- 1
      else
        iLeft = ceiling(runif(1, min = 0, max = length(ctable)-1 )) #at least one from border of the table
      leftSum <- sum(ctable[1:iLeft])
      # assure the resulting split contains at least minNodeWeight instances
      while (leftSum <= minNodeWeight) {
        iLeft <- iLeft + 1
        leftSum <- leftSum + ctable[iLeft]
      }
      while (allSum - leftSum <= minNodeWeight && iLeft > 1){
        leftSum <- leftSum - ctable[iLeft]
        iLeft <- iLeft - 1
      }                  
    }
    else stop("makeDensityNodeSplit encountered an illegal densitySplitMethod ",densitySplitMethod)
    
    # assign left and right values
    node$leftValues = factor(names(ctable)[1:iLeft],levels=attrLevels)
    node$rightValues = factor(names(ctable)[(iLeft+1):length(ctable)],levels=attrLevels)
    node$NAvalue <- factor(names(ctable)[which.is.max(ctable)],levels=attrLevels) # replace with modus
    if (node$NAvalue %in% node$leftValues)
      node$leftValues[length(node$leftValues)+1] <- NA
    else node$rightValues[length(node$rightValues)+1] <-NA
    node$leftInstances <-which(splitData %in% node$leftValues)
    node$rightInstances <- which(splitData %in% node$rightValues)    
  }
  else {
    #numeric attribute
    node$splitType = "numeric"
    if (densitySplitMethod %in% c("balancedSplit","maxVariance")) {
      # the split shall select the number of elements closest to the half
      leftSum <- 0
      iLeft <- 0
      halfSum  <- allSum/2
      while (leftSum < halfSum ) {
        iLeft <- iLeft +1
        leftSum <- leftSum + ctable[iLeft];
      }
    }
    else if (densitySplitMethod == "randomSplit") {
      if (length(ctable) == 2)
        iLeft <- 1
      else
        iLeft = ceiling(runif(1, min = 0, max = length(ctable)-1 )) # at least one from border of the table
      leftSum <- sum(ctable[1:iLeft])
      # assure the resulting split contains at least minimalNodeWeight instances
      while (leftSum < minNodeWeight) {
        iLeft <- iLeft + 1
        leftSum <- leftSum + ctable[iLeft]
      }
      while (allSum - leftSum < minNodeWeight && iLeft > 1){
        leftSum <- leftSum - ctable[iLeft]
        iLeft <- iLeft - 1
      }                  
    }
    else stop("makeDensityNodeSplit encountered an illegal densitySplitMethod ", densitySplitMethod)
    
    #assign values
    if (iLeft < length(ctable))
      node$splitValue = (as.numeric(names(ctable)[iLeft])+as.numeric(names(ctable)[iLeft+1]))/2.0  
    else # iLeft == 1
      node$splitValue = (as.numeric(names(ctable)[length(ctable)-1])+as.numeric(names(ctable)[length(ctable)]))/2.0  
    
    node$NAvalue <- mean(splitData,na.rm=TRUE)
    if (node$NAvalue <= node$splitValue ) {
      node$leftInstances <- which(splitData <= node$splitValue | is.na(splitData))
      node$rightInstances <- which(splitData > node$splitValue)
    }
    else {
      node$leftInstances <- which(splitData <= node$splitValue)
      node$rightInstances <- which(splitData > node$splitValue | is.na(splitData)) 
    }
  }
  if (length(node$leftInstances) < minNodeWeight || length(node$rightInstances) < minNodeWeight)
    return(createDensityLeaf(length(splitData)))
  
  return(node)
}
