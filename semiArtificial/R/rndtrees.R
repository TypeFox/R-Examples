# creates random tree
# response is expected to be in the last column
randomTree<-function(dataset, minNodeWeight, noSelectedAttr, 
                     problemType, densityData, estimator) {
  
  node <- list() # empty node
  
  noInst <- nrow(dataset)
  noAttr <- ncol(dataset)-1
  
  classIdx <- ncol(dataset) # response is expected to be in the last column
  
  if (problemType == "classification") {
    classStat <- table(dataset[,classIdx]) # count number of instance of each class
    # stopping criterions: leaf size or pure leaf
    if (noInst <= minNodeWeight || max(classStat) == sum(classStat)){
      return(createClassificationLeaf(noInst, classStat, levels(dataset[[classIdx]])))
    }
  } else {
    classStat <- mean(dataset[,classIdx])
    # stopping criterions: leaf size or pure leaf
    if (noInst <= minNodeWeight || max(dataset[,classIdx]) == min(dataset[,classIdx])){
      return(createRegressionLeaf(noInst, classStat))
    }
  }
  
  # select attributes, skip the ones with all values equal
  attrs <- c()
  uniqueVal <- list() # features with unique values only, used in tree
  available <- sample(1:noAttr,size=noAttr,replace=FALSE)
  while (length(attrs) < noSelectedAttr && length(available)>0) {    
    sel <- available[1]
    un <- unique(dataset[,sel])
    un <- un[!is.na(un)]
    if (length(un) > 1) {
      attrs <- c(attrs, sel)
    }
    else if (length(un)==1)
      uniqueVal[[1+length(uniqueVal)]] <- list(attr=sel,val=un[1])
    
    available <- available[-1]
  }
  
  if  (length(attrs)==0) {
    if (problemType == "classification") 
      return(createClassificationLeaf(noInst, classStat, levels(dataset[[classIdx]])))
    else 
      return(createRegressionLeaf(noInst, classStat))  
  }
  
  # construct formula from selected attributes
  frml <- paste(names(dataset)[classIdx], "~",paste(names(dataset)[attrs], sep="+",collapse="+"),sep="") 
  newFormula <- as.formula(frml)   
  
  # evaluate selected attributes
  ae <- attrEval(newFormula, dataset, estimator=estimator, outputNumericSplits=TRUE)
  
  selAttr <- attrs[which.is.max(ae$attrEval)]
  selName <- names(ae$attrEval)[which.is.max(ae$attrEval)]
  numSplit <- NA
  if (! is.factor(dataset[[selAttr]])) # numeric feature is selected
    numSplit <- ae$splitPointNum[selName]
  
  if (problemType == "classification") 
    node <-  makeClassificationNodeSplit(dataset[,c(selAttr, classIdx)], classStat, minNodeWeight, numSplit)
  else 
    node <-  makeRegressionNodeSplit(dataset[,c(selAttr, classIdx)], classStat, minNodeWeight, numSplit)
  
  
  if (node$nodeType=="leaf") {
    return(node)
  }
  else { 
    node$splitAttr <- selAttr
    
    if (densityData == "topDown") {
      node$unique <- uniqueVal
    }
    else if (densityData == "bottomUp") {
      node$unique <- uniqueVal
      node$pLeft <- length(node$leftInstances) / nrow(dataset)
    }
    
    # recursive splitting
    node$leftTree <- randomTree(dataset[node$leftInstances,], minNodeWeight, noSelectedAttr, problemType, densityData, estimator)
    node$leftInstances <- NULL
    node$rightTree <- randomTree(dataset[node$rightInstances,], minNodeWeight, noSelectedAttr, problemType, densityData, estimator)
    node$rightInstances <- NULL
    
    return(node)
  }
}


createClassificationLeaf <- function(noInst, classStat, classLevels) {
  list(nodeType = "leaf", noInst = noInst, classProb = classStat/noInst, majorityClass = factor(classLevels[which.is.max(classStat)], levels=classLevels))
}

createRegressionLeaf <- function(noInst, meanResponse) {
  list(nodeType = "leaf", noInst = noInst, meanResponse = meanResponse)
}


# create a node of a decision tree and split the instances
makeClassificationNodeSplit <- function(attrClassData, classStat, minNodeWeight, numSplit) {
  node <- list(nodeType = "interior")
  # splitData contains two columns,first contains the split attribute, the second contains class values
  
  splitData <- attrClassData # to keep original
  classLevels <- levels(attrClassData[[2]])
  #try to dichotomize between majority class and all the others
  majorityClass <- levels(attrClassData[[2]])[which.is.max(classStat)]
  
  if (length(levels(splitData[[2]])) > 2) { # number of classes > 2                                 
    # merge non-majority classes into one class
    classData <- integer(length=nrow(splitData))
    classData[splitData[,2] == majorityClass] <- 1
    classData[splitData[,2] != majorityClass] <- 2
    splitData[[2]] <- factor(classData)
  }                    
  
  ctable <- table(splitData)
  if ( nrow(ctable) <= 1 ) {
    # a single non NA value, attribute not suitable for splitting
    return(createClassificationLeaf(nrow(splitData), classStat, classLevels ))
  }
  colSum <- apply(ctable, 2, sum, na.rm=T)
  allSum <- sum(colSum)
  if (is.factor(splitData[[1]])) {
    # discrete attribute
    node$splitType <- "discrete"
    attrLevels <- levels(splitData[[1]])
    if (nrow(ctable)==2) {
      # binary split
      node$leftValues <- factor(row.names(ctable)[1], levels=attrLevels)
      node$rightValues <- factor(row.names(ctable)[2], levels=attrLevels)
      rowSum <- apply(ctable, 1, sum) # needed below
      if (min(rowSum) <= minNodeWeight) {
        return(createClassificationLeaf(nrow(splitData), classStat, classLevels))        
      }
    }
    else {
      # find optimal split
      rowSum <- apply(ctable, 1, sum)
      ptable <- ctable / rowSum
      valueOrder <- order(ptable[,1])
      gini <- 0 ;
      minGini <- maxValue
      left <- c(0,0) # class value distribution in left split
      right <- colSum # class value distribution in right split
      for (i in 1:(nrow(ctable)-1)){
        left <- left + ctable[valueOrder[i],]
        sumleft <- sum(left)
        pleft <- left / sumleft
        right <- right - ctable[valueOrder[i],]
        sumright <- sum(right)
        pright <- right / sum(right)
        gini <- sumleft/allSum * ( 1 - (sqr(pleft[1]) + sqr(pleft[2]))) + sumright / allSum * (1 - (sqr(pright[1])+sqr(pright[2])))
        if (gini < minGini && min(c(sumleft, sumright)) > minNodeWeight) {
          minGini = gini
          node$leftValues = factor(row.names(ctable)[valueOrder[1:i]],levels=attrLevels)
          node$rightValues = factor(row.names(ctable)[valueOrder[(i+1):nrow(ctable)]],levels=attrLevels)
        }
      }
      if (minGini == maxValue) { # no suitable split
        return(createClassificationLeaf(nrow(splitData), classStat, classLevels))        
      }        
    }
    node$NAvalue <- factor(row.names(ctable)[which.is.max(rowSum)],levels=attrLevels) # replace with modus
    if (node$NAvalue %in% node$leftValues)
      node$leftValues[length(node$leftValues)+1] <- NA
    else node$rightValues[length(node$rightValues)+1] <-NA
    node$leftInstances <-which(splitData[,1] %in% node$leftValues)
    node$rightInstances <- which(splitData[,1] %in% node$rightValues)    
  }
  else {
    #numeric attribute
    node$splitType <- "numeric"
    if ( is.na(numSplit) ) { # invalid value was returned by attrEval or inappropriate estimator
      minGini <- maxValue
      left <- c(0,0) # class value distribution in left split
      right <- colSum # class value distribution in right split
      for (i in 1:(nrow(ctable)-1)){
        left <- left + ctable[i,]
        sumleft <- sum(left)
        pleft <- left / sumleft
        right <- right - ctable[i,]
        sumright <- sum(right)
        pright <- right / sum(right)
        gini <- sumleft/allSum * ( 1 - (sqr(pleft[1]) + sqr(pleft[2]))) + sumright / allSum * (1 - (sqr(pright[1])+sqr(pright[2])))
        if (gini < minGini && min(c(sumleft, sumright)) > minNodeWeight) {
          minGini <- gini
          node$splitValue <- (as.numeric(row.names(ctable)[i])+as.numeric(row.names(ctable)[i+1]))/2.0
        }
      }
      if (minGini == maxValue) { # no suitable split
        return(createClassificationLeaf(nrow(splitData), classStat, classLevels))
      } 
    } 
    else {  # a valid split was provided by attrEval 
      node$splitValue <- numSplit 
    }
    node$NAvalue <- mean(splitData[[1]],na.rm=TRUE)
    if (node$NAvalue <= node$splitValue ) {
      node$leftInstances <- which(splitData[,1] <= node$splitValue | is.na(splitData[,1]))
      node$rightInstances <- which(splitData[,1] > node$splitValue)
    }
    else {
      node$leftInstances <-which(splitData[,1] <= node$splitValue)
      node$rightInstances <- which(splitData[,1] > node$splitValue | is.na(splitData[,1])) 
    }		
  }
  if (length(node$leftInstances) < minNodeWeight || length(node$rightInstances) < minNodeWeight)
    return(createClassificationLeaf(nrow(attrClassData), classStat, levels(attrClassData[[2]])))
  
  return(node)
}


# create a node of a regression tree and split the instances
makeRegressionNodeSplit <- function(attrClassData, meanResponse, minNodeWeight, numSplit) {
  node <- list(nodeType = "interior")
  
  # splitData contains two columns,first contains the split attribute, the second contains the response
  splitData <- attrClassData # to keep original
  
  if (is.factor(splitData[[1]])) { 
    # discrete attribute
    node$splitType <- "discrete"
    attrLevels <- levels(splitData[[1]])
    #ctable <- table(splitData[,c(2,1)])
    ctable <- table(splitData)
    rowSum <- apply(ctable, 1, sum) # needed below
    if (nrow(ctable)==2) {
      # binary split
      node$leftValues <- factor(row.names(ctable)[1], levels=attrLevels)
      node$rightValues <- factor(row.names(ctable)[2], levels=attrLevels)
      if (min(rowSum) <= minNodeWeight) {
        return(createRegressionLeaf(nrow(splitData), meanResponse))        
      }
    }
    else {
      # find optimal split
      valueMean <- vector(mode="numeric", length=nrow(ctable))
      valueWeight <- vector(mode="numeric", length=nrow(ctable))
      valueSum <- vector(mode="numeric", length=nrow(ctable))
      valueMean <- vector(mode="numeric", length=nrow(ctable))
      valueSquares <- vector(mode="numeric", length=nrow(ctable))
      for (i in 1:nrow(ctable)) {
        nonZero <- ctable[i, ctable[i,] > 0 ]
        values <- as.numeric(names(nonZero))
        weightedValues <- nonZero * values
        valueSum[i] <- sum(weightedValues)
        valueWeight[i] <- sum(nonZero)
        valueMean[i] <- valueSum[i] / valueWeight[i]
        valueSquares[i] <- sum(weightedValues * values)
      }
      valueOrder <- order(valueMean)
      weightAll <- sum(valueWeight)
      minMSE <- maxValue
      leftValues <- 0 
      leftWeight <- 0
      leftSquares <- 0
      rightValues <- sum(valueSum) 
      rightWeight <- sum(valueWeight)
      rightSquares <- sum(valueSquares)
      for (i in 1:(nrow(ctable)-1)){
        leftValues <- leftValues + valueSum[valueOrder[i]]
        rightValues <- rightValues  - valueSum[valueOrder[i]]
        leftWeight <- leftWeight + valueWeight[valueOrder[i]]
        rightWeight <- rightWeight  - valueWeight[valueOrder[i]]
        leftSquares <- leftSquares + valueSquares[valueOrder[i]]
        rightSquares <- rightSquares  - valueSquares[valueOrder[i]]
        pLeft <- leftWeight / weightAll  
        mse <-     pLeft * (leftSquares/leftWeight - (leftValues/leftWeight)^2 ) + 
          (1-pLeft) * (rightSquares/rightWeight - (rightValues/rightWeight)^2)
        if (mse < minMSE && min(c(leftWeight, rightWeight)) > minNodeWeight) {
          minMSE = mse
          node$leftValues = factor(row.names(ctable)[valueOrder[1:i]],levels=attrLevels)
          node$rightValues = factor(row.names(ctable)[valueOrder[(i+1):nrow(ctable)]],levels=attrLevels)
        }
      }
      if (minMSE == maxValue) { # no suitable split
        return(createRegressionLeaf(nrow(splitData), meanResponse))        
      }        
    }
    node$NAvalue <- factor(row.names(ctable)[which.is.max(rowSum)],levels=attrLevels) # replace with modus
    if (node$NAvalue %in% node$leftValues)
      node$leftValues[length(node$leftValues)+1] <- NA
    else node$rightValues[length(node$rightValues)+1] <-NA
    node$leftInstances <-which(splitData[,1] %in% node$leftValues)
    node$rightInstances <- which(splitData[,1] %in% node$rightValues)    
  }
  else {
    #numeric attribute
    node$splitType <- "numeric"
    node$NAvalue <- mean(splitData[[1]],na.rm=TRUE)
    if ( is.na(numSplit) ) { # invalid value was returned by attrEval or inappropriate estimator
      minMSE <- maxValue
      splitData <- splitData[!is.na(splitData[,1]), ]  
      splitData <- splitData[do.call(order, splitData),]
      leftValues <- 0 
      leftWeight <- 0
      leftSquares <- 0
      rightValues <- sum(splitData[[2]])
      rightWeight <- nrow(splitData) 
      rightSquares <- sum(splitData[[2]]^2) 
      weightAll <- nrow(splitData)
      i <- 1
      while ( i < nrow(splitData)) {
        lower <- i
        i <- i + 1
        while(i < nrow(splitData)-minNodeWeight && splitData[lower, 1] == splitData[i,1])
          i <- i + 1
        leftValues <- leftValues + sum(splitData[lower:i,2])
        leftWeight <- leftWeight + (i-lower)
        leftSquares <- leftSquares + sum(splitData[lower:i,2]^2)
        rightValues <- rightValues - sum(splitData[lower:i,2])
        rightWeight <- rightWeight - (i - lower)
        rightSquares <- rightSquares - sum(splitData[lower:i,2]^2)
        pLeft <- leftWeight / weightAll  
        mse <-     pLeft * (leftSquares/leftWeight - (leftValues/leftWeight)^2 ) + 
          (1-pLeft) * (rightSquares/rightWeight - (rightValues/rightWeight)^2)
        if (mse < minMSE && min(c(leftWeight, rightWeight)) > minNodeWeight) {
          minMSE = mse
          node$splitValue <- (splitData[lower,1] + splitData[i,1])/2.0
        }
      }
      if (minMSE == maxValue) { # no suitable split
        return(createRegressionLeaf(nrow(splitData), meanResponse))
      } 
    } 
    else {  # a valid split was provided by attrEval 
      node$splitValue <- numSplit 
    }
    
    if (node$NAvalue <= node$splitValue ) {
      node$leftInstances <- which(splitData[,1] <= node$splitValue | is.na(splitData[,1]))
      node$rightInstances <- which(splitData[,1] > node$splitValue)
    }
    else {
      node$leftInstances <-which(splitData[,1] <= node$splitValue)
      node$rightInstances <- which(splitData[,1] > node$splitValue | is.na(splitData[,1])) 
    }		
  }
  if (length(node$leftInstances) < minNodeWeight || length(node$rightInstances) < minNodeWeight)
    return(createRegressionLeaf(nrow(attrClassData), meanResponse))
  
  return(node)
}


predictWithTree<-function(node, instance) {
  if (node$nodeType=="leaf") {
    if (! is.null(node$majorityClass))
      return(list(majorityClass=node$majorityClass, classProb = node$classProb))
    else 
      return(node$meanResponse) 
  }
  else {
    if (node$splitType == "discrete"){
      if (instance[1,node$splitAttr] %in% node$leftValues)
        return(predictWithTree(node$leftTree, instance))
      else return(predictWithTree(node$rightTree, instance))
    }
    else{ # numeric split
      if (is.na(instance[1,node$splitAttr])){
        if (node$NAvalue <= node$splitValue) 
          return(predictWithTree(node$leftTree, instance))
        else return(predictWithTree(node$rightTree, instance))
      }
      else if (instance[1, node$splitAttr] <= node$splitValue) {
        return(predictWithTree(node$leftTree, instance))
      }
      else return(predictWithTree(node$rightTree, instance))       
    }
  }
}


fillDataWithTreeTopDown <- function(node, instance) {
  if (node$nodeType=="leaf") {
    # if class exist, check if it is already assigned
    if (!is.null(node$classProb) && is.na(instance[length(instance)])) { 
      instance[length(instance)] <- quant(runif(1,0,1),node$classProb)
    }
    else if (!is.null(node$meanResponse) && is.na(instance[length(instance)]) && 
              any(class(node$cdf) %in% c("ecdf","logspline","kde"))) 
    {
      if (inherits(node$cdf, "ecdf") )
        instance[length(instance)] <- quantile(node$cdf, probs=runif(1, 0, 1), type=8)
      else if (inherits(node$cdf, "logspline"))
        instance[length(instance)] <- rlogspline(1, node$cdf) 
      else if (inherits(node$cdf,"kde"))
        instance[length(instance)] <- rkde(1, node$cdf)
    }      
    return(instance)
  }
  else {
    # assign value if not already set
    if (is.na(instance[node$splitAttr]) &&  any(class(node$cdf) %in% c("ecdf","logspline","kde"))) { 
      rnd <- runif(1, 0, 1)
      if (node$splitType=="discrete"){
        instance[node$splitAttr] <- quantile(node$cdf, probs = rnd, type = 3) # discontinuous sample quantile, SAS definintion
      }
      else {
        if (inherits(node$cdf, "ecdf") )
          instance[node$splitAttr] <- quantile(node$cdf, probs=runif(1, 0, 1), type=8)
        else if (inherits(node$cdf, "logspline"))
          instance[node$splitAttr] <- rlogspline(1, node$cdf) 
        else if (inherits(node$cdf,"kde"))
          instance[node$splitAttr] <- rkde(1, node$cdf)
      }
    }
    # set values which were selected and unique
    for (i in seq(along.with=node$unique))  {
      if (is.na(instance[node$unique[[i]]$attr]))
        instance[node$unique[[i]]$attr] <- as.numeric(node$unique[[i]]$val)
    }
    if (node$splitType == "discrete"){
      value <- factor(instance[node$splitAttr], levels = 1:length(levels(node$leftValues)), labels=levels(node$leftValues))
      if (value %in% node$leftValues)
        return(fillDataWithTreeTopDown(node$leftTree, instance))
      else if (value %in% node$rightValues)
        return(fillDataWithTreeTopDown(node$rightTree, instance))
      else stop("fillDataWithTreeTopDown encountered an unrecognized attribute value in tree node comparison",value)
    }
    else{ # numeric split
      if (instance[node$splitAttr] <= node$splitValue) 
        return(fillDataWithTreeTopDown(node$leftTree, instance))
      else return(fillDataWithTreeTopDown(node$rightTree, instance))       
    }
  }
}

fillDataWithTreeBottomUp <- function(node, instance) {
  if (node$nodeType=="leaf") {
    # if class exists, check if it is already assigned
    if (!is.null(node$classProb) && is.na(instance[length(instance)])){
      instance[length(instance)] <- quant(runif(1,0,1),node$classProb)
    } 
    else if (!is.null(node$meanResponse) && is.na(instance[length(instance)]) && any(class(node$cdf) %in% c("ecdf","logspline","kde"))){
      if (inherits(node$cdf, "ecdf") )
        instance[length(instance)] <- quantile(node$cdf, probs=runif(1, 0, 1), type=8)
      else if (inherits(node$cdf, "logspline"))
        instance[length(instance)] <- rlogspline(1, node$cdf) 
      else if (inherits(node$cdf,"kde"))
        instance[length(instance)] <- rkde(1, node$cdf)
    }      	
    return(instance)
  }
  else {
    
    # set values which were selected and unique
    for (i in seq(along.with=node$unique))  {
      if (is.na(instance[node$unique[[i]]$attr]))
        instance[node$unique[[i]]$attr] <- as.numeric(node$unique[[i]]$val)
    }
    
    if (is.na(instance[node$splitAttr])) {
      # value missing 
      if (runif(1,0,1) <= node$pLeft) { 
        #going left
        instReturned <- fillDataWithTreeBottomUp(node$leftTree, instance)
        if (is.na(instReturned[node$splitAttr]) &&  any(class(node$cdfLeft) %in% c("ecdf","logspline","kde"))) { # value still missing on return
          if (node$splitType == "discrete") {              
            instReturned[node$splitAttr] <- quantile(node$cdfLeft, probs = runif(1, 0, 1), type = 3)
          }
          else  {
            if (inherits(node$cdf, "ecdf") )
              instReturned[node$splitAttr] <- quantile(node$cdfLeft, probs=runif(1, 0, 1), type=8)
            else if ( inherits(node$cdf, "logspline"))
              instReturned[node$splitAttr] <- rlogspline(1, node$cdfLeft) 
            else if ( inherits(node$cdf, "kde"))
              instReturned[node$splitAttr] <- rkde(1, node$cdfLeft)
          }
        }
      }
      else { #going right
        instReturned <- fillDataWithTreeBottomUp(node$rightTree, instance)
        if (is.na(instReturned[node$splitAttr]) &&  any(class(node$cdfRight) %in% c("ecdf","logspline","kde"))) {  # value still missing on return
          if (node$splitType == "discrete")  {             
            instReturned[node$splitAttr] <- quantile(node$cdfRight, probs = runif(1, 0, 1), type = 3)
          }
          else  {
            if ( inherits(node$cdf,"ecdf") )
              instReturned[node$splitAttr] <- quantile(node$cdfRight, probs=runif(1, 0, 1), type=8)
            else if ( inherits(node$cdf,"logspline"))
              instReturned[node$splitAttr] <- rlogspline(1, node$cdfRight) 
            else if ( inherits(node$cdf, "kde") )
              instReturned[node$splitAttr] <- rkde(1, node$cdfRight)
          }
        }
      }
      return(instReturned) # return further up
    }
    else { # value not missing
      if (node$splitType == "discrete")  {        
        value <- factor(instance[node$splitAttr], levels = 1:length(levels(node$leftValues)), labels=levels(node$leftValues))
        if (value %in% node$leftValues) 
          return(fillDataWithTreeBottomUp(node$leftTree, instance))
        else if (value %in% node$rightValues)
          return(fillDataWithTreeBottomUp(node$rightTree, instance))
        else stop("fillDataWithTreeBottomUp encountered an unrecognized attribute value in tree node comparison ",value)
      }
      else { # numeric attribute
        if (instance[node$splitAttr] <= node$splitValue) 
          return(fillDataWithTreeBottomUp(node$leftTree, instance))
        else return(fillDataWithTreeBottomUp(node$rightTree, instance))       
      }
    }
  }
}

fillWithInstances <- function(tree, dat, instIdx) {
  for (i in 1:nrow(dat)) {
    tree <- fillWithInstance(tree, dat[i,], instIdx[i])
  }
  tree
} 

fillWithInstance <- function(node, instance, instIdx){
  
  if (is.null(node$instances))
    node$instances <- c(instIdx)
  else node$instances <- c(node$instances, instIdx)
  
  if (node$nodeType != "leaf") {
    
    if (node$splitType == "discrete"){
      if (instance[1,node$splitAttr] %in% node$leftValues)
        node$leftTree <- fillWithInstance(node$leftTree, instance, instIdx)
      else 
        node$rightTree <- fillWithInstance(node$rightTree, instance, instIdx)
    }
    else{ # numeric split
      if (is.na(instance[1,node$splitAttr])){
        if (node$NAvalue <= node$splitValue) 
          node$leftTree <- fillWithInstance(node$leftTree, instance, instIdx)
        else 
          node$rightTree <- fillWithInstance(node$rightTree, instance, instIdx)
      }
      else if (instance[1,node$splitAttr] <= node$splitValue) {
        node$leftTree <- fillWithInstance(node$leftTree, instance, instIdx)
      }
      else 
        node$rightTree <- fillWithInstance(node$rightTree, instance, instIdx)       
    }
  }  
  return(node)
}

genFromLeaves<-function(tree, dat, cdfEstimation, ...) {
  leaf <- list()
  leaf <- collectLeaves(tree, leaf) 
  generator <- list(leaves=list(),weights=c())
  for (i in 1:length(leaf)) {
    generator$weights <- c(generator$weights, length(leaf[[i]]$instances))
    cdfs <- list()
    for (a in 1:ncol(dat)) {
      vals <- dat[leaf[[i]]$instances, a]
      if (all(is.na(vals)))
        cdfs[[a]] <- NA
      else if ( cdfEstimation == "ecdf" || is.factor(dat[[a]]) )
        cdfs[[a]] <- ecdf(vals)
      else if (cdfEstimation == "logspline")                             
        cdfs[[a]] <- robustLogspline(vals, ...)
      else if (cdfEstimation == "kde")
        cdfs[[a]] <- robustKde(vals, ...)
    }
    generator$leaves[[i]] <- list(cdfs = cdfs)
  }
  generator$weights <- generator$weights / sum(generator$weights)
  
  return(generator)
}

generatorFromTree<-function(node, dat, densityData, cdfEstimation, ...) {
  if (node$nodeType == "leaf") {    
    vals <- dat[node$instances, ncol(dat)]
    if (all(is.na(vals)))
      node$cdf <- NA
    else if (cdfEstimation == "ecdf") 
      node$cdf <- ecdf(vals)
    else if (cdfEstimation == "logspline")                              
      node$cdf <- robustLogspline(vals, ...)
    else if (cdfEstimation == "kde") 
      node$cdf <- robustKde(vals, ...)
    
    return(node)
  }
  else {
    if (densityData == "topDown") {
      vals <- dat[node$instances, node$splitAttr]
      if (all(is.na(vals)))
        node$cdf <- NA
      else if (cdfEstimation == "ecdf" || is.factor(vals) ) 
        node$cdf <- ecdf(vals)
      else if (cdfEstimation == "logspline")                            
        node$cdf <- robustLogspline(vals, ...)
      else if (cdfEstimation == "kde") 
        node$cdf <- robustKde(vals, ...)
    }
    else if (densityData == "bottomUp") {
      
      leftVals <- dat[node$leftTree$instances, node$splitAttr]
      if (all(is.na(leftVals))) 
        node$cdfLeft <- NA
      else if (cdfEstimation == "ecdf" || is.factor(leftVals) ) 
        node$cdfLeft <- ecdf(leftVals)
      else if (cdfEstimation == "logspline")                              
        node$cdfLeft <- robustLogspline(leftVals, ...)
      else if (cdfEstimation == "kde") 
        node$cdfLeft <- robustKde(leftVals, ...)
      
      rightVals <- dat[node$rightTree$instances, node$splitAttr]
      if (all(is.na(rightVals))) 
        node$cdfRight <- NA
      else if (cdfEstimation == "ecdf" || is.factor(rightVals) ) 
        node$cdfRight <- ecdf(rightVals)
      else if (cdfEstimation == "logspline")                           
        node$cdfRight <- robustLogspline(rightVals, ...)
      else if (cdfEstimation == "kde") 
        node$cdfRight <- robustKde(rightVals, ...)  
    }  
    
    node$leftTree <- generatorFromTree(node$leftTree, dat, densityData, cdfEstimation, ...)
    node$rightTree <- generatorFromTree(node$rightTree, dat, densityData, cdfEstimation, ...)
    
    return(node)
  }
}


collectLeaves <- function(node, leaf) {
  if (node$nodeType=="leaf") {
    leaf[[length(leaf)+1]] <- node
  }
  else {
    leaf <- collectLeaves(node$leftTree, leaf)
    leaf <- collectLeaves(node$rightTree, leaf)
  }
  return(leaf)
  
}



tree2dot <- function(tree, attrNames, filename=NULL, name="dottree",digits=2) {
  if (!is.null(filename))
    sink(file=filename)
  
  cat(sprintf("digraph \"%s\" {\n", name))
  
  printTree2dot(tree, attrNames, 1, digits)
  
  cat(sprintf("}\n"))
  if (!is.null(filename))
    sink()
} 

printTree2dot <- function(node, attrNames, nodeIdx, digits){
  
  if (node$nodeType == "leaf") {
    if (! is.null(node$classProb)) {
      cat(sprintf("\tn%d [shape = box, label = \"noInst=%d, classProb=%s\"]\n", nodeIdx, node$noInst, probs2str(node$classProb,digits)))
    } else if (! is.null(node$meanResponse)) {       
      cat(sprintf("\tn%d [shape = box, label = \"noInst=%d, meanResponse=%.2f\"]\n", nodeIdx, node$noInst, node$meanResponse))	
    }
    else{
      cat(sprintf("\tn%d [shape = box, label = \"noInst=%d\"]\n", nodeIdx, node$noInst))
    }
    return(nodeIdx)    
  }
  else {
    
    if (node$splitType == "discrete"){
      if (any(is.na(node$leftValues))) 
        leftVals <- node$leftValues[-match(NA,node$leftValues)]
      else leftVals <- node$leftValues
      buf <- sprintf("%s = %s", attrNames[node$splitAttr], factors2str(leftVals))
    }
    else { 
      buf <- sprintf("%s <= %s", attrNames[node$splitAttr], format(node$splitValue,digits=digits))
    }
    
    cat(sprintf("\tn%d [label = \"%s\"]\n", nodeIdx, buf))
    
    cat(sprintf("\tn%d -> n%d [label = \"yes\"]\n", nodeIdx, nodeIdx+1))
    lastIdx <- printTree2dot(node$leftTree, attrNames, nodeIdx+1, digits)
    cat(sprintf("\tn%d -> n%d [label = \"no\"]\n", nodeIdx, lastIdx+1))
    lastUsed <- printTree2dot(node$rightTree, attrNames, lastIdx +1, digits)   
    return(lastUsed)
  }
}

