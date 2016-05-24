tree <-function(formula, dataset, noTrees = 1, minNodeWeight=5, problemType="byResponse") {
  treeEnsemble(formula, dataset, noTrees=noTrees, noSelectedAttr=ncol(dataset)-1, 
               minNodeWeight=minNodeWeight, problemType=problemType, densityData="no") 
}
bagging <-function(formula, dataset, noTrees = 50, minNodeWeight=2, problemType="byResponse") {
  treeEnsemble(formula, dataset, noTrees=noTrees, noSelectedAttr=ncol(dataset)-1, 
               minNodeWeight=minNodeWeight, problemType=problemType, densityData="no") 
}
rf <-function(formula, dataset, noTrees = 50, minNodeWeight=2, noSelectedAttr=max(1,integer(sqrt(ncol(dataset)-1))),
              problemType="byResponse",densityData="no") {
  treeEnsemble(formula, dataset, noTrees=noTrees, noSelectedAttr=noSelectedAttr, 
               minNodeWeight=minNodeWeight, problemType=problemType) 
}
rfDensity <-function(formula, dataset, noTrees = 50, minNodeWeight=2, noSelectedAttr=max(1,integer(sqrt(ncol(dataset)-1))), 
                     problemType="byResponse", densityData="leaf", cdfEstimation="ecdf", ...) {
  treeEnsemble(formula, dataset, noTrees=noTrees, noSelectedAttr=noSelectedAttr, 
               minNodeWeight=minNodeWeight, problemType=problemType, densityData, cdfEstimation, ...)
}
densityEnsemble <-function(formula, dataset, noTrees = 100, minNodeWeight=2, noSelectedAttr=2,
                           problemType="density", densitySplitMethod="balancedSplit", densityData="leaf",
                           cdfEstimation = "ecdf", ...) {
  treeEnsemble(formula, dataset, noTrees=noTrees, minNodeWeight=minNodeWeight,
               noSelectedAttr=noSelectedAttr, problemType=problemType, 
               densitySplitMethod=densitySplitMethod, densityData=densityData, 
               cdfEstimation=cdfEstimation, ...)				   
}
indAttrGen <-function(formula, dataset, cdfEstimation = c("ecdf","logspline","kde"), problemType="byResponse") {
  treeEnsemble(formula, dataset, cdfEstimation=cdfEstimation, noTrees=1, minNodeWeight=nrow(dataset), 
               problemType=problemType, densityData="leaf") 
}

treeEnsemble <-function(formula, dataset, noTrees = 100, minNodeWeight=2, noSelectedAttr=2, 
                        problemType=c("byResponse","classification", "regression","density"),
                        densityData=c("leaf", "topDown", "bottomUp","no"),
                        cdfEstimation = c("ecdf","logspline","kde"), 
                        densitySplitMethod=c("balancedSplit","randomSplit","maxVariance"),					
                        estimator=NULL, ...) 
{
  if (!inherits(formula,"formula")) 
    stop("First argument must be a formula.");
  problemType <- match.arg(problemType)	
  densitySplitMethod <- match.arg(densitySplitMethod)
  densityData <- match.arg(densityData)
  cdfEstimation <- match.arg(cdfEstimation)
  #prepare data according to the supplied formula
  dat <- model.frame(formula, data=dataset, na.action=na.pass);
  terms <- attr(dat,"terms")
  # check and make target variable coherent with problemType
  if (attr(terms,"response")==0) {
    if (problemType %in% c("byResponse", "density"))
      problemType <- "density"
    else
      stop("Formula shall provide target variable in classification and regression problems")
  }
  else {
    if (problemType == "density")
      stop("Formula shall not define target variable for density ensembles.")
  }
  if (problemType == "classification" && ! is(dat[[1]],"factor")) {
    dat[[1]] <- factor(dat[[1]]);
    cat("Changing dependent variable to factor with levels:",levels(dat[[1]]),"\n");
  }
  if (problemType == "regression" && ! is(dat[[1]],"numeric"))
    stop("Prediction variable type shall be numeric for regression problems.")
  if (problemType == "byResponse") {
    if (is(dat[[1]],"factor"))
      problemType <- "classification"
    else if (is(dat[[1]],"numeric"))
      problemType <- "regression"
    else stop("Prediction variable type shall be either factor for classification or numeric for regression problems.")
  }
  
  # set and check estimator
  if ( is.null(estimator) ) {
    if (problemType == "classification")
      estimator <- "Gini"
    else estimator <-"MSEofMean"
  }
  else {
    if (problemType == "classification"){
      if (! estimator %in% infoCore(what="attrEval"))
        stop("Estimator ", estimator, " is not a valid option for classification problem.")
    } else if (problemType == "regression"){
      if (! estimator %in% infoCore(what="attrEvalReg"))
        stop("Estimator ", estimator, " is not a valid option for regression problem.")
    }
  }
  
  if (problemType != "density") {
    dat <- dat[,c(2:ncol(dat),1)] # reorder, class to the end
    classIdx <- ncol(dat)
  }
  if (problemType == "classification") {
    class.lev <- levels(dat[[classIdx]]);
    noClasses <- length(class.lev);
  } 
  else {
    class.lev <- NULL
    noClasses <- 0
  }
  noInst <- nrow(dat)
  noAttr = ncol(dat)
  
  if (densityData != "no") { # collect additional data needed for generation of artificial data
    attrClasses<-list()
    attrLevels <-list()
    attrOrdered<-logical()
    for (i in 1:noAttr) {
      attrClasses[[i]] <- class(dat[[i]])
      attrOrdered[i] <- is.ordered(dat[[i]])
      if (is.factor(dat[[i]])){
        attrLevels[[i]] <- levels(dat[[i]])
      }
      else attrLevels[[i]] <- NULL    
    }
  }
  
  #create trees
  splitVector = vector("logical",length=noInst)
  trees<-list()
  for (t in 1:noTrees) {
    splitVector[] <-TRUE
    if (noTrees==1) { #ordinary decision or regression tree
      ib <- 1:noInst
      oob <-c()
    }
    else {
      ib <- sample(1:noInst, size = noInst, replace = TRUE)
      splitVector[ib] <- FALSE
      oob <- (1:noInst)[splitVector]
    }
    if (problemType == "density")
      rt <- densityTree(dat[ib,], minNodeWeight,  noSelectedAttr, densitySplitMethod, densityData)
    else 
      rt <- randomTree(dat[ib,], minNodeWeight, noSelectedAttr=noSelectedAttr, problemType, densityData, estimator)
    
    trees[[t]] <- list(inBag = ib, outOfBag = oob, tree = rt)
  }
  
  treeEnsemble <- list(formula=formula, terms = terms, class.lev=class.lev, noClasses = noClasses,                        
                       noTrees = noTrees, trees = trees, problemType=problemType, densityData = densityData, 
                       cdfEstimation = cdfEstimation)
  if (densityData != "no") {
    treeEnsemble$noAttr = noAttr
    treeEnsemble$attrClasses = attrClasses
    treeEnsemble$attrLevels = attrLevels
    treeEnsemble$attrNames = names(dat)
    treeEnsemble$originalNames <- names(dataset)[sort(match(names(dat), names(dataset)))]
    treeEnsemble$attrOrdered = attrOrdered
    
    treeEnsemble <- treeDataGenerator(treeEnsemble, dat, ...)   
  }
  
  class(treeEnsemble) <- "TreeEnsemble"
  return(treeEnsemble)
}

predict.TreeEnsemble <- function(object, newdata, type=c("both","class","probability")) {
  type <- match.arg(type)
  if (object$problemType=="classification") {
    class.lev <- object$class.lev;
    noClasses <- length(class.lev);
  }
  noInst = nrow(newdata)
  terms <- delete.response(object$terms);
  newdat <- as.data.frame(newdata)
  dat <- model.frame(terms, data=newdat, na.action=na.pass);
  
  if (object$problemType=="classification") {
    ePred <- factor(c(), levels = class.lev)
    eProb <- matrix(NA, nrow = noInst, ncol=noClasses)
    prob <- matrix(NA, nrow = object$noTrees, ncol=noClasses)
    for (i in 1:noInst) {
      pred <- factor(c(), levels = class.lev)
      prob[,] <- NA
      for (t in 1:object$noTrees) {
        prediction <- predictWithTree(object$trees[[t]]$tree, dat[i,])
        pred[t] <- prediction$majorityClass
        prob[t,] <- prediction$classProb
      }
      #ensemble predictions
      ePred[i] <- factor(class.lev[which.is.max(table(pred))], levels=class.lev)
      eProb[i,] <- apply(prob, 2, sum) / object$noTrees
    }	  
    if (type == "both")
      returnList <- list(class=ePred,probabilities=eProb)
    else if (type=="class")
      returnList <- ePred
    else if (type == "probability")
      returnList <- eProb  
  }
  else {
    ePred <- vector(mode="numeric", length = noInst)  
    pred <- vector(mode="numeric", length = object$noTrees)  
    for (i in 1:noInst) {
      for (t in 1:object$noTrees) {
        pred[t] <- predictWithTree(object$trees[[t]]$tree, dat[i,])
      }
      #ensemble predictions
      ePred[i] <- sum(pred) / object$noTrees
    }	  
    returnList <- ePred
  }
  returnList
}

# newdata <- function(object, ...) UseMethod("newdata", object)

newdata.TreeEnsemble <- function(object, size=1, ...) {
  if (object$densityData == "no")
    stop("newdata needs proper information stored in the model to generate data")
  if (! is(object, "TreeEnsemble"))
    stop("newdata cannot generate data for objects of class ",class(object))
  
  noInst = size
  dat <- matrix(NA, nrow = size, ncol = object$noAttr )
  
  if (object$densityData == "leaf") {
    treeOrder <- sample(1:object$noTrees, noInst, replace=T)
    for (i in 1:noInst) {
      t <- treeOrder[i]
      leafIdx <- sample(1:length(object$dataGenerator[[t]]$leaves), size = 1, prob = object$dataGenerator[[t]]$weights, replace=T)
      for (a in 1:object$noAttr) {
        if (! any(class(object$dataGenerator[[t]]$leaves[[leafIdx]]$cdfs[[a]]) %in% c("ecdf","logspline","kde"))) {
          dat[i,a] <- NA
        }
        else if ("factor" %in% object$attrClasses[[a]]) {
          dat[i, a] <- quantile(object$dataGenerator[[t]]$leaves[[leafIdx]]$cdfs[[a]], probs=runif(1, 0, 1), type=3)
        }
        else if ( inherits(object$dataGenerator[[t]]$leaves[[leafIdx]]$cdfs[[a]], "ecdf") ) {
          dat[i, a] <- quantile(object$dataGenerator[[t]]$leaves[[leafIdx]]$cdfs[[a]], probs=runif(1, 0, 1), type=8, ...)
        }
        else if ( inherits(object$dataGenerator[[t]]$leaves[[leafIdx]]$cdfs[[a]], "logspline") ){
          dat[i,a] <- rlogspline(1, object$dataGenerator[[t]]$leaves[[leafIdx]]$cdfs[[a]]) 
        }
        else if ( inherits(object$dataGenerator[[t]]$leaves[[leafIdx]]$cdfs[[a]], "kde") ) {
          dat[i, a] <- rkde(1, object$dataGenerator[[t]]$leaves[[leafIdx]]$cdfs[[a]])
        }
        else stop("Invalid type of generator detected in the leaf ",class(object$dataGenerator[[t]]$leaves[[leafIdx]]$cdfs[[a]]))
      }
      if (any(is.na(dat[i,]))) {
        for (j in 1:ncol(dat))
          if (is.na(dat[i,j]))
            dat[i,j] <- useDefaultGenerator(object$defaultGenerator[[j]])
      }	  
    }
  }
  else {
    for (i in 1:noInst) {
      treeOrder <- sample(1:object$noTrees, object$noTrees)
      t <- 1
      # while not all the values are filled in
      while (any(is.na(dat[i,])) && t <= object$noTrees) {
        if (object$densityData == "topDown") {
          gen <- fillDataWithTreeTopDown(object$trees[[treeOrder[t]]]$tree, dat[i,])
        }
        else if (object$densityData == "bottomUp") {
          gen <- fillDataWithTreeBottomUp(object$trees[[treeOrder[t]]]$tree, dat[i,])     
        }
        else stop("newdata encountered unrecognized densityData type in the model ", object$densityData)
        dat[i,] <- gen 
        t <- t + 1
      }
      if (any(is.na(dat[i,]))) {
        for (j in 1:ncol(dat))
          if (is.na(dat[i,j]))
            dat[i,j] <- useDefaultGenerator(object$defaultGenerator[[j]])
      }	  
    }
  }
  newdata <- as.data.frame(dat)
  
  names(newdata) <- object$attrNames
  for (i in 1:object$noAttr){
    if ("factor" %in% object$attrClasses[[i]]) {
      newdata[[i]] <- factor(newdata[[i]],levels=1:length(object$attrLevels[[i]]), labels=object$attrLevels[[i]])
      if (object$attrOrdered[i])
        newdata[[i]] <- as.ordered(newdata[[i]])
    }
  }
  newdata[,object$originalNames]
}

treeDataGenerator<-function(ensemble, dat, ...) {
  if (ensemble$densityData == "leaf") { 
    generator <- list()
    for (t in 1:ensemble$noTrees) {
      generator[[t]] <- list()
      tree <- fillWithInstances(ensemble$trees[[t]]$tree, dat[ensemble$trees[[t]]$inBag,], ensemble$trees[[t]]$inBag)
      generator[[t]]$leaves <- list()
      generator[[t]]$weights <- c()
      generator[[t]] <- genFromLeaves(tree, dat, ensemble$cdfEstimation, ...)
    }
    ensemble$dataGenerator <- generator
  }
  else {
    for (t in 1:ensemble$noTrees) {
      
      ensemble$trees[[t]]$tree <- fillWithInstances(ensemble$trees[[t]]$tree, dat[ensemble$trees[[t]]$inBag,], ensemble$trees[[t]]$inBag)
      
      ensemble$trees[[t]]$tree <- generatorFromTree(ensemble$trees[[t]]$tree, dat, ensemble$densityData, ensemble$cdfEstimation, ...)
    }  
  }
  ensemble$defaultGenerator <- defaultGenerator(dat, ensemble$cdfEstimation, ... )
  return(ensemble)
}

robustLogspline<-function(x, ...) {
  tc <- tryCatch(model<-logspline(x, ...),error=function(e) e, warning=function(w) w)
  if (is(tc,"warning") || is(tc,"error")) {
    model <- ecdf(x)
    warning(tc, "converting to ecdf")
  }
  model
}

robustKde<-function(x, ...) {
  tc <- tryCatch(model<-kde(x, ...),error=function(e) e, warning=function(w) w)
  if (is(tc,"warning") || is(tc,"error")) {
    model <- ecdf(x)
    warning(tc, "converting to ecdf")
  }
  model
}

defaultGenerator <- function(data, cdfEstimation="ecdf", ...) {
  dgen <- list()
  for (a in 1:ncol(data)) {
    vals <- data[, a]
    if (all(is.na(vals)))
      dgen[[a]] <- NA
    else if ( cdfEstimation == "ecdf" || is.factor(data[[a]]) )
      dgen[[a]] <- ecdf(vals)
    else if (cdfEstimation == "logspline")                             
      dgen[[a]] <- robustLogspline(vals, ...)
    else if (cdfEstimation == "kde")
      dgen[[a]] <- robustKde(vals, ...)
  }
  names(dgen) <- names(data)
  return(dgen)
}

useDefaultGenerator<-function(dataGenerator) {
  if (! any(class(dataGenerator) %in% c("ecdf","logspline","kde")))
    return(NA)
  else if ( is(dataGenerator, "ecdf") )
    return (quantile(dataGenerator, probs=runif(1, 0, 1), type=8))
  else if ( is(dataGenerator, "logspline") )
    return(rlogspline(1, dataGenerator)) 
  else if ( is(dataGenerator, "kde") )
    return(rkde(1, dataGenerator))
  else stop("Default data generator is of unknown type.")
}