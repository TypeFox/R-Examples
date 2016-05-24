#########################################################################
# Categorical Network Class Methods
# Sample Categorization

.categorizeSample <- function(data, perturbations = NULL, object = NULL, nodeCats = NULL, ask=TRUE) {
  ## make sure data and object's nodes are the same
  
  if(is.matrix(data)) {
    numnodes <- nrow(data)
    numsamples <- ncol(data)

    ## save row/column names
    samplenames <- dimnames(data)
    ##if(is.numeric(data)) {
    ##  data <- matrix(as.integer(data), nrow=dim(data)[1])
    ##}
    
    maxcats <- 1
    lcats <- vector("list", numnodes)
    for(i in 1:numnodes) {

      ids <- !is.na(data[i,])
      if(sum(ids)<1)
        next
        
      if(is.numeric(data[i,ids])) {
        data[i,ids] <- as.integer(data[i,ids])
        cats <- min(data[i,ids]):max(data[i,ids])
      }
      else {
        cats <- levels(as.factor(data[i,ids]))
        ##cat(i," factors: ", cats, "\n")
      }
      if(!is.null(cats))
        cats <- sort(cats)
      else
        stop("Wrong cats for node ", i)
      
      if(!is.null(object)) {
        if(is.integer(cats)) {
          catcheck <- sapply(cats, function(cc) sum(c(1:length(object@cats[[i]]))==cc) > 0)
          cats <- 1:length(object@cats[[i]])
        }
        else {
          catcheck <- sapply(cats, function(cc) sum(object@cats[[i]]==cc)>0)
          cats <- object@cats[[i]]
        }
        if(prod(catcheck) != 1)
           stop("Incompatible object/data cats for node ", i)
      }
      else if(!is.null(nodeCats) && !is.null(nodeCats[[i]])) {
        catcheck <- sapply(cats, function(cc) length(which(nodeCats[[i]]==cc)))
        if(prod(catcheck) != 1)
           stop("Incompatible nodeCats/data cats for node ", i)
        cats <- nodeCats[[i]]
      }
      
      lencat <- length(cats)
      
      if(is.numeric(data[i,ids]) && is.null(object)) {
        data[i,ids] <- as.integer(data[i,ids]) - min(data[i,ids], cats) + 1
        data[i,data[i,ids]<1] <- 1
        data[i,data[i,ids]>lencat] <- lencat
      }
      else {
        data[i,ids] <- as.integer(sapply(data[i,ids], function(x) which(cats==x)))
      }
      
      lcats[[i]] <- cats
      
      if(maxcats < lencat)
        maxcats <- lencat
    }
    
    data <- matrix(as.integer(data), nrow=numnodes)
    dimnames(data) <- samplenames
    
    if(!is.null(perturbations) && !is.matrix(perturbations))
      stop("Perturbations should be a matrix")
    
  } ## is.matrix
  else {
    ## data is data.frame format
    numnodes <- ncol(data)
    numsamples <- nrow(data)
    
    fdata <- data
    data <- matrix(rep(NA, numnodes*numsamples), nrow=numnodes)
    
    maxcats <- 1
    lcats <- vector("list", numnodes)
    for(i in 1:numnodes) {

      ids <- !is.na(fdata[,i])
      if(sum(ids)<1)
        next
      
      if(is.numeric(fdata[ids,i])) {
        fdata[ids,i] <- as.integer(fdata[ids,i])
        cats <- min(fdata[ids,i]):max(fdata[ids,i])
      }
      else {
        cats <- levels(as.factor(fdata[ids,i]))
      }
      if(!is.null(cats))
        cats <- sort(cats)
      else
        stop("Wrong cats for node ", i)

      if(!is.null(object)) {
        if(is.integer(cats)) {
          catcheck <- sapply(cats, function(cc) sum(c(1:length(object@cats[[i]]))==cc) > 0)
          cats <- 1:length(object@cats[[i]])
        }
        else {
          catcheck <- sapply(cats, function(cc) sum(object@cats[[i]]==cc)>0)
          cats <- object@cats[[i]]
        }
        if(prod(catcheck) != 1)
           stop("Incompatible object/data cats for node ", i)
      }
      else if(!is.null(nodeCats) && !is.null(nodeCats[[i]])) {
        catcheck <- sapply(cats, function(cc) length(which(nodeCats[[i]]==cc)))
        if(prod(catcheck) != 1)
           stop("Incompatible nodeCats/data cats for node ", i)
        cats <- nodeCats[[i]]
      }
        
      lencat <- length(cats)

      if(is.numeric(fdata[ids,i]) && is.null(object)) {
        ## fdata is actually indices
        data[i,ids] <- as.integer(fdata[ids,i]) - min(fdata[ids,i], cats) + 1
        data[i, data[i,ids]<1] <- 1
        data[i, data[i,ids]>lencat] <- lencat
      }
      else {
        data[i,ids] <- as.integer(sapply(fdata[ids,i], function(x) which(cats==x)))
      }

      lcats[[i]] <- cats
      
      if(maxcats < lencat)
        maxcats <- lencat
    } ## i

    rownames(data) <- colnames(fdata)
    colnames(data) <- rownames(fdata)

    if(!is.null(perturbations)) {
      if(!is.data.frame(perturbations))
        stop("Perturbations should be a data frame")
      perturbations <- as.matrix(t(perturbations))
    }
  }

  if(length(rownames(data)) < numnodes)
    rownames(data) <- 1:numnodes
  
  if(ask && maxcats*maxcats > numsamples) {
    #cat("The sample size is too small. Continue? ('y' or 'n')\n")
    #if(scan("", what="character", nmax=1, quiet=TRUE) != "y" ) 
    #  stop("Operation canceled")
    warning("Small sample size")
  }

  if(ask && maxcats > 16) {
    cat("The data seems to have too many cats. The processing is expected to be very slow and memory consuming. Continue? ('y' or 'n')\n")
    if(scan("", what="character", nmax=1, quiet=TRUE) != "y" ) 
      stop("Operation canceled")
  }

  if(!is.null(perturbations)) {
    dims <- dim(data)
    pertdims <- dim(perturbations)
    if(dims[1] != pertdims[1] ||
       dims[2] != pertdims[2])
      stop("Incompatible perturbation dimensions.\n")
  }
  
  return(list(data=data, cats=lcats, maxcats=as.integer(maxcats), perturbations=perturbations))
}

