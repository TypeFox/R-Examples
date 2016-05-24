germ <- function(newData, refData, 
                 parameters = list("Max forward error" = 25,
                                   "Max backward error" = 25,
                                   "Max sum error" = 100,
                                   "Lower measurement limit" = 100), 
                 method = "joint", na.rm = TRUE){
  if(!is.na(pmatch(method, "joint"))) method <- "joint"
  METHODS <- c("joint", "forward", "backward", "sum")
  method <- pmatch(method, METHODS)
  if(is.na(method)) stop("invalid ranking method")
  if (method == -1) stop("ambiguous ranking method")
  
  newNames <- unique(newData$Sample)
  refNames <- unique(refData$Sample)
  
  res <- array(NA, dim = c(length(refNames), 4, length(newNames)))
  rownames(res) <- refNames
  colnames(res) <- c("Forward Max", "Backward Max", "Sum of Bands", "Joint")
  dimnames(res)[[3]] <- newNames
  
  if(na.rm){
    newData <- newData[!is.na(newData$MW),]
    refData <- refData[!is.na(refData$MW),]
  }
  
  for(i in 1:length(newNames)){
    for(j in 1:length(refNames)){
      New <- newData[newData$Sample == newNames[i],]
      Ref <- refData[refData$Sample == refNames[j],]
      nrEnz <- length(unique(Ref$Enzyme))
      res.enz.sum <- res.enz.bw <- res.enz.fw <- numeric(nrEnz)
      refenzNames <- unique(Ref$Enzyme)
      newenzNames <- unique(New$Enzyme)
      
      if(length(refenzNames) != length(newenzNames))
        stop("Number of enzymes is different for new data and reference data!")
      if(!all(refenzNames == newenzNames))
        stop("Names of enzymes are different for new data and reference data!")
      
      for(k in 1:nrEnz){
        Newk <- New[New$Enzyme == newenzNames[k],]
        Refk <- Ref[Ref$Enzyme == newenzNames[k],]
        ## forward
        ind.fw <- Newk$MW > (parameters$"Max forward error" +
                               parameters$"Lower measurement limit")
        res.enz.fw[k] <- .forward(abs(Newk$MW[ind.fw]), abs(Refk$MW))
        ## backward
        ind.bw <- Refk$MW > (parameters$"Max backward error" + 
                             parameters$"Lower measurement limit")
        res.enz.bw[k] <- .backward(abs(Newk$MW), abs(Refk$MW[ind.bw]))
        ## sum
        res.enz.sum[k] <- abs(sum(abs(Newk$MW)) - sum(abs(Refk$MW)))
      }
#      cat(Ref$Sample[1], ":\t", res.enz.sum, "\n")
      res[j,1,i] <- max(res.enz.fw, na.rm = TRUE)
      res[j,2,i] <- max(res.enz.bw, na.rm = TRUE)
      res[j,3,i] <- max(res.enz.sum, na.rm = TRUE)
      res[j,4,i] <- res[j,1,i]+res[j,2,i]
    }
  }
  res.ord <- vector("list", length = length(newNames))
  names(res.ord) <- newNames
  if(method == 1){
    for(i in 1:length(newNames)){
      res.ord[[i]] <- res[order(res[,4,i]),,i]
    }
  }
  if(method == 2){
    for(i in 1:length(newNames)){
      res.ord[[i]] <- res[order(res[,1,i]),,i]
    }
  }
  if(method == 3){
    for(i in 1:length(newNames)){
      res.ord[[i]] <- res[order(res[,2,i]),,i]
    }
  }
  if(method == 4){
    for(i in 1:length(newNames)){
      res.ord[[i]] <- res[order(res[,3,i]),,i]
    }
  }
  res.ord
}
.backward <- function(x, y){
  res <- numeric(length(y))
  ## zero padding of x
  if(length(x) < length(y))
    x <- c(x, rep(0, length(y)-length(x)))
  for(i in 1:length(y)){
    res[i] <- min(abs(y[i]-x))
  }
  max(res)
}
.forward <- function(x, y){ .backward(x = y, y = x) }
