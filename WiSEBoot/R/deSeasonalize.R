deSeasonalize <-
function(dates, X, type="daily", method="deMean"){

  ##Check function arguments do not have incorrect length
  if(length(type)>1 || length(method)>1){
    stop("The type and method conditions are length 1.")
  }

  ##Check options have the correct value.
  if(!(type %in% c("daily","monthly","custom"))){
    stop("Invalid value for type.")
  }
  if(!(method %in% c("deMean","standardize"))){
    stop("Invalid value for method.")
  }

  ##Check X is a matrix or vector.
  if(is.matrix(X)==FALSE && ((length(X)==0) || is.atomic(X)!=TRUE ) ){
    stop("X should be a matrix or vector of positive dimension or length.")
  }
  ##Check X is numeric and non missing
  if(mode(X)!="numeric" || anyNA(X)==TRUE){
    stop("X should be numeric and contain no missing values.")
  }

  ##Check dates is non missing and a vector
  if(anyNA(dates)==TRUE){
    stop("dates should contain no missing values.")
  }
  if(is.atomic(dates)!=TRUE || length(dates)==0 || !is.null(dim(dates))){
    stop("dates should be a vector of positive length.")
  }


  ##Check dates are date format for any daily or monthly deseasonalizing
  if(type!="custom"){
    if(inherits(dates, "Date")==FALSE){
      stop("The dates option must be formatted as a Date within R to use daily or monthly deseasonalizing.")
    }
  }
  
  ##Check the length of dates matches one of the dimensions in X
  if(dim(X)[1]!=length(dates) && length(X)!=length(dates)){
    stop("X and dates should have the same number of observations in the time series.")
  }



  ##Set up IDs 
  if(type=="custom"){
    ID <- dates
  }else if(type=="monthly"){
    ID <- as.character(format(dates, format="%m"))
  }else{
    ID1 <- as.character(format(dates, format="%m"))
    ID2 <- as.character(format(dates, format="%d"))
    ID <- paste(ID1, ID2, sep='')
  }
  dataID <- data.frame(ID, Order=seq(1,length(ID)))

  ##Get a data frame of the means and standard deviations by ID
  meanID <- data.frame(aggregate(X~ID, FUN=mean))
  sdID <- data.frame(aggregate(X~ID, FUN=sd))

  ##Merge this information, re-order in the same order as original X matrix
  longMeanID <- merge(dataID, meanID, all.x=TRUE, sort=FALSE)
  longSdID <- merge(dataID, sdID, all.x=TRUE, sort=FALSE)
  longMeanID <- longMeanID[order(longMeanID$Order),]
  longSdID <- longSdID[order(longSdID$Order),]
  
  ##Remove the ID and Order columns from these new long aggregate data
  longMeanID <- longMeanID[ ,!(colnames(longMeanID) %in% c("ID","Order"))]
  longSdID <- longSdID[ ,!(colnames(longSdID) %in% c("ID","Order"))]

  ##Make long aggregate data as a vector/matrix
  if(is.matrix(X)==TRUE){
    longMeanID <- as.matrix(longMeanID)
    longSdID <- as.matrix(longSdID)
  }else{ 
    longMeanID <- as.vector(longMeanID)
    longSdID <- as.vector(longSdID)
  }

  ##Warn user if any of the ID have only 1 replicate.  De-seasonalized will just return a 0 there.
  if(anyNA(longSdID)){
    warning("One or more ID have only 1 replicate in the data.  
       The de-seasonalized data for that ID will be 0.")
    longSdID[is.na(longSdID)] <- 1 
  }

  ##Deseasonalize here
  if(method=="deMean"){
    newX <- X - longMeanID
  }else{ 
    newX <- (X - longMeanID) /longSdID 
  }
  return(newX)
}
