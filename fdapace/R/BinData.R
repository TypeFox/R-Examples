BinData = function(y,t,optns){
  
  # Bin the data 'y'
  # y : n-by-1 list of vectors
  # t : n-by-1 list of vectors 
  # dataType : indicator about structure of the data 
  #   (dense (2), or  dataType data with missing values (1) or sparse (0))
  # verbose gives warning messages
  # numBins: number of bins (if set)

  BinDataOutput <- list( newy <- NULL, newt <- NULL);

  dataType = optns$dataType;
  verbose = optns$verbose;  

  n = length(t);
  ni = sapply(FUN= length,t);
  
  if (dataType == 'Sparse'){
    m = median(ni)
  } else {
    m = max(ni);
  }
  
  # Check the number of bins automatically
  if (is.null(numBins)){
    numBins = GetBinNum(n,m,dataType,verbose)
  }else if( numBins <1){
    warning("number of bins must be positive integer! We reset to the default number of bins!\n")
    numBins = GetBinNum(n,m,dataType,verbose)
  }
  
  # If it is determined to be NULL return the unbinned data
  if (is.null(numBins)){
    BinDataOutput$newt = t;
    BinDataOutput$newy = y;
    return( BinDataOutput )
  }
     
  numBins = ceiling(numBins);
  
  tt = unlist(t);  
  a0 = min(tt);
  b0 = max(tt);
  
  for (i in 1:n){
    res = GetBinnedCurve(t[[i]], y[[i]], numBins, TRUE, TRUE, c(a0, b0));
    BinDataOutput$newt[i] = res$midpoint;   
    BinDataOutput$newy[i] = res$newy;      
  }
  
  return( BinDataOutput )
}
