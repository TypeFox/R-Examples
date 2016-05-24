GetBinnedDataset <- function (y, t, optns){
   
  # Bin the data 'y'
  # y : n-by-1 list of vectors
  # t : n-by-1 list of vectors 

  BinDataOutput <- list( newy=NULL, newt=NULL);

  dataType = optns$dataType;
  verbose = optns$verbose;  
  numBins = optns$numBins;
  tt = unlist(t);  
  a0 = min(tt);
  b0 = max(tt);

  n = length(t);
  ni = sapply(FUN= length,t);
  
  if (dataType == 'Sparse'){
    m = median(ni)
  } else {
    m = max(ni);
  }
  
  # Determine the number of bins automatically if numBins is null
  if (is.null(numBins) && optns$useBinnedData =='AUTO'){
    numBins = GetBinNum(n,m,dataType,verbose)
    # and if it is still NULL return the unbinned data
    if (is.null(numBins)){
      BinDataOutput$newt = t;
      BinDataOutput$newy = y;
      return( BinDataOutput )
    }
  }
  # otherwise use the one provided by the user (ceiled)
  numBins = ceiling(numBins);

  resList <- lapply(1:n, function(i) 
    GetBinnedCurve(t[[i]], y[[i]], numBins, TRUE, TRUE, c(a0, b0)))
  BinDataOutput[['newt']] <- lapply(resList, `[[`, 'midpoint')
  BinDataOutput[['newy']] <- lapply(resList, `[[`, 'newy')
  
  # for (i in 1:n){
    # res = GetBinnedCurve(t[[i]], y[[i]], numBins, TRUE, TRUE, c(a0, b0));
    # BinDataOutput$newt[[i]] = res$midpoint;   
    # BinDataOutput$newy[[i]] = res$newy;      
  # }
     
  result <- list( 'newt' = BinDataOutput$newt, 'newy' = BinDataOutput$newy);
  # Garbage Collection
  gc()
  return(result)
}


