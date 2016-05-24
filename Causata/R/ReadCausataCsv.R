
require(doMC)

# reads csv data from causata
ReadCausataCsv <- function(causataR, include=c(), exclude=c(), 
  maxMb=1000, colFilterFunc=NA, rowIndex=NA, nrows=NA, 
  metadata=FALSE, debug=FALSE, ...) {
  # keep: a list of patterns to match against column names, matches are kept
  # exclude: a list of patterns to match against column names, matches are excluded
  tStart <- proc.time()
  
  # the number of rows to read is pulled from the causataR object unless
  # it's provided as a parameter to ReadCausataR
  if (is.na(nrows)){
    # nrows was not passed in, use the value from the causataR object
    nrows <- causataR$nrows
  }
   
  # ensure that user is not running R from an interactive app, which may have problems with doMC and foreach
  if (interactive() & getDoParWorkers()>1){stop("ReadCausataCsv may not be called from an interactive session when\nthe doMC package is used.  Run from a terminal instead.  Stopping execution.\n")}
    
  # initialize list of columns to load
  cat("Beginning data load:\n")
  cat("  ", causataR$fileData, "\n", sep="")
  colsToLoad <- causataR$col.names
  if ((causataR$nrows==-1) & is.na(nrows)) {
    nrowsMessage <- "unknown"
  } else if (is.na(nrows)) {
    nrowsMessage <- causataR$nrows
  } else {
    nrowsMessage <- nrows
  }
  cat("  Data set contains", nrowsMessage, "rows and", length(colsToLoad), "columns\n")
  # test if an exclude list was provided
  if (length(exclude)>0){
    # exclude unwanted variables
    colsToLoad <- colsToLoad[-GrepLoop(exclude, colsToLoad)]
  }
  # test if a "keep" list was provided
  if (length(include)>0){
    colsToLoad <- colsToLoad[GrepLoop(include, colsToLoad)]
  }
  
  # estimate the memory required and number of passes
  bytesPerCell <- 12 # estimated memory required per cell of a data frame in bytes
  if (nrows > 0){
    estimatedMb <- nrows * 12 * length(colsToLoad) / 1e6
    numPasses <- ceiling(estimatedMb / maxMb)
  } else {
    # rows weren't provided, default to 1 pass
    numPasses <- 1
    estimatedMb <- NA
  }
  cat("  Loading", length(colsToLoad), "columns in", numPasses,"passes.  Estimating", format(estimatedMb, digits=5), "megabytes before filtering.\n")
  
  #
  # Loop over passes starts here
  #
  i1 <- 1 # initialize for first iteration
  # columns per pass
  colsPerPass <- ceiling(length(colsToLoad) / numPasses)
  # loop
  for (NP in 1:numPasses){
    # set index for columns
    i2 <- min(i1 + colsPerPass - 1, length(colsToLoad))
    # set columns to load on this pass
    colsToLoadPass <- colsToLoad[i1:i2]
  
    # load data from csv
    idx <- causataR$col.names %in% colsToLoadPass # index of columns to load before pass filtering
    colClasses <- causataR$colClasses
    colClasses[!idx] <- "NULL" # set class to NULL for columns we don't want to load
    dfp <- read.csv(causataR$fileData, sep=",", comment.char="", na.strings=c("MV"), skip=1, 
      colClasses=colClasses, col.names=causataR$col.names, nrows=nrows)
  
    cat("  Pass ", NP, "/", numPasses, " load complete, ", nrow(dfp), " rows and ", ncol(dfp) , " columns.  Columns ", i1, " to ", i2, "\n", sep="")
    
    # check if an index of rows was provided on first pass, if not then set index to all rows
    if (NP==1){
      if (length(rowIndex)==1){ 
        rowIndex <- rep(TRUE,nrow(dfp)) 
      } else {
        # require that the input index is a logical
        stopifnot(class(rowIndex)=="logical")
      }
    }
    # filter out unwanted rows, catch the case where all rows are kept
    if (sum(rowIndex) < nrow(dfp)) {
      dfp <- dfp[rowIndex,]
    }
    
    # apply the column filter function to keep or discard columns
    colFilterOutput <- ApplyColumnFilterReadCausataCsv(dfp, colFilterFunc, debug, ...)
    # copy an index of TRUE/FALSE values indicating which columns to keep
    colFilterIdx <- colFilterOutput$colFilterIdx
    # remove unwanted columns, catch the case where all columns are kept
    if (sum(colFilterIdx) < ncol(dfp)){
      dfp <- dfp[,colFilterIdx]
    }
    
    # if this is the first pass then copy data frame
    if (NP == 1){
      # first pass, copy the data frame from the first iteration
      df <- dfp
      metadataList <- colFilterOutput$foreachOutList
    } else {
      # after first pass, join new data frane to existing data frame
      df <- cbind(df, dfp)
      metadataList <- c(metadataList, colFilterOutput$foreachOutList)
    }
    cat("  Returning", sum(rowIndex), "rows and", sum(colFilterIdx), "columns.\n\n")
    # update index for next iteration
    i1 <- i2 + 1
    # if this is the end of the columns then exit the loop
    if (i1 > length(colsToLoad)){break}
  }
  #
  # Loop over passes end
  #
  
  # write processing time
  cat("  Passes complete, returning", sum(rowIndex), "rows and", ncol(df), "columns.\n")
  cat(sprintf("  Elapsed time (minutes): %6.2f \n", (proc.time()-tStart)[3]/60 ))
  
  # return output data
  if (metadata){
    # assign class of the metadata list
    class(metadataList) <- class(colFilterOutput$foreachOutList)
    # save to a list
    outlist <- list(df=df, metadata=metadataList)
    # initializes list names variable
    vnames <- c()
    # extract variable names from list elements
    for (i in 1:length(outlist$metadata)){
      vnames[i] <- outlist$metadata[[i]]$name
    }
    # apply names to list
    names(outlist$metadata) <- vnames
    return(outlist)
  } else {
    return(df)
  }
}


ApplyColumnFilterReadCausataCsv <- function(df, colFilterFunc, debug, ...){
  # check if a column filter function was provided
  if (class(colFilterFunc) != "function"){
    # no function provided, pass all columns with a TRUE value
    colFilterIdx <- rep(TRUE, ncol(df))
    outList <- list(colFilterIdx = colFilterIdx)
    return(outList)
  } else {
    # a filter function was provided, check if it returns a list or a logical by applying to the first column
    cat("  Applying column filter...\n")
    colFilterTestOutput <- colFilterFunc(df[,1], ...)
    if (class(colFilterTestOutput) != "logical"){
      # function returns non-logical, may be list or object, ensure it has an element called "keep"
      stopifnot("keep" %in% names(colFilterTestOutput))
      # set flag to false, it's a list not a logical
      colFilterLogical <- FALSE 
    } else if (class(colFilterTestOutput) == "logical"){
      # returned a logical
      colFilterLogical <- TRUE
    } else {
      # didn't return list or logical so it's an error
      stop("Error, return value must be list or logical")
    }
  }
  
  # the user provided a column function, and we know if it returns a list or logical,
  # now run the function on all columns
  if (debug){
    # debug mode, iterate manually
    foreachOutList <- list()
    for (i in 1:ncol(df)){
      cat("  ", names(df)[i], "\n")
      foreachOutList[[i]] <- colFilterFunc(df[,i], ...)
    }
  } else {  
    # not debug mode, use foreach (and maybe parellelization)
    foreachOutList <- foreach( i=1:ncol(df) ) %dopar% colFilterFunc(df[,i], ...)
  }
  
  # process the output from the filtering function
  colFilterIdx <- rep(NA, length(foreachOutList) )
  if (colFilterLogical){
    # the filter function returns a logical, extract it
    for (i in 1:length(foreachOutList)){
      colFilterIdx[i] <- foreachOutList[[i]]
    }
  } else {
    # the filter function returns a list, look in the 'keep' element
    for (i in 1:length(foreachOutList)){
      colFilterIdx[i] <- foreachOutList[[i]]$keep
    }
  }
  
  # apply column names to the output list
  varnames <- colnames(df)
  for (i in 1:length(foreachOutList)){
    foreachOutList[[i]]$name <- varnames[i]
  }
  
  # set the class of the output
  class(foreachOutList) <- foreachOutList[[1]]$containerClass
  
  # return the output in a list
  outList <- list(colFilterIdx=colFilterIdx, foreachOutList=foreachOutList)
}


GrepLoop <- function(patternVec, x, ignore.case=TRUE, boolean=FALSE){
  # given a vector of patterns and a vector of strings x, this returns the indices of strings in x that match
  # one or more patterns
  idx <- c()
  for (pat in patternVec) {
    # index of which columns match this variable
    idx <- c( idx, grep(pat, x, ignore.case=ignore.case) )
  }
  if (boolean) {
    # return boolean index
    b.idx <- rep(FALSE, length(x))
    if (length(idx)>0){
      # found matches
      b.idx[sort(unique(idx))] <- TRUE
    }
    return(b.idx)
  } else {
    # return numeric index
    return(sort(unique(idx))) 
  }
}

