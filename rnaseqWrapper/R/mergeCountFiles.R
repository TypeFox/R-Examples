## Function to merge read/fpkm files

mergeCountFiles <-
  function(fileDir, # name of the directory to read in
           fileID = "*.genes.results$", # a character to use to limit which files are imported
           #                              regular expressions allowed
           #                              this is removed from the file names to name samples
           fileSep = "\t", # The column delimiter used in the file (e.g. "," or "\t")
           seqIDcol = 1, # Which column has the reference position?
           #                                    can be numeric or character
           colsToKeep = c("expected_count","FPKM"), ## Which columns of info should be kept?
           idCols = NULL, # Which columns of additional information should be KEPT?
           #           # note: row.names will match refPosCol.
           #           # note: format (i.e. numeric vs. character) must match seqIDcol IF merging
           minMatchToMerge = 0.5 # What portion of gene names must be in common to proceed with merge?
           #                        This is a protection step to avoid accidentally merging datasets
           #                          from different references.
           ){

    ## Returns a data.frame with the mergedCounts

    
    ## Make sure that fileDir has trailing slash
    if(length(grep("/$",fileDir))==0){
      fileDir <- paste(fileDir,"/",sep="")
    }
    
    ## Get Files
    files <- list.files(fileDir,pattern=fileID)
    sampName <- gsub(fileID,"",files)
    indFiles <- list()
    for(k in 1:length(files)){
      indFiles[[sampName[k]]] <- read.table(paste(fileDir,files[k],sep=""), header = TRUE, sep = fileSep)
    }
    
    ## Confirm that rows are in the same order:
    namesTemp <- lapply(indFiles,function(x) {as.character(x[,seqIDcol])})
    names <- do.call(cbind,namesTemp)
    
    testName <- apply(names,1,function(x) {sum(x == x[1]) })
#     table(testName)
    nNotMatchAll <- sum(testName != length(indFiles))
    
    if( nNotMatchAll > 0){
      # If not, Check that the names match, just not in order:
      checkingSet <- levels(as.factor(namesTemp[[1]]))
      namesMatched <- lapply(namesTemp[-1],function(x) {sum(x %in% checkingSet)/length(x)})
      
      if(min(unlist(namesMatched)) < minMatchToMerge){
        stop("Fewer names matched across some samples than your cut off. ",
             "There appears to be an error in your input data.")
      } else {
        warning("Gene names are not in the same order in all files. ",
                "Confirm that they were produced against the same reference",
                "before proceeding.",
                "\n\nAnalyis proceeding by merging files.")
        minPropMatched <- min(unlist(namesMatched))
      }
      
      keepOrder <- FALSE
      
    } else{ ## To tell it not to bother merging
      keepOrder <- TRUE
    }
    
    ## Limit down to just columns to keep 
    colsInt <- lapply(indFiles,function(x) { 
      out <- x[,colsToKeep]
      ## Rescue when only one col kept
      if(length(colsToKeep == 1)){
        out <- as.data.frame(out)
        if(class(colsToKeep) == "numeric"){
          names(out) <- names(x)[colsToKeep]
        } else {
          names(out) <- colsToKeep
        }
      }
      row.names(out) <- x[,seqIDcol]
      return(out)
      })
    # it appears that it is unnecessary to confirm that the columns exist,
    #  as the system with throw an "undefined columns selected" error if not
    
    colNames <- do.call(rbind,lapply(colsInt,names))
    if(any(apply(colNames,2,function(x) {sum(x == x[1]) != length(x) }) ) ){
      warning("Not all column names match -- confirm that all file types are the same.")
    }
        
    ## keep separate set with additional info
    if( !is.null(idCols) ){
      if(keepOrder){
        addInfo <- indFiles[[1]][,idCols]
        row.names(addInfo) <- indFiles[[1]][,seqIDcol]
      } else if(minPropMatched == 1){
        addInfo <- indFiles[[1]][,idCols]
        row.names(addInfo) <- indFiles[[1]][,seqIDcol]
      } else {
        tempInfo <- do.call(rbind,lapply(indFiles,function(x) { x[,c(seqIDcol,idCols)]} ))
        tempFull <- tempInfo[!duplicated(tempInfo[,1]),] ## Keep only unique info
        row.names(tempFull) <- tempFull[,1]
        addInfo <- tempFull[,-1]
      }
      
    } else{
      addInfo <- NULL
    }
        
    
    ## Combine the files
    if(keepOrder){
      tempMerge <- do.call(cbind,colsInt)
    } else{
#       colsIntTemp <- lapply(colsInt,head) ## For testing
      colsIntTemp <- colsInt
      for(k in 1:length(colsIntTemp)){
        names(colsIntTemp[[k]]) <- paste(names(colsIntTemp)[k],names(colsIntTemp[[k]]),sep=".")
      }
      
      tempMerge <- merge(colsIntTemp[[1]],colsIntTemp[[2]],by=0,all=TRUE)
      if(length(colsIntTemp) > 2){
        for(k in 3:length(colsIntTemp)){
          tempMerge <- merge(tempMerge,colsIntTemp[[k]],by.x = 1,by.y=0,all=TRUE)
        }
      } 
      row.names(tempMerge) <- tempMerge[,1]
      tempMerge <- tempMerge[,-1]
    }
    
    ## Make sure that the columns are renamed
    if(names(tempMerge)[1] == colsToKeep[1]){
      names(tempMerge) <- paste(rep(sampName,each=length(colsToKeep)),names(tempMerge),sep=".")
    }
    
    ## Add the additional info
    if(!is.null(addInfo)){
      if(keepOrder){
        mergeOut <- cbind(addInfo,tempMerge)
      } else {
        mergeOut <- cbind(addInfo[row.names(tempMerge),],tempMerge)
      }
      
    } else {
      mergeOut <- tempMerge
    }
    
    return(mergeOut)
  }
