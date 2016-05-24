readVariantFiles <-
function(fileDir, # name of the directory to read in
                             sepSymbol="_", # Which symbol separates the sample names from other info?
                             #                  set to "" if the full file name should be used
                             fileID = "*_variants.txt", # a character to use to limit which files are imported
                             #                              regular expressions allowed
                             firstColName = "SEQ_ID", # What should the first column be renamed to
                             #                            set to NULL or "" to leave the column as is.
                             fileSep = "\t", # The column delimiter used in the file (e.g. "," or "\t")
                             idCols = 5, # How many columns of position information are there?
                             refPosCol = "Reference.Position", # Which column has the reference position?
                             #                                    can be numeric or character
                             colToSort = "Coverage", # Which column should be used to keep one line per position?
                             #                          can be numeric or character
                             removeDups = TRUE, # Should duplicates at a position be removed?
                             #                    this is necessary to avoid massive over merging
                             returnMerged = TRUE, # Should the merged variants be returned?
                             returnSing = FALSE, # Should each of the separate variant files be returned?
         limitGenes = NULL, # A character vector of the genes to be analyzed (omits others)
         omitRefMatches = TRUE, # should 'variants' which match the reference be excluded?
         refAlleleCol = "Reference$", # Which column has the reference allele?
         varAlleleCol = "Allele" # Which column has the variable alleles?
){
  ## Based on returnMerged & returnSing returns:
  ##    if both TRUE:    a list with mergedVariants, and a list of the singVariants (cleaned if removeDups=TRUE)
  ##    if returnMerged: a data.frame with the mergedVariants
  ##    if returnSing:   a list of the singVariants (cleaned if removeDups=TRUE)
  
  
  ## Add slash to fileDir, if needed:
  if(length(grep("/$",fileDir))==0){
    fileDir <- paste(fileDir,"/",sep="")
  }
  
  

  
  
  ## get the file names from that directory:
  filesToAnal <- list.files(fileDir,pattern=fileID)
  sampNames <- do.call(rbind,strsplit(filesToAnal,sepSymbol))[,1]
  
  message('Reading in variants','\n')
  cohortVariants <- list()
  for(k in 1:length(sampNames)){
    message('\t reading variants for ',sampNames[k],'\n')
    cohortVariants[[k]] <- read.table(paste(fileDir,filesToAnal[k],sep=""),sep=fileSep,header=TRUE)
    if(!is.null(limitGenes)){
      cohortVariants[[k]] <-cohortVariants[[k]][cohortVariants[[k]][,1] %in% limitGenes,]
      
      ## Remove excess levels
      for(l in 1:length(colnames(cohortVariants[[k]]))){
        if(is.factor(cohortVariants[[k]][,l])){
          cohortVariants[[k]][,l] <- factor(cohortVariants[[k]][,l])
        }
      }
      
    }
    if(omitRefMatches){
      
      if(k==1){ ## set up the column IDs the first time
        if(is.character(refAlleleCol)){
          refAlleleCol <- grep(refAlleleCol,names(cohortVariants[[k]]))
          if(length(refAlleleCol) != 1){
            stop("ERROR: refAlleleCol did not match exactly one column name")
          }
        } else if(is.numeric(refAlleleCol)){
          if(length(refAlleleCol) != 1){
            stop("ERROR: refAlleleCol did not match exactly one column name")
          }
        } else{
          stop("ERROR: refAlleleCol must be numeric or character")
        }
        
        if(is.character(varAlleleCol)){
          varAlleleCol <- grep(varAlleleCol,names(cohortVariants[[k]]))
          if(length(varAlleleCol) != 1){
            stop("ERROR: varAlleleCol did not match exactly one column name")
          }
        } else if(is.numeric(varAlleleCol)){
          if(length(varAlleleCol) != 1){
            stop("ERROR: varAlleleCol did not match exactly one column name")
          }
        } else{
          stop("ERROR: varAlleleCol must be numeric or character")
        }
      } ## Done setting allele columns
      
      ## Find out which variants match the reference, then exclude thems
      isVarRef <- as.character(cohortVariants[[k]][,varAlleleCol]) == as.character(cohortVariants[[k]][,refAlleleCol])
      cohortVariants[[k]] <- cohortVariants[[k]][!isVarRef,]
    } ## finish omitting ref matches
    
    names(cohortVariants)[k] <- sampNames[k]
    ## Rename the columns with the sample names for merging
    if(returnMerged){
      names(cohortVariants[[k]])[(idCols+1):length(names(cohortVariants[[k]]))] <- paste(sampNames[k],"_",names(cohortVariants[[k]])[(idCols+1):length(names(cohortVariants[[k]]))],sep="")  
    }
    
    ## Rename the first column to match requested or first one given
    if(!is.null(firstColName) & firstColName!="" & firstColName!=0 & is.character(firstColName)){
      names(cohortVariants[[k]])[1] <- firstColName[1]
    } else {
      ## This sets all of the files to the name of the first one, if firstColName is not set above
      firstColName <- names(cohortVariants[[k]])[1]
    }
  }
  
  
  if(removeDups){
    ## Remove duplicate positions:
    
    ## Get the right column for sorting:
    if (is.character(colToSort)){
      colToSort <- grep(colToSort,names(cohortVariants[[1]]))
    }
    if(colToSort > length(names(cohortVariants[[1]])) | length(colToSort) == 0){
      stop("ERROR: argument colToSort does not select valid columns")
    }
    
    cleanVars<-list()
    message('Removing Duplicate Variants \n')
    for (k in 1:length(cohortVariants)){
      test <- split(cohortVariants[[k]],cohortVariants[[k]][,firstColName])
      ls <- unlist(lapply(test,function(x) {length(x[,1])}))
      vs <- unlist(lapply(test,function(x) {length(levels(factor(x[,refPosCol])))}))
      diff <- ls - vs #how many extra levels
      for(j in 1:length(diff))
        if(diff[j]>0){
          ## Split the file by the first idCols to ensure that both SNVs and Insertions are kept
          ##    this way, only duplicates true duplicates are removed
          tmpFull <- test[[j]]
          ## Remove excess levels
          for(l in 1:idCols){
            if(is.factor(tmpFull[,l])){
              tmpFull[,l] <- factor(tmpFull[,l])
            }
          }
          
          splittingVar <- interaction(tmpFull[,1:idCols],drop=TRUE)
          tmp <- split(tmpFull,splittingVar)
          
          tmp2 <- lapply(tmp,function(x) {x[which.max(x[,colToSort]),]})
          tmp3 <- do.call(rbind,tmp2)
          test[[j]] <- tmp3
        }
      cleanVars[[k]] <- do.call(rbind,test)
      names(cleanVars)[k]<- names(cohortVariants)[k]
    } # end sample loop
    
  } # End if(removeDups)
  
  
  if(returnMerged){
    # Merge the variant files
    #   NB: will automatically merge by all shared column names (e.g. the first idCols, which are not renamed)
    message('Merging Variants together \n\t\t')
    for(k in 1:length(cleanVars)){
      message(names(cleanVars[k]),', ')
      if(k==1){
        merged <- cleanVars[[k]]
      } else {
        merged <- merge(merged,cleanVars[[k]],all=TRUE)
      }
    }  # End loop
  } # end merging
  
  if(returnMerged & returnSing){
    out <- list(mergedVariants = merged, singVariants = cleanVars)
  } else if(returnMerged){
    out <- merged
  } else if(returnSing){
    out <- cleanVars
  }
  
  return(out)
  
}
