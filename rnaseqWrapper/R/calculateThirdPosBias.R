calculateThirdPosBias <-
function(varTable, # data.frame with rows=genes, cols=inds
                                  seqIDcol = 1, # which column is the sequence ID in?
                                  #                  can be numeric or character
                                  refPosCol = "Reference.Position", # Which column has the reference position?
                                  #                                    can be numeric or character
                                  readCutoffs = 0, # how many variable positions need to be present to calculate
                                  #                      set to 1 (or 0 or NULL) to include all
                                  #                      without a reference, small numbers will be almost meaningless
                                  colprepend = "nVar_", # what should the columns be prepended with?
                                  codonStartPos = NULL # If known, a vector giving which position each codon starts at
                                  #                            currently not handled.
                                  #                            if "cds" assumes all start at position 1.
){
  
  if(is.null(readCutoffs)){
    readCutoffs <- 1
  } else if(!is.numeric(readCutoffs)){
    stop("ERROR: readCutoffs must be numeric or NULL")
  }

  if(is.character(seqIDcol)){
    seqIDcolOut <- (1:length(colnames(varTable)) )[colnames(varTable) == seqIDcol]
    if(length(seqIDcolOut) != 1){
      stop("ERROR: seqIDcol did not match exactly one column name")
    }
  } else if(is.numeric(seqIDcol)){
    seqIDcolOut <- seqIDcol
  } else{
    stop("ERROR: seqIDcol must be numeric or character")
  }
  
  if(is.character(refPosCol)){
    refPosColOut <- (1:length(colnames(varTable)) )[colnames(varTable) == refPosCol]
    if(length(refPosColOut) != 1){
      stop("ERROR: refPosCol did not match exactly one column name")
    }
  } else if(is.numeric(refPosCol)){
    refPosColOut <- refPosCol
  } else{
    stop("ERROR: refPosCol must be numeric or character")
  }
  ## Set ref position to numeric:
  varTable[,refPosColOut] <- as.numeric(as.character(varTable[,refPosColOut]))
  
  ## Split the table by seqID, keep only the reference positions
  splitVars <- split(varTable[,refPosCol],varTable[,seqIDcol])
  nVarPos <- unlist(lapply(splitVars,length))
  limVars <- splitVars[nVarPos >= readCutoffs]

  
  if(is.null(codonStartPos)){
    
    ## calculate the module position of the variants; set max as third pos
    tempPosOfVar <- lapply(limVars,function(x){ tempFactor <- factor(x%%3,levels=0:2)
                                                tempTable <- table(tempFactor)
                                                maxPos <- names(tempTable)[which.max(tempTable)]
                                                
                                                ## change the position names so that max is 2
                                                #     base is c(0,1,2) which is what it should be if 2 is max
                                                if       (maxPos == 0){
                                                  levels(tempFactor) <- c("thirdPos","firstPos","secondPos") 
                                                } else if(maxPos == 1){
                                                  levels(tempFactor) <- c("secondPos","thirdPos","firstPos")
                                                } else if(maxPos == 2){
                                                  levels(tempFactor) <- c("firstPos","secondPos","thirdPos")
                                                }
                                                ## Set a new factor with the levels in the right order
                                                finalFactor <- factor(tempFactor,levels=c("firstPos","secondPos","thirdPos"))
                                                levels(finalFactor) <- paste(colprepend,c("firstPos","secondPos","thirdPos"),sep="")
                                                finalTable <- table(finalFactor)
                                                return(finalTable)
                                                })  
    
    finalPosOfVars <- as.data.frame(do.call(rbind,tempPosOfVar))
    
    warning("Most variable position was assumed to be 3rd position because no reference is present")
  } else if(codonStartPos == "cds"){
    ## "cds" assumes that everything is in frame
    ## calculate the module position of the variants; set max as third pos
    tempPosOfVar <- lapply(limVars,function(x){ finalFactor <- factor(x%%3,levels=0:2)
                                                levels(finalFactor) <- paste(colprepend,c("thirdPos","firstPos","secondPos"),sep="")
                                                
                                                ## Set a new factor with the levels in the right order
                                                finalFactor <- factor(finalFactor,levels=paste(colprepend,c("firstPos","secondPos","thirdPos"),sep=""))
                                                finalTable <- table(finalFactor)
                                                return(finalTable)
    })  
    
    finalPosOfVars <- do.call(rbind,tempPosOfVar)
    
  }  else {
    stop("Sorry, there currently isn't handling for a reference. I will fix this once I know what a ref looks like")
  }
  
  ## Calculate the proportion in each position
  varProp <- t(apply(finalPosOfVars,1,function(x) {x/sum(x)}))
  colnames(varProp) <- paste(colnames(varProp),"_prop",sep="")

  ## combine the files for return
  varPosOut <- cbind(finalPosOfVars,varProp)
  
  return(varPosOut)
}
