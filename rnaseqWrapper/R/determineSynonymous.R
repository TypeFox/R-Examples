determineSynonymous <-
function(varTable, # data.frame with rows=genes, cols=inds
                                seqIDcol = 1, # which column is the sequence ID in?
                                #                  can be numeric or character
                                refPosCol = "Reference.Position", # Which column has the reference position?
                                #                                    can be numeric or character
                                refAlleleCol = "Reference", # Which column has the reference allele?
                                #                                    can be numeric or character
                                varAlleleCol = "Allele", # Which column has the variable alleles?
                                #                                    will be grepped, so make it unique enough
                                readCutoffs = 1, # how many variable positions need to be present to calculate
                                #                      set to 1 (or 0 or NULL) to include all
                                #                      without a reference, small numbers will be almost meaningless
                                colprepend = "snvs_", # what should the columns be prepended with?
                                codonStartPos = "cds", # File with reference position information
                                referenceSeqs # FASTA format to pull out the sequences
                                ){

  if(!require(seqinr)){
    stop("The package 'seqinr' is required to run this function. ",
         "Please install it now using:\n\n",
         'install.packages("seqinr")\n')
  }
  
  
  if(is.null(readCutoffs)){
    readCutoffs <- 1
  } else if(!is.numeric(readCutoffs)){
    stop("ERROR: readCutoffs must be numeric or NULL")
  }
  
    
  ## Make all column inputs numeric
  alleleCols <- grep(varAlleleCol,colnames(varTable))
  
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
  
  
  if(is.character(refAlleleCol)){
    refAlleleColOut <- (1:length(colnames(varTable)) )[colnames(varTable) == refAlleleCol]
    if(length(refAlleleColOut) != 1){
      stop("ERROR: refAlleleCol did not match exactly one column name")
    }
  } else if(is.numeric(refAlleleCol)){
    refAlleleColOut <- refAlleleCol
  } else{
    stop("ERROR: refAlleleCol must be numeric or character")
  }
  
  varTableLim <- varTable[,c(seqIDcolOut,refPosColOut,refAlleleColOut,alleleCols)]
  
  ## Split the table by seqID, keep only the reference positions
  splitVars <- split(varTableLim,varTableLim[,1])
  
  nVarPos <- unlist(lapply(splitVars,function(x) length(x[,1])))
  limVars <- splitVars[nVarPos >= readCutoffs]
  
  ## for testing:
  #limVars <- limVars[1:100]
  
  codon <- c("CTT","ATG","ACA","ACG","ATC","AAC","ATA","AGG",
             "CCT","ACT","AGC","AAG","AGA","CAT","AAT","ATT",
             "CTG","CTA","CTC","CAC","AAA","CCG","AGT","CCA",
             "CAA","CCC","TAT","GGT","TGT","CGA","CAG","TCT",
             "GAT","CGG","TTT","TGC","GGG","TAG","GGA","TAA",
             "GGC","TAC","TTC","TCG","TTA","TTG","TCC","ACC",
             "TCA","GCA","GTA","GCC","GTC","GCG","GTG","GAG",
             "GTT","GCT","TGA","GAC","CGT","TGG","GAA","CGC")
  aminoAcid <- c("L","M","T","T","I","N","I","R","P","T","S",
                 "K","R","H","N","I","L","L","L","H","K","P",
                 "S","P","Q","P","Y","G","C","R","Q","S","D",
                 "R","F","C","G","*","G","*","G","Y","F","S",
                 "L","L","S","T","S","A","V","A","V","A","V",
                 "E","V","A","*","D","R","W","E","R")
  codonCode <- data.frame(codon,aminoAcid,stringsAsFactors=FALSE)
  rownames(codonCode) <- codonCode$codon
  
  if(is.null(codonStartPos)){
    stop("when codonStartPos is unknown, this analysis is fatally flawed")
  } else if(codonStartPos == "cds"){
    ## "cds" assumes that everything is in frame
    
    byGeneFullInfo <- lapply(limVars,
                             function(x){ whichSeq <- as.character(x[1,1])
                                          seq <- referenceSeqs[[whichSeq]] 
                                          # cat(whichSeq,";") ## For diagnosing problems
                                          # grabs the sequence name, and pulls the matching fasta
                                          
                                          # Get the variants, including mulitple, if they exist
                                          varAllelesTemp <- apply(x[4:length(colnames(x))],1,function(y) {levels(as.factor(y))} )
                                          
                                          varAllele <-  character()
                                          refPos <- integer()
                                          for(k in 1:length(varAllelesTemp)){
                                            ## This catches multiple alleles and separates them
                                            varAllele <- c(varAllele,varAllelesTemp[[k]])
                                            refPos <- c(refPos,rep(x[k,2],length(varAllelesTemp[[k]])))
                                          }
                                          toMerge <- data.frame(refPos,varAllele)
                                          

                                          ## Add the reference position and alleles
                                          temp <- merge(x[,2:3],toMerge,by.x=1,by.y=1)
                                          rownames(temp) <- paste(temp[,1],temp[,3],sep="")
                                          
                                          ## Get the codon position, change to proper
                                          temp$codonPosition <- temp[,1]%%3 
                                          temp$codonPosition[temp$codonPosition == 0] <- 3
                                          
                                          ## Get the codon
                                          codonStart <- temp[,1]-(temp$codonPosition - 1)
                                          tempRefCodon <- list() 
                                          for(k in 1:length(temp[,1])){
                                            tempRefCodon[[k]] <- seq[codonStart[k]:(codonStart[k]+2)]
                                          }
                                          
                                          tempVarCodon <-  do.call(rbind,tempRefCodon)
                                          # Make the variant codons
                                          for(k in 1:length(tempVarCodon[,1])){
                                            tempVarCodon[k,temp$codonPosition[k]] <- as.character(temp$varAllele[k])
                                          }
                                          
                                          temp$refCodon <- as.character(do.call(rbind,lapply(tempRefCodon,function(x) {paste(x,collapse="")})))
                                          temp$varCodon <- as.character(apply(tempVarCodon,1, function(x) {paste(x,collapse="")}))
                                          
                                          
                                          temp$refAminoAcid <- codonCode[temp$refCodon,"aminoAcid"]
                                          temp$varAminoAcid <- codonCode[temp$varCodon,"aminoAcid"]
                                          
                                          tempSynonomous <- factor(temp$refAminoAcid == temp$varAminoAcid,levels=c("TRUE","FALSE"))
                                          levels(tempSynonomous) <- c("synonymous","nonSynonymous")
                                          temp$varType <- tempSynonomous
                                          
                                          tempOut <- data.frame(SEQ_ID=whichSeq,temp)
                                          return(tempOut)
    }) ## close the lapply loop
    
    byGeneResults <- lapply(byGeneFullInfo,
                            function(x){whichSeq <- as.character(x[1,1])
                                        
                                        nSynonymous <- length(x$varType[x$varType=="synonymous"])
                                        nNonSynonymous <- length(x$varType[x$varType=="nonSynonymous"])
                                        dNdS <- nNonSynonymous/nSynonymous
                                        
                                        tempOut <- data.frame(SEQ_ID=whichSeq,nSynonymous,nNonSynonymous,dNdS)
                                        return(tempOut)
                                        })
    
    geneInfo <- do.call(rbind,byGeneResults) 
    variantInfo <- do.call(rbind,byGeneFullInfo)
    
    out <- list(variantInfo=variantInfo,geneInfo=geneInfo)
       
    return(out)
    #write.table(variantInfo,file='quickTestFullInfoCORRECTED.txt',sep="\t",row.names=FALSE,quote=FALSE)
    #write.table(synOut[[1]],file='quickTestFullInfoCORRECTED_allVariants.txt',sep="\t",row.names=FALSE,quote=FALSE)
    #write.table(geneInfo,file='quickTestGeneInfoCORRECTED.txt',sep="\t",row.names=FALSE,quote=FALSE)
    #write.table(synOut[[2]],file='quickTestGeneInfoCORRECTED_allVariants.txt',sep="\t",row.names=FALSE,quote=FALSE)
    
  }  else {
    stop("Sorry, there currently isn't handling for variable codon start sites. I will fix this once I know what a ref looks like")
  }
  
}
