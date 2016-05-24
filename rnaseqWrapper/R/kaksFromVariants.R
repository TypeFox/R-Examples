kaksFromVariants <-
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
    
    if(is.null(codonStartPos)){
      stop("when codonStartPos is unknown, this analysis is fatally flawed")
    } else if(codonStartPos == "cds"){
      ## "cds" assumes that everything is in frame
      
      ## For testing
      # limVars <- limVars[1:100]
      # x <- limVars[[1]] ## For testing
      # x <- limVars[["Beilschmiedia_comp124531_c0_seq10"]]
      
      byGeneFullInfo <- lapply(
        limVars,
        function(x){ 
          whichSeq <- as.character(x[1,1])
          seq <- referenceSeqs[[whichSeq]] 
          # cat(whichSeq,";") ## For diagnosing problems
          # grabs the sequence name, and pulls the matching fasta
          
          mults <- grep(",",x$Var)
          splitVar <- strsplit(x$Var,",")
          nVar <- unlist(lapply(splitVar,length))
          
          synthSeqs <- list(reference = paste(seq,collapse=""))
          for(k in 1:max(nVar)){
            tempSeq <- seq
            tempSeq[x$Position] <- unlist(lapply(splitVar,function(y) {y[min(k,length(y))]}))
            
            synthSeqs[[paste("variant",k,sep="")]] <- paste(tempSeq,collapse="")
          }
          
          forKaKS <- as.alignment(nb = length(names(synthSeqs)),
                          nam = names(synthSeqs),
                          seq = unlist(synthSeqs),
                          com= NA)

          kaTemp <- kaks(forKaKS,verbose=TRUE)
          toGrab <- lapply(kaTemp,function(w) {as.matrix(w)})
          
          ka <- mean(toGrab$ka[-1,1])
          ks <- mean(toGrab$ks[-1,1])
          KaKs <- ka/ks
          KaKs[ka < 0 | ks < 0] <- NA
          
          
          nSynSites <- mean(toGrab$l4[-1,1])+mean(toGrab$l2[-1,1])/3
          nNonSynSites <- mean(toGrab$l0[-1,1])+mean(toGrab$l2[-1,1])*2/3
          
          nVariants <- sum(nVar)
          
          tempOut <- data.frame(ka, ks, KaKs, nSynSites, nNonSynSites,nVariants)
          return(tempOut)
        }) ## close the lapply loop
      
      out <- do.call(rbind,byGeneFullInfo)
#       out <- list(variantInfo=variantInfo,geneInfo=geneInfo)
      
      return(out)
      
    }  else {
      stop("Sorry, there currently isn't handling for variable codon start sites. I will fix this once I know what a ref looks like")
    }
    
  }
