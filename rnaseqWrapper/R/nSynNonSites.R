nSynNonSites <-
  function(geneNames, # character vector
           codonStartPos = "cds", # File with reference position information
           referenceSeqs # FASTA format to pull out the sequences
  ){

    if(!require(seqinr)){
      stop("The package 'seqinr' is required to run this function. ",
           "Please install it now using:\n\n",
           'install.packages("seqinr")\n')
    }
    
    if(is.null(codonStartPos)){
      stop("when codonStartPos is unknown, this analysis is fatally flawed")
    } else if(codonStartPos == "cds"){
      ## "cds" assumes that everything is in frame
      
      ## For testing
      # geneNames <- geneNames[1:100]
      # x <- geneNames[[1]] ## For testing
      # x <- limVars[["Beilschmiedia_comp124531_c0_seq10"]]
      
      byGeneFullInfo <- lapply(
        geneNames,
        function(x){ 
          whichSeq <- as.character(x)
          seq <- referenceSeqs[[whichSeq]] 
          # cat(whichSeq,";") ## For diagnosing problems
          # grabs the sequence name, and pulls the matching fasta
          
          synthSeqs <- list(reference = paste(seq,collapse=""),reference2 = paste(seq,collapse=""))
          
          
          forKaKS <- as.alignment(nb = length(names(synthSeqs)),
                                  nam = names(synthSeqs),
                                  seq = unlist(synthSeqs),
                                  com= NA)
          
          kaTemp <- kaks(forKaKS,verbose=TRUE)
          toGrab <- lapply(kaTemp,function(w) {as.matrix(w)})
          
#           ka <- mean(toGrab$ka[-1,1])
#           ks <- mean(toGrab$ks[-1,1])
#           KaKs <- ka/ks
          
          nSynSites <- mean(toGrab$l4[-1,1])+mean(toGrab$l2[-1,1])/3
          nNonSynSites <- mean(toGrab$l0[-1,1])+mean(toGrab$l2[-1,1])*2/3
          
#           nVariants <- sum(nVar)
          
          tempOut <- data.frame(nSynSites, nNonSynSites,row.names=whichSeq)
          return(tempOut)
        }) ## close the lapply loop
      
      out <- do.call(rbind,byGeneFullInfo)
      #       out <- list(variantInfo=variantInfo,geneInfo=geneInfo)
      
      return(out)
      
    }  else {
      stop("Sorry, there currently isn't handling for variable codon start sites. I will fix this once I know what a ref looks like")
    }
    
  }
