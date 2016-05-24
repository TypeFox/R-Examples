
##  Function removeMismatches
##
## Remove mismatching cases between mothers and their progeny.
##
## \code{removeMismatches} calls \code{\link{getProgenyStatusGenot}}
## or \code{\link{getProgenyStatusPhenot}} (as appropriate) to
## identify cases where a progeny's mother is incapable of providing
## a gamete compatible with the progeny's allele set.
##
## Where \code{matMismatches} or fewer loci in a progeny-mother pair
## are mismatches, the offending loci are set to have no alleles in
## the progeny.  This locus will therefore be ignored in subsequent
## analyses.
##
## Where more than \code{matMismatches} loci are mismatches, the
## progeny is removed entirely from the dataset.
##
## \code{removeMismatches} will report details of such mismatches to
## the R console, to help the user to track down and check their
## origins, and also allow them to see whether the problem may lie
## with the mother, rather than in her progeny.
##
## @title Remove mismatches between mother and progeny
## @param adata data frame: an allele dataset.
## @param matMismatches an integer between 0 and \code{numLoci}-1,
## being the maximum number of mismatching alleles between mother and
## offspring that are allowed before the offspring is removed from the
## dataset.
## @return An allele dataset as a data frame, with identified
##              mismatches removed (see details)
## @author Alexander Zwart (alec.zwart at csiro.au)
##
removeMismatches <- function(adata,matMismatches) {
  ##
  checkForValidPPEDataset(adata)
  numLoci <- attr(adata,"numLoci")
  ploidy <- attr(adata,"ploidy")
  dataType <- attr(adata,"dataType")
  dioecious <- attr(adata,"dioecious")
  selfCompatible <- attr(adata,"selfCompatible")
  ##
  if (dataType=="genotype") {
    progenyStatusTable <- getProgenyStatusGenot(adata)$progenyStatusTable
  } else {  ## Phenotypic data
    progenyStatusTable <- getProgenyStatusPhenot(adata)$progenyStatusTable
  }
  ##
  locusMismatches <- progenyStatusTable == "MP.noMatch"
  if (matMismatches==0)
    {
      fewerLocusMismatches <- FALSE
    } else {
      fewerLocusMismatches <- apply(locusMismatches,1,sum) %in% (1:matMismatches)
      ## '>' below returns a named vector, '%in%' does not (shrug).
      names(fewerLocusMismatches) <- rownames(locusMismatches)
    }
  moreLocusMismatches <- apply(locusMismatches,1,sum) > matMismatches
  ##
  if (any(moreLocusMismatches)) {
    multiMismatchProgeny <- names(moreLocusMismatches)[moreLocusMismatches]
    cat("\n Progeny having greater than matMismatches = ",matMismatches,
        " mismatches between",sep="")
    cat("\n mother and progeny have been found.  The relevant progeny will be")
    cat("\n removed from the dataset.  The progeny and their mismatching")
    cat("\n loci are: \n\n")
    ML.Ptab <-  cbind(ProgenyId = multiMismatchProgeny,
                      Loci = apply(locusMismatches[moreLocusMismatches,,drop=FALSE],
                        1,
                        function(vv){
                          paste(which(vv),collapse=", ")
                        })
                      )
    rownames(ML.Ptab) <- NULL
    print(ML.Ptab,quote=FALSE)
    ##Tally the progeny from each mother that are to be removed.
    t1 <- with(adata,table(mother[id %in% multiMismatchProgeny]))
    ##The total number of progeny for each each mother
    t2 <- table(stripNAs(adata$mother))
    ML.Mtab <- data.frame(Mother=names(t1),
                          num.ML.mismatches=as.vector(t1),
                          total.progeny = as.vector(t2[match(names(t1),
                            names(t2))]))
    cat("\n The following table provides a summary for the mothers of the")
    cat("\n progeny reported above.  The table shows, for each affected")
    cat("\n mother, the number of progeny in which mismatches occur, and")
    cat("\n the total number of progeny for that mother. This may help")
    cat("\n identify possible problems with a mother's allele set.\n\n")
    print(ML.Mtab,quote=FALSE)
    ##
    ##Remove progeny with multi-locus mismatches from the dataset.
    adata <- adata[!(adata$id %in% multiMismatchProgeny),]
    ##
  }
  ##
  if (any(fewerLocusMismatches)) {
    ##
    fewerMismatchProgeny <- names(fewerLocusMismatches)[fewerLocusMismatches]
    cat("\n\n Progeny having matMismatches = ",matMismatches,
        " or fewer mismatches between",sep="")
    cat("\n mother and progeny have been found.  The relevant progeny-locus")
    cat("\n combinations will be set to contain no alleles. The progeny and")
    cat("\n their mismatching loci are: \n\n")
    cc <-  cbind(ProgenyId = fewerMismatchProgeny,
                 Loci = apply(locusMismatches[fewerLocusMismatches,,drop=FALSE],
                   1,
                   function(vv){
                     paste(which(vv),collapse=", ")
                     })
                 )
    rownames(cc) <- NULL
    cc <- as.data.frame(cc)
    cc$Loci <- as.numeric(as.character(cc$Loci))
    print(cc,quote=FALSE)
    cat("\n")
    ##
    ##Tally the progeny from each mother that are to be 'repaired'.
    t1 <- with(adata,table(mother[id %in% fewerMismatchProgeny]))
    ##The total number of progeny for each each mother
    t2 <- table(stripNAs(adata$mother))
    SL.Mtab <- data.frame(Mother=names(t1),
                          num.SL.mismatches=as.vector(t1),
                            total.progeny = as.vector(t2[match(names(t1),
                              names(t2))]))
    cat("\n The following table summarises the mothers of the progeny")
    cat("\n reported above, the number of progeny per mother in")
    cat("\n which mismatches occur, and the total number of progeny")
    cat("\n per mother.  This may help identify possible problems with")
    cat("\n a mother's allele set.\n\n")
    print(SL.Mtab,quote=FALSE)
    cat("\n")
    ##
    ##Set the affected 'fewerLocusMismatches' loci to missing.
    cc$Loci <- as.character(cc$Loci) ## To ensure strsplit() works.
    for (thisProgeny in as.character(cc$ProgenyId))
      {
        affectedLoci <- as.numeric(unlist(strsplit(cc$Loci[cc$ProgenyId==thisProgeny],", ")))
        for (thisLocus in affectedLoci)
          {
            locusRange <- (3 + dioecious) + 1:ploidy +
              (thisLocus-1)*ploidy
            is.na(adata[thisProgeny, locusRange]) <- TRUE
          }
      }
    cat("\n Note that the affected loci will subsequently be")
    cat("\n recorded as 'P.missing', since such loci have had")
    cat("\n all progeny alleles removed from the dataset\n\n")
  }
  ##
  return(adata)
}
