
##  Function getProgenyStatusGenot
##
## Obtain the status of the mother-progeny relationship at each
## locus for genotypic data.
##
## \code{getProgenyStatusGenot} determines the nature of the
## relationship between mother and progeny at each locus.  This is
## both for the benefit of subsequent processing by PolyPatEx
## routines, and for the benefit of the package user.
##
## Note that \code{getProgenyStatusGenot} assumes that the genotype
## allele sets in the input data frame are either \sQuote{full} (all
## alleles present) or \sQuote{missing} (no alleles present).
## Incomplete allele sets in the original dataset should have
## previously been set to be missing by \code{\link{preprocessData}}.
##
## The possible status codes are:
## \describe{
## \item{\code{M.missing}}{Maternal alleles set is missing.}
## \item{\code{P.missing}}{Progeny allele set is missing.}
## \item{\code{MP.missing}}{Both Maternal and Progeny
##       alleles sets are missing.}
## \item{\code{MP.noMatch}}{maternal allele set cannot provide a
##            gamete compatible with the progeny allele set.}
## \item{\code{MAO}}{\sQuote{Maternal alleles only} - Progeny
## contains only alleles that are present in the mother.}
## \item{\code{NMA}}{\sQuote{Non-Maternal alleles} - Progeny contains
## one or more alleles that are not present in the mother.}
## }
##
## @title Determine mother progeny relationship at each locus
## @param adata A data frame:
## @return A list, containing three data frames:
## \code{progenyStatusTable}, \code{MP.alleleTable} and
## \code{PNotM.alleleTable}
## @author Alexander Zwart (alec.zwart at csiro.au)
## @examples
## \dontrun{
##
## ## Assuming 'adata' is a genotype allele data frame previously
## ## returned by inputData() or preprocessData():
##
## getProgenyStatusGenot(adata)
##
## }
##
getProgenyStatusGenot <- function(adata) {
  ##Note: getProgenyStatusGenot() assumes either complete, or wholly
  ## missing allelesets only.  Insufficient allelesets should have
  ## already been set to missing by processData().
  ##
  numLoci <- attr(adata,"numLoci")
  ploidy <- attr(adata,"ploidy")
  dioecious <- attr(adata,"dioecious")
  dataType <- attr(adata,"dataType")
  selfCompatible <- attr(adata,"selfCompatible")
  ##
  progeny <- with(adata,id[!is.na(mother)])
  ##progenyMothers matches up with progeny - allows slightly
  ## faster indexing than 'adata$mother[adata$id==thisProgeny]',
  ## although a hash would be better still - see the 'hash' package?
  progenyMothers <- adata[progeny,"mother"]
  names(progenyMothers) <- progeny  ##We now index by name
  progenyStatusTable <- matrix(NA_character_,
                                nrow=length(progeny),
                                ncol=numLoci,
                                dimnames=list(progeny,
                                paste("Locus",1:numLoci,sep="")))
  ##Alleles shared by _progeny_ and _mother_
  MP.alleleTable <- matrix(NA_character_,
                           nrow=length(progeny),
                           ncol=numLoci,
                           dimnames=list(progeny,
                           paste("Locus",1:numLoci,sep="")))
  ##Alleles present in _progeny_ but not in _mother_
  PNotM.alleleTable <- matrix(NA_character_,
                              nrow=length(progeny),
                              ncol=numLoci,
                              dimnames=list(progeny,
                              paste("Locus",1:numLoci,sep="")))
  ##
  for (locus in 1:numLoci) {
    locusRange <- (3 + dioecious) + (locus-1)*ploidy + 1:ploidy
    for (thisProgeny in progeny) {
      pAlleles <- stripNAs(adata[thisProgeny,locusRange])
      mAlleles <- stripNAs(adata[progenyMothers[thisProgeny],locusRange])
      noPAlleles <- length(pAlleles)==0
      noMAlleles <- length(mAlleles)==0
      if (noPAlleles && noMAlleles) {
        progenyStatusTable[thisProgeny,locus]  <- "MP.missing"
        ## MP.alleleTable, PNotM.alleleTable left as NA
      } else if (noPAlleles) {
        progenyStatusTable[thisProgeny,locus]  <- "P.missing"
        ## MP.alleleTable, PNotM.alleleTable left as NA
      } else if (noMAlleles) {
        progenyStatusTable[thisProgeny,locus]  <- "M.missing"
        ## MP.alleleTable, PNotM.alleleTable left as NA
      } else {  ##Allele sets both present - compare 'em:
        MP.elements <- make.unique(pAlleles) %in% make.unique(mAlleles)
        PNotM.alleles <- pAlleles[!MP.elements]
        MP.alleles <- pAlleles[MP.elements]
        if (length(MP.alleles) < ploidy/2) {  ## Progeny-mother mismatch
          progenyStatusTable[thisProgeny,locus] <- "MP.noMatch"
          ##PNotM.alleleTable & MP.alleleTable entries left as NA
        } else if (length(PNotM.alleles) == 0) { ##Maternal alleles only
          progenyStatusTable[thisProgeny,locus] <- "MAO"
          MP.alleleTable[thisProgeny,locus] <- paste(MP.alleles,
                                                      collapse=" ")
          ##PNotM.alleleTable entry left as NA (as is appropriate)
        } else {  ## Non-maternal alleles present
          progenyStatusTable[thisProgeny,locus] <- "NMA"
          MP.alleleTable[thisProgeny,locus] <- paste(MP.alleles,
                                                      collapse=" ")
          PNotM.alleleTable[thisProgeny,locus] <- paste(PNotM.alleles,
                                                         collapse=" ")
        }
      }
    }  ## End thisProgeny loop
  }  ## End locus loop
  attr(progenyStatusTable,"progenyMothers") <- progenyMothers
  return(list(progenyStatusTable=as.data.frame(progenyStatusTable,
              stringsAsFactors=FALSE),
              MP.alleleTable=as.data.frame(MP.alleleTable,
              stringsAsFactors=FALSE),
              PNotM.alleleTable=as.data.frame(PNotM.alleleTable,
              stringsAsFactors=FALSE)))
}
