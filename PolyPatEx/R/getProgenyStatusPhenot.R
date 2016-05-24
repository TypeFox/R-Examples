
##  Function getProgenyStatusPhenot
##
## Obtain the status of the mother-progeny relationship at each
## locus for phenotypic data.
##
## \code{getProgenyStatusPhenot} determines the nature of the
## relationship between mother and progeny at each locus.  This is
## both for the benefit of subsequent processing by PolyPatEx
## routines, and for the benefit of the package user.
##
## Noting that a 'missing allele set' means that no alleles are
## present at the affected locus, the possibilities are:
## \describe{
## \item{\code{M.missing}}{Maternal alleles set is missing.}
## \item{\code{P.missing}}{Progeny allele set is missing.}
## \item{\code{MP.missing}}{Both Maternal and Progeny
##       alleles sets are missing.}
## \item{\code{MP.noMatch}}{maternal allele set cannot provide a
##            gamete compatible with the progeny allele set.}
## \item{\code{MAO}}{\sQuote{Maternal alleles only} - Progeny
##       contains only alleles that are present in the mother.}
## \item{\code{NMA}}{\sQuote{Non-Maternal alleles} - Progeny
## contains one or more alleles  that are not present in the mother.}
## }
##
## @title Determine mother progeny relationship at each locus
## @param adata data frame:
## @return A list, containing three data frames:
## \code{progenyStatusTable}, \code{MP.alleleTable} and
## \code{PNotM.alleleTable}
## @author Alexander Zwart (alec.zwart at csiro.au)
## @examples
## \dontrun{
##
## ## Assuming 'adata' is a phenotype allele data frame previously
## ## returned by inputData() or preprocessData():
##
## getProgenyStatusPhenot(adata)
##
## }
##
getProgenyStatusPhenot <- function(adata) {
  ##
  numLoci <- attr(adata,"numLoci")
  ploidy <- attr(adata,"ploidy")
  dioecious <- attr(adata,"dioecious")
  dataType <- attr(adata,"dataType")
  selfCompatible <- attr(adata,"selfCompatible")
  ##
  progeny <- with(adata,id[!is.na(mother)])
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
  ## Get the MP comparison lookup table:
  progenyStatusLookup <- Phenot.getProgenyStatusLookup(ploidy)
  ## Columns nM, nP, nMP, MPStatus
  progenyStatusLookup$nM.nP.nMP <- with(progenyStatusLookup,
                                        paste(nM,"_",nP,"_",
                                              nMP,sep=""))
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
        MP.elements <- pAlleles %in% mAlleles
        PNotM.alleles <- pAlleles[!MP.elements]
        MP.alleles <- pAlleles[MP.elements]
        nM.nP.nMP <- paste(length(mAlleles),"_",
                           length(pAlleles),"_",
                           length(MP.alleles),sep="")
        caseNum <- match(nM.nP.nMP,progenyStatusLookup$nM.nP.nMP)
        if(is.na(caseNum)) {
          stop("\n Error in Phenot.getProgenyStatusLookup() - caseNum=NA.  Please notify author of software!\n\n")
        } else {
          progenyStatusTable[thisProgeny,locus] <- progenyStatusLookup$MPStatus[caseNum]
        }
        if (progenyStatusTable[thisProgeny,locus] == "MAO") { ##Maternal alleles only
          MP.alleleTable[thisProgeny,locus] <- paste(MP.alleles,
                                                      collapse=" ")
          ##PNotM.alleleTable entry left as NA (the default)
        } else if (progenyStatusTable[thisProgeny,locus] == "NMA") {  ## Non-maternal alleles present
          MP.alleleTable[thisProgeny,locus] <- paste(MP.alleles,
                                                      collapse=" ")
          PNotM.alleleTable[thisProgeny,locus] <- paste(PNotM.alleles,
                                                         collapse=" ")
        }
        ## ...& if "MP.noMatch", both MP and PNotM allele table
        ## entries are left as NA (the default)
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
