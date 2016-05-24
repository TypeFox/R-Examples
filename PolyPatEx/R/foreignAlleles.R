
##  Function foreignAlleles
##
##' Identify 'foreign' alleles (alleles present in progeny, but not in
##' parents or other non-parental adults).  Note that
##' \code{foreignAlleles} does \emph{not} distinguish between
##' populations as indicated by the \code{popn} column of the allele
##' data frame.
##'
##' @title Identify foreign alleles
##' @param adata data frame: the preprocessed allele data set
##' returned by either \code{\link{inputData}} or
##' \code{\link{preprocessData}}.
##' @return A list, containing two data frames.  Data frame
##' \code{byLocus} summarises the foreign alleles found at each locus.
##' Data frame \code{byProgenyLocus} summarises the same alleles by
##' progeny and locus.  In this latter data frame, the code
##' \code{P.missing} indicates no alleles were present in the progeny
##' at this locus.
##'
##' @author Alexander Zwart (alec.zwart at csiro.au)
##' @export
##' @examples
##'
##' ## Using the example dataset 'FR_Genotype':
##' data(FR_Genotype)
##'
##' ## Since we did not load this dataset using inputData(), we must
##' ## first process it with preprocessData() before doing anything
##' ## else:
##' gData <- preprocessData(FR_Genotype,
##'                         numLoci=7,
##'                         ploidy=4,
##'                         dataType="genotype",
##'                         dioecious=TRUE,
##'                         mothersOnly=TRUE)
##'
##' head(gData)  ## Checked and Cleaned version of FR_Genotype
##'
##' fAlleles <- foreignAlleles(gData)
##'
##' ## View foreign alleles detected at each locus:
##' fAlleles$byLocus
##'
##' ## View foreign alleles detected in each progeny, at each locus:
##' fAlleles$byProgenyLocus
##'
##' ## Both of these objects are data frames, hence can be written to file
##' ## via, e.g., write.csv().
##'
foreignAlleles <- function(adata) {
  ##
  ## Defined as an explicit function for easier debugging...
  FAByPLChecker <- function(vv) {
    vv <- stripNAs(vv)
    if (length(vv)==0) {
      ##Missing progeny alleles
      return("P.missing")
    } else {
      vv <- vv[!(vv %in% popnAlleles)]
      if (length(vv) == 0) {
        ##No external alleles
        return(NA)
      } else {
        ##External alleles
        return(paste(vv,collapse=" "))
      }
    }
  }
  ##
  numLoci <- attr(adata,"numLoci")
  ploidy <- attr(adata,"ploidy")
  dioecious <- attr(adata,"dioecious")
  selfCompatible <- attr(adata,"selfCompatible")
  ##
  progeny <- with(adata,id[!is.na(mother)])
  allAdults <- with(adata,id[is.na(mother)])
  byProgenyLocus <- data.frame(row.names=progeny)
  byLocus <- vector()
  ##
  for (locus in 1:numLoci) {
    locusRange <- (3+dioecious) + (locus-1)*ploidy + 1:ploidy
    popnAlleles <- sort(unique(stripNAs(unlist(adata[allAdults,locusRange]))))
    progenyTypes <- adata[progeny,locusRange]  ##Pheno- or geno-
    byProgenyLocus[[locus]] <- as.character(apply(progenyTypes,1,FAByPLChecker))
    uniqueAllelesAtLocus <- unique(stripNAs(unlist(progenyTypes)))
    FAL <- paste(sort(uniqueAllelesAtLocus[!(uniqueAllelesAtLocus %in%
                                        popnAlleles)]),collapse=" ")
    if (length(FAL) == 0) {
      is.na(byLocus[locus]) <- TRUE
    } else {
      byLocus[locus] <- FAL
    }
  }
  loci <- paste("Locus",1:numLoci,sep="")
  names(byProgenyLocus) <- loci
  ##
  return(list(byLocus = data.frame(Locus = loci,foreignAlleles = byLocus),
              byProgenyLocus = byProgenyLocus))
}


