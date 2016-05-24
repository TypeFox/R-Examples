
##  Function multilocusTypes
##
##' Return summaries of individual- and multi- locus genotypes for
##' adults and progeny.
##'
##' Function \code{multilocusTypes} summarises the different genotypes
##' present at each locus in the dataset (separately for progeny and
##' adults), and across the loci (again, separately for progeny and
##' adults). \code{multilocusTypes} returns a list structure with
##' several elements.
##'
##' @title Genotype summaries
##' @param adata data frame: the checked and preprocessed dataset
##' returned by \code{\link{inputData}}.
##' @return A list structure, with the following components:
##'
##' \describe{
##'
##' \item{\code{uniqueProgenyTypes}}{A data frame containing, for each
##' locus, the distinct genotypes that are present in the progeny in
##' the dataset, and the numbers of progeny containing each genotype
##' at that locus.}
##'
##' \item{\code{numUniqueProgenyTypes}}{The number of unique genotypes
##' at each locus in the progeny in the dataset.}
##'
##' \item{\code{uniqueAdultTypes}}{A data frame containing, for each
##' locus, the distinct genotypes that are present in the adults in
##' the dataset, and the numbers of adults containing each genotype at
##' that locus.}
##'
##' \item{\code{numUniqueAdultTypes}}{The number of unique genotypes
##' at each locus in the adult set.}
##'
##' \item{\code{uniqueProgenyMLTypes}}{A data frame containing the
##' distinct genotypes \emph{across all loci} that are present in the
##' progeny in the dataset, and the numbers of progeny containing each
##' multilocus genotype.}
##'
##' \item{\code{numUniqueProgenyMLTypes}}{The total number of progeny
##' multilocus genotypes.}
##'
##' \item{\code{uniqueAdultMLTypes}}{A data frame containing, the
##' distinct genotypes \emph{across all loci} that are present in the
##' adults in the dataset, and the numbers of adults containing each
##' multilocus genotype.}
##'
##' \item{\code{numUniqueAdultMLTypes}}{The total number of adult
##' multilocus genotypes.}
##'
##' }
##'
##' @author Alexander Zwart (alec.zwart at csiro.au)
##' @importFrom utils stack
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
##' mTypes <- multilocusTypes(gData)
##'
##' ## mTypes is a list structure - individual components can be
##' ## printed to the screen, or saved to file via, e.g. read.csv().
##'
##' mTypes$numUniqueProgenyTypes
##'
##' ## Components of mTypes
##' names(mTypes)
##'
multilocusTypes <- function(adata) {
  ##
  checkForValidPPEDataset(adata)
  numLoci <- attr(adata,"numLoci")
  ploidy <- attr(adata,"ploidy")
  dioecious <- attr(adata,"dioecious")
  selfCompatible <- attr(adata,"selfCompatible")
  progeny <- with(adata,id[!is.na(mother)])
  allAdults <- with(adata,id[is.na(mother)])
  progenyTypes <- list(); adultTypes <- list()
  uniqueProgenyTypes <- list() ; uniqueAdultTypes <- list()
  for (locus in 1:numLoci) {
    affectedLocus <- paste("Locus",locus,sep="")
    locusRange <- (3 + dioecious) + (locus-1)*ploidy + 1:ploidy
    sTypes <- apply(adata[progeny,locusRange],
                    1,
                    function(vv) {
                      paste(stripNAs(vv),collapse=" ")})
    aTypes <- apply(adata[allAdults,locusRange],
                    1,
                    function(vv) {
                      paste(stripNAs(vv),collapse=" ")})
    progenyTypes[[affectedLocus]] <- sTypes ##For later use...
    adultTypes[[affectedLocus]] <- aTypes ##For later use...
    uniqueProgenyTypes[[affectedLocus]] <- table(sTypes[sTypes!=""])
    uniqueAdultTypes[[affectedLocus]] <- table(aTypes[aTypes!=""])
  } ##end Loci Loop
  numUniqueProgenyTypes <- sapply(uniqueProgenyTypes,length)
  numUniqueAdultTypes <- sapply(uniqueAdultTypes,length)
  ##
  progenyMLTypes <- do.call("paste",c(progenyTypes,sep=" | "))
  names(progenyMLTypes) <- progeny
  adultMLTypes <- do.call("paste",c(adultTypes,sep=" | "))
  names(adultMLTypes) <- allAdults
  uniqueProgenyMLTypes <- table(progenyMLTypes)
  uniqueAdultMLTypes <- table(adultMLTypes)
  numUniqueProgenyMLTypes <- length(uniqueProgenyMLTypes)
  numUniqueAdultMLTypes <- length(uniqueAdultMLTypes)
  ##Restructure uniqueProgenyTypes for output
  cc <- utils::stack(lapply(uniqueProgenyTypes,as.vector))
  nn <- utils::stack(lapply(uniqueProgenyTypes,names))
  rm(uniqueProgenyTypes)
  uniqueProgenyTypes <- data.frame(Locus=cc$ind,
                                    progenyType=nn$values,
                                    nIndividuals=cc$values)
  ##Restructure uniqueAdultTypes for output
  cc <- utils::stack(lapply(uniqueAdultTypes,as.vector))
  nn <- utils::stack(lapply(uniqueAdultTypes,names))
  rm(uniqueAdultTypes)
  uniqueAdultTypes <- data.frame(Locus=cc$ind,
                                 adultType=nn$values,
                                 nIndividuals=cc$values)
  ##Restructure uniqueProgenyMLTypes
  tt <- strsplit(names(uniqueProgenyMLTypes),split=" | ",fixed=TRUE)
  tt <- lapply(tt,
               function(vv){
                 if (length(vv) < numLoci) {
                   return(c(vv,""))
                 } else {
                   return(vv)
                 }
               })  ##Correct for splits at end (see ?strplit)
  tt <- as.data.frame(do.call(rbind,tt))
  uniqueProgenyMLTypes <- cbind(tt,uniqueProgenyMLTypes)
  names(uniqueProgenyMLTypes) <- c(paste("Locus",1:numLoci,sep=""),
                                    "progenyMLType","Freq")
  ##Restructure uniqueAdultMLTypes
  tt <- strsplit(names(uniqueAdultMLTypes),split=" | ",fixed=TRUE)
  tt <- lapply(tt,
               function(vv){
                 if (length(vv) < numLoci) {
                   return(c(vv,""))
                 } else {
                   return(vv)
                 }
               })  ##Correct for splits at end (see ?strplit)
  tt <- as.data.frame(do.call(rbind,tt))
  uniqueAdultMLTypes <- cbind(tt,uniqueAdultMLTypes)
  names(uniqueAdultMLTypes) <- c(paste("Locus",1:numLoci,sep=""),
                                 "adultMLType","Freq")
  ####################################################################
  ## T and C want to know the ID's of the individuals in each type.  I
  ## think they can most easily do this in Excel, with a bit of
  ## sorting, provided I dump out the original data, processed by
  ## preprocessData() and with alleleSets concatenated into single
  ## strings...
  alleleSets <- data.frame(row.names=rownames(adata))
  ##Concatenate the alleleSets into single strings:
  for (locus in 1:numLoci) {
    locusRange <- (3 + dioecious) + (locus-1)*ploidy + 1:ploidy
    alleleSets[[locus]] <- apply(adata[,locusRange],
                                     1,
                                     function(vv) {
                                       paste(stripNAs(vv),collapse=" ")})
  }
  names(alleleSets) <- c(paste("Locus",1:numLoci,sep=""))
  ## Add a Multilocus type column
  alleleSets$allMLTypes <- apply(alleleSets,1,
                                 function(vv){
                                   paste(
                                         paste(
                                               paste("Locus",1:numLoci,sep=""),
                                               vv),
                                         collapse=" - ")
                                 })
  ##Add the preliminary columns of adata
  alleleSets <- cbind(adata[,1:(3+dioecious)],alleleSets)
  ##
  return(list(uniqueProgenyTypes=uniqueProgenyTypes,
              numUniqueProgenyTypes=numUniqueProgenyTypes,
              uniqueAdultTypes=uniqueAdultTypes,
              numUniqueAdultTypes=numUniqueAdultTypes,
              uniqueProgenyMLTypes=uniqueProgenyMLTypes,
              numUniqueProgenyMLTypes=numUniqueProgenyMLTypes,
              uniqueAdultMLTypes=uniqueAdultMLTypes,
              numUniqueAdultMLTypes=numUniqueAdultMLTypes))
}
