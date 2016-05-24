
##  Function genotPPE
##
##' Conduct a paternity exclusion analysis on a genotype dataset.
##'
##' \code{genotPPE} conducts a basic paternity exclusion analysis on
##' a genotype dataset.
##'
##' For the purposes of the PolyPatEx package, the term
##' \sQuote{genotype} refers to forms of marker data where the allele
##' copy numbers (or multiplicities) are known - hence for a polyploid
##' species of ploidy \emph{p}, there should be exactly \emph{p}
##' alleles detected at each locus, some of which may be repeats of
##' the same allele state.  In PolyPatEx, no allowance is made for
##' undetected alleles in genotype data - allele sets having fewer
##' than \emph{p} alleles present should have been reset to contain no
##' alleles by \code{\link{preprocessData}}.
##'
##' For the above and other reasons, \code{genotPPE} should \bold{NOT}
##' be applied to a dataset that has not been preprocessed by
##' \code{\link{preprocessData}} (either by calling directly
##' \code{\link{preprocessData}} on the data frame directly, or by
##' using \code{\link{inputData}} to load the data from file).
##'
##' The genotype-based paternity analysis is based on simple
##' comparison of genotype allele sets between mother, progeny, and
##' candidate father.  The mother-progeny relationship is assumed to
##' be known.  For genotype data, should a progeny contain only
##' alleles also present in its mother, then a potential father is any
##' candidate that can provide a gamete compatible with the progeny's
##' genotype, given the mother's genotype.
##'
##' @title Simple paternity exclusion for genotype allele data
##' @param adata data frame: the preprocessed allele data set returned
##' by either \code{\link{inputData}} or \code{\link{preprocessData}}.
##' @return  A list whose components are described below.  The
##' components that are probably of primary interest to the user are
##' \code{adultTables$FLCount} and \code{adultTables$VLTotal}.  These
##' are likely to be large tables, so note that the functions
##' \code{\link{potentialFatherCounts}} and
##' \code{\link{potentialFatherIDs}} are available to usefully
##' summarise their content.
##'
##' The list returned by \code{genotPPE} contains two elements,
##' \code{progenyTables} and \code{adultTables}, each of which are
##' themselves lists.
##'
##' Element \code{adultTables} contains the following components:
##'
##' \describe{
##'
##' \item{\code{FLCount}}{Father Loci Count - a matrix, showing for
##' each progeny-candidate combination, the number of loci at which
##' the candidate matches (i.e., could have fathered) the progeny}
##'
##' \item{\code{VLTotal}}{Valid Loci Total - a matrix, showing for
##' each progeny-candidate combination, the total number of loci at
##' which a valid comparison between progeny and candidate could be
##' made.  (Missing allele sets, whether in the original data, or due
##' to progeny-mother mismatches found by
##' \code{\link{preprocessData}}, can result in fewer loci at which
##' progeny-candidate (father) comparisons are possible.)}
##'
##' \item{\code{fatherSummaryTable}}{A matrix, combining the results
##' of \code{FLCount} and \code{VLTotal} (see above) for each
##' progeny-candidate combination in one table.  This is purely for
##' ease of viewing purposes, but note also the functions
##' \code{\link{potentialFatherCounts}} and
##' \code{\link{potentialFatherIDs}} which may provide more useful
##' summary output.}
##'
##' \item{\code{CPNotM.alleleArray}}{A 3D array containing the
##' alleles present in both candidate (father) and progeny, but not in
##' the progeny's mother (for each progeny/candidate/locus
##' combination)}
##'
##' \item{\code{CMP.alleleArray}}{A 3D array containing the alleles
##' present in candidate, progeny and progeny's mother (for each
##' progeny/candidate/locus combination)}
##'
##' \item{\code{simpleFatherArray}}{A 3D array indicating whether each
##' candidate is compatible with each progeny, for each locus}
##'
##' }
##'
##' \code{progenyTables} contains the following components:
##'
##' \describe{
##'
##' \item{\code{progenyStatusTable}}{A data frame, indicating the
##' status of the progeny / mother allele set comparison (for each
##' progeny, at each locus).}
##'
##' \item{\code{MP.alleleTable}}{A data frame containing the alleles
##' that are found in both mother's and progeny's allele sets (for
##' each progeny, at each locus)}
##'
##' \item{\code{PNotM.alleleTable}}{A data frame, containing the
##' alleles in the progeny's allele set, that are \emph{not} present
##' in the mother's allele set(for each progeny, at each locus)}
##'
##' }
##'
##' The status codes in \code{progenyTables$progenyStatusTable} are:
##'
##' \describe{
##'
##' \item{\code{"MAO"}}{Mother Alleles Only - the progeny contains
##' only alleles found also in the mother}
##'
##' \item{\code{"NMA"}}{Non-Mother Alleles - the progeny contains
##' alleles that are not found in the mother}
##'
##' \item{\code{"P.missing"}}{No comparison was possible at this locus
##' because the progeny's allele set was missing}
##'
##' \item{\code{"P.missing"}}{No comparison was possible at this locus
##' because the mother's allele set was missing}
##'
##' \item{\code{"PM.missing"}}{No comparison was possible at this
##' locus because both progeny's and mother's allele sets were
##' missing}
##'
##' }
##'
##' Note that some of the \code{"P.missing"} or \code{"PM.missing"}
##' codes may have arisen due to progeny / mother mismatches found
##' (and corresponding progeny allele sets removed) by
##' \code{\link{preprocessData}}.
##'
##' @author Alexander Zwart (alec.zwart at csiro.au)
##' @seealso \code{\link{phenotPPE}}
##' @importFrom utils flush.console
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
##' gPPE <- genotPPE(gData)  ## Perform the exclusion analyses
##'
##' ## gPPE is a large (and rather ugly!) data structure - see
##' ## functions potentialFatherCounts() and potentialFatherIDs() for
##' ## more useful output from the gPPE object.
##'
genotPPE <- function(adata) {
  ##
  checkForValidPPEDataset(adata)
  dataType <- attr(adata,"dataType")
  if (dataType == "phenotype") {
    stop("\n You appear to have passed a *phenotype* dataset to genotPPE() - stopping...\n\n")
  }
  numLoci <- attr(adata,"numLoci")
  ploidy <- attr(adata,"ploidy")
  dioecious <- attr(adata,"dioecious")
  selfCompatible <- attr(adata,"selfCompatible")
  ##
  progenyStatusInfo <- getProgenyStatusGenot(adata)
  progenyStatusTable <- progenyStatusInfo$progenyStatusTable
  MP.alleleTable <- progenyStatusInfo$MP.alleleTable
  PNotM.alleleTable <- progenyStatusInfo$PNotM.alleleTable
  rm(progenyStatusInfo)
  ######################################################################
  ## Error check - remove once confident of content of progenyStatusTable!
  if (any(is.na(unlist(progenyStatusTable)))) {
    stop("\n NA's in progenyStatusTable - please notify PPE developer!\n\n")
  } else if (any(!(unlist(progenyStatusTable) %in%
                   c("MAO","NMA","MP.noMatch","M.missing","P.missing","MP.missing")))) {
    stop("\n Incorrect flags in progenyStatusTable - please notify PolyPatEx developer!  Stopping...\n\n")
  }
  ######################################################################
  ## Mismatching allele sets should no longer be possible at this
  ## point...
  progeny <- with(adata,id[!is.na(mother)])
  ##progenyMothers matches up with progeny
  progenyMothers <- adata[progeny,"mother"]
  ##
  ##Define the list of potential fathers ('candidates')
  if (dioecious) {
    ##Adult MALES only
    potentialFathers <- with(adata,id[is.na(mother) & gender=="M"])
    ##Technically, gender=="M" should be sufficient, assuming they
    ## have entered the data correctly - but there is no harm in the
    ## above...
  } else {
    ##Any adult
    potentialFathers <- with(adata,id[is.na(mother)])
    ##Self-compatible vs self-incompatible is dealt with later...
  }
  ##
  ##Progeny alleles not present in mother, but present in candidate
  CPNotM.alleleArray <- array(NA_character_,
                              dim=c(numLoci,
                                length(progeny),
                                length(potentialFathers)),
                              dimnames=list(1:numLoci,
                                progeny,
                                potentialFathers))
  ##Progeny alleles present in both mother and candidate
  CMP.alleleArray <- array(NA_character_,
                           dim=c(numLoci,
                             length(progeny),
                             length(potentialFathers)),
                           dimnames=list(1:numLoci,
                             progeny,
                             potentialFathers))
  ##Status of each candidate at each locus, for each progeny
  simpleFatherArray <- array(NA_character_,
                             dim=c(numLoci,
                               length(progeny),
                               length(potentialFathers)),
                             dimnames=list(1:numLoci,
                               progeny,
                               potentialFathers))
  ## Numbers of possibly paternal allele sets out of total number of
  ## 'valid' allele sets (allele sets at which a valid comparison was
  ## possible) for each progeny/candidate combination
  fatherSummaryTable <- matrix(NA_character_,nrow=length(progeny),
                               ncol=length(potentialFathers),
                               dimnames=list(progeny,potentialFathers))
  ##Number of possibly paternal loci - numerators in fatherSummaryTable
  FLCount <- matrix(-9999,nrow=length(progeny),
                    ncol=length(potentialFathers),
                    dimnames=list(progeny,potentialFathers))
  ##Number of 'valid' loci - denominators in fatherSummaryTable
  VLTotal <- matrix(-9999,nrow=length(progeny),
                    ncol=length(potentialFathers),
                    dimnames=list(progeny,potentialFathers))
  ##
  cat("\n Comparing progeny and candidate fathers...\n")
  utils::flush.console()
  ##
  for (locus in 1:numLoci) {
    cat("\n Processing locus",locus,"of",numLoci)
    utils::flush.console()
    locusRange <- (3+dioecious) + (locus-1)*ploidy + 1:ploidy
    ##
    PNotM.alleles <- strsplit(PNotM.alleleTable[,locus]," ")
    names(PNotM.alleles) <- progeny  ## For indexing purposes
    MP.alleles <- strsplit(MP.alleleTable[,locus]," ")
    names(MP.alleles) <- progeny  ## For indexing purposes
    ##
    for (candidate in potentialFathers) {
      CAlleles <- unlist(adata[candidate,locusRange])
      if (any(is.na(CAlleles))) { ## Genotype data
        ## Array indexing order:  Locus, progeny, potential father.
        simpleFatherArray[locus,,candidate] <- "C.missing"
        ##Note: assigned across ALL progeny...
        ## Leave CPNotM and CMP as NA (the default)
        next ## Skip to next candidate
      }
      ##
      for (thisProgeny in progeny) {
        PNotM.alleleSet <- PNotM.alleles[[thisProgeny]]
        MP.alleleSet <- MP.alleles[[thisProgeny]]
        if (progenyStatusTable[thisProgeny,locus]=="MAO") {
          ## Progeny contains maternal alleles only. Candidate must
          ## be able to account for at least ploidy/2 of these alleles
          matchingMPAlleles <- make.unique(CAlleles,sep=" ") %in%
                              make.unique(MP.alleleSet,sep=" ")
          CMP.alleleArray[locus,thisProgeny,candidate] <-
            paste(CAlleles[matchingMPAlleles],collapse=" ")
          ## Leave CPNotM as NA (the default)
          if (sum(matchingMPAlleles) >= ploidy/2) {
            ## Candidate IS a potential father
            simpleFatherArray[locus,thisProgeny,candidate] <- "CP.match"
          } else {
            ## Candidate is NOT a potential father
            simpleFatherArray[locus,thisProgeny,candidate] <- "CP.noMatch"
          }
        } else if (progenyStatusTable[thisProgeny,locus]=="NMA") {
          ## Progeny contains some non-maternal alleles.  The
          ## candidate must match ALL of the non-maternal alleles,
          ## plus as many maternal alleles as needed to make up the
          ## required ploidy/2 total.
          ##
          ## Alleles in candidate that match the PNotM set
          matchingPNotMAlleles <- make.unique(CAlleles,sep=" ") %in%
                                   make.unique(PNotM.alleleSet,sep=" ")
          CPNotM.alleleArray[locus,thisProgeny,candidate] <-
            paste(CAlleles[matchingPNotMAlleles],collapse=" ")
          ## Alleles in candidate, EXCLUDING those already matched to
          ## PNotM set, that match the MP set (recall that this is
          ## genotypic data, so alleles can be repeated).
          matchingMPAlleles <-
            make.unique(CAlleles[!matchingPNotMAlleles],sep=" ") %in%
              make.unique(MP.alleleSet,sep=" ")
          CMP.alleleArray[locus,thisProgeny,candidate] <-
            paste(CAlleles[!matchingPNotMAlleles][matchingMPAlleles],
                  collapse=" ")
          ##
          if (sum(matchingPNotMAlleles)<length(PNotM.alleleSet)) {
            ## Candidate cannot match all non-maternal alleles in
            ## progeny
            simpleFatherArray[locus,thisProgeny,candidate] <- "CP.noMatch"
          } else if (sum(matchingPNotMAlleles) +
                     sum(matchingMPAlleles) >= ploidy/2) {
            simpleFatherArray[locus,thisProgeny,candidate] <- "CP.match"
            ## Note: The case of sum(matchingPNotMAlleles) > ploidy/2
            ## has already been accounted for as an "MP.noMatch"
          } else {
            ## Cannot account for enough alleles
            simpleFatherArray[locus,thisProgeny,candidate] <- "CP.noMatch"
          }
        } else {
          ## Progeny status = "MP.noMatch","M.missing","P.missing" or
          ## "MP.missing"
          ##        } else if (progenyStatusTable[thisProgeny,locus] %in%
          ##                   c("MP.noMatch","M.missing","P.missing","MP.missing")) {
          ##Leave CMP, CPNotM as NA
          simpleFatherArray[locus,thisProgeny,candidate] <-
            progenyStatusTable[thisProgeny,locus]
        }
        ##  |---------------------+-------------------+---------+---------|
        ##  | ProgenyStatusTable  | SimpleFatherArray | CPNotM  | CMP     |
        ##  |---------------------+-------------------+---------+---------|
        ##  | MAO                 | Match             | NA      | alleles |
        ##  | MAO                 | noMatch           | NA      | alleles |
        ##  | NMA                 | Match             | alleles | alleles |
        ##  | NMA                 | noMatch           | alleles | alleles |
        ##  | MP.noMatch          | MP.noMatch        | NA      | NA      |
        ##  | (any code)          | C.missing         | NA      | NA      |
        ##  | M.missing           | M.missing         | NA      | NA      |
        ##  | P.missing           | P.missing         | NA      | NA      |
        ##  | MP.missing          | MP.missing        | NA      | NA      |
        ##  |---------------------+-------------------+---------+---------|
        ##
        ##  For the MAO-Match case, number of alleles is >= ploidy/2
        ##  For the MAO-noMatch case, number of alleles is < ploidy/2
        ##     Need to watch what happens when number of alleles = 0!
        ##  Empty allele sets can occur in CPNotM and CMP in the NMA cases as
        ##     well - check the behaviour...
        ##
      } ## End progeny loop
    } ## End candidate loop
  } ## End locus loop
  names(progenyStatusTable) <- paste("Locus",1:numLoci,sep="")
  names(MP.alleleTable) <- paste("Locus",1:numLoci,sep="")
  names(PNotM.alleleTable) <- paste("Locus",1:numLoci,sep="")
  ##
  ######################################################################
  ## Error check - remove once confident NA's cannot occur in
  ## simpleFatherArray!
  if (any(is.na(as.vector(simpleFatherArray)))) {
    stop("\n NA's in simpleFatherArray - please notify PolyPatEx developer!  Stopping...\n\n")
  }
  ######################################################################
  ## Correct for mothers in self-incompatible species
  if (!dioecious && !selfCompatible) {
    for (k in 1:length(progeny)) {
      is.na(CPNotM.alleleArray[,progeny[k],progenyMothers[k]]) <- TRUE
      is.na(CMP.alleleArray[,progeny[k],progenyMothers[k]]) <- TRUE
      simpleFatherArray[,progeny[k],progenyMothers[k]] <- "Mother!"
    }
  }
  for (thisProgeny in progeny) {
    for (candidate in potentialFathers) {
      ##
      totalValid <- sum(simpleFatherArray[,thisProgeny,candidate] %in%
                        c("CP.match","CP.noMatch"),na.rm=TRUE)
      ##na.rm=TRUE redundant here, since %in% does not produce NA
      ## results upon comparison with NA - but just in case this
      ## behaviour changes in future versions of R...
      totalFather <- sum(simpleFatherArray[,thisProgeny,candidate]=="CP.match",
                         na.rm = TRUE)
      fatherSummaryTable[thisProgeny,candidate] <- paste(totalFather,
                                                      totalValid,
                                                      sep=" / ")
      FLCount[thisProgeny,candidate] <- totalFather
      VLTotal[thisProgeny,candidate] <- totalValid
      ##
    } ## End candidate loop
  } ## End progeny loop
  ##
  ##The mother for each of the progeny are stored as an attribute to
  ## progenyStatusTable, for later convenience...
  attr(progenyStatusTable,"progenyMothers") <- progenyMothers
  cat("\n Done \n")
  return(list(
              progenyTables = list(
                progenyStatusTable = progenyStatusTable,
                MP.alleleTable = MP.alleleTable,
                PNotM.alleleTable = PNotM.alleleTable),
              adultTables = list(
                CPNotM.alleleArray = CPNotM.alleleArray,
                CMP.alleleArray = CMP.alleleArray,
                simpleFatherArray = simpleFatherArray,
                fatherSummaryTable = fatherSummaryTable,
                FLCount = FLCount,
                VLTotal = VLTotal)))
}
