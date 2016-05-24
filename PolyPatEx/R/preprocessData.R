
##  Function preprocessData
##
##' Check and preprocess the input allele data frame prior to
##' subsequent analysis.
##'
##' If \code{\link{inputData}} is used to load the allele data set
##' into R, then \code{preprocessData} will be called automatically on
##' the data frame before it is returned by \code{\link{inputData}}.
##' However, if the user loads their data into R by some means other
##' than \code{\link{inputData}}, then \code{preprocessData} MUST be
##' called on the data frame prior to using any other PolyPatEx
##' functions to analyse the allele data---\code{preprocessData}
##' performs a series of critical checks and preprocessing steps on
##' the data frame, without which other analysis functions in
##' PolyPatEx will fail.
##'
##' Note that \code{\link{inputData}} strips leading or trailing
##' spaces (whitespace) from each entry in the allele dataset as it is
##' read in.  If you load your data by a means other than
##' \code{\link{inputData}}, you should ensure that you perform this
##' step yourself, as \code{preprocessData} will not carry out this
##' necessary step.
##'
##' Note also that you should not use spaces in any of your allele
##' codes - PolyPatEx functions use spaces to separate allele codes as
##' they process the data - if allele codes already contains spaces,
##' errors will occur in this processing. If you need a separator, I
##' recommend using either \sQuote{code{.}} (a period) or
##' \sQuote{code{_}} (an underscore) rather than a space.
##'
##' \code{preprocessData} first performs a number of simple checks on
##' the format and validity of the data set.  These checks look for
##' the presence of certain required columns and correct naming and
##' content of these columns.  \code{preprocessData} will usually stop
##' with an error message should the data fail these basic checks.
##' Correct the indicated problem in the CSV file or R allele data
##' frame, then call \code{\link{inputData}} or \code{preprocessData}
##' again as appropriate.  If you use a spreadsheet to edit the CSV
##' file, don't forget that you may also need to call
##' \code{\link{fixCSV}} on the CSV file, prior to calling
##' \code{\link{inputData}} again.
##'
##' If the data is \sQuote{genotypic} data PolyPatEx requires that all
##' \eqn{p} alleles must be present in each allele set, where \eqn{p}
##' is the species' ploidy.  If an allele set contains fewer than
##' \eqn{p} alleles, then it is reset to contain no alleles and is
##' subsequently ignored by other PolyPatEx functions.  ID and locus
##' information is printed to the R terminal, to help the user locate
##' these cases in their original dataset.
##'
##' Further checks look for mismatches between progeny and their
##' mothers' allele sets at each locus---these are situations where a
##' progeny's allele set could not have arisen from any gamete that
##' the mother can provide.  When no more than \code{matMismatches}
##' mismatching loci occurs in a mother-progeny pair, the offending
##' allele sets in the progeny are reset to contain no alleles (we
##' term these \sQuote{missing} allele sets).  When mismatches occur
##' in more than \code{matMismatches} loci, the progeny is removed
##' entirely from the dataset.  Information is printed to the R
##' terminal to assist the user in identifying the affected
##' individuals and loci---in particular, note that removal of several
##' (or all) of a single mother's progeny may indicate an error in the
##' mother's allele data, rather than in her progeny.
##'
##' After the mother/progeny mismatch check above, a subsequent check
##' removes individuals from the dataset that have fewer than
##' \code{lociMin} non-missing allele sets remaining.  The default
##' value for \code{lociMin} is 1---an individual must have at least
##' one non-missing locus to remain in the dataset.  If any mothers
##' are removed from the dataset at this stage, all of her progeny are
##' removed also.  Again, information about these removals is printed
##' to the R terminal.
##'
##' Note that in the data frame that is returned by
##' \code{preprocessData}, the alleles in each allele set (i.e,
##' corresponding to each locus) will be sorted into alphanumeric
##' order---this sort order is necessary for the correct operation of
##' other PolyPatEx routines, and should not be altered.
##'
##' PolyPatEx needs to know the characteristics of the dataset being
##' analysed.  These are specified in the \code{\link{inputData}} or
##' \code{preprocessData} calls and are invisibly attached to the
##' allele data frame that is returned, for use by other PolyPatEx
##' functions. The required characteristics are:
##' \itemize{
##'  \item \code{numLoci}: the number of loci in the dataset
##'  \item \code{ploidy}: the ploidy (\eqn{p}) of the species
##'        (currently allowed to be 4, 6, or 8.  \code{ploidy} can
##'        also be 2, provided \code{dataType="genotype"})
##'  \item \code{dataType}: whether the data is genotypic (all \eqn{p}
##'        alleles at each locus are observed) or phenotypic (only
##'        the distinct allele states at a locus are observed -
##'        alleles that appear more than once in the genotype of
##'        a locus only appear once in the phenotype)
##'  \item \code{dioecious}: whether the species is dioecious
##'        or monoecious
##'  \item \code{selfCompatible}: whether a monoecious species is
##'        self compatible (i.e., whether an individual can
##'        fertilise itself)
##'  \item \code{mothersOnly}: whether a dioecious dataset should
##'        retain only adult females that are mothers of progeny
##'        in the dataset.
##' }
##'
##' @title Check and preprocess an allele dataset
##' @param adata data frame: an allele dataset.
##' @param numLoci integer: the number of loci in the allele dataset.
##' @param ploidy integer: the species' ploidy, one of \code{2},
##' \code{4}, \code{6}, or \code{8}.
##' @param dataType character: either \code{"genotype"} or
##' \code{"phenotype"}.
##' @param dioecious logical: is the species dioecious or monoecious?
##' @param selfCompatible logical: In monoecious species
##' (\code{dioecious=FALSE}), can individuals self-fertilise?  When
##' \code{dioecious=FALSE}, this argument may be left at its
##' default value of NULL - it will be set to \code{FALSE} by
##' \code{preprocessData}.
##' @param mothersOnly logical: in dioecious species, should females
##' without progeny present be removed from the dataset?  If
##' \code{dioecious=TRUE}, then \code{mothersOnly} must be set to
##' either \code{TRUE} or \code{FALSE}.  If \code{dioecious=FALSE},
##' argument \code{mothersOnly} should be left at its default value of
##' \code{NULL}.
##' @param lociMin integer: the minimum number of loci in a individual
##' that must have alleles present for the individual (and its
##' progeny, if any) to be retained in the dataset (default 1).
##' @param matMismatches an integer between 0 and \code{numLoci}-1,
##' being the maximum number of mismatching alleles between mother and
##' offspring that are allowed before the offspring is removed from the
##' dataset.  The default value is 0.  If an offspring has fewer
##' than \code{matMismatches} loci that mismatch with its mother, the
##' offending loci are set to contain no alleles.
##' @return  A data frame, containing the checked and pre-processed
##' allele data, ready for further analysis by other \pkg{PolyPatEx}
##' functions.  All columns in the returned data frame will be of mode
##' \code{character}.
##' @author Alexander Zwart (alec.zwart at csiro.au)
##' @export
##' @examples
##'
##' ## If using inputData to input the allele dataset from CSV file,
##' ## preprocessData() is applied automatically before the dataset is
##' ## returned by inputData().
##'
##' ## Otherwise, if the allele dataset is created or loaded into R
##' ## by other means, such preprocessData() must be applied before
##' ## other PolyPatEx analysis routines are applied:
##'
##' ## Using the example dataset 'GF_Phenotype':
##' data(GF_Phenotype)
##'
##' ## Since we did not load this dataset using inputData(), we must
##' ## first process it with preprocessData() before doing anything
##' ## else:
##' pData <- preprocessData(GF_Phenotype,
##'                         numLoci=7,
##'                         ploidy=6,
##'                         dataType="phenotype",
##'                         dioecious=FALSE,
##'                         selfCompatible=FALSE)
##'
##' head(pData)  ## Checked and Cleaned version of GF_Phenotype
##'
##' ## pData is now ready to be passed to other PolyPatEx analysis
##' ## functions...
##'
preprocessData <- function(adata,
                           numLoci,
                           ploidy,
                           dataType,
                           dioecious,
                           selfCompatible=NULL,
                           mothersOnly=NULL,
                           lociMin=1,
                           matMismatches=0) {
  ##
  if (!(ploidy %in% c(2,4,6,8))) {
    stop('\n Please specify the ploidy= argument as 2, 4, 6 or 8\n\n')
  }
  if (!(dioecious %in% c(TRUE, FALSE))) {
    stop('\n Please specify the dioecious= argument as TRUE or FALSE\n\n')
  }
  if (dioecious) {
    if (is.null(selfCompatible) || selfCompatible != FALSE)
        cat('\nSpecies is dioecious - setting selfCompatible = FALSE\n\n')
    selfCompatible <- FALSE
    if (is.null(mothersOnly)) {
      stop('\n Argument dioecious=TRUE - please set mothersOnly to TRUE or FALSE\n\n')
    }
  } else if (is.null(selfCompatible) ||
             !(selfCompatible %in% c(TRUE, FALSE))) {
    stop('\n Please specify the selfCompatible= argument (TRUE or FALSE)\n\n')
  }
  if (!(dataType %in% c("genotype","phenotype"))) {
    stop('\n Please specify the dataType= argument as "genotype" or "phenotype"\n\n')
  }
  if (ploidy==2 && dataType=="phenotype")
  {
    stop("\n Only genotype data can be processed when ploidy = 2\n")
  }
  ##
  if (!(matMismatches %in% 0:(numLoci-1)))
    {
      stop("Argument matMismatches= must be an integer between 0 and numLoci-1 (inclusive)")
    }
  ##
  nonCharColumns <- sapply(adata,"class") != "character"
  if (any(nonCharColumns)) {
    cat("\n\nThe following columns have been converted to character vectors:\n\n")
    print(names(adata)[nonCharColumns])
    cat("\n")
    adata <- data.frame(lapply(adata,
                                function(vv){as.character(vv)}),
                         stringsAsFactors = FALSE)
  }
  ## Check - blank cells with preprocessData should be represented by
  ## NA, not ""
  if (any(sapply(adata,function(vv) any(vv[!is.na(vv)]==""))))
  {
    msg <- paste('\n Blank cells in the data frame passed to preprocessData()\n',
                  ' should contain NA, rather than "".\n',collapse="")
    stop(msg)
  }
  ## Names of the important columns
  if (any(names(adata)[1:3]!= c("id","popn","mother"))) {
    stop("\n The first three columns of your data spreadsheet should
          be called 'id', 'popn', and 'mother', in that order.
          Note the capitalisation used...\n\n")
  }
  if (dioecious && names(adata)[4] != "gender") {
    stop("\n You have specified dioecious=TRUE: The fourth column of
          your dataset should be 'gender'...\n\n")
  }
  ##
  ## Does the number of columns match the specified number of loci and
  ## ploidy (& dioecious status)?
  if (ncol(adata)-(3+dioecious) != numLoci*ploidy) {
    messg <- paste("\n The number of columns in the dataset is inconsistent with \n",
                   "    the specified number of loci, and ploidy.  The number of \n",
                   "    columns should be",3+dioecious,
                   " + (number of loci) * ploidy\n\n")
    stop(messg)
  }
  ##
  ## All entries uniquely identified
  if (any(duplicated(adata$id))) {
    stop("\n The 'id' column in the dataset should UNIQUELY identify each
          row of the dataset.\n\n")
  }
  ##
  ## If all id's are unique, use 'em as rownames...
  rownames(adata) <- adata$id
  ##
  ## Comment on the number of populations:
  if (length(unique(adata$popn)) > 1) {
    cat("\n There appears to be more than one population in this dataset.\n")
    cat("\n Note that PolyPatEx functions ignore distinctions between populations...\n")
  }
  ##
  ## No progeny without a mother.
  progenyMothers <- with(adata,unique(stripNAs(mother)))
  if(any(!(progenyMothers %in% adata$id))) {
    stop("\n Every progeny in the dataset should have a mother present
          in the dataset\n\n")
  }
  ##
  ## No more than one generation of progeny...
  progenyIDs <- with(adata,id[!is.na(mother)])
  if(any((progenyMothers %in% progenyIDs))) {
    stop("\n At least one progeny has another progeny as its mother\n\n")
  }
  ##
  ## Dioecious checks
  if (dioecious) {
    if (any(!(stripNAs(unique(adata$gender)) %in% c("M","F")))) {
      stop("\n Non-blank entries in the gender column of a dioecious
            dataset should be 'F' (Female) or 'M' (Male) only.\n\n")
    }
    ## All mothers should be recorded as gender "F" (Female)
    motherGenders <- with(adata,gender[id %in% progenyMothers])
    if (any(is.na(motherGenders)) || any(motherGenders != "F")) {
      stop("\n All mothers should have gender 'F' (for Female)\n\n")
    }
    if (any(!is.na(adata$gender[!is.na(adata$mother)]))) {
      stop("\n Progeny genders should not be specified.\n\n")
    }
    if (any(is.na(adata$gender[is.na(adata$mother)]))) {
      stop("\n All adults should have specified gender: 'F' (Female) or 'M' (Male)\n\n")
    }
  }
  ##
  ## Matrix to store the number of alleles detected at each locus
  alleleCounts <- matrix(0, nrow=dim(adata)[1],
                         ncol=numLoci,
                         dimnames=list(adata$id,
                           paste("Locus",1:numLoci,sep="")))
  if (dataType=="phenotype") {
    uniqueAlleleCounts <- matrix(0, nrow=dim(adata)[1],
                                 ncol=numLoci,
                                 dimnames=list(adata$id,
                                   paste("Locus",1:numLoci,sep="")))
  }
  ##
  ## Sort each allele set into (alphanumeric) order, and count the
  ## number of alleles present in each set - store in alleleCounts
  for (locus in 1:numLoci) {
    locusRange <- (3+dioecious) + (locus-1)*ploidy + 1:ploidy
    ## Note: apply() returns the _transpose_ of the desired matrix
    ## below - hence the use of t()...
    adata[ , locusRange] <- t(apply(adata[ , locusRange],
                                     1,
                                     sort,
                                     na.last = TRUE))
    ## Store the numbers of alleles detected at each locus
    alleleCounts[,locus] <- apply(adata[ , locusRange],1,
                                      function(thisAlleleVector) {
                                        sum(!is.na(thisAlleleVector))
                                      })
    if (dataType=="phenotype") {
      uniqueAlleleCounts[,locus] <- apply(adata[ , locusRange],1,
                                              function(thisAlleleVector) {
                                                sum(!is.na(unique(thisAlleleVector)))
                                              })
    }
  }
  ##
  ## If genotypic data, then look for alleleSets with
  ## 0 < num(alleles) < ploidy, and convert these to all missing.
  if (dataType=="genotype") {
    for (locus in 1:numLoci) {
      locusRange <- (3+dioecious) + (locus-1)*ploidy + 1:ploidy
      dodgyIndividuals <- rownames(alleleCounts)[alleleCounts[,locus] > 0 &
                                         alleleCounts[,locus] < ploidy]
      if(length(dodgyIndividuals) > 0) {
        cat("\n\n Note: Allele sets with SOME missing alleles have been found at Locus",locus,"\n")
        cat(" These allele sets will be reset to have NO alleles...\n")
        cat(" The affected IDs are:\n")
        print(dodgyIndividuals)
        is.na(adata[dodgyIndividuals,locusRange]) <- TRUE
        alleleCounts[dodgyIndividuals,locus] <- 0
      }
    }
  }
  ##
  ## If allelic phenotype data, check that all alleles in each
  ## individual at each locus are distinct - stop if there are any
  ## problems
  if (dataType=="phenotype") {
    inds <- alleleCounts > uniqueAlleleCounts
    if(any(inds)) {
      errIDs <- rownames(alleleCounts)[row(alleleCounts)[inds]]
      errLoci <- colnames(alleleCounts)[col(alleleCounts)[inds]]
      errCombos <- data.frame(id=errIDs, Locus=errLoci)
      errCombos <- errCombos[sort(errCombos$id),]
      rownames(errCombos) <- NULL
      cat("\n Note: The following id-locus combos have repeated alleles:\n\n")
      print(errCombos)
      cat("\n")
      stop("\n Please remove allele repeats from allelic phenotype data\n\n")
    }
  }
  ##
  ## Append dataset descriptors as attributes (needed for
  ## removeMismatches)
  attr(adata,"numLoci") <- numLoci
  attr(adata,"ploidy") <- ploidy
  attr(adata,"dataType") <- dataType
  attr(adata,"dioecious") <- dioecious
  attr(adata,"selfCompatible") <- selfCompatible
  attr(adata,"mothersOnly") <- mothersOnly
  ## Deal with any mismatching allele sets between mothers and their
  ## progeny
  adata <- removeMismatches(adata,matMismatches)
  ##
  ## For each individual, which allele sets remain non-empty?
  allelesPresent <- matrix(0, nrow=dim(adata)[1],
                           ncol=numLoci,
                           dimnames=list(adata$id,
                             paste("Locus",1:numLoci,sep="")))
  for (locus in 1:numLoci) {
    locusRange <- (3+dioecious) + (locus-1)*ploidy + 1:ploidy
    allelesPresent[,locus] <- apply(adata[ , locusRange],
                                        1,
                                        function(vv) {
                                          any(!is.na(vv))
                                        })
  }
  ## We need to re-create progenyMothers, since removeMismatches() may
  ## have removed some mothers from the dataset. [Actually I may be
  ## wrong about that - removeMismatches() does not remove any
  ## mothers, only progeny. So the next line can't hurt, but may be
  ## redundant!]
  progenyMothers <- with(adata,unique(stripNAs(mother)))
  ##
  ## Check for individuals where fewer than lociMin loci are available
  ## (i.e., are non-missing). These are 'dudIndividuals'.
  dudIndividuals <- rowSums(allelesPresent) < lociMin
  ## The set of dud Individuals in general may include mothers, other
  ## (non-mother) adults, and progeny.  Dud progeny and "other
  ## adults" can simply be removed, but if any mothers are duds, then
  ## their progeny become duds as well.
  dudMothers <- (adata$id %in% progenyMothers) & (dudIndividuals)
  ## Flag all progeny of dud mothers as dud individuals
  if (any(dudMothers)) {
    for (dudMother in adata$id[which(dudMothers)]) {
      dudIndividuals[adata$mother==dudMother] <- TRUE
    }
  }
  ## Print message, and remove dud individuals from the dataset.
  if (any(dudIndividuals)) {
    cat("\n Individuals identified as having insufficient non-missing loci have
   been removed from this dataset, along with any progeny of removed
   mothers.  Some of these removals may be a result of the previous
   removal of mismatching allele sets (or, in genotypic data, allele
   sets with insufficent alleles)\n")
    cat("\n Individuals removed correspond to id's :\n")
    cat("  ",adata$id[dudIndividuals],"\n")
    adata <- adata[!dudIndividuals,]
  }
  ##
  ## If species is dioecious and mothersOnly=TRUE, remove any mothers
  ## whose progeny have all been removed.
  if(dioecious && mothersOnly) {
    femalesLeft <- stripNAs(with(adata,id[gender=="F"]))
    mothersLeft <- with(adata,unique(stripNAs(mother)))
    if(!all(femalesLeft %in% mothersLeft)) {
      femalesToBeRemoved <- femalesLeft[!(femalesLeft %in% mothersLeft)]
      adata <- adata[!(adata$id %in% femalesToBeRemoved),]
      cat("\n\n mothersOnly = TRUE : Adult females without progeny will be")
      cat("\n removed from the dataset.  Note that some of these may have")
      cat("\n arisen from the removal of progeny in previous steps.")
      cat("\n  The females removed are:\n\n")
      cat(femalesToBeRemoved,"\n\n")
    }
  }  # Otherwise non-mother females remain in the dataset
  ##
  ## Mother/progeny mismatches could result in no progeny (or their
  ## mothers) left in the dataset - check, and abort if this is the
  ## case.
  if (length(adata$id[!is.na(adata$mother)])==0) {
    msg <- paste0("\nThere appear to be no progeny in this dataset.  Perhaps",
                  "\n the data checking & cleaning above has removed all",
                  "\n progeny and their mothers from the dataset?  Otherwise,",
                  "\n check the format of your data file for possible errors.\n")
    stop(msg)
  }
  cat("\n Done \n\n")
  return(adata)
}


