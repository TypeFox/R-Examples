
##  Function potentialFatherIDs
##
##' Identify the potential fathers for each progeny.
##'
##' Given the output from \code{\link{genotPPE}} or
##' \code{\link{phenotPPE}}, \code{potentialFatherIDs} returns, for
##' each progeny, the IDs of candidates that are identified as
##' potential fathers.
##'
##' To decide whether a given candidate is a potential father to a
##' given progeny, \code{potentialFatherIDs} uses the quantities
##' FLCount (the number of loci at which a candidate can provide a
##' gamete compatible with the progeny) and VLTotal (the number of
##' loci at which a valid comparison was possible - \sQuote{valid}
##' loci) that are returned by \code{\link{genotPPE}} or
##' \code{\link{phenotPPE}}.
##'
##' For a candidate to be identified as a potential father of a
##' progeny, there are two criteria to be met:
##' \enumerate{
##' \item \code{VLTotal >= max(VLTMin,mismatches+1)},
##' \item \code{FLCount >= VLTotal-mismatches}.
##' }
##' Here, \code{VLTmin} and \code{mismatches} are user-specified
##' parameters. \code{VLTmin} allows the user to ensure that a
##' candidate is only considered for potential fatherhood if a
##' sufficient number of valid loci were available for comparison.
##' \code{mismatches} allows the user to specify a maximum number of
##' allowed mismatching loci between progeny and candidate, before
##' the candidate is rejected as a potential father.  Hence the user
##' may wish to relax the condition that ALL valid loci must match for
##' a candidate to be regarded as a potential father to a progeny.
##'
##' @title Identify potential fathers
##' @param dataset list: a list structure previously output from
##' \code{\link{genotPPE}} or \code{\link{phenotPPE}}.
##' @param mismatches integer: the maximum allowed number of
##' mismatching loci between candidate and progeny, before the
##' candidate is rejected as a potential father.  Default value is 0 -
##' i.e., no mismatches allowed.
##' @param VLTMin integer: the minimum number of \sQuote{valid} loci
##' (loci at which a valid progeny-candidate comparison was possible)
##' required for a candidate to be considered as a potential father.
##' Default value is 1.
##' @return A data frame, containing the columns \code{Progeny} (ID)
##' \code{Mother} (ID), \code{potentialFather} (ID or \code{None})
##' \code{FLCount} and \code{VLTotal} (the \code{FLCount} and
##' \code{VLTotal} values for the given potential father).
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
##' gPPE <- genotPPE(gData)  ## Perform the exclusion analyses
##'
##' ## Obtain IDs of potential fathers of each seedling, allowing a
##' ## single allele mismatch:
##' pFI <- potentialFatherIDs(gPPE,mismatches=1,VLTMin=2)
##'
##' ## pFC can be viewed or writted to file via, e.g. write.csv()
##'
potentialFatherIDs <- function(dataset,mismatches=0,VLTMin=1) {
  ##
  checkForValidPPEOutputObj(dataset)
  ##
  aa <- apply(with(dataset$adultTables,
                   VLTotal >= max(VLTMin,mismatches+1) &
                   (FLCount >= VLTotal-mismatches)),
              1,
              function(vv,fatherNames) {
                stripNAs(fatherNames[vv])},
              colnames(dataset$adultTables$FLCount))
  if (length(aa)==0)
    {
      warning("No potential fathers were found for any offspring")
      return(data.frame(Progeny=character(),
                 Mother=character(),
                 potentialFather=character(),
                 FLCount=numeric(),
                 VLTotal=numeric()))
    }
  aa[sapply(aa,length)==0] <- "None"
  ss <- utils::stack(aa)
  sMs <- attr(dataset$progenyTables$progenyStatusTable,
              "progenyMothers")
  aPs <- rownames(dataset$progenyTables$progenyStatusTable)
  progenyMothers <- sMs[match(ss$ind,aPs)]
  indMatrix <- with(dataset$adultTables,
                    cbind(match(ss$ind,rownames(FLCount)),
                          match(ss$values,colnames(FLCount))))
  return(data.frame(Progeny=ss$ind,
                    Mother=progenyMothers,
                    potentialFather=ss$values,
                    FLCount = dataset$adultTables$FLCount[indMatrix],
                    VLTotal = dataset$adultTables$VLTotal[indMatrix]
                    ))
}
