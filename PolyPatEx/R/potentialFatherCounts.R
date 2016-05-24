
##  Function potentialFatherCounts
##
##' Count the number of potential fathers detected for each progeny.
##'
##' Given the output from \code{\link{genotPPE}} or
##' \code{\link{phenotPPE}},  \code{potentialFatherCounts} returns,
##' for each progeny, the number of candidates that are identified as
##' potential fathers.
##'
##' To decide whether a given candidate is a potential father to a
##' given progeny, \code{potentialFatherCounts} uses the quantities
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
##' @title Count potential fathers
##' @param dataset list: a list structure previously output from
##' \code{\link{genotPPE}} or \code{\link{phenotPPE}}.
##' @param mismatches integer: the maximum allowed number of
##' mismatching loci between candidate and progeny, before the
##' candidate is rejected as a potential father.
##' @param VLTMin integer: the  minimum number of \sQuote{valid} loci
##' (loci at which a valid progeny-candidate comparison was possible)
##' required for a candidate to be considered as a potential father.
##' @return A data frame, containing columns \code{Progeny} (progeny
##' id), \code{Mother} (id of the progeny's mother) and
##' \code{potentialFatherCount} (the number of potential fathers found
##' for the given progeny, given the criteria described above).
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
##' gPPE <- genotPPE(gData)  ## Perform the exclusion analyses
##'
##' ## Obtain counts of potential fathers of each seedling, allowing a
##' ## single allele mismatch:
##' pFC <- potentialFatherCounts(gPPE,mismatches=1,VLTMin=2)
##'
##' ## pFC can be viewed or written to file via, e.g. write.csv()
##'
potentialFatherCounts <- function(dataset,mismatches=0,VLTMin=1) {
  ##
  checkForValidPPEOutputObj(dataset)
  ##
  progenyMothers <- attr(dataset$progenyTables$progenyStatusTable,
                          "progenyMothers")
  pFC <- apply(with(dataset$adultTables,
                    VLTotal >= max(VLTMin,mismatches+1) & ##Note the constraint...
                    (FLCount >= VLTotal-mismatches)),
               1,
               function(vv){sum(vv,na.rm=TRUE)})
  return(data.frame(Progeny=rownames(dataset$progenyTables$progenyStatusTable),
                    Mother=progenyMothers,
                    potentialFatherCount=pFC))
}
