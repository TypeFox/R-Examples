
##  Function inputData
##
##' Read in an allele dataset from file, and return a checked and
##' preprocessed data frame.
##'
##' \code{inputData} reads in an allele dataset from the specified
##' file, then calls \code{\link{preprocessData}} to perform a series
##' of data format checks and preprocessing steps before returning the
##' checked and preprocessed dataset as an R data frame.  The
##' reference information for \code{\link{preprocessData}} contains
##' further information on the checks and preprocessing - it is
##' strongly recommended you read that information in addition to the
##' information below.
##'
##' The use of \code{inputData} is optional, if you wish to create or
##' load the allele dataset into R by other means. However, it is then
##' necessary to call \code{\link{preprocessData}} on the data frame
##' prior to using any other analysis functions in this package.
##' Similarly, if you decide to change or manipulate the data frame
##' contents within R, you should call \code{\link{preprocessData}}
##' again on the data frame prior to using any of the \pkg{PolyPatEx}
##' analysis functions.  See the help for \code{\link{preprocessData}}
##' for further details.
##'
##' Note that \code{inputData} strips leading or trailing spaces
##' (whitespace) from each entry in the allele dataset as it is read
##' in.  If you load your data by a means other than \code{inputData},
##' you should ensure that you perform this step yourself, as
##' \code{\link{preprocessData}} will not carry out this necessary
##' step.
##'
##' Note also that you should not use spaces in any of your allele
##' codes - PolyPatEx functions use spaces to separate allele codes as
##' they process the data - if allele codes already contains spaces,
##' errors will occur in this processing. If you need a separator, I
##' recommend using either \sQuote{code{.}} (a period) or
##' \sQuote{code{_}} (an underscore) rather than a space.
##'
##' Neither \code{inputData} (nor \code{\link{preprocessData}}) will
##' alter the CSV file from which the data is loaded - they merely
##' return a checked and preprocessed version of your allele dataset
##' (in the form of an R data frame) within the R environment, ready
##' for use by other \pkg{PolyPatEx} functions.
##'
##' To load the allele dataset into R, \code{inputData} calls R's
##' \code{\link{read.csv}} function with certain arguments specified.
##' These arguments make \code{\link{read.csv}} more stringent about
##' the precise format of the input datafile, requiring in particular
##' that each row of the CSV-formatted data file contain the correct
##' number of commas.  This is not always guaranteed when the CSV file
##' has been exported from spreadsheet software.  Should you get
##' \sQuote{Error in scan} messages complaining about the number of
##' elements in a line of the input file, consider calling
##' \code{\link{fixCSV}} on the data file, before calling
##' \code{inputData} again.  \code{\link{fixCSV}} attempts to find and
##' correct such errors in a CSV file - see the help for this
##' function.  Note that if you specify the \code{skip} parameter in
##' a call to \code{\link{fixCSV}}, you should use the same value for
##' this parameter in \code{inputData} to avoid an error.
##'
##' The various \pkg{PolyPatEx} functions need to know the characteristics
##' of the dataset being analysed - these are specified in the
##' \code{inputData} or \code{\link{preprocessData}} calls and are
##' invisibly attached to the allele data frame that is returned, for
##' use by other \pkg{PolyPatEx} functions. The required characteristics
##' are:
##'
##' \itemize{
##' \item \code{numLoci}: the number of loci in the dataset.
##' \item \code{ploidy}: the ploidy \eqn{p} of the species (currently
##' allowed to be 4, 6, or 8.  \code{ploidy} can also be 2, provided
##' \code{dataType="genotype"}).
##' \item \code{dataType}: whether the data is genotypic (all \eqn{p}
##' alleles at each locus are observed) or phenotypic (only the
##' distinct allele states at a locus are observed - alleles that
##' appear more than once in the genotype of a locus only appear once
##' in the phenotype).
##' \item \code{dioecious}: whether the species is dioecious or
##' monoecious.
##' \item \code{selfCompatible} whether a monoecious species is self
##' compatible (i.e., whether an individual can fertilise itself).
##' \item \code{mothersOnly}: whether a dioecious dataset should
##' retain only adult females that are mothers of progeny in the
##' dataset.  If \code{dioecious=TRUE}, then \code{mothersOnly} must
##' be set to either \code{TRUE} or \code{FALSE}.
##' }
##'
##' @title Read in, check and preprocess the allele dataset
##' @param file character: the name of the allele data file.
##' @param numLoci integer: the number of loci in the allele dataset.
##' @param ploidy integer: the species' ploidy, one of \code{2},
##' \code{4}, \code{6}, or \code{8}.
##' @param dataType character: either \code{"genotype"} or
##' \code{"phenotype"}.
##' @param dioecious logical: is the species dioecious (\code{TRUE})
##' or monoecious (\code{FALSE})?
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
##' @param lociMin integer: the minimum number of loci in an
##' individual that must have alleles present for the individual (and
##' its progeny, if the individual is a mother) to be retained in the
##' dataset.  See the help for \code{\link{preprocessData}} for more
##' on this parameter.
##' @param matMismatches an integer between 0 and \code{numLoci}-1,
##' being the maximum number of mismatching loci between mother and
##' offspring that are allowed before the offspring is removed from
##' the dataset.  The default value is 0.  If an offspring has fewer
##' than \code{matMismatches} loci that mismatch with its mother, the
##' offending loci are set to contain no alleles.
##' @param skip integer: the number of lines in the CSV to skip
##' before the header row of the table.
##' @return A data frame, containing the checked and pre-processed
##' allele data, ready for further analysis by other \pkg{PolyPatEx}
##' functions.  All columns in the output data frame will be of mode
##' \code{character}.
##' @author Alexander Zwart (alec.zwart at csiro.au)
##' @importFrom utils read.csv
##' @export
##' @examples
##'
##' \dontrun{
##'
##' ## Obtain path to the example genotype data file
##' ## 'FR_Genotype.csv'
##' gDataFile <- system.file("extdata/FR_Genotype.csv",
##'                          package="PolyPatEx")
##' print(gDataFile)
##'
##' gData <- inputData(gDataFile,
##'                    numLoci=7,
##'                    ploidy=4,
##'                    dataType="genotype",
##'                    dioecious=TRUE,
##'                    mothersOnly=TRUE)
##'
##' ## ...or use 'mothersOnly=FALSE' if you wish to retain
##' ## non-maternal females in the dataset.
##'
##' ## gData now contains the checked and preprocessed allele dataset,
##' ## ready to be passed to other PolyPatEx analysis functions.
##'
##' ## In your own workflow, you would typically specify the path to
##' ## your allele dataset directly - e.g. if the dataset
##' ## myAlleleData.csv is on the Data subdirectory of the current R
##' ## working directory (see R function setwd()), then:
##' ##
##' ## gData <- inputData("Data/myAlleleData.csv",
##' ##                    numLoci= etc etc etc...,
##'
##' }
##'
inputData <- function(file,
                      numLoci,
                      ploidy,
                      dataType,
                      dioecious,
                      selfCompatible=NULL,
                      mothersOnly=NULL,
                      lociMin=1,
                      matMismatches=0,
                      skip=0) {
  ##
  if (!(ploidy %in% c(2,4,6,8))) {
    stop('\n Please specify the ploidy= argument as 2, 4, 6 or 8\n\n')
  }
  if (!(dioecious %in% c(TRUE, FALSE))) {
    stop('\n Please specify the dioecious= argument as TRUE or FALSE\n\n')
  }
  if (dioecious) {
    if (is.null(selfCompatible) || selfCompatible != FALSE)
        cat('\n Species is dioecious - setting selfCompatible = FALSE\n\n')
    selfCompatible <- FALSE
    if (is.null(mothersOnly)) {
      stop('\n Argument dioecious= is TRUE - please set mothersOnly to TRUE or FALSE\n\n')
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
  adata <- utils::read.csv(file,
                     header=TRUE,
                     colClasses = "character",
                     strip.white = TRUE,
                     na.strings=c("*","","NA"),
                     fill=FALSE,
                     skip=skip)
  ##
  return(preprocessData(adata = adata,
                        numLoci = numLoci,
                        ploidy = ploidy,
                        dataType = dataType,
                        dioecious = dioecious,
                        selfCompatible = selfCompatible,
                        mothersOnly = mothersOnly,
                        lociMin = lociMin,
                        matMismatches = matMismatches))
}



