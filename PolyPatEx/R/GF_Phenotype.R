
##  Dataset GF_Phenotype
##
##' Example \sQuote{phenotype} allele dataset - a monoecious non-selfing
##' hexaploid species, seven loci observed.
##'
##' The dataset is available in two forms - as a compressed data file
##' which can be loaded easily into R via the R \code{\link{data}}
##' function, i.e., \code{data(FR_Genotype)}, and as a CSV
##' (Comma-Separated-Value, a plain text format) file, to provide an
##' example of the required CSV format.
##'
##' Note that a technicality of R's package building process requires
##' the use of \code{data} to load the data in the reference help
##' examples, whereas the user would generally invoke the
##' \code{\link{inputData}} function to load their own data from file.
##' An example of the latter is demonstrated in the example section on
##' this page, but is not run.
##'
##' @name GF_Phenotype
##' @title Example phenotype allele dataset
##' @docType data
##' @author Alexander Zwart (alec.zwart at csiro.au)
##' @keywords data
##' @examples
##'
##' \dontrun{
##' ## To locate this dataset in your filesystem, use:
##'
##' pDataFile <- system.file("extdata/GF_Phenotype.csv",
##'                          package="PolyPatEx")
##' print(pDataFile)
##'
##' ## To load this file using PolyPatEx's 'inputData' function use:
##'
##' pData <- inputData(pDataFile,
##'                   numLoci=7,
##'                   ploidy=6,
##'                   dataType="phenotype",
##'                   dioecious=FALSE,
##'                   selfCompatible=FALSE)
##'
##' ## pData now contains the checked and preprocessed allele dataset,
##' ## ready to be passed to other PolyPatEx analysis functions.
##'
##'}
##'
NULL
