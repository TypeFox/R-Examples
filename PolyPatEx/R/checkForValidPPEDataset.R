
##  Function checkForValidPPEDataset
##
## Check for a valid PPE dataset.
##
## Allele datasets \emph{must} be checked and preprocessed by
## \code{\link{preprocessData}} prior to invocation of any of the
## PolyPatEx data analysis functions.  This is achieved either by
## loading the allele data from a CSV file using
## \code{\link{inputData}} (which calls \code{\link{preprocessData}}
## automatically), or by running \code{\link{preprocessData}}
## explicitly on the allele data frame prior to invoking any
## analysis functions.
##
## Function \code{checkForValidPPEDataset} checks that
## \code{\link{preprocessData}} has indeed been run on a dataset, by
## checking that the required attributes (\code{numLoci},
## \code{ploidy}, etc) have been attached to the data frame by
## \code{\link{preprocessData}}.
##
## Of course, a knowledgeable R user could easily attach the required
## attributes to the data frame themselves \emph{without} invoking
## \code{\link{preprocessData}}, but there's not much that I can do
## about that\ldots
##
## \code{checkForValidPPEDataset} is invoked by PolyPatEx's various
## data analysis functions to protect against attempts to analyse
## datasets that have not been properly prepared.
##
## @title Check for a valid PolyPatEx allele dataset
## @param dataset a data frame containing allele data
## @author Alexander Zwart (alec.zwart at csiro.au)
## @examples
## \dontrun{
##
## checkForValidPPEDataset(adata)
##
## }
##
checkForValidPPEDataset <- function (dataset) {
  ##
  requiredAttributes <- c("numLoci", "ploidy", "dataType", "dioecious",
                          "selfCompatible")
  ## Note: This routine does not check for mothersOnly when
  ## dioecious=TRUE...
  attributesPresent <- names(attributes(dataset))
  errorString <- paste("\n Your allele dataset does not appear to have been preprocessed",
                       "\n  by the preprocessData() function.  Please either load your",
                       "\n  allele dataset from a CSV file using the inputData() function",
                       "\n  (which calls preprocessData() automatically), or explicitely",
                       "\n  run preprocessData() on your allele dataset prior to using",
                       "\n  other PolyPatEx analysis functions.\n\n",
                       collapse=" ")
  if (is.null(attributesPresent) ||
      !all(requiredAttributes %in% attributesPresent)) {
    stop(errorString)
  }
  return(invisible(NULL))
}
