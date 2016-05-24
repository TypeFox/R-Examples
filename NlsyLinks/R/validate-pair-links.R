
##' @name ValidatePairLinks
##' @export
##' 
##' @title Validates the schema of a links for pairs of relatives
##' 
##' @description A helper function that verifies the linking dataset contains (A) the
##' essential columns exist, and (B) at least one row.  It is called by
##' \code{CreatePairLinks}.
##' 
##' Typical use of \pkg{NlsyLinks} will not require this function, since a
##' valid paired links are supplied for each supported sample (ie,
##' \code{\link{Links79Pair}}).
##' 
##' The \pkg{NlsyLinks} uses several types of linking schemas.  This function
##' validates the type where each relative subject has their own row.
##' 
##' The following four columns must be present: (1) \code{Subect1Tag}, (2)
##' \code{Subect2Tag}, (3) \code{R}, and (4) \code{MultipleBirth}.  They must
##' have a \code{numeric} mode/datatype.
##' 
##' @param linksPair The \code{data.frame} to validate.
##' @return Returns \code{TRUE} if the validation passes. Returns an error (and
##' associated descriptive message) if it false.
##' @seealso \code{\link{Links79Pair}}, \code{\link{Links79PairExpanded}},
##' @keywords validation
##' @examples
##' 
##' dsSingleLinks <- data.frame(
##'   ExtendedID=c(1, 1, 1, 2), 
##'   SubjectTag_S1=c(101, 101, 102, 201), 
##'   SubjectTag_S2=c(102, 103, 103, 202), 
##'   R=c(.5, .25, .25, .5), 
##'   RelationshipPath=rep("Gen2Siblings", 4)
##' )
##' dsSingleOutcomes <- data.frame(
##'   SubjectTag=c(101, 102, 103, 201, 202), 
##'   DV1=c(11, 12, 13, 41, 42), 
##'   DV2=c(21, 22, 23, 51, 52))
##' dsDouble <- CreatePairLinksDoubleEntered(
##'   outcomeDataset=dsSingleOutcomes, 
##'   linksPairDataset=dsSingleLinks, 
##'   outcomeNames=c("DV1", "DV2"), 
##'   validateOutcomeDataset=TRUE)
##' dsDouble #Show the 8 rows in the double-entered pair links
##' summary(dsDouble) #Summarize the variables
##' 
##' ValidatePairLinksAreSymmetric(dsDouble) #Should return TRUE.
##' 

ValidatePairLinks <- function( linksPair ) {
  if(!base::nrow(linksPair) > 0 ) stop("The linksPair data frame should have at least one row, but does not.")
  
  columnNames <- base::colnames(linksPair)
  if( !base::any(columnNames=="SubjectTag_S1") ) base::stop("The column 'SubjectTag_S1' should exist in the linksPair file, but does not.")
  if( !base::any(columnNames=="SubjectTag_S2") ) base::stop("The column 'SubjectTag_S2' should exist in the linksPair file, but does not.")
  if( !base::any(columnNames=="R") ) base::stop("The column 'R' should exist in the linksPair file, but does not.")
  # if( !base::any(columnNames=="MultipleBirth") ) stop("The column 'MultipleBirth' should exist in the linksPair file, but does not.")
  
  if( base::mode(linksPair$SubjectTag_S1) != 'numeric' ) base::stop("The column 'SubjectTag_S1' should have a 'numeric' mode, but does not.")
  if( base::mode(linksPair$SubjectTag_S2) != 'numeric' ) base::stop("The column 'SubjectTag_S2' should have a 'numeric' mode, but does not.")
  if( base::mode(linksPair$R) != 'numeric' ) base::stop("The column 'R' should have a 'numeric' mode, but does not.")
  # if( base::mode(linksPair$MultipleBirth) != 'numeric' ) base::stop("The column 'MultipleBirth' should have a 'numeric' mode, but does not.")
  
  return( TRUE )
}
