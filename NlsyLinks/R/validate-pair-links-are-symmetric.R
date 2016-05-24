#' @name ValidatePairLinksAreSymmetric
#' @export
#' @title Verifies that the pair relationships are symmetric.
#' 
#' @description For certain analyses, the pairs links (which can be considered a type of
#' sparse matrix) need to be symmetric. For instance, if there is a row for
#' Subjects 201 and 202 with \code{R}=0.5, there should be a second row for
#' Subjects 202 and 201 with \code{R}=0.5.
#' 
#' This validation function is useful to some types of DF methods and some
#' spatially-inspired methods.
#' 
#' 
#' @param linksPair The \code{\link{data.frame}} object that should be
#' symmetric
#' @return Returns \code{TRUE} if symmetric.  Throw an error with
#' \code{\link{stop}} if asymmetric.
#' @author Will Beasley
#' @seealso \code{\link{CreatePairLinksDoubleEntered}}
#' @keywords validation
#' @examples
#' 
#' dsSingleLinks <- data.frame(
#'   ExtendedID=c(1, 1, 1, 2), 
#'   SubjectTag_S1=c(101, 101, 102, 201), 
#'   SubjectTag_S2=c(102, 103, 103, 202), 
#'   R=c(.5, .25, .25, .5), 
#'   RelationshipPath=rep("Gen2Siblings", 4)
#' )
#' dsSingleOutcomes <- data.frame(
#'   SubjectTag=c(101, 102, 103, 201, 202), 
#'   DV1=c(11, 12, 13, 41, 42), 
#'   DV2=c(21, 22, 23, 51, 52))
#' dsDouble <- CreatePairLinksDoubleEntered(
#'   outcomeDataset=dsSingleOutcomes, 
#'   linksPairDataset=dsSingleLinks, 
#'   outcomeNames=c("DV1", "DV2"), 
#'   validateOutcomeDataset=TRUE)
#' dsDouble #Show the 8 rows in the double-entered pair links
#' summary(dsDouble) #Summarize the variables
#' 
#' ValidatePairLinksAreSymmetric(dsDouble) #Should return TRUE.
#' 
#' 
ValidatePairLinksAreSymmetric <- function( linksPair ) {
  ValidatePairLinks(linksPair)
  for( rowIndex in base::seq_along(base::nrow(linksPair)) ) {
#     tag1 <- linksPair$SubjectTag_S1[rowIndex]
#     tag2 <- linksPair$SubjectTag_S2[rowIndex]    
#     r <- linksPair$R[rowIndex]    
#     path <- linksPair$RelationshipPath[rowIndex]
    tag1 <- linksPair[rowIndex, 'SubjectTag_S1']
    tag2 <- linksPair[rowIndex, 'SubjectTag_S2']    
    r <- linksPair[rowIndex, 'R']    
    path <- linksPair[rowIndex, 'RelationshipPath']    
    #oppositeCount <- base::nrow(subset(linksPair, SubjectTag_S1==tag2 & SubjectTag_S2==tag1 & R==r & RelationshipPath==path))
    oppositeCount <- base::nrow(linksPair[linksPair$SubjectTag_S1==tag2 & linksPair$SubjectTag_S2==tag1 & linksPair$R==r & linksPair$RelationshipPath==path, ])
    if( oppositeCount != 1 ) {
      base::stop(paste0("The 'linksPair' dataset doesn't appear to be double-entered & symmetric.  The reciprocal of (SubjectTag_S1, SubjectTag_S2, R)=(",
                 tag1, ", ", tag2, ", ", r, ") was found ", oppositeCount, " time(s)."))
    }
  }  
  return( TRUE )
}
