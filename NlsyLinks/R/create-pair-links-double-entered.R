#' @name CreatePairLinks
#' @aliases CreatePairLinksSingleEntered CreatePairLinksDoubleEntered
#' CreatePairLinksDoubleEnteredWithNoOutcomes
#' @export CreatePairLinksSingleEntered CreatePairLinksDoubleEntered CreatePairLinksDoubleEnteredWithNoOutcomes
#' 
#' @title Creates a pairs linking file.
#' @description Creates a linking file for BG designs using this file structure (e.g., DF analysis, other ACE modeling).
#' A DF analysis requires a double-entered file that contains the \code{R}
#' value for the pair, and their two outcome variable values.
#' 
#' \code{CreatePairLinksDoubleEnteredWithNoOutcomes} is intended to be a
#' primarily a helper function for \code{\link{CreateSpatialNeighbours}}.
#' 
#' @usage CreatePairLinksDoubleEntered(outcomeDataset, linksPairDataset, outcomeNames, 
#'    linksNames = c("ExtendedID", "R", "RelationshipPath"), validateOutcomeDataset = TRUE, 
#'    subject1Qualifier = "_S1", subject2Qualifier = "_S2")
#' 
#'  CreatePairLinksSingleEntered(outcomeDataset, linksPairDataset, outcomeNames, 
#'     linksNames = c("ExtendedID", "R", "RelationshipPath"), validateOutcomeDataset = TRUE, 
#'     subject1Qualifier = "_S1", subject2Qualifier = "_S2")
#'  
#'  CreatePairLinksDoubleEnteredWithNoOutcomes(linksPairDataset, 
#'     linksNames = c("ExtendedID", "R", "RelationshipPath"))
#' 
#' @param outcomeDataset A data frame containing the outcome variable(s)
#' @param linksPairDataset A data frame containing the \code{SubjectTag}s of
#' each subject in the pair and their \code{R} coefficient.
#' @param outcomeNames The column names of the outcome variable(s)
#' @param linksNames The column names desired to be prseent in the newly
#' created data frame.  \code{SubjectTag_S1} and \code{SubjectTag_S2} are included
#' automatically.
#' @param validateOutcomeDataset Indicates if characteristics of the outcomeDataset should be validated.
#' @param subject1Qualifier Indicates how the outcome variable for the pair's first subject is distinguished from the other subject.  The default is \code{_S1}.
#' @param subject2Qualifier Indicates how the outcome variable for the pair's second subject is distinguished from the other subject.  The default is \code{_S2}.
#' @author Will Beasley
#' @references For more information about a DF analysis, see Rodgers, Joseph Lee, & Kohler, Hans-Peter (2005).
#' \href{http://www.springerlink.com/content/n3x1v1q282583366/}{Reformulating and simplifying the DF analysis model.}
#' \emph{Behavior Genetics, 35} (2), 211-217.
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
CreatePairLinksDoubleEntered <- function( 
  outcomeDataset, linksPairDataset, outcomeNames, 
  linksNames=c("ExtendedID", "R", "RelationshipPath"), validateOutcomeDataset=TRUE,
  subject1Qualifier="_S1", subject2Qualifier="_S2" ) {
  
  ValidatePairLinks(linksPairDataset)
  if(validateOutcomeDataset) ValidateOutcomeDataset(dsOutcome=outcomeDataset, outcomeNames=outcomeNames)
  
  dsLinksLeftHand <- base::subset(linksPairDataset, select=c("SubjectTag_S1","SubjectTag_S2", linksNames)) #'Lefthand' is my slang for Subjec1Tag is less than the SubjectTag_S2
  dsLinksRightHand <- base::subset(linksPairDataset, select=c("SubjectTag_S1", "SubjectTag_S2", linksNames))
  
  base::colnames(dsLinksRightHand)[base::colnames(dsLinksRightHand)=="SubjectTag_S1"] <- "SubjectTempTag"
  base::colnames(dsLinksRightHand)[base::colnames(dsLinksRightHand)=="SubjectTag_S2"] <- "SubjectTag_S1"
  base::colnames(dsLinksRightHand)[base::colnames(dsLinksRightHand)=="SubjectTempTag"] <- "SubjectTag_S2"
  
  dsOutcomeSubject1 <- base::subset(outcomeDataset, select=c("SubjectTag", outcomeNames))
  dsOutcomeSubject2 <- base::subset(outcomeDataset, select=c("SubjectTag", outcomeNames))
  
  for( j in 1:ncol(dsOutcomeSubject1) ) {
    columnName <- base::colnames(dsOutcomeSubject1)[j]
    if( columnName %in% outcomeNames ) {
      base::colnames(dsOutcomeSubject1)[colnames(dsOutcomeSubject1) == columnName] <- base::paste0(columnName, subject1Qualifier)
      base::colnames(dsOutcomeSubject2)[colnames(dsOutcomeSubject2) == columnName] <- base::paste0(columnName, subject2Qualifier)      
    }
  }
  
  dsLeftHand <- base::merge(x=dsLinksLeftHand, y=dsOutcomeSubject1, by.x="SubjectTag_S1", by.y="SubjectTag", all.x=TRUE)
  dsLeftHand <- base::merge(x=dsLeftHand, y=dsOutcomeSubject2, by.x="SubjectTag_S2", by.y="SubjectTag", all.x=TRUE)
  
  dsRightHand <- base::merge(x=dsLinksRightHand, y=dsOutcomeSubject2, by.x="SubjectTag_S2", by.y="SubjectTag", all.x=TRUE)
  dsRightHand <- base::merge(x=dsRightHand, y=dsOutcomeSubject1, by.x="SubjectTag_S1", by.y="SubjectTag", all.x=TRUE)
  
  base::rm(dsLinksLeftHand, dsLinksRightHand, dsOutcomeSubject1, dsOutcomeSubject2)
  
  ds <- base::rbind(dsLeftHand, dsRightHand) #'RowBind' the two datasets
  
  firstTwoNames <- c("SubjectTag_S1", "SubjectTag_S2")
  remaining <- base::setdiff(colnames(ds), firstTwoNames)
  ds <- ds[, c(firstTwoNames, remaining)]
  
  return( ds )
}
