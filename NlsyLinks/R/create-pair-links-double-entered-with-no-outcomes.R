#roxygen2 documentation in CreatePairLinksDoubleEntered.R
CreatePairLinksDoubleEnteredWithNoOutcomes <- function( linksPairDataset, linksNames=c("ExtendedID", "R", "RelationshipPath") ) {
  ValidatePairLinks(linksPairDataset)
  
  dsLinksLeftHand <- base::subset(linksPairDataset, select=c("SubjectTag_S1","SubjectTag_S2", linksNames)) #'Lefthand' is my slang for Subjec1Tag is less than the SubjectTag_S2
  dsLinksRightHand <- base::subset(linksPairDataset, select=c("SubjectTag_S2", "SubjectTag_S1", linksNames))
  
  base::colnames(dsLinksRightHand)[colnames(dsLinksRightHand)=="SubjectTag_S1"] <- "SubjectTempTag"
  base::colnames(dsLinksRightHand)[colnames(dsLinksRightHand)=="SubjectTag_S2"] <- "SubjectTag_S1"
  base::colnames(dsLinksRightHand)[colnames(dsLinksRightHand)=="SubjectTempTag"] <- "SubjectTag_S2"
  
  ds <- base::rbind(dsLinksLeftHand, dsLinksRightHand) #'RowBind' the two datasets
  ds <- ds[, c("SubjectTag_S1", "SubjectTag_S2", linksNames)]
  base::rm(dsLinksLeftHand, dsLinksRightHand)
  return( ds )
}
