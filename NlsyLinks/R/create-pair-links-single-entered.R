#roxygen2 documentation in CreatePairLinksDoubleEntered.R

CreatePairLinksSingleEntered <- function( outcomeDataset, linksPairDataset, outcomeNames, 
  linksNames=c("ExtendedID", "R", "RelationshipPath"), validateOutcomeDataset=TRUE,
  subject1Qualifier="_S1", subject2Qualifier="_S2") {
  
  ValidatePairLinks(linksPairDataset)
  if(validateOutcomeDataset) ValidateOutcomeDataset(dsOutcome=outcomeDataset, outcomeNames=outcomeNames)
  
  dsLinksLeftHand <- base::subset(linksPairDataset, select=c("SubjectTag_S1","SubjectTag_S2", linksNames)) #'Lefthand' is my slang for Subjec1Tag is less than the SubjectTag_S2
  
  dsOutcomeSubject1 <- base::subset(outcomeDataset, select=c("SubjectTag", outcomeNames))
  dsOutcomeSubject2 <- base::subset(outcomeDataset, select=c("SubjectTag", outcomeNames))
  
  for( j in 1:ncol(dsOutcomeSubject1) ) {
    columnName <- base::colnames(dsOutcomeSubject1)[j]
    if( columnName %in% outcomeNames ) {
      colnames(dsOutcomeSubject1)[colnames(dsOutcomeSubject1) == columnName] <- base::paste0(columnName, subject1Qualifier)
      colnames(dsOutcomeSubject2)[colnames(dsOutcomeSubject2) == columnName] <- base::paste0(columnName, subject2Qualifier)      
    }
  }
 
  ds <- base::merge(x=dsLinksLeftHand, y=dsOutcomeSubject1, by.x="SubjectTag_S1", by.y="SubjectTag", all.x=TRUE)
  ds <- base::merge(x=ds, y=dsOutcomeSubject2, by.x="SubjectTag_S2", by.y="SubjectTag", all.x=TRUE)
 
  base::rm(dsLinksLeftHand, dsOutcomeSubject1, dsOutcomeSubject2)
  
  firstTwoNames <- c("SubjectTag_S1", "SubjectTag_S2")
  remaining <- base::setdiff(colnames(ds), firstTwoNames)
  ds <- ds[, c(firstTwoNames, remaining)]
    
  return( ds )
}
