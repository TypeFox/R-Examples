GetHierID <- function(MyHier){
  if(class(MyHier) == "list")
    conceptID <- gsub("^\\D+|\\D+$", "", names(MyHier))
  if(class(MyHier) == "character")
    conceptID <- gsub("^\\D+|\\D+$", "", MyHier)
  if(any(conceptID == "") || any(is.na(conceptID)))
    conceptID[which(conceptID == "")] <- NA
  return(conceptID)
}