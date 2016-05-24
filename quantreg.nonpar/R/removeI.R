removeI <-
function(inString){
  
  tempstr<-""
  newstr<-inString
  while(tempstr != newstr){
    tempstr <- newstr
    newstr <- gsub("I\\(([^()]*)\\)","\\1",newstr,perl=TRUE)
    newstr <- gsub("^\\(([^()]*)\\)","SPAREN\\1EPAREN",newstr,perl=TRUE)
    newstr <- gsub("([^I])\\(([^()]*)\\)","\\1SPAREN\\2EPAREN",newstr,perl=TRUE)
  }  
  # get rid of words we introduced
  newstr <- gsub("SPAREN","(",newstr,perl=TRUE)
  newstr <- gsub("EPAREN",")",newstr,perl=TRUE)
  
  return (newstr)
}
