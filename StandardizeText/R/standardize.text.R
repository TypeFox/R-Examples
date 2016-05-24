standardize.text <-
function(input,input.column=NULL,standard,standard.column=NULL,regex=NULL,codes=NULL,match=FALSE,only.names=FALSE,na.rm=FALSE,suggest="prompt",print.changes=TRUE,verbose=FALSE) {
  return(standardize.names(input,input.column,standard,standard.column,regex,codes,match,only.names,na.rm,suggest,print.changes,verbose))
}
