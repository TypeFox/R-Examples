pvsRequest <-
function (request,inputs) {
  pvs.url <- paste("http://api.votesmart.org/",request,"key=",get('pvs.key',envir=.GlobalEnv),inputs,sep="") #generate url for request
  
  httpresp <- GET(url=pvs.url)
  xmltext <- content(x=httpresp, as="text")
  errors <-  getXMLErrors(xmltext) # check if xml can be parsed properly
  
  if (length(errors) != 0) {
    
    if (names(errors[[1]]$code) == "XML_ERR_CDATA_NOT_FINISHED") { # if not, try to fix 
      
      xmltext <- gsub(pattern="\003", replacement="", x=xmltext, fixed=TRUE)
      
    }
  }
  
  output.base <- xmlRoot(xmlTreeParse(xmltext, useInternalNodes=TRUE))
  output <- t(xmlSApply(removeChildren(output.base,kids=1), function(x) xmlSApply(x, xmlValue)))
  
  output.df <- data.frame(output, row.names=NULL)
  output.df

}
