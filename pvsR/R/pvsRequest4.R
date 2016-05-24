pvsRequest4 <-
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
  
  
  if (names(output.base)[1]=="errorMessage") {
    
    # if the requested data is not available, return an empty (NA) data frame and give a warning
    warning(gsub(pattern="&", replacement=" ", x=paste("No data available for: ", inputs,". The corresponding rows in the data frame are filled with NAs.", sep=""), fixed=TRUE), call.=FALSE)
    
    
    output.df <- data.frame(matrix(nrow=1,ncol=0))
    output.df
    

    } else {
  
  output <- t(xmlSApply(removeChildren(output.base,kids=1), function(x) xmlSApply(x, xmlValue)))
  output.df <- data.frame(output, row.names=NULL)
  output.df

  }
}
