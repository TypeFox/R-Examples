pvsRequest3 <-
  function (request,inputs) {
    
    pvs.url <- paste("http://api.votesmart.org/",request,"key=",get('pvs.key',envir=.GlobalEnv),inputs,sep="") #generate url for request
    
    # Handle slow connection (try max. 3 times to get data from PVS, otherwhise move on)
    httpresp <- try(GET(url=pvs.url, timeout(5)), silent=TRUE)
    timedout <- any(grepl("OPERATION_TIMEDOUT", x=class(attributes(httpresp)$condition)))
    t <- 0
    while (timedout & t<3) {
      httpresp <- try(GET(url=pvs.url, timeout(5)), silent=TRUE)
      timedout <- any(grepl("OPERATION_TIMEDOUT", x=class(attributes(httpresp)$condition)))
      t <- t +1
    }
    
    
    xmltext <- content(x=httpresp, as="text")
    errors <-  getXMLErrors(xmltext) # check if xml can be parsed properly
    
    if (length(errors) != 0) {
    
    if (names(errors[[1]]$code) == "XML_ERR_CDATA_NOT_FINISHED") { # if not, try to fix 
      
      xmltext <- gsub(pattern="\003", replacement="", x=xmltext, fixed=TRUE)
      
    }
    }
    
    # in case of (still) malformed XML, force parsing problems into an empty function (does not
    # break parsing process. parser might still work.)
    output.base <- xmlRoot(xmlTreeParse(xmltext, useInternalNodes=TRUE, error=function(...){}))
    
    output <- xmlSApply(output.base, function(x) data.frame(t(xmlSApply(x, xmlValue))))
    output  
  }