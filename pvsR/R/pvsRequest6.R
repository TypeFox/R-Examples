pvsRequest6 <-
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
    
    
    
    output <- output.base
    
    nodenames <- unique(names(output)) # get names of nodes
    
    freq.names <- summary(as.factor(names(output))) # check frequency, 
    
    
    # if there are some nodes with one entry, and others with many, the ones with one contain data that belongs to every single entry of the one
    # with several entries (like invariant variable over different observations, --> makrodata)
    
    # hence, the following approach: 
    
    
    # check if there are several nodes with the same name
    
    if (sum(freq.names)>length(nodenames)) { 
      
      
      freq0 <- unique(names(freq.names[freq.names==0]))
      freq1 <- unique(names(freq.names[freq.names==1]))
      freqh <- unique(names(freq.names[freq.names>1]))
      
      # ignore generalInfo:
      freq1b <- freq1[freq1!="generalInfo"]
      
      
      # scrap all data from nodes that come up once, and cbind the resulting dfs
      freq1.list <- lapply(freq1b, FUN=function(i) {
        
        x <- output[[i]]
        
        
        data.frame(t(xmlSApply(x, xmlValue)))
        
        
      })
      
      freq1.df <- do.call("cbind", freq1.list)
      
      
      # remove freq1-nodes from output.
      
      for (i in freq1)   output <- removeChildren(output,i)
      
      
      
      # now scrap the remaining output
      
      output.list <- lapply(1:length(names(output)), FUN=function(i) {
        
        x <- output[[i]]
        
        data.frame(t(xmlSApply(x, xmlValue)))
        
        
      })
      
      
      output2 <- dfList(output.list)
      
      output2 <- cbind(freq1.df,output2)
      output2
      
    } else {
      
      output <- t(xmlSApply(removeChildren(output.base,kids=1), function(x) xmlSApply(x, xmlValue)))
      output.df <- data.frame(output, row.names=NULL)
      output.df
      
    }
    
    
  }
}
