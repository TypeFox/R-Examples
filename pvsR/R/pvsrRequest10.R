pvsRequest10 <-
  function (request,inputs, separate=NULL) {
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
      
      if (length(separate)==0) {  #extract all xml-values as usual (result: one df with one row, many columns )
                  
      output <- output.base # useInternal=TRUE allows for this simple processing approach, but does not allow splitting up the output! Hence, in the following else{} it is the other way around!
      
      nodenames <- names(output) # get names of nodes
      
      # remove unnecessary child-nodes
      
      
      if (nodenames[1]=="generalInfo") {
        
        
        output <- removeChildren(output,kids="generalInfo")
        
      } else { if (nodenames[1]=="generalinfo") {
        
        output <- removeChildren(output,kids="generalinfo")
          }
      }

      
      
      # extract values of remaining nodes
      
      
      outputode.list <- lapply(1:length(names(output)), function(i) {
        
        if (length(unlist(xmlToList(output[[i]])))>1) {
          
          data.frame(t(unlist(xmlToList(output[[i]]))))
          
          
        } else {
          
          val <-  xmlValue(output[[i]])
          val.df <- data.frame(i=t(val), row.names=NULL)
          names(val.df) <- xmlName(output[[i]])
          val.df
          
        }
      }
      )
      
      
      outputode.df <- do.call("cbind",outputode.list)
      outputode.df
      
 }else{ # if some subnodes have to be processed and returned separately.       
      # the xml-doc has to be parsed several times because useInternalNodes=TRUE is necessary for the processing but doesn't allow the separation of the parsed xml-doc.
      # approach: generate a nested list of the separate-list (with the arguments) and a list of all nodenames but the ones from the separate list, then lapply  
  
  
  
    
 
  # download the xml-doc twice, process separately, result is returned in a list:
  
  #--------- Fist the nodes that were defined in the separate-input -----------------
  # these nodes contain several subnodes of the same type with several values, such as sponsors or candidates, hence they are to be processed normally:
  
  
  output <- output.base 
  nodenames <- names(output) # get names of nodes
  
  # remove generally unnecessary child-nodes:
  
  if (nodenames[1]=="generalInfo") {
    
    
    output <- removeChildren(output,kids="generalInfo")
    
  } 
  
  nonsepl <- names(output) %in% separate
  nonseparate <- names(output)[!nonsepl]
  
    
  # remove specified nodes for separate processing (here remove the ones in nonseparate, because the nodes of separate should be processed)
             
      for (i in nonseparate)   output <- removeChildren(output,i)
  
  # extract values for each remaining mainnode
  
  seplist <- lapply(separate, function(s) {
  
  data.frame(t(xmlSApply(output[[s]], function(x) xmlSApply(x, xmlValue))),row.names=NULL)
  
  })
  
  
 # give each list entry the coresponding name
  for (i in 1:length(seplist)) names(seplist)[i] <- separate[i]
  
 
  
  # ------- Second, the remaining nodes ---------
  
  
  output <- xmlRoot(xmlTreeParse(xmltext, useInternalNodes=TRUE))

  
  nodenames <- names(output) # get names of nodes
  
  # remove generally unnecessary child-nodes:
  
  if (nodenames[1]=="generalInfo") {
    
    output <- removeChildren(output,kids="generalInfo")
    
  } 
  
  # remove specified nodes for separate processing (here remove the ones in separate, because all but the nodes in separate should be processed )
  
  for (i in separate)   output <- removeChildren(output,i)
  
  
  
  # now process as usual:      

  outputode.list <- lapply(1:length(names(output)), function(i) {
    
    if (length(unlist(xmlToList(output[[i]])))>1) {
      
      data.frame(t(unlist(xmlToList(output[[i]]))))
      
      
    } else {
      
      val <-  xmlValue(output[[i]])
      val.df <- data.frame(i=t(val), row.names=NULL)
      names(val.df) <- xmlName(output[[i]])
      val.df
      
    }
  }
  )
  
  
  nonsep.df <- do.call("cbind",outputode.list)

  seplist$main <- nonsep.df
  
  output.list <- seplist[c("main",subset(names(seplist),!names(seplist) %in% "main"))]
  
  output.list
}
  
    }
  
}
