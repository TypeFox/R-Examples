
pvsRequest9 <- 
  function (request,inputs) {
  pvs.url <- paste("http://api.votesmart.org/",request,"key=",get('pvs.key',envir=.GlobalEnv),inputs,sep="") #generate url for request
  
  if (names(xmlRoot(xmlTreeParse(pvs.url,useInternalNodes=TRUE)))[1]=="errorMessage") {
    
    # if the requested data is not available, return an empty (NA) data frame and give a warning
    warning(gsub(pattern="&", replacement=" ", x=paste("No data available for: ", inputs,". The corresponding rows in the data frame are filled with NAs.", sep=""), fixed=TRUE), call.=FALSE)
    
    
    output.df <- data.frame(matrix(nrow=1,ncol=0))
    output.df
    
    
  } else {
    
    
    
    
    output <- xmlRoot(xmlTreeParse(pvs.url,useInternalNodes=TRUE))
    
    nodenames <- names(output) # get names of nodes
    
    # remove unnecessary child-nodes
    
    
    if (nodenames[1]=="generalInfo") {
    
    
    output <- removeChildren(output,kids="generalInfo")
    
    } else { if (nodenames[1]=="generalinfo") {
    
    output <- removeChildren(output,kids="generalinfo")
    
    }
    }
    
    nodenames <- names(output) # get names of nodes
    
  
    # process each relevant subnode of output (e.g., without generalinfo) individually
    # specifically check if subnode has again subnodes with several entries such as stage in the Election.getElection above.Then process each individually
    
    subnodes.list <- lapply(1:length(nodenames), FUN= function(x) {
      
      subn <- output[[x]]
      
      freq.names <- summary(as.factor(names(subn))) # check frequency
      n.subnames <- sapply(names(freq.names), function(y) { length(names(subn[[y]]))}) # check number of subnodes in each element
      
      
      if (sum(n.subnames)==length(n.subnames)) {
        
        subnode.df <- data.frame(t(xmlSApply(subn, xmlValue)))
        
        
      } else {
        
        
        # first, extract xml values of subnodes of subnodes (would otherwise generate problematic columns)
        
        severalsubn <- which(freq.names>1)
        
        sevsubs <- lapply (1:length(severalsubn), FUN=function(i) {
          
          
          .sevsub <- which(names(subn)==names(severalsubn)[i]) # which subnodes (of subn) have several entries?
          
          sevsub.list <- lapply(1:length(.sevsub), function(j){ 
            
            df <-   data.frame(t(xmlSApply(subn[[.sevsub[[j]] ]], xmlValue)))
            
            df.names <- sapply(1:length(names(df)), FUN=function(z) {paste(names(.sevsub[j]),j,".",names(df)[z], sep="")})
            
            names(df) <- df.names
            
            df
            
          }
          )
          
          
          sevsub.df <- do.call("cbind", sevsub.list)
          
          
        }
        )
        
        sevsubs.extracted <-  do.call("cbind", sevsubs)
        
        
        # second, extract all xml values as usual (generates some problematic columns due to subnodes in subnodes)
        # therefore only keep the normal ones:
        
        df.ok <- data.frame(t(xmlSApply(subn, xmlValue)))
        ok.names <- names(subn)[(names(subn)!=names(severalsubn)[1])]
        df.ok <- df.ok[,ok.names]
        
        
        # third, cbind the two data frames
        
        subnode.df <- cbind(df.ok,sevsubs.extracted)
        
      }
      
      subnode.df
      
    })
    
    output.df <- dfList(subnodes.list)
    
    output.df
    
  }
  
  
}

