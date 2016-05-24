#' @import RJSONIO
#' @import RCurl
#' @import plyr

executeURL<-
  function(fullURL, type){
    
    #print(fullURL)
    resp <- getURL(fullURL)
    out <- fromJSON(resp)
    
    if(length(out) != 0){
      if(class(out[[1]]) == 'list'){
        if(type == 'm'){ # Use a special, fast, reformatter
          # Unique is added as there is some duplication in the database
          # This fix removes duplicate entries in the data and 
          # removes duplicate Valuetype entries. This would not
          # be needed if the database did not contain these issues
          out <- YoutheriaToDF(out)
          out1 <- unique(out)
          tab1 <- table(out$MeasurementSetID)
          tab2 <- table(out1$MeasurementSetID)
          dups <- as.numeric(names(tab1[!tab1==tab2]))
          for(i in dups){            
            tempOut <- out1[out1$MeasurementSetID==i,]
            tabTemp <- names(table(tempOut$ValueType)[table(tempOut$ValueType)>1])  
            for(j in tabTemp){              
              rnam <- row.names(out1[out1$MeasurementSetID==i & out1$ValueType==j,])
              out1 <- out1[!row.names(out1) %in% rnam[2:length(rnam)],]
            }            
          }
          out <- out1
          
        } else {
          out <- ldply(out, data.frame, stringsAsFactors=FALSE)
        } 
      } else {
        out <- as.data.frame(out)
      }
      
    } else {
      out <- NULL
    }
    return(out)
  }