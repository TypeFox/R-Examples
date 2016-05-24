select <-
function(fila, vardep.summary, ...) {
    
    if(length(which(fila==max(fila)))>1)
    {predclass <-names(vardep.summary[which(fila==max(fila))])[
      order(vardep.summary[which(fila==max(fila))],decreasing=TRUE)[1]]
    }
    else{predclass<- as.character(names(vardep.summary)[(order(fila,decreasing=TRUE)[1])])} 
    
   # ans<-predclass
   predclass
  }
