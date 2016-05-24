`write.Data` <-
function(x,FileNameExtension="Data"){
      if (length(x)==3){
            fillmat <- matrix("",nrow=nrow(x[[2]]),ncol=ncol(x[[3]]))
            
            upperpart <- cbind(fillmat,as.matrix(x[[2]]))
            upperpart[1:nrow(upperpart),ncol(x[[3]])] <- rownames(x[[2]])
            
            lowerpart <- cbind(as.matrix(x[[3]]),x[[1]])
            lowerpart <- rbind(colnames(lowerpart),lowerpart)
            
            alldat <- rbind(upperpart,lowerpart)
            colnames(alldat) <- colnames(lowerpart)
                   
            write.table(alldat,file=paste(labels(x)[1],".txt",sep="")
                        ,sep="\t",,row.names=F,col.names=F)
       }     
            if (length(x)==4){
            
                  fillmat <- matrix("",nrow=nrow(x[[3]]),ncol=ncol(x[[4]]))
            
                  upperpart <- cbind(fillmat,as.matrix(x[[3]]))
                  upperpart[1:nrow(upperpart),ncol(x[[4]])] <- rownames(x[[3]])
                  
                  lowerpart <- cbind(as.matrix(x[[4]]),x[[1]])
                  lowerpart <- rbind(colnames(lowerpart),lowerpart)
                  lowerpartII <- cbind(as.matrix(x[[4]]),x[[2]])
                  lowerpartII <- rbind(colnames(lowerpart),lowerpartII)
                  
                  alldat <- rbind(upperpart,lowerpart)
                  alldat[1:nrow(x[[3]]),ncol(x[[4]])] <- rownames(x[[3]])
                  alldatII <- rbind(upperpart,lowerpartII)
                  alldatII[1:nrow(x[[3]]),ncol(x[[4]])] <- rownames(x[[3]])
                  
                  colnames(alldat) <- colnames(lowerpart)
                  colnames(alldatII) <- colnames(lowerpart)
              
                  write.table(alldat,file=paste(FileNameExtension,labels(x)[1],".txt",sep="")
                              ,sep="\t",row.names=F,col.names=F)
                  write.table(alldatII,file=paste(FileNameExtension,labels(x)[2],".txt",sep="")
                              ,sep="\t",row.names=F,col.names=F)
            }      
}

