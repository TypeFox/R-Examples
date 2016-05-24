savepajek<-function(pajekList,filename,twomode="default",asMatrix=FALSE,symetric=NULL){
    if(length(grep(patt="w32",x=version["os"]))){
        eol<-"\n"
    }else{eol<-"\r\n"}
    cat("",file = filename,append=FALSE)
    pajekListElements<-names(pajekList)
    if("Networks" %in% pajekListElements){
        tmpList<-pajekList[["Networks"]]
        n<-length(tmpList)
        tmpNames<-names(tmpList)
        if(is.null(tmpNames)) tmpNames<-1:n
        for(i in 1:n){
            if(asMatrix){
                cat(paste("*Matrix", tmpNames[i]),eol, file = filename,append=TRUE)            
                savematrix(n=tmpList[[i]],filename,twomode=1,cont=TRUE)
            } else {
                cat(paste("*Network", tmpNames[i]),eol, file = filename,append=TRUE)
                savenetwork(n=tmpList[[i]],filename,twomode="default",symetric=NULL,cont=TRUE)
            }
            cat(eol, file = filename,append=TRUE)            
        }
    }
    if("Partitions" %in% pajekListElements){
        tmpList<-pajekList[["Partitions"]]
        n<-length(tmpList)
        tmpNames<-names(tmpList)
        if(is.null(tmpNames)) tmpNames<-1:n
        for(i in 1:n){
            cat(paste("*Partition", tmpNames[i]),eol, file = filename,append=TRUE)
            savevector(v=tmpList[[i]],filename,cont=TRUE)
            cat(eol, file = filename,append=TRUE)            
        }
    }    
    if("Vectors" %in% pajekListElements){
        tmpList<-pajekList[["Vectors"]]
        n<-length(tmpList)
        tmpNames<-names(tmpList)
        if(is.null(tmpNames)) tmpNames<-1:n
        for(i in 1:n){
            cat(paste("*Vector", tmpNames[i]),eol, file = filename,append=TRUE)
            savevector(v=tmpList[[i]],filename,cont=TRUE)
            cat(eol, file = filename,append=TRUE)            
        }
    }    
    
    if("Permutations" %in% pajekListElements){
        tmpList<-pajekList[["Permutations"]]
        n<-length(tmpList)
        tmpNames<-names(tmpList)
        if(is.null(tmpNames)) tmpNames<-1:n
        for(i in 1:n){
            cat(paste("*Permutation", tmpNames[i]),eol, file = filename,append=TRUE)
            savevector(v=tmpList[[i]],filename,cont=TRUE)
            cat(eol, file = filename,append=TRUE)            
        }
    }    
    if("Clusters" %in% pajekListElements){
        tmpList<-pajekList[["Clusters"]]
        n<-length(tmpList)
        tmpNames<-names(tmpList)
        if(is.null(tmpNames)) tmpNames<-1:n
        for(i in 1:n){
            cat(paste("*Cluster", tmpNames[i]),eol, file = filename,append=TRUE)
            cat(paste(tmpList[[i]],collapse=eol),file = filename,append=TRUE)
            cat(eol, file = filename,append=TRUE)            
        }
    }    
        
}
