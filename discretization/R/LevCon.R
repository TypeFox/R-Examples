LevCon <-
function(data){
        p <- ncol(data)
        n <- nrow(data)
        data.ex <- data[,-p] #data except class
        if(p==2)
        {data.ex <- as.matrix(data.ex,ncol=1)}else{data.ex <- data[,-p]}    
        dtf <- duplicated(data.ex)
        if(all(dtf)==FALSE)return(LevelConsis=1)
        dup <- data.ex[dtf,] #removed data
        if(is.numeric(dup)==TRUE) dup=as.matrix(dup) 
        duplic <- unique(dup) # duplicated data except class
        n.dup <- nrow(duplic) 
        num <- 0; difC <- 0
          for(i in 1:n.dup){
                same <- apply(as.numeric(duplic[i,])==t(data.ex),2,prod) #except class, only same data
                Sdata <- data[same==1,] #same data include class
                Stable <- table(Sdata[,p])
                difC[i] <- sum(Stable)-max(Stable) #number of differenct class, same data
          }                                     #e.g. : miss classified
        SdifC <- sum(difC)
        Lapp <- (n-SdifC) #classified with certainty:#lower approximation 
        LevelConsis <- Lapp/n
        return(LevelConsis)
}
