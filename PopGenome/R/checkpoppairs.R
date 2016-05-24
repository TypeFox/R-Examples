checkpoppairs <- function(npops,popmissing,pairs,nn){

  if(length(popmissing)!=0){respop <- (1:npops)[-popmissing]}else{respop <- 1:npops}

      if(length(popmissing)!=0)
       {
           id <- NULL
           for(zz in 1:(length(popmissing))){
           id2        <-  which(pairs==popmissing[zz],arr.ind=TRUE)[,2]
           id         <-  c(id,id2)
           }
           respairpop <- (1:length(nn))[-unique(id)]


       }else{respairpop <- 1:length(nn)}
       
 return(list(respop=respop,respairpop=respairpop))

}