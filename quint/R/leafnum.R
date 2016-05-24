leafnum <-
function(nnum){
    #compute number of layers of the final tree
     nlayer<-floor(log2(max(nnum)))
    #order the leafs
    oleaf<-numeric(length(nnum))
    for (i in 1:length(nnum)){
      oleaf[i]<-nnum[i]*2^(nlayer-floor(log2(nnum[i])))
      }
     index<-sapply(1:length(oleaf),function(kk,vec){which(sort(vec)== vec[kk])},vec=oleaf )
      return(order(index))}
