mbrd <-
function(l2levels=c(2,2), fraction=0, gen=NULL, fnames1=NULL, fnames2=NULL){

    nfac <- length(l2levels)
    
    if(is.null(fnames1)){
      fnames1 <- character()
      for(i in 1:nfac){
        fnames1 <- c(fnames1,paste(letters[i],1:l2levels[i],sep=""))
      }
    }
    if(is.null(fnames2)){
      fnames2 <- LETTERS[1:nfac]      
    }

    D1 <- FrF2(nfactors=sum(l2levels),randomize=FALSE,nruns=2^(sum(l2levels)-fraction),generators=gen,factor.names=fnames1)
    D2 <- ifelse(D1==-1,0,1)
    
    D3 <- t(apply(D2,1,.bits2int,l2levels=l2levels))
    colnames(D3) <- fnames2
    res <- list(BitDesign = D1, Design = D3+1)
    
    return(res)
}
