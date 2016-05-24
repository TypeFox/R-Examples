"fun.by.blocks.default" <-
function(x=M, M=x, clu, ignore.diag = identical(ss(diag(M)), 0)&&!is.list(clu), FUN = "mean",sortNames=TRUE,...)
{
    if(is.list(clu)){
      nmode<-length(clu)
      if(nmode>2){
      clu<-unlist(clu)
      clu<-list(clu,clu)
      }
    } else {
      clu<-list(clu,clu)
      nmode<-1
    }
    
    if(sortNames) {
    	k <- lapply(clu,function(x)sort(as.character(unique(x))))
    }else {
    	k <- lapply(clu,function(x)as.character(unique(x)))
    }
    IM.V <- matrix(NA, nrow = length(k[[1]]), ncol = length(k[[2]]))
    dimnames(IM.V)<-k
    for (i in k[[1]]) {
        for (j in k[[2]]) {
            B <- M[clu[[1]] == i, clu[[2]] == j, drop = FALSE]
            if (nmode==1 && i == j && dim(B)[1] > 1 && ignore.diag)
                diag(B) <- NA
            IM.V[i, j] <- do.call(FUN, list(x = B, na.rm = TRUE))
        }
    }
    return(IM.V)
}

