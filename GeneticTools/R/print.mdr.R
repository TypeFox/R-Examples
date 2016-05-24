# Version: 30-11-2012, Daniel Fischer

`print.mdr` <- function(x,...){
  X <- list()
  tempFold <- x$fold
  fix <- x$fix
  precPos <- c(3,7,11,16)
  for(i in 1:tempFold){
    temp <- unlist(x$mdr[precPos[i]])[1:x$top]
    if(i==1){
      if(fix!=0) temp <- temp[x$top]
      temp2 <- colnames(x$X)[unlist(x$mdr[precPos[i]+1])[1:x$top]]
      if(fix!=0) temp2 <- temp2[1]
    } else if(i==2){
      tempPos1 <- unlist(x$mdr[precPos[i]+1])[1:x$top]
      tempPos2 <- unlist(x$mdr[precPos[i]+2])[1:x$top]
      temp2 <- paste(colnames(x$X)[tempPos1],colnames(x$X)[tempPos2],sep=",")
    } else if(i==3){
      tempPos1 <- unlist(x$mdr[precPos[i]+1])[1:x$top]
      tempPos2 <- unlist(x$mdr[precPos[i]+2])[1:x$top]
      tempPos3 <- unlist(x$mdr[precPos[i]+3])[1:x$top]
      temp2 <- paste(colnames(x$X)[tempPos1],colnames(x$X)[tempPos2],colnames(x$X)[tempPos3],sep=",")
    } else if(i==4){
      tempPos1 <- unlist(x$mdr[precPos[i]+1])[1:x$top]
      tempPos2 <- unlist(x$mdr[precPos[i]+2])[1:x$top]
      tempPos3 <- unlist(x$mdr[precPos[i]+3])[1:x$top]
      tempPos4 <- unlist(x$mdr[precPos[i]+4])[1:x$top]
      temp2 <- paste(colnames(x$X)[tempPos1],colnames(x$X)[tempPos2],colnames(x$X)[tempPos3],colnames(x$X)[tempPos4],sep=",")
    }
    X[[i]] <- apply(t(data.frame(Prec=temp,SNP=temp2)),1,rev) 
    if((i==1) && (fix!=0)){

    }else {
      rownames(X[[i]]) <- 1:x$top
    }
  }
  X$cvRes <- x$cvRes
  print(X,...)
} 
