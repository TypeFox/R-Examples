`chao2` <-
function(x,taxa.row=TRUE) {
  if (taxa.row==FALSE) x<-t(x)
  if (ncol(as.matrix(x))==1) ch<-(int.chao(x))
  else if (nrow(as.matrix(x))==1) ch<-(int.chao(x))
    else {mat.vectorized <- numeric(length(x[,1]))
     for (i in 1:length(x[,1])) mat.vectorized[i] <- length(x[i,][x[i,]>0])
     ch<-int.chao(mat.vectorized)
       }
  return(ch)
}

