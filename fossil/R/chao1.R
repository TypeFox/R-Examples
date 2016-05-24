`chao1` <-
function(x,taxa.row=TRUE) {
  if (taxa.row==FALSE) x<-t(x)
  if (ncol(as.matrix(x))==1 || nrow(as.matrix(x))==1) ch<-(int.chao(x))
  else {ch<-int.chao(rowSums(x)) 
  if (sum(x[x>1])==0) warning("This data appears to be presence/absence based, but this estimator is for abundance data only")
    }
return(ch)
}

