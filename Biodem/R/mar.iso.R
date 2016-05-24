mar.iso=function(x){
  Pt=rep(0,dim(x)[3])
  Pr=rep(0,dim(x)[3])
  for(i in 1:dim(x)[3]){
    Pt[i]=(sum(diag(x[,,i])))/sum(x[,,i])
    Pr[i]=(sum(rowSums(x[,,i])* colSums(x[,,i])))/
      ((sum(rowSums(x[,,i])))*(sum(colSums(x[,,i]))))
  }
  pop=dimnames(x)[[3]]
  data.frame(pop,Pt,Pr)
}
