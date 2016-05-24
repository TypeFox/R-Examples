r.pairs = function(x){
  RP = rep(0,dim(x)[3])
  RPr = rep(0,dim(x)[3])
  perc.diff = rep(0,dim(x)[3])
  for (i in 1:dim(x)[3]){
    RP[i] = (sum(x[,,i]*(x[,,i]-1)))/(sum(x[,,i])*(sum(x[,,i])-1))
    RPr[i] = (((1/(sum(x[,,i])*(sum(x[,,i])-1)))
               *sum((rowSums(x[,,i]))^2))-(1/(sum(x[,,i])-1)))*
                 (((1/(sum(x[,,i])*(sum(x[,,i])-1)))
                   *sum((colSums(x[,,i]))^2))-(1/(sum(x[,,i])-1)))
    perc.diff[i] = ((RP[i]-RPr[i])/RPr[i])
  }
 pop = dimnames(x)[[3]]
  data.frame(pop,RP,RPr,perc.diff)
}
