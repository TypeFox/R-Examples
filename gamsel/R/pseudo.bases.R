pseudo.bases=function(x,degree=8,df=6,parallel=FALSE,...){
  p=ncol(x)
  degree=rep(degree,length=p)
  df=rep(df,length=p)
  out=as.list(1:p)
  
#  if (parallel && require(foreach)) {
  if (parallel){
    out = foreach(i = seq(p), .packages = c("gamsel")) %dopar% 
    {
      basis.gen(x[,i],degree=degree[i],df=df[i],...)
    }
  }
  else {
    for(i in 1:p){
      out[[i]]=basis.gen(x[,i],degree=degree[i],df=df[i],...)
    }
  }
  out
}



