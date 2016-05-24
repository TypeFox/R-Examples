plapply<-function(X, FUN, ...){
if (exists("parLapply") & exists(".GlobalEnv$cl"))
  {o<-parallel::parLapply(.GlobalEnv$cl,X, FUN,...)}
else
  {o<-lapply(X, FUN,...)}
o
}
