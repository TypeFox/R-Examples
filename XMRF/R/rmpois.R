rmpois <-
function(n,lambda.v){
  tmp=c();
  # for(i in lambda.v){
  # tmp=rbind(tmp,rpois(n,i))
  # }
  ## require(multicore)
  tmp =  do.call(rbind, parallel::mclapply(lambda.v,function(i) { rpois(n,i)}))
  return(tmp)
}
