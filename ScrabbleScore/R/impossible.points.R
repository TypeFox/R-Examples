data(sysdata,envir=environment())
impossible.points <-
function(cl){
 
  ct <- lapply(cl,table)
  extra.points <- lapply(ct,function(t){
    sapply(names(t),function(l){
      ifelse(t[l] > letter.dists[l],(t[l]-letter.dists[l])*sls(l),0)
    })
  })
  sapply(extra.points,sum)
}
