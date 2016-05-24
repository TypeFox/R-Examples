simu.new <- function(simutry, harv)
{
  if(length(rownames(harv))>0){
  nuniq <- unique(harv$rownn)
  sim.try<-simutry
  for (m in 1:length(sim.try)) 
  {
    for (k in 1:length(nuniq))
    { 
      nnode <- harv[harv$rownn==nuniq[k],]
      if (sim.try[[m]]$label==nuniq[k]) 
      {
        tot <- simutry[[m]]$total
        acti <- simutry[[m]]$active
        sim.try[[m]]$total <- tot-nrow(nnode)
        sim.try[[m]]$active <- acti-sum(nnode$y)
      }  
    }
  }
  logi <- unlist(lapply(sim.try, function(x) {x$total==0}))
  sim.try <- sim.try[!logi]
  }
  else{
    sim.try<-simutry
    logi <- unlist(lapply(sim.try, function(x) {x$total==0}))
    sim.try <- sim.try[!logi]
  }
  return(sim.try)
}