delta.log <- function(bnode, noden, simutry)
{
  bnode1 <- bnode[bnode$rownn!=noden$label,]
  abc <- order(bnode1$rownn)
  bnode1 <- bnode1[abc,]  
  deltalog <- likeli(nrow(bnode), sum(bnode$y)) - likeli(noden$total,noden$active)
  nuniq <- unique(bnode1$rownn)
  if (length(nuniq)>=1)
  {        # tot=acti=node_cha=acti_cha=numeric(length(simutry)*length(nuniq))
          #  dim(tot)=dim(acti)=dim(acti_cha)=dim(node_cha)=c(length(simutry),length(nuniq))
    for (k in 1:length(nuniq))
    {
      nnode <- bnode1[bnode1$rownn==nuniq[k],]
      for (m in 1:length(simutry)) 
      { 
        if (simutry[[m]]$label==nuniq[k]) 
        {
          tot=simutry[[m]]$total
          acti=simutry[[m]]$active
          
         # tot[m,k] <- simutry[[m]]$total
        #  acti[m,k] <- simutry[[m]]$active
        #  node_cha[m,k] <- nrow(nnode)
        #  acti_cha[m,k] <- sum(nnode$y)
        }  
      }
      deltalog <- likeli(tot-nrow(nnode), acti-sum(nnode$y)) - likeli(tot, acti) + deltalog
    }
  }
  else deltalog <- 0
  return(deltalog)
}
