'coi' <- 
function(mst, groups) {
  coitab<-table(groups)
  grpnms<-names(coitab)
  for (i in 1:length(coitab)) coitab[i] <- (sum(mst[groups==i, groups==i])/2)/(coitab[i]-1)
  return(coitab)
}
