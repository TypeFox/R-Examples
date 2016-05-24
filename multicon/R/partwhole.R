partwhole <-
function(x, nitems=1, nomiss=.8) {
  
  Utarg <- composite(data.frame(x), nomiss=nomiss)
  Ftarg <- principal(x, scores=T)$scores
  
  totItems <- ncol(x)
  
  comps <- matrix(NA, nrow=nrow(x), ncol=choose(totItems,nitems))
  combs <- combn(totItems,nitems)
  combnames <- apply(combs, 2, paste, collapse="_")
  for(i in 1:ncol(comps)) {
    comps[,i] <- composite(x[,combs[,i]], nomiss=nomiss)
  }
  Umatch <- cor(Utarg, comps, use="pair")
  Fmatch <- cor(Ftarg, comps, use="pair")
  out <- rbind(Umatch, Fmatch)
  colnames(out) <- combnames
  rownames(out) <- c("UnitWgt", "Component")
  return(out)
}
