`weightsGivenSize` <-
function(vec, modlist){
  modsizes <- sapply(modlist, num.Terms)
  ttab <- table(modsizes)
  if(length(vec)!=length(names(ttab)))stop("vec is wrong length")
  sizes <- as.numeric(names(ttab))
  wts <- rep(0, length(modlist))
  count = 0
  for(j in unique(sizes)){
    count = count+1
    for(k in 1:length(modlist)){
      if(modsizes[k]==j){
        wts[k] <- vec[count]/as.numeric(ttab[count])}}}
  wts
}

