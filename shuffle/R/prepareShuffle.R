prepareShuffle <-
function(des, perm){
  prep = getAveraging(des)
  prep$des = des
  prep$perm = perm
  normali = getNormalizer(prep, perm)
  prep$facA = normali$facA
  if(abs(prep$facA-1)>0.00001) {
    cat("Warning, original data factor < 1")
  }
  prep$facB = normali$facB
  prep$norm = normali$norm
  return(prep)
}
