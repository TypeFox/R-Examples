meanif <-
function(set, nomiss=.8, tr=0) {
  comp <- (length(set) - sum(is.na(set))) / length(set)
  ifelse(comp >= nomiss, out <- mean(set, na.rm=T, tr=tr), out <- NA)
  return(out)
}
