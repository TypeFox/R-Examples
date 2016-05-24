ipsatize <-
function(set) {
  out <- data.frame(t(apply(set, 1, scale2)))
  colnames(out) <- paste(names(set), ".ip", sep="")
  return(out)
}
