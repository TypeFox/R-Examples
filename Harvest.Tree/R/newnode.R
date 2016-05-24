newnode <- function(ruleset, j, data){
  nnn <- ncol(ruleset)
  lmatr <- ruleset[,-j]
  if (nnn > 2) logivec <- apply(lmatr, 1, all)
  else  logivec <- lmatr       
  bnode <- data[logivec==T,]
  return(bnode)
}
