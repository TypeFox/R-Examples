ind.twoway <-
function(formula, data, sig.level=.05, digits=3){
  depvar  <- unlist(strsplit(as.character(formula), " ")[2])
  indvar1 <- unlist(strsplit(as.character(formula), " ")[3])[1]
  indvar2 <- unlist(strsplit(as.character(formula), " ")[3])[3]  
  
  m.mat  <- tapply(data[,depvar], list(data[,indvar1], data[,indvar2]), mean)
  sd.mat <- tapply(data[,depvar], list(data[,indvar1], data[,indvar2]), sd)
  n.mat  <- tapply(data[,depvar], list(data[,indvar1], data[,indvar2]), length)
  
  output <- ind.twoway.second(m=m.mat, sd=sd.mat, n=n.mat, sig.level=sig.level,digits=digits)
  return(output)
}

