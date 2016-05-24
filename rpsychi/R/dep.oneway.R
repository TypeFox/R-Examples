dep.oneway <-
function(formula, data, block, contr=NULL, sig.level=.05, digits=3){
  xx <- aggregate(formula, data=data, FUN=sumfun)[[2]]
  n <- nlevels(factor(data[, block]))
  indvar <- unlist(strsplit(as.character(formula), " ")[3])
  depvar <- unlist(strsplit(as.character(formula), " ")[2])
  datawide <- reshape(data, direction="wide", idvar=block, timevar=indvar)  
  corr <- cor(datawide[, paste(depvar, ".", levels(data[, indvar]), sep="")])
  output <- dep.oneway.second(m=xx[,1], sd=xx[,2], n=n, corr=corr, contr=contr, sig.level=sig.level,digits=digits)
  return(output)
}

