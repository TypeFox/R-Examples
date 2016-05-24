dep.t.test <-
function(formula, data, block, sig.level=.05, digits=3){
  xx <- aggregate(formula, data=data, FUN=sumfun)[[2]]
    
  indvar <- unlist(strsplit(as.character(formula), " ")[3])
  depvar <- unlist(strsplit(as.character(formula), " ")[2])
  datawide <- reshape(data, direction="wide", idvar=block, timevar=indvar)
    
  corr <- cor(datawide[, paste(depvar, ".", levels(data[, indvar]), sep="")])[1,2]
  output <- dep.t.test.second(m=xx[,1], sd=xx[,2], n=xx[1,3], corr=corr, sig.level=sig.level,digits=digits)
  return(output) 
}

