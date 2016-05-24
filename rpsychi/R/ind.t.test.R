ind.t.test <-
function(formula, data, correct=TRUE, sig.level=.05, digits=3){
  x <- aggregate(formula, data=data, FUN=sumfun)[[2]]
  output <- ind.t.test.second(m=x[,1], sd=x[,2], n=x[,3], correct=correct, sig.level=sig.level,digits=digits)
  return(output)
}

