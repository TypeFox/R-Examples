`mu.kruskal.test` <-
function(y,groups,blocks,...) {
  funCall <- match.call()
  funCall[[1]] <- as.name("prentice.test")
  return(eval(funCall, sys.parent()))
}
