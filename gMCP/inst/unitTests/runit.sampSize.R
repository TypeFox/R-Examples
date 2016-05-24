test.sampSize <- function() {
  target <- 0.8
  
  graph <- BonferroniHolm(2)  
  powerReqFunc <- function(x) { x[1] || x[2] }
  result <- sampSize(graph, esf=c(1,1), effSize=c(1,1), corr.sim=diag(2), powerReqFunc=powerReqFunc, target=target, alpha=0.05)
  
  result2 <- calcPower(graph=graph, alpha=0.05, mean = rep(sqrt(result$sampSize[1]+1), 2), f=powerReqFunc)
  checkTrue(result2[[5]]>target)
  
  if (gMCP:::tests("extended")) {
    graph <- BonferroniHolm(4)
    
    powerReqFunc=list('all(x[c(1,2)]'=function(x) {all(x[c(1,2)])},
                      'any(x[c(0,1)]'=function(x) {any(x[c(0,1)])})
    
    sampSize(graph=graph, 
             effSize=list("Scenario 1"=c(1, 0.1, 0.1, 0.1), 
                          "Scenario 2"=c(0.1, 2, 0.1, 0.1)), 
             esf=c(0.5, 0.707, 0.5, 0.707),
             powerReqFunc=powerReqFunc, 
             corr.sim=diag(4), target=c(0.8, 0.8), alpha=0.025, n.sim=100)
  }
}
