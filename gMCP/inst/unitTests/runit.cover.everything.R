test.everything.analysis <- function() {
  graph <- HungEtWang2010()
  graphAnalysis(graph)
  checkEquals(gMCP:::accessible(graph, "H_{1,NI}"), c(FALSE, TRUE, TRUE, TRUE))
  gMCP:::checkOptimal(graph)
  checkEquals(gMCP:::as.num(c(1, NA, NaN, Inf, "1")), c(1, NA, NaN, Inf, 1))
}

test.everything.calcPower <- function() {
  
  # calcPower calls (and therefore tests in some way) extractPower and graphTest
  
  weights <- c(1/2, 1/2, 0, 0)  
  G <- simpleSuccessiveII()@m
  ## alternative (mvn distribution)
  corMat <- rbind(c(1, 0.5, 0.5, 0.5/2),
                  c(0.5,1,0.5/2,0.5),
                  c(0.5,0.5/2,1,0.5),
                  c(0.5/2,0.5,0.5,1))
  theta <- c(3, 0, 0, 0)
  calcPower(weights=weights, alpha=0.025, G=G, mean=theta, corr.sim=corMat, n.sim = 100000)
  
  G <-BonferroniHolm(7)
  C<-t(matrix(c(0, 1, 0,   -1,    0,
                0, 0, 1,   -1,    0,
                1, 0, 0,   -1,    0,
                0, 1, 0, -0.5, -0.5,
                0, 0, 1, -0.5, -0.5,
                0, 1, 0,    0,   -1,
                0, 0, 1,    0,   -1), ncol=7));
  cov<-C%*%t(C)
  dcov<-diag(cov)
  corr<-cov/sqrt(dcov%*%t(dcov))
  result <- calcPower(weights=G@weights, alpha=0.025, G=G@m, mean=rep(1,7), corr.sim=corr, n.sim = 100000)
  checkEquals(unname(round(result$LocalPower, 2)), rep(0.05, 7)) # No theory behind this result - just the current calculation result (which should not change without reason).
  
  # The next command calls also resultL2Text.
  
  gMCP:::calcMultiPower(weights=G@weights, alpha=0.025, G=G@m, muL=list(rep(1,7), rep(0,7)), sigmaL=list(rep(1,7), rep(2,7)), nL=list(rep(10,7), rep(20,7)),
                             # sigma = diag(length(muL[[1]])), cr = NULL,
                             n.sim = 10000, type = "quasirandom",
                             f=list(), digits=4, variables=NULL)
  
}

# pretty minor TODO but for completeness: test updateGraphToNewClassDefinition

