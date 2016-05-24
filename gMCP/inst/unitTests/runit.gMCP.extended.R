test.gMCP.ext <- function() {
  graph <- BonferroniHolm(4)
  pvalues <- c(0.01, 0.05, 0.03, 0.02)
  alpha <- 0.05
  
  # Bonferroni test:
  r1 <- gMCP.extended(graph=graph, pvalues=pvalues, test=bonferroni.test, verbose=TRUE)
  r2 <- gMCP(graph=graph, pvalues=pvalues, verbose=TRUE)
  checkTrue(all.equal(r1@adjPValues, r2@adjPValues, check.attributes=FALSE)) # No names for gMCP.extended
  
  # Simes test:
  r1 <- gMCP.extended(graph=graph, pvalues=pvalues, test=simes.test, verbose=TRUE)
  r2 <- gMCP(graph=graph, pvalues=pvalues, test="Simes", verbose=TRUE)
  checkTrue(all.equal(r1@adjPValues, r2@adjPValues)) # TODO Both don't have names
  
  # Parametric test:
  rho1<-(-0.9)
  rho2<-0.9
  cr1<-matrix(c(1,rho1,rho1,1),nrow=2)
  cr2<-matrix(c(1,rho2,rho2,1),nrow=2)
  correlation <- bdiagNA(cr1, cr2)
  
  r1 <- gMCP.extended(graph=graph, pvalues=pvalues, test=parametric.test, verbose=TRUE, correlation=correlation)
  r2 <- gMCP(graph=graph, pvalues=pvalues, verbose=TRUE, correlation=correlation)
  checkTrue(all.equal(r1@adjPValues, r2@adjPValues, check.attributes=FALSE)) # No names for gMCP.extended
  
  # Trimmed Simes test:
  r1 <- gMCP.extended(graph=graph, pvalues=pvalues, test=bonferroni.trimmed.simes.test, verbose=TRUE, alpha=0.05, adjPValues = FALSE)
  r2 <- gMCP.extended(graph=graph, pvalues=pvalues, test=bonferroni.test, verbose=TRUE)
  checkTrue(all(as.numeric(r1@rejected) >= as.numeric(r2@rejected)))
  
  # Simes with partition:
  r1 <- gMCP.extended(graph=graph, pvalues=pvalues, test=simes.on.subsets.test, verbose=TRUE, alpha=0.05, subsets=list(1:2, 3:4))
  #TODO r1 <- gMCP.extended(graph=graph, pvalues=pvalues, test=simes.on.subsets.test, verbose=TRUE, alpha=0.05, subsets=list(c("H1", "H2"), c("H3", "H4")))
  r2 <- gMCP.extended(graph=graph, pvalues=pvalues, test=bonferroni.test, verbose=TRUE)
  checkTrue(all(as.numeric(r1@rejected) >= as.numeric(r2@rejected)))
}

test.gMCP.ext.warnings.and.errors <- function() {
  graph <- BonferroniHolm(4)
  pvalues <- c(0.01, 0.05, 0.03, 0.02)
  alpha <- 0.05
  
  old.warn.level <- getOption("warn")
  options(warn=2)
  
  checkException(gMCP.extended(graph=graph, pvalues=pvalues, test=bonferroni.test, subsets=list(1:2, 3:4)))
  checkException(gMCP.extended(graph=graph, pvalues=pvalues, test=bonferroni.test, correlation=diag(4)))
  checkException(gMCP.extended(graph=graph, pvalues=pvalues, test=simes.on.subsets.test, verbose=TRUE))
  checkException(gMCP.extended(graph=graph, pvalues=pvalues, test=parametric.test, verbose=TRUE))
  
  options(warn=old.warn.level)
}