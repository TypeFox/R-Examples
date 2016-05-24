GT.wrapper <-
function( TestStatistic, alpha=0.05, eta=alpha, pi1.ini=0.7, pi2.1.ini=0.4, L=2, muL.ini=c(-1,1), sigmaL.ini=c(1,1), cL.ini=c(0.5,0.5), DELTA=0.001, sigma.KNOWN=FALSE )
  {
    em.esti <- GT.em( TestStatistic, pi1.ini, pi2.1.ini, L, muL.ini, sigmaL.ini, cL.ini, DELTA, sigma.KNOWN ) ## Calculate the estimator based on EM algorithm
    TestStatistic <- GT.localfdr(TestStatistic, em.esti$pi1, em.esti$pi2.1, L, em.esti$muL, em.esti$sigmaL, em.esti$cL) ## Calculate the local fdr scores for all the hypothesis
    GroupTest <- GT.decision(TestStatistic, alpha, eta=alpha)
    GroupTest$parameter <- em.esti
    GroupTest
  }
