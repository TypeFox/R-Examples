GT.em <-
function( TestStatistic, pi1.ini, pi2.1.ini, L, muL.ini, sigmaL.ini, cL.ini, DELTA, sigma.KNOWN=FALSE)
## em <- function( TestStatistic, L=2, pi1.ini=0.5, pi2.1.ini=0.5, muL.ini=rnorm(L), sigmaL.ini=rgamma(L, 2, 2), cL.ini=rdirichlet(1, array(1,L), DELTA=0.001 )
{
  G <- length( TestStatistic )
  mgs <- array(0, G)
  pi1.old <- pi1.ini;  pi2.1.old <- pi2.1.ini;  muL.old <- muL.ini;  sigmaL.old <- sigmaL.ini;  cL.old <- cL.ini ## Set the initial value as the old iteration
  pi1.new <- pi1.old;  pi2.1.new <- pi2.1.old;  muL.new <- muL.old;  sigmaL.new <- sigmaL.old;  cL.new <- cL.old ## Set the new iteration as the old iteration
  delta <- 1 ## Set the

  while( delta> DELTA){
    ##while( itr < 7 ){
    TestStatistic <- GT.localfdr(TestStatistic, pi1.old, pi2.1.old, L, muL.old, sigmaL.old, cL.old)
    total.fdr.g <- 0
    total.theta.g.1 <- 0
    total.theta.1.theta.j.g.1 <- 0
    total.m.j.g.L <- array(0, L)
    total.muL <- array(0, L)
    total.sigmaL <- array(0, L )
    for(g in 1: G)
      {
        mgs[g] <- TestStatistic[[g]]$mg
        total.fdr.g <- total.fdr.g + TestStatistic[[g]]$fdr.g
        total.theta.g.1 <- total.theta.g.1 + ( 1-TestStatistic[[g]]$fdr.g ) * mgs[g]
        total.theta.1.theta.j.g.1 <- total.theta.1.theta.j.g.1 + sum( TestStatistic[[g]]$prob.theta.1.theta.j.g.1 )
        total.m.j.g.L <- total.m.j.g.L + apply( TestStatistic[[g]]$m.j.g, 2, sum )
        total.muL <- total.muL + apply( (TestStatistic[[g]]$X %x% array(1, c(1,L)) ) * TestStatistic[[g]]$m.j.g, 2, sum )
      }
    ## Update all the parameters
    pi1.new <- 1 - total.fdr.g/G
    pi2.1.new <- total.theta.1.theta.j.g.1 / total.theta.g.1
    cL.new <- total.m.j.g.L/total.theta.1.theta.j.g.1
    muL.new <- total.muL/total.m.j.g.L
    if( sigma.KNOWN==FALSE ){
      for(g in 1:G)
        total.sigmaL <- total.sigmaL + apply( ( TestStatistic[[g]]$X %x% array(1, c(1,L)) - array(muL.new, c(1,L))%x% array(1, c(TestStatistic[[g]]$mg, 1) ))^2 * TestStatistic[[g]]$m.j.g, 2, sum )
      sigmaL.new <- sqrt( total.sigmaL/total.m.j.g.L )
    }else{
      sigmaL.new <- sigmaL.old
    }
    ## Calculate the squared distance of all the parameters
    delta <- max( abs(pi1.new-pi1.old), abs(pi2.1.new-pi2.1.old), abs(muL.new-muL.old), abs(sigmaL.new-sigmaL.old), abs(cL.new-cL.old ))
    ## Update the iteration by setting the current iteration as the old iteration
    pi1.old <- pi1.new; pi2.1.old <- pi2.1.new; muL.old <- muL.new; sigmaL.old <- sigmaL.new; cL.old <- cL.new;
    print(delta)
  }

  em.esti <- list( pi1=pi1.new, pi2.1=pi2.1.new, muL=muL.new, sigmaL=sigmaL.new, cL=cL.new, L )
  em.esti
}
