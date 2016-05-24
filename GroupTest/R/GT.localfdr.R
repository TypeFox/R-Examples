GT.localfdr <-
function(TestStatistic, pi1, pi2.1, L, muL, sigmaL, cL)
{
  pi0 <- 1-pi1
  pi2.0 <- 1-pi2.1
  G <- length( TestStatistic )
  for( g in 1:G )
    {
      mg <- TestStatistic[[g]]$mg
      TestStatistic[[g]]$f0x <- dnorm( TestStatistic[[g]]$X ) ## Calculate the density under the null for each x_{gj}
      TestStatistic[[g]]$f1x <- unlist( lapply( TestStatistic[[g]]$X, function(x, L, muL, sigmaL, cL){ sum( cL * dnorm(x, muL, sigmaL) ) }, L, muL, sigmaL, cL ) ) ## Calculate the density function udner the alternative for each x_{gj}
      TestStatistic[[g]]$fx <- pi2.0 * TestStatistic[[g]]$f0x + pi2.1 * TestStatistic[[g]]$f1x ## Calcualte the marginal density for each x_{gj}
      ## Temparary quantities
      f0.prod <- prod( 10* TestStatistic[[g]]$f0x )
      f.prod <- prod( 10* TestStatistic[[g]]$fx )
      f1.prod <- prod( 10* TestStatistic[[g]]$f1x )
      ## Calcualte the group-local fdr
      TestStatistic[[g]]$fdr.g <- pi0 * f0.prod/( pi0*f0.prod + pi1 * ( f.prod - pi2.0^mg * f0.prod )/( 1- pi2.0^mg) )
      ## Calculate the within group local fdr fdr_{j|g}
      den <- f.prod - pi2.0^mg * f0.prod
      num <- pi2.0 * TestStatistic[[g]]$f0x/TestStatistic[[g]]$fx * f.prod - pi2.0^mg * f0.prod
      TestStatistic[[g]]$fdr.j.g <- num/den
      ## Calculate the posterior probability taht both theta and theta_{j|g}=1
      TestStatistic[[g]]$prob.theta.1.theta.j.g.1 <- (1 - TestStatistic[[g]]$fdr.j.g)*( 1-TestStatistic[[g]]$fdr.g ) 

      ## Calculate the posterior prob that theta=1, theta_{j|g}=1, and x_{jg} is generated from the l-th component.
      TestStatistic[[g]]$m.j.g <- array(0, c(mg, L) )
      for( j in 1:mg )
        {
          TestStatistic[[g]]$m.j.g[j, ] <- cL * dnorm( TestStatistic[[g]]$X[j], muL, sigmaL )
          TestStatistic[[g]]$m.j.g[j, ] <- TestStatistic[[g]]$m.j.g[j,] / sum( TestStatistic[[g]]$m.j.g[j, ] )
          TestStatistic[[g]]$m.j.g[j, ] <- TestStatistic[[g]]$m.j.g[j,] * TestStatistic[[g]]$prob.theta.1.theta.j.g.1[j]
        }
    }
  TestStatistic
}
