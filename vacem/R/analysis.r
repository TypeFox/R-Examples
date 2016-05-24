#*  ----------------------------------------------------------------------------
#*  Copyright (C) 2011-2012 - Justin Lessler, Jessica Metcalf
#*  
#*  This program is free software; you can redistribute it and/or modify
#*  it under the terms of the GNU General Public License as published by
#*  the Free Software Foundation; version 2 of the License.
#*  
#*  This program is distributed in the hope that it will be useful,
#*  but WITHOUT ANY WARRANTY; without even the implied warranty of
#*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#*  GNU General Public License for more details.
#*  
#*  You should have received a copy of the GNU General Public License
#*  along with this program; if not, write to the Free Software
#*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#*  
#*  $Id: analysis.r 392 2012-01-31 04:59:06Z ken $
#*  ----------------------------------------------------------------------------


##' Analysis functions for the Vaccination Activities Coverage Estimation
##' Model package
##'
##' @author Justin Lessler \email{jlessler@@jhsph.edu}
##' @author Jessica Metcalf \email{cmetcalf@@princeton.edu}
##' @@maintainer Ken Cline \email{kcline@@jhsph.edu}
##' @references Lessler J, Metcalf CJE, Grais RF, Luquero FJ, Cummings DAT,
##'             et al. (2011) Measuring the Performance of Vaccination Programs
##'             Using Cross-Sectional Surveys: A Likelihood Framework and
##'             Retrospective Analysis. PLoS Med 8(10): e1001110.
##'             doi:10.1371/journal.pmed.1001110
##'             \url{ http://www.plosmedicine.org/article/info:doi/10.1371/journal.pmed.1001110 }
##' @keywords vaccine immunization likelihood models
##' 
##' @include vacem.r
##' @include coverage.r
##'
NA
# TBD: 'NA' or 'NULL' or 'roxygen()' or 'analysis <- function () NULL' or ???



##' Quick and dirty roll your own MCMC for determining the CI for rho
##' and alpha.
##'
##' @param obs        the number of observations
##' @param campaigns  TBD
##' @param keep.from  TBD
##' @param n.iter     the number of interations to do
##' @param n.chains   the number of chains to run
##' @param trans.rho.start    starting point for rho
##' @param alpha.start   starting point for alpha
##' @param trans.rho.sd  standard deviation to use for proposal rhos
##' @param alpha.sd      standard deviation to us for proposal alphas
##' @param routine.cdf        TBD
##' @param routine.par.start  TBD
##' @param routine.par.sd     TBD
##' @param z             TBD
##' @param plot.as.go    TBD
##' @param plot.freq     TBD
##' @return matrix of all of the chains for alpha and trans.rho.start
##' @export
##' 
mcmc.estimate <- function ( obs,
                            campaigns,
                            keep.from = 1,
                            n.iter = 1000,
                            n.chains = 5,
                            trans.rho.start = 0.5,
                            alpha.start = 0,
                            trans.rho.sd = .1,
                            alpha.sd = .1,
                            routine.cdf = "9m.pointmass",
                            routine.par.start = c(),
                            routine.par.sd = rep(.1, length(routine.par.start)),
                            z = NULL,
                            plot.as.go = FALSE,
                            plot.freq = 10 )
{
  require( foreach )

  # Make z-matrix
  if ( is.null(z) ) {
    z <- z.matrix( obs, campaigns )
  }


  run.chain <- function () {
    # t.rho, alpha, ll, routine.pars
    rc <- matrix( nrow=n.iter+1, ncol=3+length(routine.par.start) )

    if ( routine.cdf == "9m.pointmass" ) {
      ll <- ll.coverage( exp(trans.rho.start)/(1+exp(trans.rho.start)),
                         alpha.start, obs, campaigns, z,
                         cdf.fun = nine.month.pointmass )
    } else if ( routine.cdf == "const.post.9m" ) {
      card.wts  <- obs$samp.wt[    ! is.na(obs$age.at.vac) ]
      card.ages <- obs$age.at.vac[ ! is.na(obs$age.at.vac) ]
      card.wts  <- card.wts[  card.ages > 8.5 ]
      card.ages <- card.ages[ card.ages > 8.5 ] - 8.5 # make sure all have positive prob
      trans.lambda <- routine.par.start

      ll <- ll.coverage( exp(trans.rho.start)/(1+exp(trans.rho.start)),
                         alpha.start, obs, campaigns, z,
                         cdf.fun = constant.post.9mo, lambda=exp(trans.lambda) )

      ll <- ll + sum( card.wts * dexp( card.ages, exp(trans.lambda), log=TRUE ) )
    } else {
      stop( "Unknown routine CDF" )
    }

    rc[1,] <- c( trans.rho.start, alpha.start, ll, routine.par.start )

    for ( i in 1:n.iter ) {
      # propose new values
      t.rho.prop <- rc[i,1] + rnorm( 1, 0, trans.rho.sd )
      alpha.prop <- rc[i,2] + rnorm( 1, 0, alpha.sd     )

      # calculate relative likelihood
      # todo: update this

      if ( routine.cdf == "9m.pointmass" ) {
        ll.prop <- ll.coverage( exp(t.rho.prop)/(1+exp(t.rho.prop)),
                                alpha.prop, obs, campaigns, z,
                                cdf.fun = nine.month.pointmass )
      }  else if ( routine.cdf == "const.post.9m" ) {
        trans.lambda.prop <- rc[i,4] + rnorm( 1, 0, routine.par.sd )
        ll.prop <- ll.coverage( exp(t.rho.prop)/(1+exp(t.rho.prop)),
                                alpha.prop, obs, campaigns, z,
                                cdf.fun = constant.post.9mo, lambda = exp(trans.lambda.prop) )
        ll <- ll + sum( card.wts * dexp( card.ages, exp(trans.lambda), log=TRUE ) )
      } else {
        stop( "Unknown routine CDF" )
      }

      alpha <- exp( ll.prop - rc[i,3] )

      if ( ! is.nan(alpha) &&  runif(1) < alpha ) {
        # accept
        if ( routine.cdf == "9m.pointmass" ) {
          rc[i+1,] <- c( t.rho.prop, alpha.prop, ll.prop )
        }  else if ( routine.cdf == "const.post.9m" ) {
          rc[i+1,] <- c( t.rho.prop, alpha.prop, ll.prop, trans.lambda.prop )
        }
      } else {
        # reject
        rc[i+1,] <- rc[i,]
      }

      if ( n.chains == 1 && plot.as.go && i%%plot.freq == 0 ) {
        par( mfcol=c( 2 + length(routine.par.start), 1 ) )
        par( mar=c(2,2,2,2) )
        plot( rc[,1], type="l" )
        plot( rc[,2], type="l" )

        if ( length(routine.par.start) > 0 ) {
          for ( i in 1:length(routine.par.start) ) {
            plot( rc[,3+i], type="l" )
          }
        }
      }

    }

    if ( routine.cdf == "9m.pointmass" ) {
      rc <- rc[,1:3]
      colnames(rc) <- c( "trans.rho", "alpha", "LL" )
    } else if ( routine.cdf == "const.post.9m" ) {
      rc <- rc[,c(1,2,4,3)]
      colnames(rc) <- c( "trans.rho", "alpha", "trans.lambda", "LL" )
    }

    return ( rc )

  }  # end of function run.chain()


  rc <- foreach( 1:n.chains, .combine='cbind' ) %dopar% run.chain()

  rc <- rc[ keep.from:n.iter, ]

  if ( routine.cdf == "9m.pointmass" ) {
    rc <- list( trans.rho     = rc[,seq(1,n.chains*3,3)],
                alpha         = rc[,seq(2,n.chains*3,3)],
                LL            = rc[,seq(3,n.chains*3,3)]  )
  } else if ( routine.cdf == "const.post.9m" ) {
    rc <- list( trans.rho     = rc[,seq(1,n.chains*4,4)],
                alpha         = rc[,seq(2,n.chains*4,4)],
                trans.lambda  = rc[,seq(3,n.chains*4,4)],
                LL            = rc[,seq(4,n.chains*4,4)]  )
  }

  return ( rc )

}  # end of function mcmc.estimate(obs,campaigns,...)



##' Function to assess MCMC convergance of a chain given a matrix of MCMC
##' iterations.
##'
##' @param chains   a matrix of n.interations X n.chains with results
##' @return TBD
##' @export
##' 
mcmc.converge <- function ( chains ) {

  par( mfcol=c(2,1) )

  # first plot the chains
  plot( chains[,1], pch="+", ylim=range(chains) )
  for ( i in 2:ncol(chains) ) {
    points( chains[,i], col=i, pch="+" )
  }

  # now plot them as histograms
  plt.mx <- 0
  for ( i in 1:ncol(chains) ) {
    plt.mx <- max( plt.mx, density(chains[,i])$y )
  }

  plot( density(chains[,1]), ylim=c(0,plt.mx) )
  for ( i in 2:ncol(chains) ) {
    lines( density(chains[,i]), col=i )
  }

  # Compute R.hat using method from
  # http://warnercnr.colostate.edu/~gwhite/mark/markhelp/mcmc_diagnostics.htm
  n <- nrow( chains )
  m <- ncol( chains )

  X  <- sapply( 1:ncol(chains), function(i) { mean(chains[,i]) } )
  Xmean <- mean( X )

  SD <- sapply( 1:ncol(chains), function(i) { sd(chains[,i])   } )
  SDmean <- mean( SD )

  B <- sum( (X-Xmean)^2 ) * n / (m-1)
  W <- SDmean

  R.hat <- ( W * (n-1) / n + B/n ) / W

  cat( "R.hat=", R.hat, "\n" )

}  # end of function mcmc.converge(chains)



