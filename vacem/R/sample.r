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
#*  $Id: sample.r 392 2012-01-31 04:59:06Z ken $
#*  ----------------------------------------------------------------------------


##' Sample data generation functions for the Vaccination Activities Coverage
##' Estimation Model package
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
##'
NA
# TBD: 'NA' or 'NULL' or 'roxygen()' or 'sample <- function () NULL' or ???


## TBD: ADJUST FOR WT

##' Generates a sample survey population with ages uniformly distributed
##' in the specified range, (@@p age.low, @@p age.high), and all assigned
##' the provided observation @@p date.
##'
##' @note The lower and upper age limits will not be generated.  This
##'       function uses @@code{runif} to generate age values and @@c runif
##'       does not generate either extreme unless @@code{ low = high }
##'       or @@code{ max - min } is small.
##'
##' @note If @@p age.high < @@p age.low, all generated ages will be @@c NaN.
##'
##' @param N        the number of samples, i.e. population size
##' @param age.low  the minimum age (in months) for the sample population
##' @param age.high the maximum age (in months) for the sample population
##' @param date     the observation date to use for all samples
##' @return a @@c data.frame with @@p N rows and two columns: @@c date and
##'         @@c age
##'
##' @@see runif
##'
##' @examples
##'    make.sample.pop( N=10, age.low=0, age.high=25, date="2012-01-01" )
##'
##'    N <- 10
##'    dates <- c( "2005-05-01" ,"2007-07-02", "2009-09-03" )
##'    make.sample.pop( N, 0, 25, array(dates,dim=N) )
##' 
##' @export
##'
make.sample.pop <- function ( N, age.low, age.high, date ) {
  rc <- data.frame( date=date, age=floor( runif(N, age.low, age.high) ) )
  return ( rc )
}  # end of function make.sample.pop(N,age.low,age.high,date)



## TBD: ADJUST FOR WT

##' Vaccinates a sample population based on probabilities generated from the
##' campaigns descriptions and the @@p rho and @@p alpha values.
##'
##' TBD: Add details
##'
##' @note The inefficiency measure @@p alpha is defined as the natural log of psi,
##'       i.e. @@code{psi = e^alpha}.  @@code{psi = 0} (@@code{alpha = -inf})
##'       represents perfect efficiency, i.e. when every dose results in an
##'       additional vaccinee.  @@code{psi = 1} (@@code{alpha = 0}) represents
##'       random efficiency, i.e. when probability of any dose resulting in an
##'       additional vaccinee is independent of receiving a dose previously
##'       during the same campaign.
##'
##' @param obs    the synthetic population, i.e. a data frame with one row
##'               per observation and columns for @@c date of observation
##'               and @@c age of individual on that date;
##'               <i>Note: all ages should be in months</i>
##' @param camps  a data frame with one row per activity and columns providing
##'               the activity's @@c date, targeted age range (@@c age.low and
##'               @@c age.high), size of targeted population (@@c N) and the
##'               number of vaccine doses nominally distributed (@@c v)
##' @param rho    the proportion of the population (scalar) that can be vaccinated,
##'               i.e @@code{1 - rho} is the @@e inaccessible portion of the
##'               general population.
##' @param alpha  the vaccination efficiency parameter (scalar) for the campaign;
##'               @@code{alpha = -inf} implies perfect efficiency,
##'               @@code{alpha = 0} implies random efficiency;
##'               see note above
##' @param cdf.fun   a cumulative distribution function (CDF) used to
##'               calculate the probability of routine vaccination based
##'               upon age, e.g. @@code{ constant.post.9mo }.
##' @param ...    any additional arguments to the @@p cdf.fun, e.g. @@c lambda
##' @return a copy of @@p obs data frame with two new columns appended:
##'         @@c y (vaccination status, i.e. @@c 0 or @@c 1) and @@c p.vac
##'         (vaccination probability)
##'
##' @@see g
##' @@see z.matrix
##' @@see w.matrix
##' @@see nine.month.pointmass
##' @@see constant.post.9mo
##'
##' @examples
##'    N <- 100
##'    dates <- as.Date( c( "2005-05-01", "2007-07-02", "2009-09-03" ) )
##'    sample.obs   <- make.sample.pop( N=N, age.low=0, age.high=25,
##'                                     date=array(dates,dim=N) )
##'    sample.camps <- data.frame( date=(dates - 3*DAYS.PER.MONTH),
##'                                N=(0.5*N), v=(0.5*N),
##'                                age.low=8, age.high=20, is.SIA=1 )
##'
##'    vaccinate.sample.pop( sample.obs, sample.camps, rho=0.9, alpha=0,
##'                          cdf.fun=nine.month.pointmass )
##' 
##' @export
##'
vaccinate.sample.pop <- function ( obs, camps, rho, alpha, cdf.fun, ... ) {

  z <- z.matrix( obs, camps )
  w <- w.matrix( obs, camps, cdf.fun, ... )

  p.vac <- sapply( 1:nrow(obs), function(i) {
    g( z[i,], w[i,], camps$v, camps$N, alpha, rho )
  })

  obs$y <- as.numeric( runif(nrow(obs),0,1) < p.vac )
  obs$p.vac <- p.vac

  return ( obs )

}  # end of function vaccinate.sample.pop(obs,camps,rho,alpha,cdf.fun,...)



## TBD: ADJUST FOR WT

##' Calculates the number of vaccinations each member of a sample population
##' would receive based upon probabilities generated from the campaigns
##' descriptions and provided the @@p rho and @@p alpha values.
##'
##' TBD: Add details
##'
##' @note The inefficiency measure @@p alpha is defined as the natural log of psi,
##'       i.e. @@code{psi = e^alpha}.  @@code{psi = 0} (@@code{alpha = -inf})
##'       represents perfect efficiency, i.e. when every dose results in an
##'       additional vaccinee.  @@code{psi = 1} (@@code{alpha = 0}) represents
##'       random efficiency, i.e. when probability of any dose resulting in an
##'       additional vaccinee is independent of receiving a dose previously
##'       during the same campaign.
##'
##' @param obs    the synthetic population, i.e. a data frame with one row
##'               per observation and columns for @@c date of observation
##'               and @@c age of individual on that date;
##'               <i>Note: all ages should be in months</i>
##' @param camps  a data frame with one row per activity and columns providing
##'               the activity's @@c date, targeted age range (@@c age.low and
##'               @@c age.high), size of targeted population (@@c N) and the
##'               number of vaccine doses nominally distributed (@@c v)
##' @param rho    the proportion of the population (scalar) that can be vaccinated,
##'               i.e @@code{1 - rho} is the @@e inaccessible portion of the
##'               general population.
##' @param alpha  the vaccination efficiency parameter (scalar) for the campaign;
##'               @@code{alpha = -inf} implies perfect efficiency,
##'               @@code{alpha = 0} implies random efficiency;
##'               see note above
##' @param cond.on.accessible   a flag indicating whether the probabilities
##'               should be conditioned on the individual's inclusion in the
##'               accessible population, i.e. @@p rho; the default is @@c FALSE
##' @param z      the eligibility matrix for each individual/campaign pairing,
##'               i.e. @@code{z[i,j] = 1} if this individual @@c i is eligible
##'               for campaign @@c j, otherwise @@code{z[i,j] = 0};
##'               <i>Note: if @@c z is @@c NULL, then the eligibilty matrix
##'               is calculated from @@p obs and @@p camps parameters using
##'               @@code{ z.matrix } function.</i>
##' @param w      the weight matrix for each individual/campaign pairing, i.e.
##'               @@code{ w[i,k] } gives the probability of individual @@c i
##'               receiving a routine vaccination in year @@c k
##'               <i>Note: if @@c w is @@c NULL, then the weight matrix is
##'               calculated from @@p obs, @@p camps, @@p cdf.fun and @@p ...
##'               parameters using @@code{ w.matrix } function.</i>
##' @param cdf.fun   a cumulative distribution function (CDF) used to
##'               calculate the probability of routine vaccination based
##'               upon age, e.g. @@code{ constant.post.9mo }.
##'               <i>Note: this parameter is only required if the weight
##'               matrix, @@p w, is <u>not</u> provided.</i>
##' @param ...    any additional arguments to the @@p cdf.fun, e.g. @@c lambda
##'               <i>Note: this parameter is only required if the weight
##'               matrix, @@p w, is not provided.</i>
##' @return a vector containing the number of vaccinations each individual
##'         in observation set, @@p obs, would have received
##'
##' @@see z.matrix
##' @@see w.matrix
##' @@see nine.month.pointmass
##' @@see constant.post.9mo
##' @@see get.win.smooth.expected.mcmcres
##' @@see get.win.smooth.thresh.mcmcres
##'
##' @examples
##'    N <- 100
##'    dates <- as.Date( c( "2005-05-01", "2007-07-02", "2009-09-03" ) )
##'    sample.obs   <- make.sample.pop( N=N, age.low=0, age.high=25,
##'                                     date=array(dates,dim=N) )
##'    sample.camps <- data.frame( date=(dates - 3*DAYS.PER.MONTH),
##'                                N=(0.5*N), v=(0.5*N),
##'                                age.low=8, age.high=20, is.SIA=1 )
##'
##'    sim.n.vacc.sample.pop( sample.obs, sample.camps, rho=0.9, alpha=0,
##'                          cdf.fun=nine.month.pointmass )
##' 
##' @export
##'
sim.n.vacc.sample.pop <- function ( obs, camps, rho, alpha,
                                    cond.on.accessible = FALSE,
                                    z = NULL, w = NULL, cdf.fun = NULL, ... )
{
  if ( is.null(z) ) {
    z <- z.matrix( obs, camps )
  }

  if ( is.null(w) ) {
    w <- w.matrix( obs, camps, cdf.fun,  ... )
  }

  # TBD: Add data validity checks, e.g. nrow(obs) == nrow(z) == nrow(w), and etc

  if ( alpha == 0 ) {
    f <- ( exp( -camps$v / (camps$N * rho) ) )
  } else {
    f <- ( ( 1 - camps$v/(camps$N * rho) * (1-exp(alpha)) )^( 1/(1-exp(alpha)) ) )
  }

  c.prob <- 1 - f

  if ( ! cond.on.accessible ) {
    c.prob <- rho * c.prob
  }

  # TBD: Change 'nrow(obs)' to 'nrow(w)' because then obs is not needed at all
  #      by this method if z & w are provided

  n.vac <- sapply( 1:nrow(obs), function(i) {
    prob.pseudo <- 1 - ( sum(w[i,] * f) + ( 1 - sum(w[i,]) ) )
    if ( ! cond.on.accessible ) {
      prob.pseudo <- rho * prob.pseudo
    }

    my.prob <- c( c.prob * z[i,], prob.pseudo )
    # print( my.prob )

    rc <- sum( runif( length(my.prob) ) < my.prob )
  })

  return ( n.vac )

}  # end of function sim.n.vacc.sample.pop(obs,camps,rho,alpha,cond.on.accessible=FALSE,z=NULL,w=NULL,cdf.fun=NULL,...)




