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
#*  $Id: coverage.r 392 2012-01-31 04:59:06Z ken $
#*  ----------------------------------------------------------------------------


##' Core coverage calculation functions for the Vaccination Activities
##' Coverage Estimation Model package
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
# TBD: 'NA' or 'NULL' or 'roxygen()' or 'coverage <- function () NULL' or ???


##' The average number of days per month.
##'
##' The @@c DAYS.PER.MONTH constant is used when converting age values
##' from days to months.
##' @export
##'
DAYS.PER.MONTH <- 30.41667



##' Defines a simple logistic function, @@f$ e^x / (1 + e^x) @@f$.
##'
##' The logistic function is used repeatedly in the coverage analysis
##' so it defined here to improve the readability of those calculations.
##'
##' @param x   the variable to calculate the logistic function for
##' @return @@code{ exp(x) / (1 + exp(x)) }
##'
##' @examples
##'    x <- c( -10:-2, seq(-2,2,0.2), 2:10 )
##'    plot( x=x, y=logistic(x) )
##'
##' @references \url{ http://en.wikipedia.org/wiki/Logistic_function }
##' @export
##'
logistic <- function ( x ) {
  return ( exp(x) / (1 + exp(x)) )
}


##' Defines a simple cumulative distribution function (CDF) with
##' a point mass at nine (9) months.
##'
##' The cumulative distribution function determines an individual's
##' probabilty of receiving a routine vaccination.  It is passed
##' into the weighting calculation function (@@code{w.matrix}) by
##' various coverage analysis functions such as @@code{ ll.coverage }
##' and @@code{ mcmc.estimate }.  A user specifies, by name, which
##' CDF to use when an analysis function is invoked, e.g.
##' @@code{ mcmc.estimate( ..., routine.cdf = "9m.pointmass", ... ) }.
##' See the documentation for the respective analysis functions for
##' details about which CDFs are supported.
##'
##' @param age   the age of the individual
##' @return zero (0) if \code{ @@p age < 9 } and one (1) otherwise
##'
##' @@see w.matrix
##' @@see ll.coverage
##' @@see mcmc.estimate
##'
##' @examples
##'    x <- 0:24
##'    plot( x=x, y=nine.month.pointmass(x) )
##' 
##' @export
##'
nine.month.pointmass <- function ( age ) {
  return ( ifelse( age < 9, 0, 1 ) )
}


##' Defines a simple cumulative distribution function (CDF) with
##' a constant "hazard" rate @@p lambda of routine vaccination
##' starting at age 8.5 months.
##'
##' The cumulative distribution function determines an individual's
##' probabilty of receiving a routine vaccination.  It is passed
##' into the weighting calculation function (@@code{w.matrix}) by
##' various coverage analysis functions such as @@code{ ll.coverage }
##' and @@code{ mcmc.estimate }.  A user specifies, by name, which
##' CDF to use when an analysis function is invoked, e.g.
##' @@code{ mcmc.estimate( ..., routine.cdf = "const.post.9m", ... ) }.
##' See the documentation for the respective analysis functions for
##' details about which CDFs are supported.
##'
##' @param age      the age of the individual
##' @param lambda   the routine vaccination "hazard" rate
##' @return \code{ pexp( @@p age - 8.5, @@p lambda ) }
##'
##' @@see w.matrix
##' @@see ll.coverage
##' @@see mcmc.estimate
##'
##' @examples
##'    x <- 0:24 ; lambda <- 0.5
##'    plot( x=x, y=constant.post.9mo(x,lambda) )
##'
##'    x <- 0:24 ; lambda <- seq( 1.0, 0.1, -0.2 )
##'    plot( type="n", x=x, xlab="age", ylab="constant.post.9mo(age,lambda)",
##'          xlim=c(min(x),max(x)), ylim=c(0,1),
##'          xaxp=c(min(x),max(x),(max(x)-min(x))/2) )
##'    legend( "bottomright", "lambda", lambda, inset=0.1,
##'            fill=c(1:length(lambda)) )
##'    for ( i in 1:length(lambda) ) {
##'      lines( x=x, y=constant.post.9mo(x,lambda[i]), col=i )
##'    }
##' 
##' @export
##'
constant.post.9mo <- function ( age, lambda ) {
  return ( pexp( age - 8.5, lambda ) )
}


##' Calculates the vaccination probability for individuals.
##'
##' The @@c g function is used to calculate the probability that an individual
##' has been vaccinated during a set of vaccination activities.  The vaccination
##' activities are described by inputs @@p v and @@p N; @@p v specifies the number
##' of doses distributed and @@p N specifies the target population size for each
##' activity.  An individual's probability is determined by inputs @@p z and @@p w;
##' @@p z indicates which activities the individual is eligible for and @@p w
##' provides the pseudo-campaign weighting that specifies the probability of
##' child receiving a routine vaccination each year.  Finally, the scalar
##' parameters, @@p alpha and @@p rho, quantify campaign efficiency and the size
##' of accessible population, respectively.
##'
##' @note The inefficiency measure @@p alpha is defined as the natural log of psi,
##'       i.e. @@code{psi = e^alpha}.  @@code{psi = 0} (@@code{alpha = -inf})
##'       represents perfect efficiency, i.e. when every dose results in an
##'       additional vaccinee.  @@code{psi = 1} (@@code{alpha = 0}) represents
##'       random efficiency, i.e. when probability of any dose resulting in an
##'       additional vaccinee is independent of receiving a dose previously
##'       during the same campaign.
##'
##' @param z      the eligibility vector for the individual, i.e. @@code{z[j] = 1}
##'               if this individual is eligible for campaign @@c j, otherwise
##'               @@code{z[j] = 0}
##' @param w      the weight vector for the individual, i.e. @@code{ w[k] } gives
##'               the probability of the individual receiving the routine
##'               vaccination in year @@c k
##' @param v      the doses vector providing the number of vaccine doses nominally
##'               distributed in each campaign
##' @param N      the population size vector specifying the number of people in
##'               the target age range for each campaign
##' @param alpha  the vaccination efficiency parameter (scalar) for the campaign;
##'               @@code{alpha = -inf} implies perfect efficiency,
##'               @@code{alpha = 0} implies random efficiency; see note above
##' @param rho    the proportion of the population (scalar) that can be vaccinated,
##'               i.e @@code{1 - rho} is the @@e inaccessible portion of the
##'               general population.
##' @param log    a flag indicating whether the log probability should be
##'               returned; the default is @@c FALSE
##' @return the probability of vaccination for an individual with eligibility
##'         vector @@p z and weight vector @@p w during the campaigns described
##'         by @@p v and @@p N
##'
##' @examples
##'   g( z = c(0,1), w = c(0.5,0.5), v = c(500,500), N = c(1000,1000), alpha = 0, rho = 0.9 )
##' 
##' @export
##'
g <- function ( z, w, v, N, alpha, rho, log = FALSE ) {

  # TBD: validate inputs, e.g. values of z should be 0 or 1, length of
  #      vectors should match, rho in [0,1], etc.
  # TBD: Check: if v > rho * N, then g < 0 ?

  # Equation 3 from paper
  if ( alpha == 0 ) {
    tmp <- exp( -v/(N*rho) )
  } else {
    tmp <- ( 1 - v/(N*rho) * (1-exp(alpha)) )^( 1 / (1-exp(alpha)) )
  }

  # Equation 6 from paper
  # need 'na.rm=TRUE' (ie ignore missing values) to keep things
  # from blowing up in weird situations where weight is 0
  f.pseudo <- sum( w * tmp, na.rm=TRUE ) + ( 1 - sum( w, na.rm=TRUE ) )

  # Equation 4 from paper
  rc <- 1 - ( (1-rho) + rho * prod(tmp^z) * f.pseudo )

  if ( log ) {
    return( log(rc) )
  }
  return(rc)

  # TBD: Check that return value always in [0,1]

}  # end of function g(z,w,v,N,alpha,rho,log=FALSE)


##' Calculates the expected number of vaccinations for individuals.
##'
##' The @@c E.vacc function complements the probability calculation @@c g and
##' computes the expected vaccination count that an individual would have
##' received during a specified set of vaccination activities.  The vaccination
##' activities are described by inputs @@p v and @@p N; @@p v provides the number
##' of doses distributed and @@p N provides the target population size for each
##' activity.  An individual's probability is determined by inputs @@p z and @@p w;
##' @@p z indicates which activities the individual is eligible for and @@p w
##' provides the pseudo-campaign weighting that specifies the probability of
##' child receiving a routine vaccination each year.  Finally, the scalar
##' parameters, @@p alpha and @@p rho, quantify campaign efficiency and the size
##' of accessible population, respectively.
##'
##' @note The inefficiency measure @@p alpha is defined as the natural log of psi,
##'       i.e. @@code{psi = e^alpha}.  @@code{psi = 0} (@@code{alpha = -inf})
##'       represents perfect efficiency, i.e. when every dose results in an
##'       additional vaccinee.  @@code{psi = 1} (@@code{alpha = 0}) represents
##'       random efficiency, i.e. when probability of any dose resulting in an
##'       additional vaccinee is independent of receiving a dose previously
##'       during the same campaign.
##'
##' @param z      the eligibility vector for the individual, i.e. @@code{z[j] = 1}
##'               if this individual is eligible for campaign @@c j, otherwise
##'               @@code{z[j] = 0}
##' @param w      the weight vector for the individual, i.e. @@code{ w[k] } gives
##'               the probability of the individual receiving the routine
##'               vaccination in year @@c k
##' @param v      the doses vector providing the number of vaccine doses nominally
##'               distributed in each campaign
##' @param N      the population size vector specifying the number of people in
##'               the target age range for each campaign
##' @param alpha  the vaccination efficiency parameter (scalar) for the campaign;
##'               @@code{alpha = -inf} implies perfect efficiency,
##'               @@code{alpha = 0} implies random efficiency; see note above
##' @param rho    the proportion of the population (scalar) that can be vaccinated,
##'               i.e @@code{1 - rho} is the @@e inaccessible portion of the
##'               general population.
##' @param log    a flag indicating whether the log probability should be
##'               returned; the default is @@c FALSE
##'               <i><b>NOTE: the @@c log parameter is @@b not currently used.</b></i>
##' @param cond.on.accessible   a flag indicating whether the return value
##'               should be conditioned on the individual's inclusion in the
##'               accessible population, i.e. @@p rho; the default is @@c FALSE
##' @return the expected number of vaccinations for an individual with
##'         eligibility vector @@p z and weight vector @@p w during the
##'         campaigns described by @@p v and @@p N
##'
##' @examples
##'   E.vacc( z = c(0,1), w = c(0.5,0.5), v = c(500,500), N = c(1000,1000), alpha = 0, rho = 0.9 )
##' 
##' @export
##'
E.vacc <- function ( z, w, v, N, alpha, rho, log = FALSE,
                     cond.on.accessible = FALSE )
{
  # Equation 3 from paper
  # bit from the SIAs
  if ( alpha == 0 ) {
    tmp <- exp( -v/(N*rho) )
  } else {
    tmp <- ( 1 - v/(N*rho) * (1-exp(alpha)) )^( 1 / (1-exp(alpha)) )
  }

  # Equation 6 from paper
  # bit from the pseudo campaign
  f.pseudo <- sum( w * tmp ) + ( 1 - sum(w) )
  # TBD: NAs not removed in these sums like is done in 'g'???

  rc <- sum( (1-tmp) * z ) + ( 1 - f.pseudo )

  if ( ! cond.on.accessible ) {
    rc <- rho * rc
  }

  return ( rc )

  # TBD: Check that return value always in [0,Inf]

}  # end of function E.vacc(z,w,v,N,alpha,rho,log=FALSE,cond.on.accessible=FALSE)



##' Calculates the vaccination probability for all individuals in a group
##' of observations and under a described set of vaccination activities
##' (e.g. campaigns).
##'
##' The @@c get.prob.vacc function simply applies the @@code{ g } function
##' to all the entries in observations data frame, @@p obs.  The @@p obs
##' argument should contain the immunization information as derived, for
##' example, from a Demographic and Health Survey (DHS).  The @@c get.prob.vacc
##' function uses the observations and vaccination activity descriptions,
##' @@p camps, to calculate the eligibility matrix, @@code{ z.matrix } and
##' the weight matrix @@code{ w.matrix }, if not provided.  Then the
##' vaccination probabilities are calculated by calling @@c g with each
##' individual's eligibility and weight vector as well as the relevant
##' vaccination activity information (i.e. @@p camps$N and @@p camps$v),
##' the efficiency parameter (@@p alpha) and the accessibility multiplier
##' (@@c rho).
##'
##' @note The @@p trans.rho argument represents the @@link{logit} of @@c rho,
##'       i.e. @@code{ trans.rho = log( rho / (1-rho) ) }, where @@c rho is
##'       the proportion of population that is accessible to any vaccination
##'       activity and @@code{ 1-rho } the inaccessible portion.  It follows
##'       then that @@c rho is calculated as the inverse-logit (or logistic
##'       function) of @@p trans.rho, that is:
##'       @@code{ rho <- exp( trans.rho ) / ( 1 + exp(trans.rho) ) }
##'
##' @note The inefficiency measure @@p alpha is defined as the natural log of psi,
##'       i.e. @@code{psi = e^alpha}.  @@code{psi = 0} (@@code{alpha = -inf})
##'       represents perfect efficiency, i.e. when every dose results in an
##'       additional vaccinee.  @@code{psi = 1} (@@code{alpha = 0}) represents
##'       random efficiency, i.e. when probability of any dose resulting in an
##'       additional vaccinee is independent of receiving a dose previously
##'       during the same campaign.
##'
##' @param obs    a data frame with one row per observation and columns for
##'               @@c date of observation and @@c age of individual on that date;
##'               <i>Note: all ages should be in months</i>
##' @param camps  a data frame with one row per activity and columns providing
##'               the activity's @@c date, targeted age range (@@c age.low and
##'               @@c age.high), size of targeted population (@@c N) and the
##'               number of vaccine doses nominally distributed (@@c v)
##' @param trans.rho    the log-odds (aka logit) of @@c rho, the proportion
##'               of the population (scalar) that can be vaccinated;
##'               see note above
##' @param alpha  the vaccination efficiency parameter (scalar) for the campaign;
##'               @@code{alpha = -inf} implies perfect efficiency,
##'               @@code{alpha = 0} implies random efficiency;
##'               see note above
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
##' @return a vector containing the probability of vaccination for each
##'         individual in observations data frame, @@p obs.
##'
##' @@see g
##' @@see z.matrix
##' @@see w.matrix
##' @@see nine.month.pointmass
##' @@see constant.post.9mo
##' @export
##'
get.prob.vacc <- function ( obs, camps, trans.rho, alpha,
                            z=NULL, w=NULL, cdf.fun=NULL, ... )
{
  if ( is.null(z) ) {
    z <- z.matrix( obs, camps )
  }

  rho <- exp( trans.rho ) / ( 1 + exp(trans.rho) )

  if ( is.null(w) ) {
    w <- w.matrix( obs, camps, cdf.fun, ... )
  }

  p.vac <- sapply( 1:nrow(obs), function(i) {
    g( z[i,], w[i,], camps$v, camps$N, alpha, rho )
  } )

  return ( p.vac )

}  # end of function get.prob.vacc(obs,camps,trans.rho,alpha,z=NULL,w=NULL,cdf.fun=NULL,...)



##' Calculates the expected number of vaccinations for all individuals
##' in a group of observations and under a described set of vaccination
##' activities (e.g. campaigns).
##'
##' The @@c get.E.vacc function simply applies the @@code{ E.vacc } function
##' to all the entries in observations data frame, @@p obs.  The @@p obs
##' argument should contain the immunization information as derived, for
##' example, from a Demographic and Health Survey (DHS).  The @@c get.E.vacc
##' function uses the observations and vaccination activity descriptions,
##' @@p camps, to calculate the eligibility matrix, @@code{ z.matrix } and
##' the weight matrix @@code{ w.matrix } if not provided.  Then the expected
##' vaccination doses are calculated by calling @@c E.vacc with each
##' individual's eligibility and weight vector as well as the relevant
##' vaccination activity information (i.e. @@p camps$N and @@p camps$v),
##' the efficiency parameter (@@p alpha), the accessibility multiplier
##' (@@c rho) and conditionality flag (@@p cond.on.accessible).
##'
##' @note The @@p trans.rho argument represents the @@link{logit} of @@c rho,
##'       i.e. @@code{ trans.rho = log( rho / (1-rho) ) }, where @@c rho is
##'       the proportion of population that is accessible to any vaccination
##'       activity and @@code{ 1-rho } the inaccessible portion.  It follows
##'       then that @@c rho is calculated as the inverse-logit (or logistic
##'       function) of @@p trans.rho, that is:
##'       @@code{ rho <- exp( trans.rho ) / ( 1 + exp(trans.rho) ) }
##'
##' @note The inefficiency measure @@p alpha is defined as the natural log of psi,
##'       i.e. @@code{psi = e^alpha}.  @@code{psi = 0} (@@code{alpha = -inf})
##'       represents perfect efficiency, i.e. when every dose results in an
##'       additional vaccinee.  @@code{psi = 1} (@@code{alpha = 0}) represents
##'       random efficiency, i.e. when probability of any dose resulting in an
##'       additional vaccinee is independent of receiving a dose previously
##'       during the same campaign.
##'
##' @param obs    a data frame with one row per observation and columns for
##'               @@c date of observation and @@c age of individual on that date;
##'               <i>Note: all ages should be in months</i>
##' @param camps  a data frame with one row per activity and columns providing
##'               the activity's @@c date, targeted age range (@@c age.low and
##'               @@c age.high), size of targeted population (@@c N) and the
##'               number of vaccine doses nominally distributed (@@c v)
##' @param trans.rho    the log-odds (aka logit) of @@c rho, the proportion
##'               of the population (scalar) that can be vaccinated;
##'               see note above
##' @param alpha  the vaccination efficiency parameter (scalar) for the campaign;
##'               @@code{alpha = -inf} implies perfect efficiency,
##'               @@code{alpha = 0} implies random efficiency;
##'               see note above
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
##' @param cond.on.accessible   a flag indicating whether the return values
##'               should be conditioned on the individual's inclusion in the
##'               accessible population, i.e. @@p rho; the default is @@c FALSE
##' @param ...    any additional arguments to the @@p cdf.fun, e.g. @@c lambda
##'               <i>Note: this parameter is only required if the weight
##'               matrix, @@p w, is not provided.</i>
##' @return a vector containing the expected number of vaccinations for each
##'         individual in observations data frame, @@p obs.
##'
##' @@see E.vacc
##' @@see z.matrix
##' @@see w.matrix
##' @@see nine.month.pointmass
##' @@see constant.post.9mo
##' @export
##'
get.E.vacc <- function ( obs, camps, trans.rho, alpha, z=NULL, w=NULL,
                         cdf.fun = NULL, cond.on.accessible=FALSE, ... )
{
  if ( is.null(z) ) {
    z <- z.matrix( obs, camps )
  }

  if ( is.null(w) ) {
    w <- w.matrix( obs, camps, cdf.fun, ... )
  }

  rho <- exp( trans.rho ) / ( 1 + exp(trans.rho) )

  e.vac <- sapply( 1:nrow(obs), function(i) {
    E.vacc( z[i,], w[i,],  camps$v, camps$N, alpha, rho,
            cond.on.accessible = cond.on.accessible )
  })

  return ( e.vac )

}  # end of function get.E.vacc(obs,camps,trans.rho,alpha,z=NULL,w=NULL,cdf.fun=NULL,cond.on.accessible=FALSE,...)



##' Creates an eligibility matrix, i.e. a @@c z matrix, from immunization
##' information for a set of individuals, e.g. observations derived from
##' a Demographic and Health Survey (DHS), and a set of vaccination
##' activity descriptions (e.g. campaigns).
##'
##' The @@c z.matrix function loops through all the provided observations
##' and campaigns.  For each pairing, it is first determined if the
##' campaign occurred before or after the observation (survey) date.  If
##' the vaccination campaign was after the observation date, then the
##' individual's eligibility status is set to zero (0).  Otherwise, if
##' the campaign was completed before survey, then the eligibility status
##' is based on age.  If the individual's age, at the time of the campaign,
##' was within age range, @@code{[age.low, age.high]}, of the immunization
##' activity's target population, then eligibility is set to one (1).  If
##' the individual's age was not within the target age range, then
##' eligibility is zero (0).
##'
##' @note The @@c z.matrix function only considers campaigns, i.e. Supplemental
##'       Immunization Activities (SIAs).  For non-SIA immunization activities,
##'       e.g. routine vaccinations, the eligibility status is always set to
##'       zero (0).
##'
##' @note All vaccination activities are assumed to have been completed on
##'       the date provided (i.e. @@p campaigns$date); separate start and
##'       end dates are not supported.
##'
##' @param obs    a data frame with one row per observation and columns for
##'               @@c date of observation and @@c age of individual on that date;
##'               <i>Note: all ages should be in months</i>
##' @param campaigns  a data frame with one row per activity and columns
##'               providing the activity's @@c date and targeted age range,
##'               i.e. @@c age.low and @@c age.high
##' @return a matrix with one row per observation and one column per campaign;
##'         matrix cells contain a @@c 1 if that individual was eligible for
##'         the campaign and @@c 0 otherwise
##' @export
##'
z.matrix <- function ( obs, campaigns ) {
  rc <- matrix( nrow = nrow(obs), ncol = nrow(campaigns) )

  for ( i in 1:nrow(obs) ) {
    for( j in 1:nrow(campaigns) ) {
      if ( campaigns$is.SIA[j] ) {
        t.since <- as.numeric(obs$date[i]) - as.numeric(campaigns$date[j])
        # cat( sprintf( "[DEBUG] obs$date[%d]: %s (as.numeric: %g), campaigns$date[%d]: %s (as.numeric: %g), t.since: %f\n",
        #               i, obs$date[i], as.numeric(obs$date[i]), j, campaigns$date[j], as.numeric(campaigns$date[j]), t.since ) )
        if( t.since > 0 ) {
          age.at.c <- obs$age[i] - t.since / 30.41667
          rc[i,j] <- (age.at.c >= campaigns$age.low[j]) * (age.at.c <= campaigns$age.high[j])
          #cat( sprintf( "[DEBUG] age.at.c: %f, campaigns$age.low[%d]: %d, campaigns$age.high[%d]: %d ==> eligibility: %d\n",
          #              age.at.c, j, campaigns$age.low[j], j, campaigns$age.high[j], rc[i,j] ) )
        } else {
          rc[i,j] <- 0
        }
      } else {
        # right now just set all non-true campaigns to 0...
        # there may be a more elegant way to do this
        rc[i,j] <- 0
      }
    }
  }

  return (rc)

}  # end of function z.matrix(obs,campaigns)



##' Creates a weight matrix from (1) immunization information for a set
##' of individuals, e.g. observations derived from a Demographic and
##' Health Survey (DHS), (2) a set of vaccination activity descriptions,
##' and (3) a cumulative distribution function (CDF).
##'
##' The @@c w.matrix function loops through all the provided activities
##' and determines weights for all observations.  These weights give the
##' probability of an individual having the "opportunity" to be vaccinated
##' during a year's routine activities.
##'
##' @note For all Supplemental Immunization Activities (SIAs), i.e. campaigns,
##'       the routine vaccination opportunity weight defaults to zero (0).
##'
##' For each routine vaccination activity, the opportunity weight is based
##' on the individual's age at the start of the activity and the length of
##' exposure to that year's activity (12 months most years, but truncated
##' in the survey year).  Specifically, the weight is calculated as:
##'    @@f[ cdf.fun( xij + lj ) - cdf.fun( xij ) @@f]
##' where @@c xij is the age of individual @@c i at the start of year @@c j
##' and @@c lj is the length of exposure to routine vaccination in year @@c j.
##'
##' For length of the exposure to the routine vaccination, @@c lj, there
##' are three possible cases:
##' @@li{ a. if survey date is before activity year, then the individual
##'          is not "exposed" to this routine vaccination and hence
##'          @@f$ lj = 0 months @@f$ }
##' @@li{ b. if survey date is during activity year, then the individual
##'          is "exposed" only up until the observation date and hence
##'          @@f$ lj = (obs$date - activities$date) @@f$ }
##' @@li{ c. if survey date is after activity year, then the individual
##'          is "exposed" to this routine vaccination for the entire year
##'          and hence @@f$ lj == 12 months @@f$ }
##'
##' @note For routine vaccinations, the date provided (i.e. @@p activities$date)
##'       is used as start of the activity year and hence is typically given
##'       as January 1st of that year, e.g. @@c 2012-01-01.
##'
##' @param obs        a data frame with one row per observation and columns for
##'                   date of observation and age of individual on that date;
##'                   <i>Note: all ages should be in months</i>
##' @param activities a data frame with one row per activity and columns
##'                   providing the activity's date and the @@c is.SIA flag
##'                   indicating whether the activity is a Supplemental
##'                   Immunization Activity (SIA) or routine vaccination.
##' @param cdf.fun    a cumulative distribution function (CDF) used to
##'                   calculate the probability of routine vaccination during
##'                   each activity as a function of age, e.g.
##'                   @@code{ constant.post.9mo }.
##' @param ...        additional arguments to the @@p cdf.fun, e.g. @@c lambda
##' @return a matrix with one row per observation and one column per activity;
##'         matrix cells contain an "opportunity" weight which is 0 for SIA
##'         campaigns and calculated as using the @@p cdf.fun for routine
##'         vaccination activities
##'
##' @@see nine.month.pointmass
##' @@see constant.post.9mo
##' @export
##'
w.matrix <- function ( obs, activities, cdf.fun, ... ) {
  rc <- matrix( nrow = nrow(obs), ncol = nrow(activities) )

  rc[, activities$is.SIA == 1 ] <- 0

  for ( i in 1:nrow(activities) ) {
    if( !activities$is.SIA[i] ) {
      t.since <- as.numeric(obs$date) - as.numeric(activities$date[i])
      age.at.c <- obs$age - t.since / 30.41667

      # t.at.risk == length of exposure, lj
      t.at.risk <- sapply( t.since/30.41667, min, 12 )
      t.at.risk[ t.at.risk < 0 ] <- 0   # obs$date before activities$date
      ## print( "---" )
      ## print( age.at.c )
      ## print( t.at.risk )
      ## print( "tmp:" )
      ## print( tmp )
      ## t.at.risk[ t.at.risk < 0 ] <- 0

      ## print( (cdf.fun(age.at.c+t.at.risk,...) -
      ##       cdf.fun(age.at.c,...))*tmp)
      ## print(1-cdf.fun(age.at.c+t.at.risk,...) +
      ##       cdf.fun(age.at.c,...))
      ## print("==========")
      ## print( (cdf.fun(age.at.c+t.at.risk,...) -
      ##         cdf.fun(age.at.c,...))*tmp+
      ##       (1-cdf.fun(age.at.c+t.at.risk,...) +
      ##        cdf.fun(age.at.c,...)))

      rc[,i] <- cdf.fun( age.at.c + t.at.risk, ... ) - cdf.fun( age.at.c, ... )
    }
  }

  return (rc)

}  # end of function w.matrix(obs,activities,cdf.fun,...)



##' Calculate the expected coverage of each campaign given
##' the trans.rho and alpha.
##'
##' Corresponds to equation 5 in paper
##'
##' @param camps    TBD
##' @param trans.rho  TBD
##' @param alpha    TBD
##' @return TBD
##' @export
##'
calc.camp.coverage <- function ( camps, trans.rho, alpha ) {
  rho <- exp(trans.rho) / ( 1 + exp(trans.rho) )

  z <- diag( 1, nrow(camps) )
  w <- matrix( 1, nrow=nrow(camps), ncol=nrow(camps) )

  rc <- sapply( 1:nrow(camps), function(i) {
    g( z[i,], w[i,], camps$v, camps$N, alpha, rho )
  })

  return ( rc )

}  # end of function calc.camp.coverage(camps,trans.rho,alpha)




##' Gets the log likelihood of a set of observations.
##'
##' @param rho      the rho parameter
##' @param alpha    the alpha parameter
##' @param obs      the observations. Must have columns y, date, age
##' @param campaigns  campaign definitions. Must have columns date, v, N,
##'                 age.low, age.high
##' @param z        the z matrix. calculated if null
##' @param w        the w matrix. calculated if null
##' @param cdf.fun  the cdf function used to calculate the w.matrix,
##'                 not needed if w is provided.
##'                 This is large F in paper eq. 7,8
##' @param ...      additional arguments to cdf.fun
##' @return the log likelihood of the observations given the parameters
##' @export
##'
ll.coverage <- function ( rho, alpha, obs, campaigns, z=NULL, w=NULL,
                          cdf.fun=NULL, ... )
{

  # give sample weights if missing
  if ( is.null(obs$samp.wt) ) {
    obs$samp.wt <- 1
  }

  if ( is.null(z) ) {
    z <- z.matrix( obs, campaigns )
  }

  if ( is.null(w) ) {
    w <- w.matrix( obs, campaigns, cdf.fun,  ... )
  }

  p.vac <- sapply( 1:nrow(obs), function(i) {
    g( z[i,], w[i,], campaigns$v, campaigns$N, alpha, rho )
  })

  # paper equation 10 with weights
  ll <- sum( obs$samp.wt * log(p.vac^obs$y) ) +
        sum( obs$samp.wt * log((1-p.vac)^(1-obs$y)) )

  # ll <- sum( log(p.vac^obs$y) ) + sum( log((1-p.vac)^(1-obs$y)) )

  return ( ll )

}  # end of function ll.coverage(rho,alpha,obs,campaigns,z=NULL,w=NULL,cdf.fun=NULL,...)



####################
# Various functions for fitting  confidence intervals, etc.
####################


##' Log likelihood for using optim.
##'
##' Assumes weight matrix parameter is not being fit.
##'
##' @param theta    the parameter vector (trans.rho, alpha),
##'                 where rho = exp(trans.rho)/(1+exp(trans.rho))
##' @param obs      the observations
##' @param campaigns  the campaigns
##' @param z        the z matrix
##' @param cdf.fun  the cdf function
##' @param ...      additional arguments to cdf.fun
##' @return the log likelihood
##' @export
##'
ll.coverage.optim <- function ( theta, obs, campaigns, z=NULL,
                                cdf.fun=NULL, ... )
{
  rc <- -ll.coverage( exp(theta[1])/(1+exp(theta[1])), theta[2],
                       obs, campaigns, z, w=NULL, cdf.fun=cdf.fun, ... )
  # print( rc )
  return ( rc )
  # return ( -ll.coverage( exp(theta[1])/(1+exp(theta[1])),
  #                        theta[2], obs, campaigns, z ) )

}  # end of function ll.coverage.optim(theta,obs,campaigns,z=NULL,cdf.fun=NULL,...)




##' Log likelihood for using optim.
##'
##' Assumes that there  there is a constant vaccination rate after 9 months
##' of age
##'
##' @param theta    the parameter vector (trans.rho, alpha, lambda),
##'                 where rho = exp(trans.rho)/(1+exp(trans.rho))
##' @param obs      the observations
##' @param campaigns  the campaigns
##' @param z        the z matrix
##' @return the log likelihood
##' @export
##'
ll.coverage.optim.const.post.9m <- function ( theta, obs, campaigns, z=NULL )
{
  cdf.fun <- constant.post.9mo

  card.wts  <- obs$samp.wt[    ! is.na(obs$age.at.vac) ]
  card.ages <- obs$age.at.vac[ ! is.na(obs$age.at.vac) ]
  card.wts  <- card.wts[  card.ages > 8.5 ]
  card.ages <- card.ages[ card.ages > 8.5 ] - 8.5

  # BUG ??? -- Should lambda = e^theta[3] or just theta[3] ???
  #            Is theta[3] suppose to be "trans.lambda" as it is
  #            in the mcmc.estimate function ???
  rc <- -ll.coverage( exp(theta[1])/(1+exp(theta[1])), theta[2],
                      obs, campaigns, z, w=NULL,
                      cdf.fun=cdf.fun, lambda=exp(theta[3]) )
  # cat( theta[3], "-", exp(theta[3]), "-",
  #      sum(dexp(card.ages,exp(theta[3]),log=TRUE)), "\n" )

  rc <- rc - sum( card.wts * dexp( card.ages, exp(theta[3]), log=TRUE ) )
  return ( rc )

}  # end of function ll.coverage.optim.const.post.9m(theta,obs,campaigns,z=NULL)



