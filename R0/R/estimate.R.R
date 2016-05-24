# Name   : est.R0
# Desc   : Key function, used to trigger R0 or R(t) estimations with any supported
#          method.
# Date   : 2011/11/09
# Author : Boelle, Obadia
###############################################################################


# Function declaration

estimate.R <- function#Estimate R0 for one incidence dataset using several methods
### Estimate R0 for one incidence dataset using several methods.
##details<< Currently, supported methods are Exponential Growth (EG), Maximum Likelihood (ML), Attack Rate (AR), Time-Dependant (TD), and Sequential Bayesian (SB). See references below.
##references<< \code{\link{est.R0.EG}}: Wallinga, J., and M. Lipsitch. "How Generation Intervals Shape the Relationship Between Growth Rates and Reproductive Numbers." Proceedings of the Royal Society B: Biological Sciences 274, no. 1609 (2007): 599.
##references<< \code{\link{est.R0.ML}}: White, L.F., J. Wallinga, L. Finelli, C. Reed, S. Riley, M. Lipsitch, and M. Pagano. "Estimation of the Reproductive Number and the Serial Interval in Early Phase of the 2009 Influenza A/H1N1 Pandemic in the USA." Influenza and Other Respiratory Viruses 3, no. 6 (2009): 267-276.
##references<< \code{\link{est.R0.AR}}: Dietz, K. "The Estimation of the Basic Reproduction Number for Infectious Diseases." Statistical Methods in Medical Research 2, no. 1 (March 1, 1993): 23-41.
##references<< \code{\link{est.R0.TD}}: Wallinga, J., and P. Teunis. "Different Epidemic Curves for Severe Acute Respiratory Syndrome Reveal Similar Impacts of Control Measures." American Journal of Epidemiology 160, no. 6 (2004): 509.
##references<< \code{\link{est.R0.SB}}: Bettencourt, L.M.A., and R.M. Ribeiro. "Real Time Bayesian Estimation of the Epidemic Potential of Emerging Infectious Diseases." PLoS One 3, no. 5 (2008): e2185.

(epid=NULL, ##<< Name of epidemic dataset
 GT=NULL,  ##<< Generation Time repartition function
 t=NULL, ##<< Date vector
 begin=NULL, ##<< Begin date for estimation. Can be an integer or a date (YYYY-mm-dd or YYYY/mm/dd)
 end=NULL,  ##<< End date for estimation. Can be an integer or a date (YYYY-mm-dd or YYYY/mm/dd)
 date.first.obs=NULL, ##<< Optional date of first observation, if t not specified
 time.step=1, ##<< Optional. If date of first observation is specified, number of day between each incidence observation
 AR=NULL, ##<< Attack rate as a percentage from total population
 pop.size=NULL, ##<< Population size in which the incident cases were observed. See more details in est.R0.AR documentation
 S0=1, ##<< Initial proportion of the population considered susceptible
 methods=NULL, ##<< List of methods to be used for R0 estimation/comparison. Must be provided as c("method 1", "method 2", ...)
 checked=TRUE, ##<< Internal flag used to check whether integrity checks were ran or not.
 ... ##<< Parameters passed to inner functions
) 


# Code

{
  #Check if at least one method was provided
  if (is.null(methods)) {
    stop("Please enter at least one supported method of estimatation ('EG', 'ML', 'AR', 'TD' or 'SB').")
  }
  
  # Checks on GT, begin and end are run
  parameters <- integrity.checks(epid, t, GT, begin, end, date.first.obs, time.step, AR, S0, methods)
  begin <- parameters$begin
  end <- parameters$end
   
  
  #List of results, currently empty. We will append results when computed.
  res<-list()
  
  #If user inputs an unsupported method, stop. 
  sel.met=pmatch(methods, c('EG','TD','ML','AR', 'SB'))
  if (any(is.na(sel.met))) {
    stop("Invalid 'methods' argument. Supported methods are 'EG',' ML', 'TD', 'AR', 'SB'.")
  }
  
  #'estimates' will contain results objects in a list, associated to their name
  estimates = as.list(rep(NA,length(methods)))
  names(estimates) = methods

  #AR needs arguments very different from other methods, so it is caled separately if required
  #Optional arguments for each method must be input when calling est.R0. They will be passed
  #to their respective method.
  for (met in 1:length(methods)) {
    if (methods[met] == "AR") {
      estimates[[met]] = do.call(paste("est.R0",methods[met],sep="."), args=list(incid=epid, AR=AR, pop.size=pop.size, S0=S0, checked=checked, ...))
    }
    else {
      if(is.null(epid)) {
        stop("Argument epid must be provided.")
      }
      estimates[[met]] = do.call(paste("est.R0",methods[met],sep="."), args=list(epid=epid, GT=GT, t=t, begin=begin, end=end, date.first.obs=date.first.obs, time.step=time.step, checked=checked, ...))
    }
  }
  
  #Gets every common return argument
  res$epid = epid
  res$GT = GT
  res$begin = begin
  res$end = end
  res$t = t
  res$estimates = estimates
  
  
  return(structure(res,class="R0.sR"))
  
  ### A list with components:
  ### \item{estimates}{List containing all results from called methods.}
  ### \item{epid}{Epidemic curve.}
  ### \item{GT}{Generation Time distribution function.}
  ### \item{t}{Date vector.}
  ### \item{begin}{Begin date for estimation.}
  ### \item{end}{End date for estimation.}
}
