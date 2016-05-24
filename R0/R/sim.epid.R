# Name   : sim.epid
# Desc   : Simulate epidemic outbreaks of specified R0 and generation time distribution
# Date   : 2012/04/11
# Author : Boelle, Obadia
###############################################################################


# Function declaration

sim.epid <- function#Epidemic outbreak simulation
### Generates several epidemic curves with specified distribution and reproduction number.

(epid.nb, ##<< Number of outbreaks to be generated.
 GT, ##<< Generation time distribution for the pathogen. Must be a R0.GT-class object.
 R0, ##<< Basic reproduction number.
 epid.length, ##<< Length of the epidemic.
 family, ##<< Distribution type for the new cases, either "poisson" or "negbin".
 negbin.size=NULL, ##<< Over-dispersion parameter, if family is set to "negbin".
 peak.value=50 ##<< Threashold value for incidence before epidemics begins decreasing
)
  
  
  # Code
  
{
  ##details<< This function is only used for simulation purposes. The output is a matrix of n columns (number of outbreaks)
  ## by m rows (maximum length of an outbreak).
  
  #Various content integrity checks
  if (class(GT) != "R0.GT") {
    stop("GT object must be of class R0.GT.")
  }
  
  #If negbin.size is omitted, size parameter is set so that variance = 10*R0. See details.
  if (family=="negbin" & is.null(negbin.size)) {
    negbin.size <- R0/4
  }
  GT <- GT$GT
  #R0.orig <- R0
  
  #Each epidemic is stored as a matrix column
  epidemics <- matrix(data=0, nrow=epid.length, ncol=epid.nb)
  
  #Loop for the required number of epidemic curves
  for (n in 1:epid.nb) {
    
    #Vector of cases with one index case
    sim.epid = c(1,rep(0,epid.length-1))
    
    #Loop on epidemic duration
    for(t in 1:epid.length) { 
      
      #New cases depend on family type
      if (family=="poisson") {
        new <- rpois(sim.epid[t], R0)
      }
      else if (family=="negbin") {
        ##details<< When using rnbinom with "mean" and "size" moments, the variance is given by mean + mean^2/size (see ?rnbinom).
        ## One should determine the size accordingly to the R0 value to increase the dispersion.
        ## From the previous variance formula, if Var(X) = k*R0, size = R0/(k-1)
        new <- rnbinom(sim.epid[t], size=negbin.size, mu=R0)
      }
      newd <-  rmultinom(1, sum(new), GT)[,1]
      sim.epid[t:(t+length(GT)-1)] <- sim.epid[t:(t+length(GT)-1)]+newd
      
      # Threshold so that epidemics eventualy dies out
      if (sim.epid[t+1] > peak.value & t < (epid.length-1)) {
        
        #Changing the R value is like implementing control measure. Uncomment if want to doing so
        #R0 <- 0.7
        
        #One the other hand, we can just stop the epidemic and start another one, as we are only interested in the growth period
        sim.epid[(t+2):epid.length] <- 0
        break
      }
    }
    sim.epid <- sim.epid[!is.na(sim.epid)]
    epidemics[,n] <- sim.epid
    #R0 <- R0.orig
  }
  
  return(epidemics)
}
