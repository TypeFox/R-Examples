if(getRversion() >= "2.15.1")  utils::globalVariables(c("anglers", "trueeffort"))

SimulateBusRoute <- structure(
function # Simulate a bus route survey

  # ##############################################################################
  # File:  SimulateBusRoute.R
  ## author<< Steven H. Ranney
  ## Contact: \email{sranney@gw-env.com}
  # Created: 12/19/13  
  # Last Edited: 9/9/14 by SHR
  ##description<<This function uses the output from \code{MakeAnglers} 
  ## and \code{GetTotalValues} to  conduct a bus-route or traditional access 
  ## point creel survey of the population of anglers from \code{MakeAnglers} 
  ## and provide clerk-observed counts of anglers and their effort.
  # Returns: Estimate catch (Ehat), the catch rate calculated by the ratio of means, 
  # the true, observed catch, and the actual catch rate (meanLambda).
  #
  # TODO: add RData for example
  # TODO: add testing section
  # ##############################################################################

  (startTime, ##<<The start time of the surveyor at each site.  This can be a vector of
              ## start times to simulate a bus route or one \code{startTime} to simulate
              ## a traditional access survey.
  waitTime, ##<<The wait time of the surveyor at each site.  This can be a vector
            ## of wait times to simulate a bus route or one \code{waitTime} to 
            ## simulate a traditional access survey.
  nanglers, ##<<The number of anglers at each site, either a vector or a single
            ## number.
  nsites, ##<<How many sites are being visited?
  samplingProb = 1, ##<<What is the sampling probability for the survey?  If 
                   ## all sites will be visited during the first or second half
                   ## of the fishing day, \code{samplingProb=0.5}.  If the
                   ## survey will take the entire fishing day, then 
                   ## \code{samplingProb=1}.
  meanCatchRate, ##<< The mean catch rate for the fishery.
  ... ##<<Arguments to be passed to other subfunctions, specifically to the 
      ## \code{\link{MakeAnglers}} function, including \code{meanTripLength} and 
      ## \code{fishingDayLength}.
  ){

  #Create a dataFrame to fill with the results
  dF <- as.data.frame(matrix(data = NA, nrow = nsites, ncol = 10, byrow=TRUE))
  names(dF) <- c("nObservedTrips", "totalObservedTripEffort", 
                "nCompletedTrips", "totalCompletedTripEffort", 
                "totalCompletedTripCatch", "startTime", "waitTime", 
                "totalCatch", "trueEffort", "meanLambda")


  #Run MakeAnglers() and GetTotalValues() iteratively for however many sites are 
  # provided in the nsites argument
  for(i in 1:nrow(dF)){
    MakeAnglers(nanglers=nanglers[i], ...)
#    MakeAnglers(nanglers=nanglers[i], meanTripLength, fishingDayLength)
    dF[i,] <- GetTotalValues(ang = anglers, teffort = trueeffort, nanglers = length(anglers$starttime), 
                             startTime = startTime[i], waitTime = waitTime[i], 
                             endTime = NULL, samplingProb, meanCatchRate, ...)
  }  
  
  ##seealso<<\code{\link{MakeAnglers}}
  ##seealso<<\code{\link{GetTotalValues}}
  
  bigT <- (startTime + waitTime)[length(startTime + waitTime)]-startTime[1]

  ##details<<Effort and catch are estimated from the the Bus Route Estimator 
  ## equation in Robson and Jones (1989), Jones and Robson (1991) and Pollock 
  ## et al. 1994.
  

  #########
  #Calculate effort based upon the bigT equation
  #Incomplete Effort
  sumEffort <- apply(data.frame(dF$totalObservedTripEffort), 1, sum)
  Ehat <- bigT*sum(1/dF$waitTime * sumEffort)
 
  #Complete Effort
  sumCompletedEffort <- dF$totalCompletedTripEffort
  completedEffort <- bigT*sum(1/dF$waitTime * sumCompletedEffort)

  ########
  #Complete catch
  #Calculate Catch based on the bigT equation
  sumCompletedCatch <- dF$totalCompletedTripCatch
  completedCatch <- bigT*sum(1/dF$waitTime * sumCompletedCatch)
  
  ##details<<Catch rate is calculated from the Ratio of Means equation (see Malvestuto (1996) and
  ## Jones and Pollock (2012) for discussions).
  
  ##references<<Jones, C. M., and D. Robson. 1991. Improving precision in angler surveys:
  ## traditional access design versus bus route design. American Fisheries Society
  ## Symposium 12:177-188.
  
  ##references<<Jones, C. M., and K. H. Pollock. 2012. Recreational survey 
  ## methods: estimation of effort, harvest, and released catch. Pages 883-919 
  ## in A. V. Zale, D. L. Parrish, and T. M. Sutton, editors. Fisheries 
  ## Techniques, 3rd edition. American Fisheries Society, Bethesda, Maryland.
 
  ##references<<Malvestuto, S. P. 1996. Sampling the recreational creel. Pages 
  ## 591-623 in B. R. Murphy and D. W. Willis, editors. Fisheries techniques, 
  ## 2nd edition. American Fisheries Society, Bethesda, Maryland.
  
  ##references<<Pollock, K. H., C. M. Jones, and T. L. Brown. 1994. Angler survey 
  ## methods and their applications in fisheries management. American Fisheries 
  ## Society, Special Publication 25, Bethesda, Maryland.
 
  ##references<<Robson, D., and C. M. Jones. 1989. The theoretical basis of an access 
  ## site angler survey design. Biometrics 45:83-98.
  
  #Total ROM catchRate
  catchRateROM <- completedCatch/completedEffort
  
  #trueTotalCatch
  trueCatch <- sum(dF$totalCatch)
  
  #totalTrueEffort
  trueEffort <- sum(dF$trueEffort)

  #meanLambda
  meanLambda <- mean(dF$meanLambda)
  
  dF <<- dF
  return(cbind(Ehat, catchRateROM, trueCatch, trueEffort, meanLambda)) 

  ##details<<The Ratio of means is calculated by  
  ##\deqn{
  ##\widehat{R_1} = \frac{\sum\limits_{i=1}^n{c_i/n}}{\sum\limits_{i=1}^n{L_i/n}}
  ##}
  ##
  ##where
  ##
  ##\emph{\eqn{c_i}} is the catch for the \emph{\eqn{i^{th}}} sampling unit 
  ##
  ##and
  ##
  ##\emph{\eqn{L_i}} is the  length of the fishing trip at the time of the interview. 
  ##For incomplete surveys, \emph{\eqn{L_i}} represents in incomplete trip.
  
  ##details<<The bus route estimator is 
  ##\deqn{
  ##\widehat{E} = T\sum\limits_{i=1}^n{\frac{1}{w_{i}}}\sum\limits_{j=1}^m{\frac{e_{ij}}{\pi_{j}}}
  ##}
  ##
  ##
  ## where
  ##
  ##\emph{E} = estimated total party-hours of effort;
  ##
  ##\emph{T} = total time to complete a full circuit of the route, including travelling and waiting;
  ##
  ##\emph{\eqn{w_i}} = waiting time at the \emph{\eqn{i^{th}}} site (where \emph{i} = 1, ..., \emph{n} sites); 
  ##
  ##\emph{\eqn{e_{ij}}} = total time that the \emph{\eqn{j^{th}}} car is parked at the \emph{\eqn{i^{th}}} site while the agent is at
  ##that site (where \emph{j} = 1, ..., \emph{n} sites).
  ##

  }, ex = function() {

  # To simulate one bus route survey that takes place in the morning, these values are used
  #start time at access sites
  startTimeAM <- c(1, 2,3,4,5) 
  #wait time at access sites
  waitTimeAM <- c(.5, .5, .5, .5, 2) 
  #the number of anglers that will visit access site throughout the day
  nanglersAM <- c(10,10,10,10,50) 
  # the number of sites to be visited
  nsitesAM <- 5
  # the sampling probability.  Here it is .5 because we are only conducting this 
  # survey during the first 50% of the fishing day
  samplingProb <- .5
  # the mean catch rate.  Here it is 2.5 which equals 2.5 fish/hour
  meanCatchRate <- 2.5
  
  SimulateBusRoute(startTimeAM, waitTimeAM, nanglersAM, nsitesAM, samplingProb, 
                   meanCatchRate)

  # To simulate one traditional access point survey where the creel clerk arrives, 
  # counts anglers, and interviews anglers that have completed their trips
  startTime = 0.001 
  waitTime = 8
  #nanglers can be informed by previously-collected data
  nanglers = 1000 
  nsites = 1
  # sampling probability here is 8/12 because we are staying at the access site
  # for 8 hours of a 12-hour fishing day.  To adjust the fishing day length, an
  # additional 'fishingDayLength' argument needs to be passed to this function.
  samplingProb <- (8/12)
  # the mean catch rate.
  meanCatchRate <- 5
  
  SimulateBusRoute(startTime, waitTime, nanglers, nsites, samplingProb, meanCatchRate)
    
  })
