if(getRversion() >= "2.15.1")  utils::globalVariables(c("anglers", "trueeffort"))

GetTotalValues <- structure(
function # Conduct a creel survey of a population of anglers at a site

  # ##############################################################################
  # File:  GetTotalValues.R
  ## author<< Steven H. Ranney
  ## Contact: \email{sranney@gw-env.com}
  # Created: 12/19/13  
  # Last Edited: 4/9/14 by SHR
  ##description<<This function uses the output from \code{MakeAnglers} to 
  ## conduct a bus-route or traditional access point creel survey of the 
  ## population of anglers from \code{MakeAnglers} and provide clerk-observed
  ## counts of anglers and their effort.
  # Returns: This function returns a dataFrame called 'tmp' of creel survey 
  # metrics sampled from the MakeAnglers function.  
  #
  # TODO: add RData for example
  # TODO: add testing section
  # ##############################################################################

  (ang = anglers, ##<<This argument renames the output from the \code{MakeAnglers}
                  ## function.
  teffort = trueeffort, ##<<This argument renames the output from the 
                        ## \code{MakeAnglers} function.
  nanglers = length(anglers$starttime), ##<<Defines the size of the angler population
                                        ## based upon the length of the \code{anglers}
                                        ## from the \code{MakeAnglers} function
  startTime = NULL, ##<< The start time of the creel clerk at this site 
  endTime = NULL, ##<< The end time of the creel clerk at this site
  waitTime = NULL, ##<< The wait time of the creel clerk at this site
  samplingProb = 1, ##<<The sampling probability for the survey.  The default is 
                   ## \code{1} but will need to be changed if the survey is conducted
                   ## during only half of the fishing day (\code{.5}) or over 
                   ## longer time periods (e.g., \code{9.5/12}, if the survey is
                   ## 9.5 hours long).
  meanCatchRate = NULL, ##<< The mean catch rate for the fishery.  
  ... ##<<Arguments to be passed to other functions
  ){
  
  ##details<<Catch rates are assigned to anglers based upon the Gamma distribution
  ## with a mean of \code{meanCatchRate}.
  lambda <- rgamma(nanglers, meanCatchRate, 1)
  
  #Calculate true total catch for all anglers
  totalcatch <- sum(ang$triplength * lambda)  
  
  ang$catch <- ang$triplength * lambda
  
  # Obtain a starting time for the surveyor
  
  ##details<<If both \code{endTime=NULL} and \code{waitTime=NULL} then \code{waitTime} 
  ## will be 0.5 (one-half hour).  If a value is passed to \code{endTime}, then 
  ## \code{waitTime} becomes \code{endTime - startTime}.
  
  #Provide a 'standard' wait time of .5 hours for the clerk
  if(is.null(waitTime) & is.null(endTime)){
  waitTime <- 0.5
  }
  
  if(!is.null(endTime)){
  waitTime <- endTime - startTime
  }
  
  ##details<<If \code{startTime=NULL}, then a \code{startTime} is generated from the 
  ## uniform distribution between \code{0} and \code{11.5} hours into the fishing day.
  if(is.null(startTime)){
    startTime <- runif(1, 0, 11.5)
  }
  
  ##details<<If \code{endTime=NULL}, then \code{endTime = startTime+waitTime}
  # how long into the fishing day did the creel clerk arrive?
  if(is.null(endTime)){
    endTime <- startTime+waitTime # how long into the fishing day did the creel clerk depart?
  }
   
  ##details<<Incomplete trip effort is observed two ways: 1) by counting anglers
  ## that were at the site for the entire time that the surveyor was at the site
  ## and 2) counting anglers that arrived after the surveyor arrived at the site
  ## but remained at the site after the surveyor left.  These anglers are counted
  ## and their effort calculated based upon surveyor \code{startTime} and 
  ## \code{endTime}.
  
  ################
  #Effort of anglers that were onsite for the duration of the time that the clerk
  # was onsite
  #How many anglers were at the site the entire time the creel surveyor was there?
  nAnglersEntireTime <- length(which(ang$starttime < startTime & ang$departureTime > endTime))
  entireTime <- which(ang$starttime < startTime & ang$departureTime > endTime)

  #how long were the anglers that arrived after the creel there before the clerk left?
  if(nAnglersEntireTime > 0){
    entireTimeSumEffort <- nAnglersEntireTime * (waitTime)
  } else {
    entireTimeSumEffort <- 0
  }
  
  ################
  #Effort of anglers that arrived after the clerk arrived and stayed beyond the 
  # clerk's wait time
  #how many anglers arrived while the clerk was on site?
  anglerArrivals <- length(which(ang$starttime > startTime & ang$starttime < endTime & ang$departureTime > endTime))
  arrivals <- which(ang$starttime > startTime & ang$starttime < endTime & ang$departureTime > endTime)

  #how long were the anglers that arrived after the creel there before the clerk left?
  if(anglerArrivals > 0){
    arrivalSumEffort <- sum(endTime - ang$starttime[arrivals])
  } else {
    arrivalSumEffort <- 0
  }
  
  ##details<<Completed trip effort is observed two ways: 1) by interviewing anglers 
  ## that left while the surveyor was at the site.  The surveyor can determine
  ## effort and catch.  2) by interviewing anglers that both arrived and departed 
  ## while the surveyor was on site.  When \code{waitTime} is short, these cases are
  ## are rare; however, when \code{waitTime} is long (e.g., all day), then these 
  ## cases are much more likely.
  
  ################
  #Completed trip information; i.e., anglers that LEFT while the creel clerk 
  # was on site
  #Did any anglers depart (complete their trips?) while the creel clerk was there
  #OR did any anglers both arrive AND depart while the clerk was on site?
  anglerDepartures <- length(which(ang$starttime < startTime & (startTime < ang$departureTime) & (ang$departureTime < endTime)))
  whichAnglerDepartures <- which(ang$starttime < startTime & (startTime < ang$departureTime) & (ang$departureTime < endTime))
  arrDep <- length(which(ang$starttime > startTime & ang$departureTime < endTime))
  whichArrDep <- which(ang$starttime > startTime & ang$departureTime < endTime)
  completedTrips <- c(which(ang$starttime < startTime & ang$departureTime < endTime & ang$departureTime > startTime), 
                          which(ang$starttime > startTime & ang$departureTime < endTime))
  
  if((anglerDepartures + arrDep) > 0){
    totalCompletedTripEffort <- sum(ang$triplength[completedTrips]/samplingProb)
    totalCompletedTripCatch <- sum(ang$catch[completedTrips]/samplingProb)
  } else {
    totalCompletedTripEffort <- 0
    totalCompletedTripCatch <- 0
  }
  
  #Convert tripLengths
  ang$triplength[entireTime] <- waitTime
  ang$triplength[arrivals] <- endTime - ang$starttime[arrivals]
  ang$triplength[whichAnglerDepartures] <- ang$departureTime[whichAnglerDepartures] - startTime
  ang$triplength[whichArrDep]
  
  ##details<<Trip lengths of observed trips (both incomplete and complete) are 
  ## scaled by the \code{samplingProb} value.  The \code{samplingProb} is used to estimate
  ## effort and catch.
  
  ##references<<Pollock, K. H., C. M. Jones, and T. L. Brown. 1994. Angler survey 
  ## methods and their applications in fisheries management. American Fisheries 
  ## Society, Special Publication 25, Bethesda, Maryland. 

  #Scale triplength based upon the sampling probability
  ang$tripLengthAdj <- ang$triplength/samplingProb
  
  observedTrips <- ang$tripLengthAdj[c(entireTime,  arrivals, whichAnglerDepartures, whichArrDep)]
  nObservedTrips <- length(observedTrips)
  totalObservedTripEffort <- sum(observedTrips)
  
 
  #Create dataFrame of anglers, effort
  tmp <- as.data.frame(cbind(nObservedTrips, totalObservedTripEffort,
                             sum(anglerDepartures, arrDep), totalCompletedTripEffort, 
                             totalCompletedTripCatch, 
                             startTime, waitTime, totalcatch, trueeffort, 
                             mean(lambda)))
  
  names(tmp) <- c("nObservedTrips", "totalObservedTripEffort", 
                 "nCompletedTrips", "totalCompletedTripEffort", 
                  "totalCompletedTripCatch", "startTime", "waitTime", 
                  "totalCatch", "trueEffort", "meanLambda")

  return(tmp)

  }, ex = function() {
  
  #The MakeAnglers() function must be run before GetTotalValues(); otherwise, there
  # will be no data from which GetTotalValues() can calculate values.
  
  #MakeAnglers(100)
  
  startTime = .001 #start of fishing day
  endTime = 12 #end of fishing day
  meanCatchRate = 0.1 #this will cause VERY few fish to be caught!
  
  #GetTotalValues()

  #MakeAnglers(100)
  
  startTime = .001 #start of fishing day
  endTime = 6 #halfway through the fishing day
  samplingProb = .5 #this needs to be .5 because we are sampling only 50% of the fishing day
  meanCatchRate = 0.1 #this will cause VERY few fish to be caught!
  
  #GetTotalValues()
  
  })
