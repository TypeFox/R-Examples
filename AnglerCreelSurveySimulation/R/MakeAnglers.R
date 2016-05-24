
MakeAnglers <- structure(
function # Create a population of anglers

  # ##############################################################################
  # File:  MakeAnglers.R
  ## author<< Steven H. Ranney
  ## Contact: \email{sranney@gw-env.com}
  # Created: 12/19/13/14
  # Last Edited: 4/8/14 by SHR
  ##description<<Creates a population of \code{nanglers} with trip length and
  ##fishing day length provided by the user.
  # Returns: This function returns a list called anglers that includes variables
  # starttime, triplength, and departureTime as well as a value called 'trueeffort,
  # the sum of all of the angler efforts.
  #
  # TODO: add RData for example
  # TODO: add testing section
  # ##############################################################################

  (nanglers = 100,##<<The number of anglers in the population
  meanTripLength = 3.88, ##<<The mean trip length to be used in the function. 
                        ## \code{3.88} is the default.  The default is from data
                        ## from the 2008 Lake Roosevelt Fishing Evaluation Program.
  fishingDayLength = 12 ##<<The fishing day length to be used in the function.
                        ## Anglers are not be allowed to be fishing past this 
                        ## day length.  The default here is set to 12 hours, which 
                        ## may not be a suitable day length for fisheries at higher
                        ## latitudes (i.e., sunrise-sunset is > 12 hours) or
                        ## during shorter seasons.
  ){

  anglers <- list() # The anglers location, start time, and trip length are
                    # stored in a list, as three seperate vectors, each of equal
                    # length equal to the number of anglers (nanglers)

  # Give all anglers a start time representing 1 hour into the fishing day
  # and limit their fishing day to fishingDayLength hours long
  
  ##details<<All trip lengths will be limited so that anglers have finished
  ## their fishing trip by the end of the fishing day.  The function uses a 
  ## \code{while} loop to ensure that the number of angles = \code{nanglers} 
  ## provided in the function argument.  \code{fishingDayLength} is passed to the 
  ## argument.  The default is set to 12 hours.
  
  i=1
  startTime=tripLength=departureTime=NULL
  while(i <= nanglers){

  ##details<<\code{starttimes} are assigned by the uniform distribution
    startTime.tmp <- c(runif(1, 0, fishingDayLength - 0.25))

  ##details<<\code{triplengths} are assigned by the gamma distribution where the default mean 
  ## value comes from the 2008 Lake Roosevelt Fisheries Evaluation Program data.
  
    tripLength.tmp <- rgamma(1, meanTripLength, scale = 1.3)  
    departureTime.tmp <- startTime.tmp+tripLength.tmp
  
    if(departureTime.tmp < fishingDayLength){
      i=i+1
      startTime <- c(startTime, startTime.tmp)
      tripLength <- c(tripLength, tripLength.tmp)
      departureTime <- c(departureTime, departureTime.tmp)
    }
  }
 
  anglers$starttime <- startTime
  anglers$triplength <- tripLength
  anglers$departureTime <- departureTime
  
  trueeffort <- sum(anglers$triplength)

  anglers <<- anglers
  trueeffort <<- trueeffort
  
#  out <- list(anglers, trueeffort)
  
#  names(out) <- c("anglers", "trueeffort")
  
#  return(out)
  
  }, ex = function() {

  MakeAnglers(100, meanTripLength = 4, fishingDayLength = 10)
  MakeAnglers(10000)
  
  })
