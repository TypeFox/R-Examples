fireSeason <- function(input, fs.start = 12, fs.end = 5, method = "WF93"){
  #############################################################################
  # Description:
  #   Calculation of fire season. The intent of this function is to allow
  #   for calculation of the start and end dates of a fire season. Knowledge
  #   of the start and end dates allows for more accurate calculation of Fire
  #   Weather Indices. Only 1 method is presented here, however future 
  #   modifications will allow the inclusion of other methods.
  #
  #   The current method contained within are contained in the following 
  #   article:
  #
  #   Wotton B.M. and Flannigan, M.D. 1993. Length of the fire season in a 
  #     changing climate. Forestry Chronicle, 69: 187-192.
  
  #   **Note that snow depth/cover is not currently included due to a lack of
  #     widespread snow coverage data, however this will be a future addition.
  #
  # Args:
  #   input:    Single station weather data input stream - view fireSeason.Rd 
  #             documentation for full description.
  #   fs.start: Temperature threshold to start the fire season (deg celcius)
  #   fs.end:   Temperature threshold to end the fire season (deg celcius)
  #   method:   Fire Season calculation method. The default and only current
  #             method being used is "WF93".
  #   
  #   
  # Returns:
  #   seasonStartEnd: data.frame containing start and end dates of fire
  #                   seasons for a particular weather station.
  #
  #############################################################################
  names(input) <- tolower(names(input))
  #Ensure that required variables exist, if not, stop and give message
  # If variables exist, conver to local scope variables
  if ("yr" %in% names(input)) {
    yr <- input$yr
  } else {
    stop("Year was not provided, year is required for this function.")
  }
  if ("mon" %in% names(input)) {
    mon <- input$mon
    if(mon<1 || mon>12){
      stop("mon value is out of bounds (1-12), stopping execution.")
    }
  } else {
    stop("Month was not provided, month is required for this function.")
  }
  if ("day" %in% names(input)) {
    day <- input$day
    if(day < 1 || day > 31){
      stop("day value is out of bounds (1-31), stopping execution.")      
    }
  } else {
    stop("day was not provided, day is required for this function, stopping 
         execution..")
  }
  
  if ("tmax" %in% names(input)) {
    tmax <- input$tmax
  } else {
    stop("Maximum Daily Temperature (tmax) was not provided, daily tmax is 
         required for this function.")
  }
  method <- tolower(method)
  if(method != "wf93"){
    stop(paste("Selected method \"", method, "\" is unavailable, read 
               documentation for available methods.",sep=""))
  }
  
  #Get length of dataset
  n0<-length(tmax)
  seasonActive <- FALSE
  seasonStartEnd <- NULL
  
  #Run Wotton & Flannign 1993 fire season calculation method (without snow)
  if(method == "wf93"){
    for (k in 1:n0){
      if (k > 3){
        #If we are not already in a fire season AND
        #  the last 3 days are all greater than the starting temperature 
        #  threshold, then set start date.
        if(!seasonActive & (all(tmax[seq(k-3, k-1, 1)] > fs.start))){
          seasonActive <- TRUE #set season to active
          theday <- day[k]
          #If the data is starting us on January 4th, then we should have 
          #  started on January 1st Should be modified to work on any starting 
          #  data, but for now, just calendar year to allow for year-round 
          #  calculations.
          if(mon[k] == 1 && day[k] == 4){
            theday <- day[k-3]
          }
          #Bind the start date to the data.frame
          seasonStartEnd <- rbind(seasonStartEnd, 
                                  data.frame(yr = yr[k], 
                                             mon = mon[k],
                                             day = theday,
                                             fsdatetype = "start"))
        }
        #I we are already in a fire season AND
        #  the last 3 days are less than the end season temperature threshold
        #  thn set end date.
        if(seasonActive & all(tmax[seq(k-3, k-1, 1)] < fs.end)){
          seasonActive <- FALSE
          #Bind the end date to the data frame.
          seasonStartEnd <- rbind(seasonStartEnd, 
                                  data.frame(yr = yr[k],
                                             mon = mon[k],
                                             day = day[k],
                                             fsdatetype = "end"))
        }
      }
    }
  }
  #Return the season start end dates data.frame to caller
  return(seasonStartEnd)
}
