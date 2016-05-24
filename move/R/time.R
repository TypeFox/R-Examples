setGeneric("timeSummary", function(x, units="hours"){standardGeneric("timeSummary")})
setMethod("timeSummary", 
          signature=".MoveTrackSingle",
          definition=function(x, units){           
            date <- timestamps(x)
            TimeDiff <- timeLag(x, units=units) 
            df <- data.frame(Duration=difftime(date[length(date)], date[1], units=units)) 
            df$AverDur <- mean(TimeDiff)       #mean time difference between relocations
            df$SDDur <- sd(TimeDiff)           #standard deviation of time differences between relocations
            df$dupl <- any((TimeDiff)<(1/3600))            #check whether any two relocations are closer than a second to each other
            df$multseason <-  any(TimeDiff > (24*30))      #check whether any two subsequent relocations are more than one month apart
            return(df)
          })

setMethod("timeSummary", 
          signature=".MoveTrackStack", 
          definition=function(x, units){
            lst <- lapply(split(x), timeSummary, units=units)
            return(lst)
          })
