#### this is to calculate the difference in seconds between to timestamps ####
#### arguments can be numeric or character ####
#### negative values means time1 < time2 ####

TimeDiff <- function(time1, time2) {
  ## formalize argumaent ##
  time1 <- as.numeric(time1)
  time2 <- as.numeric(time2)
  
  ## seconds ##
  sec1 = time1 %% 100
  sec2 = time2 %% 100
  if (sec1 > 60 | sec2 > 60) {
    stop('Wrong time given.')
  } 
  
  ## minutes ##
  min1 = (time1 %/% 100) %% 100
  min2 = (time2 %/% 100) %% 100
  if (min1 > 60 | min2 > 60) {
    stop('Wrong time given.')
  } 
  
  ## hours ##
  hr1 = (time1 %/% 10000) %% 100
  hr2 = (time2 %/% 10000) %% 100
  if (hr1 > 24 | hr2 > 24) {
    stop('Wrong time given.')
  } 
  
  ## calculate difference in seconds ##
  diffSec = (hr2 - hr1) * 3600 + 
            (min2 - min1) * 60 +
            (sec2 - sec1)
  
  return (diffSec)
}