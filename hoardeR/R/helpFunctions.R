timeToSec <- function(x){
  temp <- strsplit(x,":")
  hours <- as.numeric(temp[[1]][1])
  minutes <- as.numeric(temp[[1]][2])
  seconds <- as.numeric(temp[[1]][3])
  
  totalSec <- seconds + minutes * 60 + hours * 3600
  totalSec
}

secToTime <- function(x){
  hours <- x%/%3600
  x <- x-hours*3600
  minutes <- x%/%60
  x <- x-minutes*60  
  seconds <- floor(x)
  hours <- as.character(hours)
  minutes <- as.character(minutes)
  seconds <- as.character(seconds)
  if(nchar(hours)==1) hours <- paste("0",hours,sep="")
  if(nchar(minutes)==1) minutes <- paste("0",minutes,sep="")
  if(nchar(seconds)==1) seconds <- paste("0",seconds,sep="")
  
  paste(hours,minutes,seconds,sep=":")
}

timeStat <- function(x){
  timeInSec <- c()
  for(i in 1:length(x)){
    timeInSec[i] <- timeToSec(x[i])
  }
  avgTime <- round(sum(timeInSec)/length(timeInSec))
  secToTime(avgTime)
}
