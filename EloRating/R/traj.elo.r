# library(nigradom); library(zoo)
# data(PB); data(PBpres)
# av <- elo.seq(PB$Winner, PB$Loser, PB$Date, presence=PBpres, init="average")
# eloobject <- av
# (from=min(eloobject$stability$date))
# (to=max(eloobject$stability$date))
# ID <- "LL"


traj.elo <- function(eloobject, ID, from=min(eloobject$stability$date), to=max(eloobject$stability$date)){

  # check integrity of dates
  if(as.Date(to)<as.Date(from)) stop("'from' date is later than 'to' date")
  
  # get lines that correspond to date range
  DR <- seq(from=as.Date(eloobject$misc["minDate"]), to=as.Date(eloobject$misc["maxDate"]), by="day")
  if((as.Date(from) %in% DR & as.Date(to) %in% DR) == FALSE) stop("one of the dates is out of data range")
  DR <- which(DR==as.Date(from)):which(DR==as.Date(to))
  
  # check whether IDs are among individuals
  excl <- NULL
  for(i in ID) {
    if(!i %in% eloobject$allids) excl <- c(excl, i)
  }
  
  if(length(excl)>0) warning(paste("the following IDs do not occur in the data: ", paste(excl, collapse=", "), sep=""))
  #if(ID %in% eloobject$allids == FALSE) stop(paste("\"", ID, "\"", " is not among IDs!", sep=""))
  
  # create output object
  res <- data.frame(ID=ID, fromDate=as.Date(from), toDate=as.Date(to), slope=NA, Nobs=NA)
  
  for(i in 1:length(ID)) {
    if(ID[i] %in% eloobject$allids) {
      # extract ratings for ID during the date range
      traj <- eloobject$mat[DR, ID[i]]
    
      # how many data points 
      res$Nobs[i] <- length(na.omit(traj))
      
      # calculate slope (but only if there were actual observations...)
      if(res$Nobs[i]  > 1 ) res$slope[i] <- as.numeric(lm(traj ~ DR)$coefficients["DR"])
      if(res$Nobs[i] <= 1 ) warning(paste("no (or only one) observation for", ID[i], "during specified date range\n"))
    }
  }
  
  # create output object
  #res <- data.frame(ID=ID, fromDate=as.Date(from), toDate=as.Date(to), slope=slo, Nobs=N)

  return(res)
}


#elo.traj(el, "aa")

