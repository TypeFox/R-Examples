stab.elo <-
function(eloobject, from=min(eloobject$stability$date), to=max(eloobject$stability$date), weight=TRUE) {
  
  # create errors if the date range lies outside the supplied data
  if(as.Date(from)<min(eloobject$stability$date) | as.Date(from)>max(eloobject$stability$date)) stop("start date outside date range")
  if(as.Date(to)<min(eloobject$stability$date) | as.Date(to)>max(eloobject$stability$date)) stop(paste(to, "outside date range"))
  
  # subset the stability dataframe in the eloobject according to the specified dates
  stab <- subset(eloobject$stability, date >= as.Date(from) & date <= as.Date(to))
  
  # calculate and return stability either with (default) or without the weighing factor
  if(weight==TRUE) {
    S <- round(sum(stab$rankdiffs * stab$eloweights)/(sum(floor(stab$Idspresent^2/2))), 4)
    return(1 - S)
  }else{
    S <- round(sum(stab$rankdiffs)/(sum(floor(stab$Idspresent^2/2))), 4)
    return(1 - S)
  }
}
