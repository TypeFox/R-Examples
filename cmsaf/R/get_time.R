get_time <- function(time.unit,time.step){
  
  stopifnot(is.character(time.unit))
  stopifnot(is.numeric(time.step) && !any(is.na(time.step)))
  
  # convert time.unit to POSIXct
  
  if (unlist(strsplit(time.unit," "))[2]=="since"){
    t.unit <- unlist(strsplit(time.unit," since "))[1]
    ref.date <- unlist(strsplit(time.unit," since "))[2]
  } else {stop("Time unit is not conform to CF-Convention!")}
  
  # check for Julian date
  julian <- FALSE
  if (ref.date=="-4712-01-01 12:00:00"){
    ref.date <- "1900-01-01 00:00:00"
    julian <- TRUE
  }
  ref.date <- as.POSIXct(ref.date,tz="UTC")
  
  # get factor to convert time.step to seconds
  # for months it is only an estimation for the average days per month
  
  factor <- 0
  
  # check reference time unit
  if (t.unit=="minutes"|t.unit=="Minutes"|t.unit=="Mins"|t.unit=="Min"|t.unit=="min"|t.unit=="mins")(factor <- 60)
  if (t.unit=="seconds"|t.unit=="Seconds"|t.unit=="Secs"|t.unit=="Sec"|t.unit=="sec"|t.unit=="secs")(factor <- 1)
  if (t.unit=="Hours"|t.unit=="Hour"|t.unit=="hour"|t.unit=="hours")(factor <- 60*60)
  if (t.unit=="Days"|t.unit=="Day"|t.unit=="day"|t.unit=="days")(factor <- 24*60*60)
  if (t.unit=="Weeks"|t.unit=="Week"|t.unit=="week"|t.unit=="weeks")(factor <- 7*24*60*60)
  if (t.unit=="Months"|t.unit=="Month"|t.unit=="month"|t.unit=="months")(factor <- 30.4375*24*60*60)
  if (factor==0)(stop(cat(paste("Error! Non-compliant time unit: ",t.unit,sep="")),"\n"))
  
  # calculate times
  
  if (julian){
    # 1582 conversion from julian to gregorian
    ts <- time.step-2415020.5
    if(min(ts)<0)(ts <- time.step-2414982.5)
    time.step <- ts
  }
  check <- ifelse((time.step*factor)<=.Machine$integer.max,FALSE,TRUE)
  if (sum(check)>0){
    ffactor <- factor/10
      times <- ref.date+(time.step*ffactor) 
    for (i in 1:9){
      times <- times+(time.step*ffactor)
    }
    check <- ifelse((time.step*ffactor)<=.Machine$integer.max,FALSE,TRUE) 
    if (sum(check)>0)(cat(paste("Some times exeed maximum integer value and may be wrong: ",times[check])))
  } else {
      times <- ref.date+(time.step*factor)
  }
  
  return(times)  
}