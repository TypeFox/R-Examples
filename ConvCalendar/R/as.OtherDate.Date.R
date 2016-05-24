`as.OtherDate.Date` <-
function(x,calendar,...){
  calendars<-c("gregorian","julian","hebrew","islamic","frenchrev",
               "persian","chinese","modpersian")
  calendar<-as.integer(pmatch(calendar, calendars))
  if ((length(calendar)!=1) || (is.na(calendar)))
    stop("wrong calendar")
  
  jd<-as.integer(julian(x))
  n<-as.integer(length(jd))
  tmp<-.C("to_calendar", jd, n, calendar, day=integer(n), month=integer(n), year=integer(n))
  tmp<-as.data.frame(tmp[-(1:3)])
  attr(tmp,"calendar")<-calendars[calendar]
  class(tmp)<-"OtherDate"
  tmp
}

