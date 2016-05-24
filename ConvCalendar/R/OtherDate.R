`OtherDate` <-
function(day,month,year,calendar){
  calendars<-c("gregorian","julian","hebrew","islamic","frenchrev",
               "persian","chinese","modpersian")
  calendar<-as.integer(pmatch(calendar, calendars))
  if ((length(calendar)!=1) || (is.na(calendar)))
    stop("wrong calendar")

  rval<-data.frame(day=day,month=month,year=year)
  attr(rval,"calendar")<-calendars[calendar]
  class(rval)<-"OtherDate"
  rval
}

