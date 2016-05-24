`as.Date.OtherDate` <-
function(x,...){
  calendars<-c("gregorian","julian","hebrew","islamic","frenchrev",
               "persian","chinese","modpersian")
  calendar<-as.integer(pmatch(attr(x,"calendar"), calendars))
  if ((length(calendar)!=1) || (is.na(calendar)))
    stop("This can't happen")

  n<-as.integer(length(x$day))
  tmp<-.C("from_calendar", integer(n),n,calendar,as.integer(x$day), as.integer(x$month), as.integer(x$year))
  as.Date(tmp[[1]], origin=as.Date("1970-1-1"))
}

