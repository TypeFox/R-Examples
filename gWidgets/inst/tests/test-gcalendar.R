require(RUnit)

test.gcalendar <- function() {
  w <- gwindow()
  g <- ggroup(cont = w, horizontal=FALSE)

  date <- "2009-01-20"

  ## this fails in rJava -- can't set the initial date
  
  ## vanilla
  cal <- gcalendar(date, cont = g)
  checkEquals(svalue(cal), date)

  ## coerce.with
  cal <- gcalendar(date, coerce.with="as.Date", cont = g)
  checkEquals(svalue(cal), as.Date(date))

}

  
