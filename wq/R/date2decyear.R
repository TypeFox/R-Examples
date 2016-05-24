date2decyear <- 
function(w) {
  old_ops <- options(digits = 8)
  on.exit(options(old_ops))
  posx <- as.POSIXlt(w)
  yr <- 1900 + posx$year
  dy <- posx$yday + 0.5
  yr + dy/ifelse(leapYear(yr), 366, 365)
}