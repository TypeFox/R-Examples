week.to.date <-
function(week) {
  as.Date(paste(week, "1", sep=""), "%YW%U%u") - 1
}
