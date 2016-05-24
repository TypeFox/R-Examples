date.to.week <-
function(date) {
  if (as.numeric(format(date+3, "%U")) == 0) {
    Recall(date-1)
  } else {
    format(date+3, "%YW%U")
  }
}
