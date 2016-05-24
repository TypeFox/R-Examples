MSbetAvg <-
function(dat, avgmat) {
  # Currently not suited for matrices
  sqs = ((avgmat$B-avgmat$G)%*%dat)^2
  msqs= sqs/avgmat$ns # in case unbalanced design
  MSbet = sum(msqs)/(avgmat$m-1)
  return(MSbet)
}
