#Dealing with NAs

lfnacheck <- function(lfobj){
  lfcheck(lfobj)
  total <- sum(is.na(lfobj$flow))

  percentage <- total/length(lfobj$flow)

  year <- aggregate(is.na(flow)~hyear, data = lfobj, sum)

  dur <- rle(is.na(lfobj$flow))
  duration <- table(dur$length[dur$value])

  na <- list(total,percentage, year, duration)
  names(na) <- c("total","percentage", "hydrologicalyear", "duration")
  na
}


lfnainterpolate <- function(lfobj){
  lfcheck(lfobj)
  if(sum(is.na(lfobj$flow) > 0)){
    lfobj$flow <-approx(1:length(lfobj$flow),lfobj$flow,n=length(lfobj$flow))$y
  }
  lfobj
}
