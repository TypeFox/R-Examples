lnq <-
function(x, q)
{
  if (q == 1) {
    return (log(x))
  } else {
    Log <- (x^(1-q)-1)/(1-q)
    Log[x < 0] <- NA
    return (Log)
  }
}


lnq.CommunityProfile <-
function(Profile)
{
  if (!is.CommunityProfile(Profile))
    stop("Profile must be a CommunityProfile")
  
  CP <- Profile
  CP$y <- vapply(1:length(CP$x), function(i) lnq(CP$y[i], CP$x[i]), 0)
  if (!is.null(CP$low))
    CP$low <- vapply(1:length(CP$x), function(i) lnq(CP$low[i], CP$x[i]), 0)
  if (!is.null(CP$hi))
    CP$hi <- vapply(1:length(CP$x), function(i) lnq(CP$hi[i], CP$x[i]), 0)
  
  return (CP)
}
