expq <-
function(x, q)
{
  if (q == 1) {
    return (exp(x))
  } else {
    Exponential <- (x*(1-q)+1)^(1/(1-q))
    if (q > 1) {
      Exponential[x > 1/(q-1)] <- NA
    }
    return (Exponential)
  }
}


expq.CommunityProfile <-
function(Profile)
{
  if (!is.CommunityProfile(Profile))
    stop("Profile must be a CommunityProfile")
  
  CP <- Profile
  CP$y <- sapply(1:length(CP$x), function(i) expq(CP$y[i], CP$x[i]))
  if (!is.null(CP$low))
    CP$low <- sapply(1:length(CP$x), function(i) expq(CP$low[i], CP$x[i]))
  if (!is.null(CP$high))
    CP$high <- sapply(1:length(CP$x), function(i) expq(CP$high[i], CP$x[i]))
  
  return (CP)
}
