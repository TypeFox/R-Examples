#============================================================================================
# Compute STI events
#============================================================================================
stiEvents <- function(values) {
  values <- as.numeric(values)
  sti.extremely.hot   <- sum(values >= 2.00, na.rm=T)
  sti.very.hot        <- sum(values < 2.00 & values >= 1.50, na.rm=T)
  sti.moderately.hot  <- sum(values < 1.50 & values >= 1.00, na.rm=T)
  sti.near.normal     <- sum(values < 1.00 & values > -1.00, na.rm=T)
  sti.moderately.cold <- sum(values <= -1.00 & values > -1.50, na.rm=T)
  sti.very.cold       <- sum(values <= -1.50 & values > -2.00, na.rm=T)
  sti.extremely.cold  <- sum(values <= -2.00, na.rm=T)
  r <- c(sti.extremely.hot, sti.very.hot, sti.moderately.hot, sti.near.normal, sti.moderately.cold, sti.very.cold, sti.extremely.cold)
  names(r) = c("Extremely hot","Very hot","Moderately hot","Near normal","Moderately cold","Very cold","Extremely cold") 
  if (sti.near.normal == 0) {
    return (NA)
  } else {
    return (r)
  }
}
