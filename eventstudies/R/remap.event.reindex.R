# A function which consumes a zoo object where there are lots of events
# (as columns)
# The contents are all levels
# For each column, the event date value is set to 100 and all other
# values are scaled accordingly.
remap.event.reindex <- function(z) {
  eventvals <- as.numeric(window(z, start=0, end=0))
  for (i in 1:ncol(z)) {
    z[,i] <- 100*z[,i]/eventvals[i]
  }
  z
}
