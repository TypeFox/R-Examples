## This demo considers creating a sliced functional time series object

library(rainbow)

sfts(ts(as.numeric(ElNino$y), frequency = 12), xname = "Month", yname = "Sea surface temperature")

