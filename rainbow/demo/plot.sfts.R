## This demo considers functional data visualization using sfts().

library(rainbow)

plot(sfts(ts(as.numeric(ElNino$y), frequency = 12), xname = "Month", yname = "Sea surface temperature"))
   