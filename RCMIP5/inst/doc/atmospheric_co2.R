## ----, echo=FALSE--------------------------------------------------------
df <- structure(list(Z = c(1e+05, 1e+05, 1e+05, 1e+05, 1e+05), time = c(2006, 
2007, 2008, 2009, 2010), value = c(380.717267354329, 382.683990478516, 
384.720212300618, 386.920120239258, 389.219032287598)), .Names = c("Z", 
"time", "value"), row.names = c(NA, -5L), class = "data.frame")
plot(df$time, df$value, type='b')

## ----, eval=FALSE, echo=FALSE--------------------------------------------
#  # Here's the short script that produced the outputs in this vignette.
#  # (Because it uses very large CMIP5 data files, we can't run it as
#  # part of the vignette itself.)
#  sink("sink.txt")
#  
#  library(RCMIP5)
#  
#  c5files<-getFileInfo()
#  
#  co2files <- subset(c5files, variable =="co2")
#  co2filechk <- checkTimePeriod(co2files)
#  print(co2filechk[-1])
#  
#  co2 <- loadCMIP5("co2", "CanESM2", "rcp85", verbose=T, yearRange = c(2006,2010))
#  print(summary(co2))
#  save(co2, file="co2")
#  
#  print(str(co2))
#  print(co2$Z)
#  
#  co2surf <- filterDimensions(co2, Zs=co2$Z[1], verbose=T)
#  print(co2surf)
#  
#  co2annual <- makeAnnualStat(co2surf, verbose=T)
#  co2summary <- makeGlobalStat(co2annual, verbose=T)
#  print(summary(co2summary))
#  print(as.data.frame(co2summary))
#  
#  write.csv(as.data.frame(co2summary), "co2_canesm2.csv")
#  
#  save(co2summary, file="co2summary")
#  
#  sink()

