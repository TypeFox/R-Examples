# # # hillmakeR example using "hflights" data package
 require(hflights)
 require(plyr)
# # *** Data cleanup and formating ***
# # Add leading 0's to arrival times
# hflights$ArrTime <- sprintf("%04d", hflights$ArrTime)
# hflights$DepTime <- sprintf("%04d", hflights$DepTime)
# 
# # paste together readable timestamps, and parse as POSIXlt
# hflights$timestampArr <- paste(hflights$Year, "/", hflights$Month, "/", hflights$DayofMonth, " ", hflights$ArrTime, sep = "")
# hflights$timestampDep <- paste(hflights$Year, "/", hflights$Month, "/", hflights$DayofMonth, " ", hflights$DepTime, sep = "")
# dateFormat = "%Y/%m/%d %H%M"
# Arrivals <- as.POSIXct(as.character(hflights$timestampArr), format = dateFormat)
# Departures <- as.POSIXct(as.character(hflights$timestampDep), format = dateFormat)
# 
# LOS <- as.numeric(Departures - Arrivals)
# Departures[which(LOS < 0)] <- Departures[which(LOS < 0)] + 60*60*24 #add one day
# 
# # *** hillmakeR example ***
# # Run hillmakeR occupancy pattern over data
# library(plyr)
# planeCount <- occupancyPatternFunc.processdata(startTimes=Arrivals, stopTimes=Departures, resolution="min", fillup=0.95)
# planeCount$hour <- as.POSIXlt(planeCount$timeSequence)$hour
# planeCount$wday <- as.POSIXlt(planeCount$timeSequence)$wday
# byHourOfDay <- ddply(planeCount, c("hour"), function(x) c(mean = mean(x$census), median = median(x$census), q90 = quantile(x$census, 0.9)))
# 
# byHourOfWeek <- ddply(planeCount, c("wday", "hour"), function(x) c(mean = mean(x$census), median = median(x$census), q90 = quantile(x$census, 0.9)))
# 
# byHourOfWeek$DayName <- "Sunday"
# byHourOfWeek[which(byHourOfWeek$wday == 1), "DayName"] <- "Monday"
# byHourOfWeek[which(byHourOfWeek$wday == 2), "DayName"] <- "Tuesday"
# byHourOfWeek[which(byHourOfWeek$wday == 3), "DayName"] <- "Wednesday"
# byHourOfWeek[which(byHourOfWeek$wday == 4), "DayName"] <- "Thursday"
# byHourOfWeek[which(byHourOfWeek$wday == 5), "DayName"] <- "Friday"
# byHourOfWeek[which(byHourOfWeek$wday == 6), "DayName"] <- "Saturday"
# 
# #ggplot(byHourOfWeek, aes(hour, median,colour=DayName)) + 
#   geom_line(aes(group = DayName)) + 
#   geom_point()
# 
# # Repeat by Carrier, wrapping ddply around call to hillmakeR functions
# inFrame <- data.frame(Arrivals, Departures, hflights$UniqueCarrier)
# inFrame <- subset(inFrame, )
# planeCountbyCarrier <- ddply(inFrame, "hflights.UniqueCarrier", function(x) occupancyPatternFunc.processdata(startTimes=x$Arrivals, stopTimes=x$Departures, resolution="min", fillup=0.95))
# planeCountbyCarrier$hour <- as.POSIXlt(planeCountbyCarrier$timeSequence)$hour
# planeCountbyCarrier$wday <- as.POSIXlt(planeCountbyCarrier$timeSequence)$wday
# # 
# # byHourOfWeekAndCarrier <- ddply(planeCountbyCarrier, c("hflights.UniqueCarrier","wday", "hour"), function(x) c(mean = mean(x$census), median = median(x$census), q90 = quantile(x$census, 0.9)))
# # byHourOfWeekAndCarrier$wdayhour <- paste(byHourOfWeekAndCarrier$wday, byHourOfWeekAndCarrier$hour)
# # ggplot(byHourOfWeekAndCarrier, aes(wdayhour, median,colour=hflights.UniqueCarrier)) + 
# #   geom_line(aes(group = hflights.UniqueCarrier), size =2)