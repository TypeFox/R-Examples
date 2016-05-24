summaryData <-
function(data, validCut=600, perMinuteCts=1, markingString = "w"){
	if(perMinuteCts==1){
		unit = "1 min"
		} else if(perMinuteCts==60){
			unit = "1 sec"			
			} else{
				epoch <- 60/perMinuteCts
				unit = paste(epoch, "sec")
			}
	validCut <- validCut*perMinuteCts
	data$weekend <- ifelse(data$weekday=="Saturday" | data$weekday=="Sunday", 1, 0)
	data$weekend <- factor(data$weekend, levels=c(0,1), labels=c("weekday", "weekend") )
	# total number of week and weekend days
	totalNumWeekWeekend <- tapply(data$weekend, data$weekend, length)/(1440*perMinuteCts)
	totalNumWeekWeekend[is.na(totalNumWeekWeekend)==1] <- 0
	# total number of days
	totalNumDays <- (length(data[,1])/(1440*perMinuteCts))
	totalNumDays[is.na(totalNumDays)==1] <- 0
	wearTime <- sumVct(data, markingString = markingString)
	wearTime$weekend <- ifelse(wearTime$weekday=="Saturday" | wearTime$weekday=="Sunday", 1, 0)
	wearTime$weekend <- factor(wearTime$weekend, levels=c(0,1), labels=c("weekday", "weekend") )
	wearTimeByDay <- tapply(wearTime$duration, wearTime$days, sum, na.rm=T)
	wearTimeByDay
	validWearTimeByDay <- wearTimeByDay[ifelse(wearTimeByDay > validCut, 1, 0)==1]
	validWearTimeByDay
	valid.days <- as.numeric(names(validWearTimeByDay))
	valid.days
	valid.wearTime <- wearTime[wearTime$days%in%valid.days,]
	valid.wearTime
	valid.data <- data[data$days%in%valid.days,]
	valid.data$weekend <- ifelse(valid.data$weekday=="Saturday" | valid.data$weekday=="Sunday", 1, 0)
	valid.data$weekend <- factor(valid.data$weekend, levels=c(0,1), labels=c("weekday", "weekend") )
	# total number of week and weekend days for valid days
	totalValidNumWeekWeekend <- tapply(valid.data$weekend, valid.data$weekend, length)/(1440*perMinuteCts)
	totalValidNumWeekWeekend[is.na(totalValidNumWeekWeekend)==1] <- 0
	# total number of days for valid days
	totalValidNumDays <- (length(valid.data[,1])/(1440*perMinuteCts))
	totalValidNumDays[is.na(totalValidNumDays)==1] <- 0
	meanWeartimeValidDays <- tapply(valid.wearTime$duration, valid.wearTime$weekend, sum, na.rm=T)/totalValidNumWeekWeekend
	meanWeartimeValidDays
	meanWeartimeValidDays[is.na(meanWeartimeValidDays)==1] <- 0
	meanWeartimeOverallValidDays <- sum(tapply(valid.wearTime$duration, valid.wearTime$weekend, sum, na.rm=T), na.rm=T)/totalValidNumDays
	meanWeartimeOverallValidDays
	meanWeartimeOverallValidDays[is.na(meanWeartimeOverallValidDays)==1] <- 0
	return(list(unit=unit, totalNumDays=totalNumDays,totalNumWeekWeekend=totalNumWeekWeekend, validCut=validCut, totalValidNumDays=totalValidNumDays, totalValidNumWeekWeekend=totalValidNumWeekWeekend, wearTimeByDay=wearTimeByDay, validWearTimeByDay=validWearTimeByDay, meanWeartimeValidDays=meanWeartimeValidDays, meanWeartimeOverallValidDays=meanWeartimeOverallValidDays) )
	}

