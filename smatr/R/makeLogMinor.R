#returns vector of minor tick spacings approrpaite for log 10 scaled axis with major ticks given by 'major'.
makeLogMinor<-function(major){	
	temp <- NULL
	if(length(major) > 1) 
		for (i in 1:(length(major)-1))
			temp<-c(temp, seq(major[i], major[i+1]-major[i], major[i]))
	
	temp <- c(temp, major[length(major)])
	return(temp)
} 
