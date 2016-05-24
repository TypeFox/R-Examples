bin.by.date <- function(dates, values, binwidth=7, ...){
	bindates <- structure((unclass(dates) %/% binwidth  + 0.5) * binwidth, class="Date")
	quantiles <- sapply(split(values, bindates), quantile, probs=c(0, .25, 0.5, .75, 1), na.rm=T)
	myData <- as.data.frame(t(quantiles));
	myData <- na.omit(myData)
	myData <- cbind(Date=as.Date(row.names(myData)), myData);
	row.names(myData) <- NULL;
	names(myData) <- c("Date", "Min", "Q1", "Median", "Q3", "Max");
	myData$Mean <- sapply(split(values, bindates), mean)
	return(myData);
}