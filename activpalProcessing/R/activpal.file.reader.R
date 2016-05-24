activpal.file.reader <-
function(file.name.and.path)
{
  data <- read.csv(file.name.and.path, stringsAsFactors=FALSE)
  data <- data[,(1:6)]
  names(data) <- c("time","datacount","interval","activity","cumulativesteps","methrs")
	
	data$time <- sub("#","",data$time)
	data$time <- sub("#","",data$time)
	data[,2] <- as.numeric(as.character(data[,2]))
	data[,3] <- as.numeric(as.character(data[,3]))
	data[,4] <- as.numeric(as.character(data[,4]))
	data[,5] <- as.numeric(as.character(data[,5]))*2 #event files have half the actual number of steps for some reason
	data[,6] <- as.numeric(as.character(data[,6]))
	
		t <- dim(data)[1]
		
		data <- data[!(data[,"time"] == "1899-12-30"),]
		data <- data[!(data[,"time"] == "0"),]
		n <- dim(data)[1]		

		if(is.character(data$time)==TRUE&t==n)
		{
		data$time <- as.numeric(data$time)
		data$time <- as.POSIXct(as.Date(data$time,origin="1899-12-30"))
		data$time <- as.POSIXlt(data$time,tz="GMT")
		data$time <- strptime(data$time,format="%Y-%m-%d %H:%M:%S")
			}
		
	return(data)
}

