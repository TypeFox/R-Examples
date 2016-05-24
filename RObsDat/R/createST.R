createST <- function(sp, siteID, location, timeobject, timelong, variables, thedata){ # timelong wird bisher nicht genutzt.
  
  # zur Erstellung eines leeren DataFrames wird Matrix benoetigt.
  # Spaltenanzahl ist Anzahl der unterschiedlichen Variablen, Zeilenanzahl Messzeiten * Orte
  dataValues <- as.data.frame(matrix(NA, ncol=(length(variables) + 2), nrow=length(timeobject) *length(sp)))

  # wiederhole die siteID solange bis es auf die Laenge der Time bzw. DataValues kommt
	time.col <- rep(timeobject, length(sp)) 
	site.col <- rep(siteID, each=length(timeobject))

  # the first both columns get the just created values.
  dataValues[,1] <- site.col
  dataValues[,2] <- time.col
  
  
  # generate the data.frame
  # the data-list is sort at location (this was created by loading the data from the db-file
for(colNumber in 1:length(variables)){
	for(rowNumber in 1:nrow(dataValues)){
  
		dataValues[rowNumber,colNumber+2] <- thedata[rowNumber+(nrow(dataValues)*(colNumber-1))]
	}
 }
  # the columns get a name  
  colname = c("sp", "time", variables)
  colnames(dataValues) <- colname
  
  # the data.frame will be sort at time.
		dataValues2 <- dataValues[order(dataValues$time),]
		
  # to create the actual data.frame, which will be used for the STFDF-Object the column site and time must be removed
  dataValues2 <- subset(dataValues2, select = -c(sp, time))


  # noch abfangen: if das Ding ne Endzeit hat, return (STFDF(sp, timeobject, dataValues, endtime))

  return(STFDF(sp, timeobject, dataValues2))  
} 
  
 
