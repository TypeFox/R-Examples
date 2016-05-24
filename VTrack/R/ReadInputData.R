ReadInputData <-
function(infile,iHoursToAdd=0,fAATAMS=FALSE,fVemcoDualSensor=FALSE,dateformat=NULL,sVemcoFormat='Default')
{

  header.row <- c("DATETIME","TRANSMITTERID","SENSOR1","UNITS1","RECEIVERID","STATIONNAME")
  
  #Specify time/date format
  if(is.null(dateformat)==TRUE) 
    newDate <- as.POSIXct(strptime(infile[,1],format="%Y-%m-%d %H:%M:%S"))
  if(is.null(dateformat)==FALSE)
    newDate <- as.POSIXct(strptime(infile[,1],format=dateformat))
  
  # convert GMT to local times
  sDateTime <- newDate + (iHoursToAdd * 3600)
  
  # write the relevant fields to the output file discarding irrelevant fields
  if (fAATAMS==TRUE)
  {
    #  input data is in AATAMS format; SENSOR1 & UNITS1 are absent
    newdat <- data.frame(sDateTime,as.character(infile[,6]),"","",as.character(infile[,5]),as.character(infile[,2]))
    names(newdat) <- header.row   
  }else{
    # input data is in VEMCO export format 1.0
    if (sVemcoFormat == '1.0')
    {
      newdat <- data.frame(sDateTime,as.character(infile[,3]),infile[,4],infile[,5],as.character(infile[,11]),as.character(infile[,12]))
      names(newdat) <- header.row
    }
    if (sVemcoFormat == '1.0' & fVemcoDualSensor == TRUE)
    {
      # write a duplicate row for VEMCO dual sensor data; there is a value for Sensor 2 & Units 2
      dualsensor <- data.frame(sDateTime,as.character(infile[,3]),infile[,6],infile[,7],as.character(infile[,11]),as.character(infile[,12]))
      names(dualsensor) <- header.row
      newdat <- rbind(newdat,dualsensor)
    }
    # input data is in VEMCO export format Default 
    if (sVemcoFormat == 'Default')
    {
      newdat <- data.frame(sDateTime,as.character(infile[,3]),infile[,6],iconv(infile[,7], "UTF-8", "latin1"),as.character(infile[,2]),as.character(infile[,8]))
      names(newdat) <- header.row
    }
  }
    
  newdat1 <- unique(newdat[order(as.character(newdat$DATETIME)),1:6])
  
  # Replace blank STATIONNAME cells with "Unknown"
  newdat1[,6] <- as.character(newdat1[,6])
  addval <- which(is.na(newdat1[,6]))
  if(length(addval) > 0)
    newdat1[addval,6] <- "Unknown"
  return(newdat1)
}
