NonResidenceExtractId <-
function(sResidenceEventFile,sDistanceMatrix=NULL)
{
  if(length(unique(sResidenceEventFile$TRANSMITTERID)) > 1)
    stop('length(TRANSMITTERID) > 1. Unable to extract non-residences.')  

  # Non-residence must contain at least two residences for that transmitter
  if(length(unique(sResidenceEventFile[,5])) < 2)
  {
    # Write an empty table
    newnonresidencetable <- data.frame(STARTTIME=sResidenceEventFile[1,1],ENDTIME=sResidenceEventFile[1,2],
                                       NONRESIDENCEEVENT=0,TRANSMITTERID="",RECEIVERID1="",RECEIVERID2="",
                                       DURATION=0,DISTANCE=0,ROM=0)[NULL,]
  }else{
    # Run the non-residence extraction function
    RECEIVERID1 <- sResidenceEventFile[c(1:length(sResidenceEventFile[,5])-1),5]
    RECEIVERID2 <- sResidenceEventFile[c(2:length(sResidenceEventFile[,5])),5]
    STARTTIME <- sResidenceEventFile[c(1:length(sResidenceEventFile[,2])-1),2]
    ENDTIME <- sResidenceEventFile[c(2:length(sResidenceEventFile[,2])),1]
    TRANSMITTERID <- sResidenceEventFile[c(2:length(sResidenceEventFile[,2])),4]
    DURATION <- as.numeric(difftime(as.POSIXct(ENDTIME), as.POSIXct(STARTTIME), units = "secs"))
    
    # Remove data where animal moves between the same receivers/stations
    NONRESIDENCEEVENT <- ifelse(as.character(RECEIVERID1) == as.character(RECEIVERID2),NA,1)
    nonresidencetable <- na.omit(data.frame(STARTTIME,ENDTIME,NONRESIDENCEEVENT,TRANSMITTERID,RECEIVERID1,
                                            RECEIVERID2,DURATION))
    if(is.null(sDistanceMatrix)!=TRUE)
    {
      DISTANCE <- ReturnVR2Distance(nonresidencetable,sDistanceMatrix)*1000 #Convert km into meters
      ROM <- DISTANCE/nonresidencetable$DURATION
    }else{
      DISTANCE <- 0
      ROM <- 0
    }
    
    newnonresidencetable <- data.frame(nonresidencetable,DISTANCE,ROM)
    
  }
  # Assign non-residence id to non-residence table
  if(nrow(newnonresidencetable) >= 1)  
    newnonresidencetable$NONRESIDENCEEVENT <- c(1:length(newnonresidencetable[,1]))
  
  return(newnonresidencetable)
}
