RunSensorEventExtraction <- function(sInputFile,
                                     iEventType,            
                                     sLocation,               
                                     iSensor,               
                                     rTriggerThreshold,     
                                     iTimeThresholdStart,  
                                     iTimeThreshold,        
                                     rTerminationThreshold) 
{

  i <- NULL 
  
#  if(Sys.info()["sysname"]=="Darwin"){
#    registerDoMC()
#  }else{
    registerDoParallel(cl=(cl <- makeCluster(2)))
#  }
  
  # Extract sensor data
  sInputFile2 <- ExtractData(sInputFile, sQueryDataType = iSensor)

  # Check data are in chronological order and duplicates-free
  sInputFile2 <- unique(sInputFile2[order(as.character(sInputFile2$DATETIME)),]) 
 
  # Define if user wants to extract increasing or decreasing sensor events  
  if(iEventType=="DECREASE")
      sInputFile2$SENSOR1 <- -1 *(sInputFile2$SENSOR1)
  
  # Extract list of transmitter names from sensor filtered data frame
  TransmitterNames <- ExtractUniqueValues(sInputFile2,2)
  iTransmitterCount <- length(TransmitterNames)

  # Define whether working with receivers or stations 
  if(sLocation=="STATIONNAME")
    iLocationCol <- 6
  if(sLocation=="RECEIVERID")
    iLocationCol <- 5
  
  SensorExtractId <- function(sTransmitterId)
  {
    
    # Initialise variables
    iCount <- rTotalDepth <- iTotalDepth <- iSensorEvents <- 0
    ilogevent <- 1
    
    # Create empty output dataframe with header row
    event <- data.frame(STARTTIME=as.POSIXct("0001-01-01"),ENDTIME=as.POSIXct("0001-01-01"),SENSOREVENT=0,
                        TRANSMITTERID="",RECEIVERID="",DURATION=0,STARTSENSOR=0,ENDSENSOR=0,MAXSENSOR=0,
                        ENDREASON="",NUMRECS=0,stringsAsFactors=FALSE)[NULL,]
    
    # Create empty output log dataframe with header row
    logtable <- data.frame(DATETIME=as.POSIXct("0001-01-01"),SENSOREVENT=0,RECORD=0,
                           TRANSMITTERID="",RECEIVERID="",SENSOR1=0,ELAPSED=0,stringsAsFactors=FALSE)[NULL,]
    
    # Create empty temporary log dataframe with header row
    fTLog <- logtable
    
    # Extract depth records for individual transmitter
    infile<-ExtractData(sInputFile2, sQueryTransmitterList = TransmitterNames[sTransmitterId])
  
    # initialise variables for this transmitter
    fSensorEventStarted <- FALSE
    fSensorThresholdPassed <- FALSE
    fEventComplete <- FALSE
    sEndReason <- ""
    iRecordsProcessed <- 0
    iNumberOfRecords <- 0
    sRECEIVERID <- ""
    sPreviousReceiver <- ""
    sEndReceiver <- ""
    sSTARTTIME <- ""
    sENDTIME <- ""
    sPreviousTime <- ""
    rPreviousSensor <- 0
    rSensorStart <- 0
    rSensorEnd <- 0
    rMinimumSensorValue <- 0
    iElapsedTime <- 0
    iElapsedTimeStart <- 0
    iEventRecords <- 0
    
    # For each line in source file
    while(iCount < nrow(infile))
    {
      # Increment line counter
      iCount <- iCount + 1
    
      # Decompose the input line to its component fields
      sDATETIME <- infile[iCount,1]
      rSENSOR1 <- infile[iCount,3]
      sRECEIVERID <- as.character(infile[iCount,iLocationCol])
    
      # Evaluate this record
      iRecordsProcessed <- iRecordsProcessed + 1
      if (fSensorEventStarted)
      {
        # During a sensor event
        iElapsedTime <- as.numeric(difftime(as.POSIXct(as.character(sDATETIME)), as.POSIXct(as.character(sPreviousTime)) , units = "secs"))
        iElapsedTimeStart <- as.numeric(difftime(as.POSIXct(as.character(sDATETIME)), as.POSIXct(as.character(sSTARTTIME)) , units = "secs"))      
        logtable[ilogevent,] <- data.frame(DATETIME=as.character(sDATETIME),SENSOREVENT=iSensorEvents,RECORD=iNumberOfRecords,
                                           TRANSMITTERID=as.character(TransmitterNames[sTransmitterId]),RECEIVERID=as.character(sRECEIVERID),
                                           SENSOR1=rSENSOR1,ELAPSED=iElapsedTimeStart,stringsAsFactors=FALSE)
        ilogevent <- ilogevent + 1
        iNumberOfRecords <- iNumberOfRecords + 1
                            
        # After passing predefined sensor theshold, transmitter returns to sensor value within predefined theshold from starting value
        if(fSensorEventStarted)
          if(fSensorThresholdPassed)
          {
            if (rPreviousSensor > rSensorStart + rTerminationThreshold)
              if (rSENSOR1 <= rSensorStart + rTerminationThreshold)
              {
                fSensorEventStarted <- FALSE
                fSensorThresholdPassed <- FALSE
                fEventComplete <- TRUE            
                sEndReason <- "return"
                rSensorEnd <- rSENSOR1
                sENDTIME <- sDATETIME 
                sEndReceiver <- sRECEIVERID
                iEventRecords <- iNumberOfRecords
              }
          }else{
            # Transmitter crosses predefined theshold from original reading
            if (rPreviousSensor < rSensorStart + rTerminationThreshold)
              if (rSENSOR1 >= rSensorStart + rTerminationThreshold)
                fSensorThresholdPassed <- TRUE
          }
          
        # Check if time elapsed between pings is greater than threshold 
        if (iElapsedTime > iTimeThreshold)
        {
          fSensorEventStarted <- TRUE
          fSensorThresholdPassed <- TRUE
          fEventComplete <- FALSE
          sEndReason <- "timeout"
          rSensorEnd <- rPreviousSensor
          sENDTIME <- sPreviousTime 
          sEndReceiver <- sPreviousReceiver
          iEventRecords <- iNumberOfRecords - 1
        }
      
        # If sensor value is greater than the maximum logged sensor value during that event, log this as the new maximum
        if (fSensorEventStarted)
          if (rSENSOR1 > rMaximumSensorValue)
            rMaximumSensorValue <- rSENSOR1
        
        # Sensor event has ended, write record to event table
        if (sEndReason!="")
        {
          iElapsedTimeStart <- as.numeric(difftime(as.POSIXct(as.character(sDATETIME)), as.POSIXct(as.character(sSTARTTIME)) , units = "secs"))
          event[iSensorEvents,] <- data.frame(STARTTIME=sSTARTTIME,ENDTIME=sENDTIME,SENSOREVENT=iSensorEvents,
                                              TRANSMITTERID=as.character(TransmitterNames[sTransmitterId]),RECEIVERID=as.character(sEndReceiver),
                                              DURATION=iElapsedTimeStart,STARTSENSOR=rSensorStart,ENDSENSOR=rSensorEnd,MAXSENSOR=rMaximumSensorValue,
                                              ENDREASON=sEndReason,NUMRECS=iEventRecords,stringsAsFactors=FALSE)
          # Reset counters
          fSensorEventStarted <- FALSE
          fSensorThresholdPassed <- FALSE
          fEventComplete <- FALSE
        
          # Read lines from previous ping to ensure all events recorded
          iCount <- iCount - 1
          sDATETIME <- infile[iCount,1]
          rSENSOR1 <- infile[iCount,3]
          sRECEIVERID <- as.character(infile[iCount,iLocationCol])
          sSTARTTIME <- sDATETIME      
          rSensorStart <- rSENSOR1
          iNumberOfRecords <- 0
          fTLog <- fTLog[NULL,]
          sEndReason <-""
          sPreviousTime <- infile[iCount-1,1]
          rPreviousSensor <- infile[iCount-1,3]
          sPreviousReceiver <- as.character(infile[iCount-1,iLocationCol])
        }
    
      }else{     
        # Waiting for a sensor event
        if (iRecordsProcessed > 1)
        {
          iElapsedTime <- as.numeric(difftime(as.POSIXct(as.character(sDATETIME)), as.POSIXct(as.character(sPreviousTime)) , units = "secs")) 
          # Does sensor value increase beyond threshold
          if (rSENSOR1 > rSensorStart + rTriggerThreshold)
          {
           # Does the sensor threshold change within the given time window
           iElapsedTimeStart <- as.numeric(difftime(as.POSIXct(as.character(sDATETIME)), as.POSIXct(as.character(sSTARTTIME)) , units = "secs"))
           if (iElapsedTimeStart <= iTimeThresholdStart & iElapsedTime <= iTimeThreshold)
            {
              # Threshold reached therefore sensor event starting
              iSensorEvents <- iSensorEvents + 1
              fSensorEventStarted <- TRUE
              sEndReason <- ""
              rMaximumSensorValue <- rSENSOR1
              # Patch temporary sensor data frame to logtable and clear temporary log table
              # fTLog needs to be present for starting values to work and logtable to process
              if(nrow(fTLog) == 0)
              {
                fTLog <- data.frame(DATETIME=sSTARTTIME,SENSOREVENT=iSensorEvents,RECORD=0,TRANSMITTERID=as.character(TransmitterNames[sTransmitterId]),
                                    RECEIVERID=as.character(sPreviousReceiver),SENSOR1=rSensorStart,ELAPSED=0,stringsAsFactors=FALSE) 
                iNumberOfRecords <- iNumberOfRecords+1
              }
            
              logtable[ilogevent:(ilogevent + nrow(fTLog)-1),] <- fTLog
              ilogevent <- ilogevent + nrow(fTLog)
              fTLog <- fTLog[NULL,]
              
              # Add new data to logtable
              logtable[ilogevent,] <- data.frame(DATETIME=sDATETIME,SENSOREVENT=iSensorEvents,RECORD=iNumberOfRecords,
                                                TRANSMITTERID=as.character(TransmitterNames[sTransmitterId]),RECEIVERID=as.character(sRECEIVERID),
                                                SENSOR1=rSENSOR1,ELAPSED=iElapsedTime,stringsAsFactors=FALSE) 
              iNumberOfRecords <- iNumberOfRecords + 1
              ilogevent <- ilogevent+1
            
              # Transmitter crosses predefined theshold from original reading
              if (rPreviousSensor < rSensorStart + rTerminationThreshold)
                if (rSENSOR1 >= rSensorStart + rTerminationThreshold)
                 fSensorThresholdPassed <- TRUE
            
            }else{
              # If sensor value has increased sufficiently but time threshold exceeded
              # Reset sensor and date/time starting values
              sSTARTTIME <- sDATETIME      
              rSensorStart <- rSENSOR1
              fSensorEventStarted <- FALSE
              iNumberOfRecords <- 0
            
             # Reset temporary log file
             fTLog <- fTLog[NULL,]
             fTLog[iNumberOfRecords+1,] <- data.frame(DATETIME=sSTARTTIME,SENSOREVENT=iSensorEvents+1,RECORD=iNumberOfRecords,TRANSMITTERID=as.character(TransmitterNames[sTransmitterId]),
                                                      RECEIVERID=as.character(sRECEIVERID),SENSOR1=rSensorStart,ELAPSED=0,stringsAsFactors=FALSE)
             iNumberOfRecords <- iNumberOfRecords + 1
            }
          }else{
            # If sensor value is greater than the previous sensor value and is within the predefined time thresold
            # Add to the current temporary log file
            if(rSENSOR1 > rPreviousSensor & iElapsedTime <= iTimeThreshold)
            {
              fSensorEventStarted <- FALSE
              fTLog[iNumberOfRecords+1,] <- data.frame(DATETIME=sDATETIME,SENSOREVENT=iSensorEvents+1,RECORD=iNumberOfRecords,TRANSMITTERID=as.character(TransmitterNames[sTransmitterId]),
                                                       RECEIVERID=as.character(sRECEIVERID),SENSOR1=rSENSOR1,ELAPSED=iElapsedTime,stringsAsFactors=FALSE)
              iNumberOfRecords <- iNumberOfRecords + 1
            }            
            if(rSENSOR1 <= rPreviousSensor | iElapsedTime > iTimeThreshold)
            {
              # Clear the temporary logfile and reset the counters 
              sSTARTTIME <- sDATETIME      
              rSensorStart <- rSENSOR1
              fSensorEventStarted <- FALSE
              iNumberOfRecords <- 0
              fTLog <- fTLog[NULL,]
              fTLog[iNumberOfRecords+1,] <- data.frame(DATETIME=sSTARTTIME,SENSOREVENT=iSensorEvents+1,RECORD=0,TRANSMITTERID=as.character(TransmitterNames[sTransmitterId]),
                                                       RECEIVERID=as.character(sRECEIVERID),SENSOR1=rSensorStart,ELAPSED=0,stringsAsFactors=FALSE)
              iNumberOfRecords <- iNumberOfRecords + 1
            }
          }
        }else{
          sSTARTTIME <- sDATETIME
          rSensorStart <- rSENSOR1     
        }
      }
      # If No. events < 2, then log previous date/time, sensor and receiver values
      sPreviousTime <- sDATETIME
      rPreviousSensor <- rSENSOR1
      sPreviousReceiver <- as.character(sRECEIVERID)
    }
    # Reorganise function results to ensure all transmitters occur in relevent list fields (i.e. event, logfile)  
    ilist <- new.env()
    ilist$event <- event
    ilist$logtable <- logtable
  
    return(as.list(ilist))
  }
  ###
  
  # Run Sensor Event Extraction funtions on all 1:n transmitters
#  results <- lapply(1:iTransmitterCount,SensorExtractId)
 
  results <- foreach(i=1:iTransmitterCount, .packages='VTrack') %dopar% SensorExtractId(i)
  
  
  # Reorganise function results to ensure all transmitters occur in relevent list fields (i.e. event, logfile, nonresidences)
  # Ensure data headers are correct and the relevent fields are represented (i.e. Station and Receiver)
  ilist2 <- new.env()
  ilist2$event <- do.call(rbind,(do.call(rbind,results)[,2]))
  names(ilist2$event)[5] <- sLocation
  ilist2$logtable <- do.call(rbind,(do.call(rbind,results)[,1]))
  names(ilist2$logtable)[5] <- sLocation
  
  # Return Sensor values to positive if a decreasing sensor event 
  if(iEventType=="DECREASE")
  {
    ilist2$event[,7:9]<- -1*(ilist2$event[,7:9])
    ilist2$logtable[,6]<- -1*(ilist2$logtable[,6])
    names(ilist2$event)[9]<- "MINSENSOR"
   }
  return(as.list(ilist2))

if (.Platform$OS == "windows") stopCluster()
}
