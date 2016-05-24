RunResidenceExtraction <-
function(sInputFile,sLocation,iResidenceThreshold,iTimeThreshold,sDistanceMatrix=NULL)
{
  i <- NULL 
#  if(Sys.info()["sysname"]=="Darwin"){
#    registerDoMC()
#  }else{
  registerDoParallel(cl=(cl <- makeCluster(2)))
#  }
    
  # Check data are in chronological order and duplicates-free
  sInputFile <- unique(sInputFile[order(as.character(sInputFile$DATETIME)),]) 
  
  # Extract list of transmitter names
  TransmitterNames <- ExtractUniqueValues(sInputFile,2)
  iTransmitterCount <- length(TransmitterNames)
  
  # Define whether working with receivers or stations 
  if(sLocation=="STATIONNAME")
    iLocationCol <- 6
  if(sLocation=="RECEIVERID")
    iLocationCol <- 5

  # Function to extract residence events, residence logfile, plus nonresidence events
  ResidenceExtractId <- function(sTransmitterId)
  {
    # Initialise variables
    iCount <- 1
    iResidenceEvents <- 0
    ilogevent <- 1
    
    # Create empty output event and log dataframe with header row
    event <- data.frame(STARTTIME=sInputFile[1,1],ENDTIME=sInputFile[1,1],RESIDENCEEVENT=0,
                        TRANSMITTERID=TransmitterNames[1],RECEIVERID=as.character(sInputFile[1,iLocationCol]),DURATION=0,
                        ENDREASON=NA,NUMRECS=0,stringsAsFactors=FALSE)[NULL,]
    logtable <- data.frame(DATETIME=sInputFile[1,1],RESIDENCEEVENT=0,RECORD=0,TRANSMITTERID=TransmitterNames[1],
                           RECEIVERID=as.character(sInputFile[1,iLocationCol]),ELAPSED=0,stringsAsFactors=FALSE)[NULL,]
    
    # For each transmitter, extract records for only this transmitter
    infile<-ExtractData(sInputFile, sQueryTransmitterList = TransmitterNames[sTransmitterId])
      
    # Initialise variables for this transmitter
    fResidenceEventStarted <- FALSE
    sEndReason <- ""
    iRecordsProcessed <- 0
    iNumberOfRecords <- 0
    sRECEIVERID <- ""
    sPreviousReceiver <- ""
    sSTARTTIME <- ""
    sPreviousTime <- ""
    
    # For each line in source file
    while (iCount <= nrow(infile))
    {
    
      # Decompose the input line to its component fields and choose which Location (RECEIVER or STATION) to work with
      sDATETIME <- infile[iCount,1]
      sRECEIVERID <- as.character(infile[iCount,iLocationCol])
      
      # Evaluate this record
      iRecordsProcessed <- iRecordsProcessed + 1
      if (fResidenceEventStarted)
      {
        # During a residence event
        iNumberOfRecords <- iNumberOfRecords + 1
        iElapsedTime <- as.numeric(difftime(as.POSIXct(as.character(sDATETIME)), as.POSIXct(as.character(sPreviousTime)) , units = "secs"))
        logtable[ilogevent,] <- data.frame(DATETIME=as.character(sDATETIME),RESIDENCEEVENT=iResidenceEvents,RECORD=iNumberOfRecords,
                                           TRANSMITTERID=as.character(TransmitterNames[sTransmitterId]),RECEIVERID=as.character(sRECEIVERID),
                                           ELAPSED=iElapsedTime,stringsAsFactors=FALSE)
        ilogevent <- ilogevent + 1
        # Timeout period
        if (iElapsedTime > iTimeThreshold)
        {
          fResidenceEventStarted <- FALSE
          sEndReason <- "timeout"
        }
        # New receiver/station
        if (sPreviousReceiver != sRECEIVERID)
        {
          fResidenceEventStarted <- FALSE
          sEndReason <- "receiver"
        }
        # Lost signal
        if (iCount == nrow(infile))
        {
          fResidenceEventStarted <- FALSE
          sEndReason <- "signal lost"
        }
        if (fResidenceEventStarted)
        {
          # Residence is continuing
        }else{
          if (iNumberOfRecords >= iResidenceThreshold)
          {
            # Residence has ended, write record to event table
            iElapsedTime <- as.numeric(difftime(as.POSIXct(as.character(sPreviousTime)), as.POSIXct(as.character(sSTARTTIME)) , units = "secs"))
            event[iResidenceEvents,] <- data.frame(STARTTIME=as.character(sSTARTTIME),ENDTIME=as.character(sPreviousTime),
                                                  RESIDENCEEVENT=iResidenceEvents,TRANSMITTERID=as.character(TransmitterNames[sTransmitterId]),
                                                  RECEIVERID=as.character(sPreviousReceiver),DURATION=iElapsedTime,ENDREASON=sEndReason,
                                                  NUMRECS=iNumberOfRecords,stringsAsFactors=FALSE)
#            if(sEndReason=="signal lost")
#              iResidenceEvents <- iResidenceEvents + 1
          }
        }
      }else{
        # Waiting for a residence event
        if (iRecordsProcessed > 1)
        {
          iElapsedTime <- as.numeric(difftime(as.POSIXct(as.character(sDATETIME)), as.POSIXct(as.character(sPreviousTime)) , units = "secs"))
          if (iElapsedTime < iTimeThreshold)
          {
            # Residence event starting write records to log table
            iResidenceEvents <- iResidenceEvents + 1
            iNumberOfRecords <- 1
            sSTARTTIME <- sPreviousTime
            sEndReason <- ""
            logtable[ilogevent,] <- data.frame(DATETIME=sPreviousTime,RESIDENCEEVENT=iResidenceEvents,RECORD=0,
                                               TRANSMITTERID=as.character(TransmitterNames[sTransmitterId]),
                                               RECEIVERID=as.character(sPreviousReceiver),ELAPSED=0,stringsAsFactors=FALSE)
            logtable[ilogevent+1,] <- data.frame(DATETIME=sDATETIME,RESIDENCEEVENT=iResidenceEvents,RECORD=iNumberOfRecords,
                                                 TRANSMITTERID=as.character(TransmitterNames[sTransmitterId]),
                                                 RECEIVERID=as.character(sRECEIVERID),ELAPSED=iElapsedTime,stringsAsFactors=FALSE)
            ilogevent <- ilogevent+2
            
            if(as.character(sRECEIVERID)==as.character(sPreviousReceiver))
            {
              fResidenceEventStarted <- TRUE
            }else{
              fResidenceEventStarted <- FALSE
            }
          }
        }
      }
    
      sPreviousTime <- sDATETIME
      sPreviousReceiver <- as.character(sRECEIVERID)
      iCount <- iCount + 1
    }
    # Collect all residence events and logs for that transmitter and write dataframe into a list object
    ilist <- new.env()
    ilist<- list(na.omit(logtable),
                  na.omit(event),
                  NonResidenceExtractId(na.omit(event),sDistanceMatrix))

    return(ilist)
  }

  # Function to extract residence events, residence logfile, plus nonresidence events
  ResidenceExtractId1 <- function(sTransmitterId)
  { 
    # Initialise variables
    iCount <- 1
    iResidenceEvents <- 0
    ilogevent <- 1
  
    # Create empty output event and log dataframe with header row
    event <- data.frame(STARTTIME=sInputFile[1,1],ENDTIME=sInputFile[1,1],RESIDENCEEVENT=0,
                        TRANSMITTERID=TransmitterNames[1],RECEIVERID=as.character(sInputFile[1,iLocationCol]),DURATION=0,
                        ENDREASON=NA,NUMRECS=0,stringsAsFactors=FALSE)[NULL,]
    logtable <- data.frame(DATETIME=sInputFile[1,1],RESIDENCEEVENT=0,RECORD=0,TRANSMITTERID=TransmitterNames[1],
                           RECEIVERID=as.character(sInputFile[1,iLocationCol]),ELAPSED=0,stringsAsFactors=FALSE)[NULL,]
  
    # For each transmitter, extract records for only this transmitter
    infile <- ExtractData(sInputFile, sQueryTransmitterList = TransmitterNames[sTransmitterId])
   
    # Initialise variables for this transmitter
    fResidenceEventStarted <- FALSE
    sEndReason <- ""
    iRecordsProcessed <- 0
    iNumberOfRecords <- 0
    sRECEIVERID <- ""
    sPreviousReceiver <- ""
    sSTARTTIME <- ""
    sPreviousTime <- ""
    iElapsedTime <- 0 #new
  
    # For each line in source file
    while (iCount <= nrow(infile))
    {  
      # Set data entry for this row
      sDATETIME <- infile[iCount,1]
      sRECEIVERID <- as.character(infile[iCount,iLocationCol])
      
      # Evaluate this record
      iRecordsProcessed <- iRecordsProcessed + 1
    
      # Third Record in a Residence event
      if(iNumberOfRecords > 1){
        if(fResidenceEventStarted ==TRUE){
          # During a residence event
          iElapsedTime <- as.numeric(difftime(as.POSIXct(as.character(sDATETIME)), as.POSIXct(as.character(sPreviousTime)) , units = "secs"))
          ilogevent <- ilogevent + 1  
          logtable[ilogevent,] <- data.frame(DATETIME=as.character(sDATETIME),RESIDENCEEVENT=iResidenceEvents,RECORD=iNumberOfRecords,
                                             TRANSMITTERID=as.character(TransmitterNames[sTransmitterId]),RECEIVERID=as.character(sRECEIVERID),
                                             ELAPSED=iElapsedTime,stringsAsFactors=FALSE)
          iNumberOfRecords <- iNumberOfRecords +1
        
          # Timeout period
          if(iElapsedTime > iTimeThreshold)
          {
            fResidenceEventStarted <- FALSE
            sEndReason <- "timeout"
          }
          # New receiver/station
          if(sPreviousReceiver != sRECEIVERID)
          {
            fResidenceEventStarted <- FALSE
            sEndReason <- "receiver"
          }
          # Lost signal
          if(iCount == nrow(infile))
          {
            fResidenceEventStarted <- FALSE
            sEndReason <- "signal lost"
          }
        
          if(fResidenceEventStarted==FALSE)
          {
            if(sEndReason != "signal lost")
            {
              iElapsedTimeTotal <- as.numeric(difftime(as.POSIXct(as.character(sPreviousTime)), as.POSIXct(as.character(sSTARTTIME)) , units = "secs"))
              event[iResidenceEvents,] <- data.frame(STARTTIME=as.character(sSTARTTIME),ENDTIME=as.character(sPreviousTime),
                                                     RESIDENCEEVENT=iResidenceEvents,TRANSMITTERID=as.character(TransmitterNames[sTransmitterId]),
                                                     RECEIVERID=as.character(sPreviousReceiver),DURATION=iElapsedTimeTotal,ENDREASON=sEndReason,
                                                     NUMRECS=iNumberOfRecords-1,stringsAsFactors=FALSE)
              iNumberOfRecords <- 0
              ilogevent <- ilogevent + 1
            }else{
              iElapsedTimeTotal <- as.numeric(difftime(as.POSIXct(as.character(sDATETIME)), as.POSIXct(as.character(sSTARTTIME)) , units = "secs"))
              event[iResidenceEvents,] <- data.frame(STARTTIME=as.character(sSTARTTIME),ENDTIME=as.character(sDATETIME),
                                                     RESIDENCEEVENT=iResidenceEvents,TRANSMITTERID=as.character(TransmitterNames[sTransmitterId]),
                                                     RECEIVERID=as.character(sPreviousReceiver),DURATION=iElapsedTimeTotal,ENDREASON=sEndReason,
                                                     NUMRECS=iNumberOfRecords,stringsAsFactors=FALSE)        
            }
          }    
        }
      }
    
      # Second Record in Residence event
      if(iNumberOfRecords == 1){   
        iElapsedTime <- as.numeric(difftime(as.POSIXct(as.character(sDATETIME)), as.POSIXct(as.character(sPreviousTime)) , units = "secs"))
        # If residence ends after 1 record
        if(iElapsedTime >= iTimeThreshold | sPreviousReceiver != sRECEIVERID){
        # Timeout period
          if(iElapsedTime > iTimeThreshold)
            sEndReason <- "timeout"
        # New receiver/station
          if(sPreviousReceiver != sRECEIVERID)
            sEndReason <- "receiver"
        # Lost signal
          if(iCount == nrow(infile))
            sEndReason <- "signal lost"
        if(sEndReason == "receiver")
          event[iResidenceEvents,] <- data.frame(STARTTIME=as.character(sSTARTTIME),ENDTIME=as.character(sPreviousTime),
                                                 RESIDENCEEVENT=iResidenceEvents,TRANSMITTERID=as.character(TransmitterNames[sTransmitterId]),
                                                 RECEIVERID=as.character(sPreviousReceiver),DURATION=0,ENDREASON=sEndReason,
                                                 NUMRECS=iNumberOfRecords,stringsAsFactors=FALSE)  
        
        if(sEndReason != "receiver" & sEndReason != "")
          event[iResidenceEvents,] <- data.frame(STARTTIME=as.character(sSTARTTIME),ENDTIME=as.character(sPreviousTime),
                                                 RESIDENCEEVENT=iResidenceEvents,TRANSMITTERID=as.character(TransmitterNames[sTransmitterId]),
                                                 RECEIVERID=as.character(sPreviousReceiver),DURATION=iElapsedTime,ENDREASON=sEndReason,
                                                 NUMRECS=iNumberOfRecords,stringsAsFactors=FALSE)  
        
        iNumberOfRecords <- 0
        #iResidenceEvents <- iResidenceEvents+1
        fResidenceEventStarted <- FALSE
        
        ilogevent <- ilogevent+1
        logtable[ilogevent,] <- data.frame(DATETIME=sDATETIME,RESIDENCEEVENT=iResidenceEvents,RECORD=1,
                                           TRANSMITTERID=as.character(TransmitterNames[sTransmitterId]),
                                           RECEIVERID=as.character(sRECEIVERID),ELAPSED=iElapsedTime,stringsAsFactors=FALSE) 
        } 
      
        # If residence continues 
        if(iElapsedTime < iTimeThreshold & sPreviousReceiver == sRECEIVERID){      
          iNumberOfRecords <- iNumberOfRecords +1
          ilogevent <- ilogevent+1
          logtable[ilogevent,] <- data.frame(DATETIME=sDATETIME,RESIDENCEEVENT=iResidenceEvents,RECORD=1,
                                             TRANSMITTERID=as.character(TransmitterNames[sTransmitterId]),
                                             RECEIVERID=as.character(sRECEIVERID),ELAPSED=iElapsedTime,stringsAsFactors=FALSE) 
          fResidenceEventStarted <- TRUE
        }
      }
      
      # First Record in Residence event
      if(iNumberOfRecords == 0){      
      sSTARTTIME <- sDATETIME
      iResidenceEvents <- iResidenceEvents + 1
      iNumberOfRecords <- iNumberOfRecords + 1
      logtable[ilogevent,] <- data.frame(DATETIME=sSTARTTIME,RESIDENCEEVENT=iResidenceEvents,RECORD=0,
                                         TRANSMITTERID=as.character(TransmitterNames[sTransmitterId]),
                                         RECEIVERID=as.character(sRECEIVERID),ELAPSED=0,stringsAsFactors=FALSE)
    }
    
    iCount <- iCount + 1 # Move to next row 
    sPreviousTime <- sDATETIME
    sPreviousReceiver <- as.character(sRECEIVERID)
  }
  
  # Collect all residence events and logs for that transmitter and write dataframe into a list object
  ilist <- new.env()
  ilist<- list(na.omit(logtable),
               na.omit(event),
               NonResidenceExtractId(na.omit(event),sDistanceMatrix))
  
  return(ilist)
  }

  # Run Residence and Non Residence Extraction funtions on all 1:n transmitters 

  # Runs old function for >1 detections at a receiver
  if(iResidenceThreshold > 1)
    results <- foreach(i=1:iTransmitterCount, .packages='VTrack') %dopar% ResidenceExtractId(i)
  # Runs updated function to include single detections at a receiver 
  if(iResidenceThreshold == 1)
    results <- foreach(i=1:iTransmitterCount, .packages='VTrack') %dopar% ResidenceExtractId1(i)


  # Reorganise function results to ensure all transmitters occur in relevent list fields (i.e. event, logfile, nonresidences)
  # Ensure data headers are correct and the relevent fields are represented (i.e. Station and Receiver)
  ilist2 <- new.env()
  ilist2$residences <- do.call(rbind,(do.call(rbind,results)[,2]))
  names(ilist2$residences)[5] <- sLocation
  ilist2$residenceslog <- do.call(rbind,(do.call(rbind,results)[,1]))
  names(ilist2$residenceslog)[5] <- sLocation
  ilist2$nonresidences <- do.call(rbind,(do.call(rbind,results)[,3]))
  names(ilist2$nonresidences)[5] <- paste(sLocation,"1",sep="")
  names(ilist2$nonresidences)[6] <- paste(sLocation,"2",sep="")
  return(as.list(ilist2))

if (.Platform$OS == "windows") stopCluster(cl)
}
