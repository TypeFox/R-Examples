blastSeq <- function(seq, n_blast=20, delay_req=3, delay_rid=60, email=NULL, xmlFolder=NULL, logFolder=NULL, keepInMemory=TRUE, database="chromosome", verbose=TRUE, createLog=FALSE){

  startTime <- Sys.time()
  
# Polite system sleeps as requested from NCBI  
  if(delay_req<3) stop("Sending more requests than once every 3 seconds is considered to be rude from NCBI!")  
  if(delay_rid<60) stop("Polling RID results more often than once per minute is considered to be rude from NCBI!")  
  if(is.null(email)) stop("NCBI requires to provide an email address, please give one in the function call!")

# Getting some needed variables
  totalSeq <- length(seq)
  if(is.null(names(seq))) names(seq) <- 1:totalSeq

# Check about the xml folder settings
  writeXML <- FALSE
  if(!is.null(xmlFolder)){
  # check if the path is provided properly
    
    writeXML <- TRUE
    dir.create(xmlFolder, showWarnings = FALSE)
  }
  
# Set up the log folder
  if(createLog){
    if(is.null(logFolder)){
      if(is.null(xmlFolder)){
        stop("No log/xml path given. For small projects, please use the option 'createLog=FALSE'")
      } else {
        logFolder <- strsplit(xmlFolder,"/")[[1]]
        logFolder <- logFolder[nchar(logFolder)>0]
        logFolder <- logFolder[1:(length(logFolder)-1)]
        logFolder <- paste(c(logFolder,"hoardeRlogs/"),collapse="/")
        logFolder <- paste("/",logFolder,sep="")
      }
    }
    dir.create(logFolder, showWarnings = FALSE)    
  }

# Create the RID/sequence info table
  seqInfo <- data.frame(seqNames=names(seq),
                        seqRID=rep("0",length(seq)),
                        seqFinished=rep(FALSE,length(seq)),
                        seqRuntime=rep("00:00:00",length(seq)),
                        stringsAsFactors=FALSE)

# Read/Write the RID/sequence info table
  if(createLog){
    if(file.exists(paste(logFolder,"seqRID-info.csv",sep=""))){
      seqInfoImported <- read.table(paste(logFolder,"seqRID-info.csv",sep=""), header=TRUE, stringsAsFactors=FALSE)
      if(sum(seqInfo$seqNames==seqInfoImported$seqNames)!=nrow(seqInfo)) stop("Log-file mismatch! Please provide the right log file for the existing project or change the logFolder option.")
      seqInfo <- seqInfoImported
    } else {
      write.table(seqInfo,file=paste(logFolder,"seqRID-info.csv",sep=""), quote=FALSE, row.names=FALSE)      
    }
  }

# Write the project Log
  if(createLog){
    if(file.exists(paste(logFolder,"hoardeR.log",sep=""))){
      fileConn<-file(paste(logFolder,"hoardeR.log",sep=""))
      logLines <- readLines(fileConn)
      close(fileConn)      
      startTime <- as.POSIXct(logLines[1])
    } else {  
      fileConn<-file(paste(logFolder,"hoardeR.log",sep=""))
      writeLines(c(as.character(startTime),
                   "---------------------------------",
                   "Settings:",
                   paste("n_blast",n_blast)),fileConn)
      close(fileConn)      
      
    }

  }

# Store here the blast RIDs
  RID <- rep(0,totalSeq)

# Store here the results (In case no log file is created, otherwise get the last stored result)
  if(createLog){
  # Adjust this here still according to 'keepInMemory' settings
    res <- list()
    sendThis <- min(which(seqInfo$seqRID==0))
    curRunning <- sum((seqInfo$seqRID!=0) & seqInfo$seqFinished==FALSE)
    ready <- sum(seqInfo$seqFinished==TRUE)
    active <- which(((seqInfo$seqRID!=0) & seqInfo$seqFinished==FALSE)==TRUE)
    RID <- seqInfo$seqRID
    timeAvg <- c("00:00:00")  
  } else {
    res <- list()
    sendThis <- 1
    curRunning <- 0
    ready <- 0
    active <- NULL
    timeAvg <- c("00:00:00")  
  }
# This is very optimistic, maybe I should take also a time break, in case one result doesn't get ready
  while(ready < totalSeq){
    if((curRunning < n_blast) & (sendThis <= totalSeq)){
      Sys.sleep(delay_req)
      RID[sendThis] <- sendFA(seq[sendThis],email=email, database=database)
      curRunning <- curRunning + 1    
      active <- c(active,sendThis)
      seqInfo$seqRID[sendThis] <- RID[sendThis]
      sendThis <- sendThis + 1
      write.table(seqInfo,file=paste(logFolder,"seqRID-info.csv",sep=""), quote=FALSE, row.names=FALSE) 
    } else {
      Sys.sleep(delay_rid) 
      for(i in active){
        res[[i]] <- getBlastResult(RID[i]) 
        Sys.sleep(delay_req)
        if(res[[i]]$ready==TRUE){
           ready <- ready + 1
           curRunning <- curRunning - 1
           active <- active[-which(active==i)]
           seqInfo$seqFinished[seqInfo$seqNames==names(seq)[i]] <- TRUE
           timings <- seqInfo$seqRuntime[(seqInfo$seqFinished==TRUE)]
           timeAvg <- timeStat(timings[timings!="0"])
        # Write here then the XML file to the folder
           if(writeXML){
             file.create(paste(xmlFolder,names(seq)[i],".xml",sep=""))
             fileConn <- file(paste(xmlFolder,names(seq)[i],".xml",sep=""))
             writeLines(res[[i]]$blastRes, fileConn)
             close(fileConn)  
             write.table(seqInfo,file=paste(logFolder,"seqRID-info.csv",sep=""), quote=FALSE, row.names=FALSE)  
           }
           if(!keepInMemory){
             res[[i]] <- NULL
           }
        } else {
          seqInfo$seqRuntime[seqInfo$seqNames==names(seq)[i]] <- res[[i]]$time
          write.table(seqInfo,file=paste(logFolder,"seqRID-info.csv",sep=""), quote=FALSE, row.names=FALSE)  
        }
      }
    }
    if(verbose){
      cat("Missing:",totalSeq-ready,"\n")
      cat("Running:",length(active),"\n")
      cat("Finished:",ready,"\n")
      cat("Avg. Blast Time:",timeAvg,"\n")
      cat("Total running time:",secToTime(as.numeric(Sys.time() - startTime, units="secs")),"\n")
      cat("---------------------------------------------------------------\n")
    }
  }
 result <- list(RID=RID, res=res, database=database)
 class(result) <- "blastRes"
 result
}
