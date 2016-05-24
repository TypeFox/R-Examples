
writeMonthlyStations <- function(filename = MASTER.STATION.LIST, outfile = MONTHLY.STATION.LIST ){
  
    LookUpInfo <- read.csv(filename,stringsAsFactors=FALSE)
    # first row is bogosity
    LookUpInfo <- LookUpInfo[-1,]
     
     
    Data       <- data.frame(Name    = LookUpInfo$name,
                             Year    = substr(LookUpInfo$mlyRange,1,4),
                             WebId   = LookUpInfo$webid)

    dex <-  which(Data$Year == "N/A")
    Data <- Data[-dex,]
    nex  <- which(!is.na(Data$WebId))
    Data <- Data[nex,]
    newId <-seq(from = STARTING.STATION.ID,length =nrow(Data))
    Data <- data.frame(Id=newId,Data)
    
    if (!is.null(outfile)) write.csv(Data,outfile) 
  
    return(Data)
  
}