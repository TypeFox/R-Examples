ExtractReliabilityForTransmittersFromArchive <-
function(sInputFile,sTransmitterFile,sOutputFile)
{
    # load the list of transmitter names from a file
    TransmitterNames <- read.csv(sTransmitterFile)
    iTransmitterCount <- dim(TransmitterNames)[1]

    # create blank matrix
    ReliabilityMatrix <- matrix(nrow = iTransmitterCount, ncol = 24)
    for (i in 1:iTransmitterCount)
        for (j in 1:24)
          ReliabilityMatrix[i,j] <- 0
    sPreviousDate <- ""
              
    # create output file 
    outfile <- file(sOutputFile,"wt")
    writeLines("DATETIME,TRANSMITTERID,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,TOTAL",outfile)

    # open source file
    infile <- file(sInputFile,"rt") 
    sHeader <- readLines(infile,1)
    iLength <- length(sHeader)
    
    # for each line in source file
    while (iLength > 0)
    {
        sLine <- readLines(infile,1)
        iLength <- length(sLine)
        if (iLength > 0)
        {
            # increment line counter
            iCount <- iCount + 1

            # decompose the input line to its component fields
            aValuesVector <- unlist(strsplit(sLine,","))

            # when the date changes, write the reliability matrix to file and reinit the reliability matrix
            sDate <- substring(aValuesVector[1],1,10)
            iHour <- as.integer(substring(aValuesVector[1],12,13))
            if (sPreviousDate != sDate)
            {
                for (i in 1:iTransmitterCount)
                {
                    sLine <- paste(sPreviousDate,TransmitterNames[i,1],sep=",")
                    iAmount <- 0
                    for (j in 1:24)
                    {
                        sLine <- paste(sLine,ReliabilityMatrix[i,j],sep=",")
                        iAmount <- iAmount + ReliabilityMatrix[i,j]
                    }                           
                    
                    if (iAmount > 0)
                        writeLines(paste(sLine,iAmount,sep=","),outfile)
                }
                
                for (i in 1:iTransmitterCount)
                    for (j in 1:24)
                        ReliabilityMatrix[i,j] <- 0
                sPreviousDate <- sDate
            }
            
            # find transmitter index
            iTransmitterIndex <- 0
            for (i in 1:iTransmitterCount)
                if (aValuesVector[2] == TransmitterNames[i,1])
                    iTransmitterIndex <- i

            # update reliability matrix for this record
            ReliabilityMatrix[iTransmitterIndex,iHour+1] <- ReliabilityMatrix[iTransmitterIndex,iHour+1] + 1
        }
    }    

    # write reliability matrix, a row for each transmitter
    for (i in 1:iTransmitterCount)
    {
        sLine <- paste(sDate,TransmitterNames[i,1],sep=",")
        iAmount <- 0
        for (j in 1:24)
        {
            sLine <- paste(sLine,ReliabilityMatrix[i,j],sep=",")
            iAmount <- iAmount + ReliabilityMatrix[i,j]
        }                           
                    
        if (iAmount > 0)
           writeLines(paste(sLine,iAmount,sep=","),outfile)
    }            

    # close the files
    close(infile)
    close(outfile)
}
