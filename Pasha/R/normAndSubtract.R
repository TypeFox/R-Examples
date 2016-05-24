
# #########################################################################################
# Function that permits to rescale a WIG file
# dataWIG = the data structure of the WIG to rescale
# mulFactor = the factor by which the WIG data will be multiplied
# addFactor = the factor that will be added to the WIG data (after multiplication
#              by the previous factor
# #########################################################################################

.rescaleWIG=function(dataWIG, mulFactor, addFactor)
{
    for(chr in names(dataWIG))
    {
        dataWIG[[chr]][dataWIG[[chr]]>=0]=(dataWIG[[chr]][dataWIG[[chr]]>=0]*mulFactor)+addFactor
    }
    return(dataWIG)
}


# #########################################################################################
# Function that permits to rescale provided WIG files and subtract Input WIG file if asked 
#
# inputFile = Path of the WIG file that must be used as input
# wigFileList = Table of pathes of WIG files to rescale (and input subtract if asked)
# rescaleInput = TRUE if input WIG must be rescaled, FALSE if not
# subtractInput = TRUE if the input WIG must be subtracted from other WIG files, FALSE if not
# binSize = size of the binsof the WIG files
# #########################################################################################

normAndSubtractWIG=function(wigFileList, inputFile=NA, rescaleInput=FALSE, meanBGLevelRescale=c(0,10), subtractInput=TRUE, binSize=50)
{
    
    if( subtractInput){
        if(rescaleInput)
        {
            cat("\nReading Input WIG file to rescale\n")
            inputData <- readWIG(inputFile)
            
            cat("|-- Computing INPUT average value\n")
            unlisted <- unlist(inputData)
            inputMean <- mean(unlisted[unlisted>=0])
            
            cat("|-- Rescaling INPUT\n")
            inputRescaled <- .rescaleWIG(inputData, (meanBGLevelRescale[2]-meanBGLevelRescale[1])/inputMean, meanBGLevelRescale[1])
            
            cat("|-- Saving INPUT\n")
            writeWIG(inputRescaled, paste(substr(inputFile,1,nchar(inputFile)-4), "_Scaled",sep=""), folder=NULL, fixedStep=binSize)
            
        }
        else
        {
            cat("\nReading rescaled Input file\n")
            inputRescaled <- readWIG(inputFile)
        }
    }
    
    
    cat("\nRescaling all IP\n")
    # Rescaling wig data mean (considered as a good background level estimation)
    for(currentWIG in wigFileList[!is.na(wigFileList)])
    {
        cat("|-- Considering IP = ", currentWIG, "\n")
        cat("|---|-- Reading WIG file\n")
        dataWIG <- readWIG(currentWIG)
        
        cat("|---|-- Computing average value\n")
        unlisted <- unlist(dataWIG)
        dataWIGMean <- mean(unlisted[unlisted>=0])
        
        cat("|---|-- Rescaling\n")
        dataWIGRescaled <- .rescaleWIG(dataWIG, (meanBGLevelRescale[2]-meanBGLevelRescale[1])/dataWIGMean, meanBGLevelRescale[1])
        
        # Subtracting scaled input data to ip (if required)
        if(subtractInput)
        {
            cat("|---|-- Subtracting Input\n")
            for(chr in intersect(names(dataWIGRescaled), names(inputRescaled)))
            {
                maxSize <- max(c(length(dataWIGRescaled[[chr]]), length(inputRescaled[[chr]])))
                length(dataWIGRescaled[[chr]]) <- maxSize
                length(inputRescaled[[chr]]) <- maxSize
                dataWIGRescaled[[chr]][is.na(dataWIGRescaled[[chr]])] <- 0
                inputRescaled[[chr]][is.na(inputRescaled[[chr]])] <- 0
                
                dataWIGRescaled[[chr]][dataWIGRescaled[[chr]]>=0] <- dataWIGRescaled[[chr]][dataWIGRescaled[[chr]]>=0]-inputRescaled[[chr]][dataWIGRescaled[[chr]]>=0]
            }
        }
        
        # Saving rescaled and subtracted file
        output_name  <-  paste(substr(currentWIG,1,nchar(currentWIG)-4), "_Scaled",sep="")
        if( subtractInput)
        {
            if( !is.na(inputFile)){
                output_name  <-  paste( output_name, "_BGSub",sep="")
            }
        }
        cat("|---|-- Saving result :", output_name, "\n")
        writeWIG(dataWIGRescaled, output_name, folder=NULL, fixedStep=binSize)
    }
    
}




