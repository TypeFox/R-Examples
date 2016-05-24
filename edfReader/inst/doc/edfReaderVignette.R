## ------------------------------------------------------------------------
libDir <- system.file ("extdata", package="edfReader")
AFile <- paste (libDir, '/edfAnnonC.edf', sep='') # a file with 2 annotation signals
BFile <- paste (libDir, '/bdfPlusC.bdf' , sep='') # a continuously recorded BDF file
CFile <- paste (libDir, '/edfPlusC.edf' , sep='') # a continuously recorded EDF file
DFile <- paste (libDir, '/edfPlusD.edf' , sep='') # a discontinuously recorded EDF file

## ------------------------------------------------------------------------
require (edfReader)
AHdr  <- readEdfHeader (AFile)
BHdr  <- readEdfHeader (BFile)
CHdr  <- readEdfHeader (CFile)            
DHdr  <- readEdfHeader (DFile)                  

## ------------------------------------------------------------------------
BHdr
summary (AHdr)

## ------------------------------------------------------------------------
AHdr$sHeader
summary (CHdr$sHeader)

## ------------------------------------------------------------------------
ASignals <- readEdfSignals (AHdr)
ASignals

## ------------------------------------------------------------------------
DSignals <- readEdfSignals (DHdr)
DSignals

## ------------------------------------------------------------------------
DSignalsF <- readEdfSignals (DHdr, fragments = TRUE)
DSignalsF

## ------------------------------------------------------------------------
CSignals8 <- readEdfSignals (CHdr, signals=c(8, "8", "sine 8.5 Hz"))
CSignals8

## ------------------------------------------------------------------------
ASignalsPeriod    <- readEdfSignals (AHdr, from=0.7, till=1.8)
ASignalsPeriod

## ------------------------------------------------------------------------
CSignals <- readEdfSignals (CHdr)
summary (CSignals$pulse)         # edfReader names signals after their label

## ------------------------------------------------------------------------
CDSignals <- readEdfSignals (DHdr, from=5.1, till=18)
FDSignals <- readEdfSignals (DHdr, fragments=TRUE)

## ------------------------------------------------------------------------
summary (FDSignals$`sine 8.5 Hz`)         # note the "`" quotes for a name with spaces.

## ------------------------------------------------------------------------
CSignals$`EDF Annotations`
summary(ASignalsPeriod$`EDF Annotations`)

## ------------------------------------------------------------------------
format (AHdr$startTime, format="%Y-%m-%d %H:%M:%OS3",  usetz = FALSE)
ASignalsPlusStartTimes <- readEdfSignals(AHdr, signals='EDF Annotations-1', recordStarts=TRUE)
annots <- ASignalsPlusStartTimes$annotations
annots[annots$isRecordStart==TRUE,'onset'][1]

## ------------------------------------------------------------------------
str (CHdr,  max.level=1)

## ------------------------------------------------------------------------
str (CHdr$sHeader, max.level=1)

## ------------------------------------------------------------------------
str(CSignals$pulse, max.level=1) 

## ------------------------------------------------------------------------
str(FDSignals$`sine 8.5 Hz`, max.level=1) 

## ------------------------------------------------------------------------
str(FDSignals$`sine 8.5 Hz`$fragments[[1]], max.level=1) 

## ------------------------------------------------------------------------
str(ASignals$`EDF Annotations`, max.level=1) 

## ------------------------------------------------------------------------
str(ASignals$`EDF Annotations`$annotations, max.level=1) 

## ---- fig.width=7.2, fig.height=2.5--------------------------------------
plotEdfSignals <- function (signals,labels, from=0, till=Inf) {
    nLabels <- length (labels)
    sRate   <- numeric (length = nLabels)
    fromS   <- integer (length = nLabels)
    tillS   <- integer (length = nLabels)
    sLength <- integer (length = nLabels)
    for (i in 1:nLabels) {
        sRate[i]    <- signals[[labels[i]]]$sRate
        fromS[i]    <- ceiling (sRate[i] * max (from, 0)) +1
        tillS[i]    <- ceiling (sRate[i] * till)
        tillS[i]    <- min (tillS[i], length(signals[[labels[i]]]$signal))
        sLength[i]  <- tillS[i] - fromS[i] + 1 
    }
    totLength  <- sum (sLength)
    cat (" totLength=",  totLength)
    time    <- numeric   (length = totLength)
    signal  <- numeric   (length = totLength)
    label   <- character (length = totLength)
    from <- 1
    for (i in 1:nLabels) {
        till <- from + sLength[i] - 1
        time  [from:till]   <- seq (from=fromS[i]-1, to=(tillS[i]-1)) / sRate[i]
        signal[from:till]   <- signals[[labels[i]]]$signal[fromS[i]:tillS[i]]
        label [from:till]   <- rep(labels[i], sLength[i])
        from <- till + 1
    }
    cat (" | from-1=", from-1,'\n')
    
    ggplotDF <- data.frame (time=time, signal=signal, label=label)
    ggplot (ggplotDF, aes(x=time, y=signal, colour=label)) + geom_line()
}

if (require(ggplot2)) {
    CSignals <- readEdfSignals (CHdr)
    plotEdfSignals (CSignals, labels=c('sine 8.5 Hz', 'sine 50 Hz'), from=.2, till=0.5)
}

