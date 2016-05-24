peakDetection <-
function(IPdata, controlData, chrName, 
                        windowSize = 150, windowOverlap = 50,
                        outputName = "bPeaks_results", 
                        baseLineIP = NULL, baseLineControl = NULL,
                        IPthreshold = 6, controlThreshold = 4,
                        ratioThreshold = 2, 
                        averageThreshold = 0.7, peakDrawing = TRUE){


#######
# written by G. LELANDAIS <gaelle.lelandais@univ-paris-diderot.fr>    
#######


    # argument names were changed for user convenience
    INPUTdata = controlData
    INPUTthreshold = controlThreshold

    baseLineINPUT = baseLineControl

    if(peakDrawing == TRUE){
        genomeLandscape = 1
    }else{
        genomeLandscape = 0
    }

    # WARNING ! IP and INPUT data can start and end with NA values (because of the
    # smoothing procedure

    # number of sliding windows
    windowNumber = floor((length(IPdata) - windowSize) / windowOverlap)     

    startPos = 1
    endPos   = windowSize

    # table to store mean values
    meanValues = NULL

    # to print percentages during the analysis
    pourValues = c(1, 5, 10, 20, 50, 60, 70, 80, 90, 100)
    windowVal = floor((pourValues * windowNumber)/100)

    for(i in 1:(windowNumber + 1)){
         
        if(i %in% windowVal){
           print(paste(pourValues[windowVal == i], "% of windows were analyzed", sep = "")) 
        }

        # calculate mean values (coverage per base)
        meanIP    = mean(IPdata[startPos:endPos], na.rm = T)
        meanINPUT = mean(INPUTdata[startPos:endPos], na.rm = T)

        # ratio between IP and INPUT values
        meanRatio = log2((meanIP + 1)/(meanINPUT + 1))
		
        # criterion to evaluate signal intensities (IP and INPUT)
        meanAverage = (log2(meanIP + 1) + log2(meanINPUT + 1))/2

        meanValues = rbind(meanValues, c(startPos, endPos, 
                                        meanIP, meanINPUT, meanRatio, meanAverage)) 
 
        # change genomic positions
        startPos = startPos + windowOverlap
        endPos   = startPos + windowSize
    }

    colnames(meanValues) = c("Start", "End", "IP", "Control", "log2FC", "averageLog2")

    # calculate threshold values
    thresholdINPUT = INPUTthreshold * baseLineINPUT
    thresholdIP    = IPthreshold * baseLineIP

    if(averageThreshold > 1){
        print("problem with averageThreshold parameter (should be [0,1]), value is set to 0.7")
        thresholdAverage = quantile(meanValues[,6], probs = 0.7)

    }else{
        thresholdAverage = quantile(meanValues[,6], probs = averageThreshold)
    }

    # number of detected significant windows
    winNumber = sum((meanValues[,4] < thresholdINPUT) & 
           (meanValues[,3] > thresholdIP) & 
           (meanValues[,5] > ratioThreshold) &
           (meanValues[,6] > thresholdAverage))

        # in case of detected windows
        if(winNumber > 0){

        print(paste(winNumber, "significant window(s) were detected..."))

        seedPositions = meanValues[(meanValues[,4] < thresholdINPUT) & 
                                   (meanValues[,3] > thresholdIP) &
                                   (meanValues[,5] > ratioThreshold) &
                                   (meanValues[,6] > thresholdAverage),]

        # merging procedure if at least two windows were detected
        if(winNumber > 1){

            print("... starting merging procedure")
            finalSeeds = NULL

            startP = seedPositions[1,1]
            endP   = seedPositions[1,2]
            startL = 1
            endL   = 1

            for(i in 2:nrow(seedPositions)){
            
                if(seedPositions[i,1] < endP){
                    endP = seedPositions[i,2]
                    endL = i
                }else{
                    finalSeeds = rbind(finalSeeds, 
                                    c(startP, endP, mean(seedPositions[startL:endL, 3]), 
                                                    mean(seedPositions[startL:endL, 4]),
                                                    mean(seedPositions[startL:endL, 5]),
                                                    mean(seedPositions[startL:endL, 6]))) 
                    startP = seedPositions[i,1]
                    endP   = seedPositions[i,2]

                    startL = i
                    endL   = i
                # end of else()
                }
            # end of for()
            }
    
            # last peak has to be written
            finalSeeds = rbind(finalSeeds, 
                               c(startP, endP, mean(seedPositions[startL:endL, 3]), 
                                               mean(seedPositions[startL:endL, 4]),
                                               mean(seedPositions[startL:endL, 5]),
                                               mean(seedPositions[startL:endL, 6]))) 
     
        # check if the merging procedure was OK 
        if(!is.null(finalSeeds)){

          if(is.null(nrow(finalSeeds))){
                finalSeeds = t(as.matrix(finalSeeds))
          }

        }
        else{

            finalSeeds = seedPositions
        }

    # end of if(winNumber)
    }else{
        # in case of one peak
        seedPositions = t(as.matrix(seedPositions))
        finalSeeds = seedPositions
    }
    
    colnames(finalSeeds) = c("Start", "End", "IP", "Control", "log2FC", "averageLog2")    
    print("")
    print(paste("# of detected basic peaks (bPeaks) : ", nrow(finalSeeds)))
    print("")

    print("** Saving chromosome information in PDF file:") 
    print(paste(outputName, "_dataSummary.pdf", sep = ""))
    print("")
    pdf(paste(outputName, "_dataSummary.pdf", sep = ""))
 
    # graphical representations
    par(mfrow = c(2,2))

    # parameter summary
    vecInfo = c(paste("Chromosome name:\n", chrName),
                paste("Window size:", windowSize),
                paste("Window overlap:", windowOverlap),
                paste("IP threshold (T1): > x", IPthreshold),
                paste("Control threshold (T2): < x", INPUTthreshold),
                paste("Log2FC threshold (T3): > ", ratioThreshold),
                paste("Average log2 signals ([0,1]) (T4): >", averageThreshold))

    plot(1, axes = F, xlab = "", ylab = "", main = "Parameter summary", 
            col = "white", ylim = c(0,10), xlim = c(0, 30))
    
    text(rep(15,7), c(9, 7, 6, 5, 4, 3, 2), labels = vecInfo)

    # result summary
    vecInfo = c(paste("Chromosome name:\n", chrName),
                paste("# of significant windows \n(before the merging procedure):", nrow(seedPositions)),
                paste("# of detected bPeaks \n(after the merging procedure):", nrow(finalSeeds)),
                paste("Bed file name:\n", paste(outputName, ".bed", sep = "")))
 
    plot(1, axes = F, xlab = "", ylab = "", main = "Result summary", 
            col = "white", ylim = c(0,10), xlim = c(0, 40))
    
    text(rep(20,4), c(9, 7, 5, 3), labels = vecInfo)


    # IP and INPUT signals
 #   boxplot(meanValues[,3:4], ylab = "signal intensity", main = paste("Chromosome:\n", chrName),
 #           col = c("red", "blue"))

 #   hist(meanValues[,5], main = "Log2FC", xlab = "log2FC\n(average values, sliding windows)",
 #       col = "orange")
    
 #   hist(meanValues[,6], main = "Average log2 signals", 
 #                       xlab = "(log2(IP) + log2(control)) / 2\n(average values, sliding windows)",
 #       col = "purple")

    plot(meanValues[,4], meanValues[,3], pch = 20, main = paste("Chromosome:\n", chrName), 
         ylab = "T1: IP signal", xlab = "T2: Control signal")

    # to be selected a "good" window must have : 
    # 1) low INPUT signal (< baseLineINPUT)
    # 2) high IP signal (> IPthreshold * baseLineIP)
    # 3) log2(IP/INPUT) > ratioThreshold
    # 4) Average log2 > 

    abline(v = thresholdINPUT, col = "red")
    abline(h = thresholdIP, col = "red")

    points(seedPositions[,4], seedPositions[,3], col = "red", pch = 20)

    # MA plot representation
    plot(meanValues[,6], meanValues[,5], xlab = "T4: (log2(IP) + log2(control)) / 2", ylab = "T3: log2FC",
        main = paste("Chromosome:\n", chrName), pch = 20)  

    abline(v = thresholdAverage, col = "red")
	abline(h = ratioThreshold, col = "red")
	
    points(seedPositions[,6], seedPositions[,5], col = "red", pch = 20)
    
    # end of PDF dataSummary
    dev.off()

    # create bed file
    print("** Bed file saving in:")
    print(paste(outputName, ".bed", sep = ""))
    print("")
    bedData = finalSeeds[order(finalSeeds[,5], decreasing = T),]

    # in case of one peak
    if(is.null(nrow(bedData))){

        bedData = c(chrName, bedData[1:2], 
                    paste(chrName, "_bPeaks_", nrow(finalSeeds), sep = ""),
                    paste(bedData[3], "|", bedData[4], "|", bedData[5], "|", bedData[6],
                    sep = ""))

        bedData = t(as.matrix(bedData))

    }else{

        bedData = cbind(rep(chrName, nrow(finalSeeds)), bedData[,1:2], 
                        paste(chrName, "_bPeak_", 1:nrow(finalSeeds), sep = ""),
                        paste(bedData[,3], "|", bedData[,4], "|", bedData[,5], "|", bedData[,6],
                        sep = ""))

    }

    write.table(bedData, file = paste(outputName, ".bed", sep = ""), quote = F,
                sep = "\t", row.names = F, col.names = F)  

    # result visualization
    
    if(genomeLandscape == 1){

    print("** Peak drawing in PDF file:")
    print(paste(outputName, "_bPeaksDrawing.pdf", sep = ""))

    pdf(paste(outputName, "_bPeaksDrawing.pdf", sep = ""))

    for(i in 1:nrow(bedData)){

        peakDrawing(vecIP = IPdata, vecControl = INPUTdata, 
                      lineIP = thresholdIP, lineControl = thresholdINPUT, 
                      lineFC = ratioThreshold,
                      lineAverage = thresholdAverage,
                      posInf = as.numeric(bedData[i,2]), 
                      posSup = as.numeric(bedData[i,3]), add = 20, title = bedData[i,4])

    }

    # end of PDF bPeaksDrawing()
    dev.off()
    
    # end of if() genome landscape
    }

    # end of if()
    }else{
 
    ################### Duplicated CODE #############################"   
    print("** Saving chromosome information in PDF file:") 
    print(paste(outputName, "_dataSummary.pdf", sep = ""))
    print("")
    pdf(paste(outputName, "_dataSummary.pdf", sep = ""))
 
    # graphical representations
    par(mfrow = c(2,2))

    # parameter summary
    vecInfo = c(paste("Chromosome name:\n", chrName),
                paste("Window size:", windowSize),
                paste("Window overlap:", windowOverlap),
                paste("IP threshold (T1): > x", IPthreshold),
                paste("Control threshold (T2): < x", INPUTthreshold),
                paste("Log2FC threshold (T3): > ", ratioThreshold),
                paste("Average log2 signals ([0,1]) (T4): >", averageThreshold))

    plot(1, axes = F, xlab = "", ylab = "", main = "Parameter summary", 
            col = "white", ylim = c(0,10), xlim = c(0, 30))
    
    text(rep(15,7), c(9, 7, 6, 5, 4, 3, 2), labels = vecInfo)

    # IP and INPUT signals
#    boxplot(meanValues[,3:4], ylab = "signal intensity", main = paste("Chromosome:\n", chrName),
#            col = c("red", "blue"))

#    hist(meanValues[,5], main = "Log2FC", xlab = "log2(IP / control)\n(average values, sliding windows)",
#        col = "orange")
    
#    hist(meanValues[,6], main = "Average log2 signals", 
#                        xlab = "(log2(IP) + log2(control)) / 2\n(average values, sliding windows)",
#        col = "purple")

    # parameter summary
    vecInfo = c(paste("Chromosome name:\n", chrName),
                paste("# of significant windows \n(before the merging procedure):", 0),
                paste("# of detected bPeaks \n(after the merging procedure):", 0))
                
    plot(1, axes = F, xlab = "", ylab = "", main = "Result summary", 
            col = "white", ylim = c(0,10), xlim = c(0, 40))
    
    text(rep(20,4), c(9, 7, 5), labels = vecInfo)


    plot(meanValues[,4], meanValues[,3], pch = 20, main = paste("Chromosome:\n", chrName), 
         ylab = "T1: IP signal", xlab = "T2: Control signal")

    # to be selected a "good" window must have : 
    # 1) low INPUT signal (< baseLineINPUT)
    # 2) high IP signal (> IPthreshold * baseLineIP)
    # 3) log2(IP/INPUT) > ratioThreshold
    # 4) Average log2 > 

    abline(v = thresholdINPUT, col = "red")
    abline(h = thresholdIP, col = "red")

    # MA plot representation
    plot(meanValues[,6], meanValues[,5], xlab = "T4: (log2(IP) + log2(control)) / 2", ylab = "T3: log2(IP / control)",
        main = paste("Chromosome:\n", chrName), pch = 20)  

    abline(v = thresholdAverage, col = "red")
	abline(h = ratioThreshold, col = "red")
    
    # end of PDF dataSummary
    dev.off()

   ################### Duplicated CODE #############################"   
         
        print("Sorry, no bPeak was detected")

        finalSeeds = NULL
    }

    return(finalSeeds)

# end of function peakDetection()
}
