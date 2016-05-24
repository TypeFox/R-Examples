peakLocation <-
function(bedFile, cdsPositions, withoutOverlap = FALSE,
           outputName = "bPeaksLocation", promSize = 800){

#######
# written by G. LELANDAIS <gaelle.lelandais@univ-paris-diderot.fr>    
#######


    # argument names were changed for user convenience
    ORFinfo = cdsPositions

    print("**********************************************")
    # data reading
    print("Opening BED file with peak information:")
    print(bedFile) 
    peakRes = as.matrix(read.table(bedFile, sep = "\t")) 
    print("**********************************************")

    # to store all the results    
    locResProm  = NULL
    locResCDS   = NULL
 
    # statistics
    numPeaks = nrow(peakRes)
    peakInProm  = NULL
    peakInCDS   = NULL

    print("")
    print("Starting peak location regarding ORF/CDS positions...")
    print("")
 
    # for each detected peak....
    for(i in 1:nrow(peakRes)){

        currentPeakRes = NULL

        # chromosome
        chrmInfo = peakRes[i,1]
        peakStart = as.numeric(peakRes[i,2])
        peakEnd   = as.numeric(peakRes[i,3])        
        # the middle position will be used to map the peak on the genome
        middlePos = mean(c(peakEnd, peakStart))

	# in CDS ?
	inCDS = 0

        # get ORF on the analysed chromosome
        subORFinfo = ORFinfo[ORFinfo[,1] == chrmInfo,]
    
        if(nrow(subORFinfo) > 0){
            
            promDirect = (middlePos < as.numeric(subORFinfo[,2])) & 
                  (middlePos > (as.numeric(subORFinfo[,2]) - promSize)) &
		  (subORFinfo[,4] == "W")

            promIndirect = (middlePos > as.numeric(subORFinfo[,3])) & 
                  (middlePos < (as.numeric(subORFinfo[,3]) + promSize)) &
		  (subORFinfo[,4] == "C")
                
            ORFloc = (middlePos > as.numeric(subORFinfo[,2])) &
                     (middlePos < as.numeric(subORFinfo[,3]))  

            numElements = sum(promDirect) + sum(promIndirect) + sum(ORFloc)

       #     print(paste("number of detected elements for", 
                       # peakRes[i,4], "=", numElements))    

            if(sum(ORFloc) > 0){
                locResCDS = c(locResCDS, paste(peakRes[i,4], subORFinfo[ORFloc,5], sep = "\tin CDS\t"))
		# the detected peak is stored in a vector
		peakInCDS = c(peakInCDS, peakRes[i,4]) 
		# peak is in CDS...
	        if(withoutOverlap == TRUE){
			inCDS = 1
		 # end of if()
		 }
            }

	    if(inCDS == 0){
	
	            if(sum(promDirect) > 0){
	                locResProm = c(locResProm, paste(peakRes[i,4], subORFinfo[promDirect,5], 
						sep = "\tpromoter (W strand)\t"))
			# the detected peak is stored in a vector
			peakInProm = c(peakInProm, peakRes[i,4])
	            }
    
	            if(sum(promIndirect) > 0){
	                locResProm = c(locResProm, paste(peakRes[i,4], subORFinfo[promIndirect,5], 
						sep = "\tpromoter (C strand)\t"))
			# the detected peak is stored in a vector
			peakInProm = c(peakInProm, peakRes[i,4])		
	            }
            # end of if()      
	    }
	    
	    # before analyzing a new peak...		
	    inCDS = 0

        }else{
            print(paste("WARNING : chromosome", chrmInfo, "is unknown..."))
        }

    # end of for()
    }
   
	numProm = length(unique(peakInProm))
	numCDS  = length(unique(peakInCDS))

    # statistics printing
    print(paste("# of analyzed peaks: ", numPeaks))
    if(withoutOverlap == FALSE){
	    print(paste("# of peaks UPSTREAM annotated CDS : ", numProm))
    }
    if(withoutOverlap == TRUE){
	    print(paste("# of peaks UPSTREAM annotated CDS (without overlap) : ", numProm))
    }

    print(paste("# of peaks IN annotated CDS : ", numCDS)) 

    print("")

    pdf(paste(outputName, "_peakLocation.pdf", sep = ""))
    barplot(c(numPeaks, numProm, numCDS), 
		names = c("all peaks", paste("in promoters \n(", promSize, " bp)", sep = ""), "in CDS"),
            col = c("grey", "yellow", "green"),
            main = paste("Peak location regarding CDS positions\n# of analyzed peaks:", numPeaks,
			 "\nwithout overlap (prom/CDS) =", withoutOverlap))    
    dev.off()

    print("Saving the results in:")
    print(paste(outputName, "_peakLocation_inPromoters.txt", sep = ""))
 
    # writing the results
    write.table(locResProm, file = paste(outputName, "_peakLocation_inPromoters.txt", sep = ""), quote = F, 
                row.names = F, col.names = F)

    print(paste(outputName, "_peakLocation_inCDS.txt", sep = ""))
 
    # writing the results
    write.table(locResCDS, file = paste(outputName, "_peakLocation_inCDS.txt", sep = ""), quote = F, 
                row.names = F, col.names = F)

    statsRes = list(numPeaks =  numPeaks, upFeatures = numProm, 
                inFeatures = numCDS)   

   print("**********************************************")

   return(statsRes)


# end of function peakLocation()
}
