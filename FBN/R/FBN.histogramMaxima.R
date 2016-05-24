FBN.histogramMaxima <-
function(inputData, minSpan = .2, breaksData = NULL){
    if(is.null(inputData) ){
        cat("WARNING: hist.profile -> Please input a valid inputData\n") 
        return(NULL)
    }
    maximas = c()
    idxMaximas = c()
    if(is.null(breaksData) ){
	    histData = hist(inputData, plot = FALSE, breaks = "FD")
	} else
	    histData = hist(inputData, plot = FALSE, breaks = breaksData)
    if(minSpan < 0){
	    filteredHistData = medianFilter(inputData = histData$counts, windowSize = 5)
    	filteredHistData = meanFilter(inputData = filteredHistData, windowSize = 9)
    	histData$counts = filteredHistData
    	rm(filteredHistData)
    }
    inputValues = histData$mids
    deltaHist = histData$counts[2:length(histData$counts)] - histData$counts[1:(length(histData$counts)-1)]
    deltaHist[deltaHist < 0] = -1
    deltaHist[deltaHist > 0] = 1
    i = 1
    while(i < length(deltaHist) ){
        if(deltaHist[i] != 1){ #constant or descending region
            i = i+1
            next
        }
        j = i+1
        if(deltaHist[j] == 1){ #an ascending region
            i = j
            next
        }
        #if here than it might be a clear maximum or a flat maximum region...
        while( (deltaHist[j] == 0)&&(j < length(deltaHist) ) ){
            j = j+1
        }
        if(deltaHist[j] == -1){ #clear maximum found /\
            maximas = c(maximas, inputValues[j])
            idxMaximas = c(idxMaximas, j)
            #maximas = c(maximas, inputValues[i+floor( (j-i)/2)])
            #idxMaximas = c(idxMaximas, i+floor( (j-i)/2)+2)
        }
        i = j
    }
    outlier = 1
    while((length(maximas) >= 2)&&(outlier != 0)){
    	outlier = 0
		deltaMaximas = maximas[2:length(maximas)] - maximas[1:(length(maximas)-1)]
		idx = which(deltaMaximas < abs(minSpan))
		if(length(idx) > 0){
			if(histData$counts[idxMaximas[idx[1]]] < histData$counts[idxMaximas[idx[1]+1]])
				outlier = idx[1]
			else
				outlier = idx[1]+1
			if(outlier == 1){
				maximas = maximas[2:length(maximas)]
				idxMaximas = idxMaximas[2:length(idxMaximas)]
			}
			else{
				if(outlier == length(maximas)){
					maximas = maximas[1:(length(maximas)-1)]
					idxMaximas = idxMaximas[1:(length(idxMaximas)-1)]
				}
				else{
					maximas = c(maximas[1:(outlier-1)], maximas[(outlier+1):length(maximas)])
					idxMaximas = c(idxMaximas[1:(outlier-1)], idxMaximas[(outlier+1):length(idxMaximas)])
				}
			}
		}
	}
	backMaximas = maximas	
	#Check if they are too many maximas -> apply the smoothing filter to reduce the noise...
	if((length(maximas) > 8)&&(minSpan >= 0)){
	    maximas = FBN.histogramMaxima(inputData, minSpan = -0.1, breaksData = breaksData)
	}
	if(length(maximas)==0)
		return(backMaximas)
	#for(i in 1:length(maximas)){
	#	tmp = abs(inputData - maximas[i])
	#	idx = which(tmp == min(tmp))
	#	maximas[i] = inputData[idx[1]]
	#}
    return(maximas)
}

