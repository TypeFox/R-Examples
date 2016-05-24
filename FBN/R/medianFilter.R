medianFilter <-
function(inputData=NULL, windowSize=3){
    if(is.null(inputData) ){
        cat("WARNING: medianFilter -> Please input a valid inputData\n") 
        return(NULL)
    }
    if(windowSize <= 1){
        cat("WARNING: medianFilter -> no filtering performed: outData = inputData\n")
        return(inputData)
    }
    outData = vector(mode = mode(inputData), length = length(inputData) )
    for(i in 1:length(inputData) ){
        if(windowSize/2 == floor(windowSize/2) ){
            startTmp = max(1, i-windowSize/2)
            endTmp = min(length(inputData), i+windowSize/2-1)
        } else{
            startTmp = max(1, i-(windowSize-1)/2)
            endTmp = min(length(inputData), i+(windowSize-1)/2)
        }
        outData[i] = median(inputData[startTmp:endTmp])
    }
    return(outData)
}

