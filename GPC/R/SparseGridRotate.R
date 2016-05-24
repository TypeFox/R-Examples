SparseGridRotate <- function(inputCoord, inputJList, dStart){
    d       <- length(inputCoord)
    counter <- 0
    for (ii in dStart:d){
        if (inputCoord[ii] !=0){
            I <- diag(1,d)
            I[ii,ii] <- -1
            if (counter == 0){
                rotatedNodes <- inputCoord%*%I
                rotatedJList <- inputJList%*%I
            } else {
                rotatedNodes <- rbind(rotatedNodes,inputCoord%*%I)
                rotatedJList <- rbind(rotatedJList,inputJList%*%I)                
            }
            counter <- counter + 1;            
            if (ii < d){               
                RotatedGrid <- SparseGridRotate(rotatedNodes[counter,],rotatedJList[counter,],ii+1);                
                if (RotatedGrid$counter != 0){
                    rotatedNodes <- rbind(rotatedNodes,RotatedGrid$RotatedNodes)
                    rotatedJList <- rbind(rotatedJList,RotatedGrid$RotatedJList)
                    counter <- counter + RotatedGrid$counter
                } 
            }   
        } 
    } 
    if (counter==0){ res <- list(counter=counter) } 
    else { res <- list(RotatedNodes=rotatedNodes,RotatedJList=rotatedJList,counter=counter) }
    return(res)
}