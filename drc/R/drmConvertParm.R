"drmConvertParm" <- 
function(startVec, startMat, factor1, colList)
{
#print(startMat)
    startMat2 <- startMat
#    print(factor1)
    if (length(unique(factor1)) == 1) {return(startVec)}
    
    mmat <- model.matrix(~factor(factor1) - 1)
#    print(dim(mmat))
        
    pm <- list()
    for (i in 1:length(colList))
    {
        clElt <- colList[[i]]
        ncclElt <- dim(clElt)[2]   
            
        indVec <- !is.na(startMat2[, i, drop = FALSE])
        indVal <- min(c(sum(indVec), dim(clElt)[2]))
#        print(indVec)
#        print(indVal)
        
        indVec2 <- (1:ncclElt)[indVec]
        if (length(indVec2) > ncclElt) {indVec2 <- 1:ncclElt}
#        print(indVec2)
#        pm[[i]] <- (ginv(t(clElt)%*%clElt)%*%t(clElt))[indVec2, ,drop=FALSE]%*%mmat[,indVec]%*%startMat2[indVec, i, drop=FALSE]
# then cabanne works
        
        
#        print( ((ginv(t(clElt)%*%clElt)%*%t(clElt))[, ,drop=FALSE]%*%mmat[,indVec]%*%startMat2[indVec, i, drop=FALSE])[1:length(indVec2)] )

         pm[[i]] <- (ginv(t(clElt)%*%clElt)%*%t(clElt))[1:indVal, ,drop = FALSE]%*%mmat[,indVec]%*%startMat2[indVec, i, drop = FALSE]
#            print( (ginv(t(clElt)%*%clElt)%*%t(clElt))[1:indVal, ,drop = FALSE]%*%mmat[,indVec]%*%startMat2[indVec, i, drop = FALSE] )
#        print(pm[[i]])
            
                        
#        if ((length(posVec[[i]])>0) && !(upperPos==i)) {pm[[i]] <- pm[[i]][-pos]}
    }  
    tempVec <- unlist(pm)
    tempVec <- tempVec[!is.na(tempVec)]
#    print(tempVec)

    
    ## Checking whether the intercept column has been removed
    indVec3 <- ( abs(tempVec)<1e-10 )
    if (any(indVec3)) 
    {
#        tempVec[indVec3] <- startVec[indVec3]
        tempVec <- startVec
    }

    return(tempVec)
}
