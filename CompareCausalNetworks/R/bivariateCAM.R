
bivariateCAM <- function(X, parentsOf = dim(X)[2], variableSelMat = NULL, silent = TRUE)
{
   
    
    p <- ncol(X)
    causalParents <- matrix(FALSE, p, length(parentsOf))
    scoreMat <- matrix(NA, p, length(parentsOf))
    for(t in 1:length(parentsOf))
    {
        currentTarget <- parentsOf[t]
        if( !is.null(variableSelMat) )
        {
            possibleParents <- which(variableSelMat[,currentTarget])
        } else
        {
            possibleParents <- (1:p)[-currentTarget]
        }
        
        for(i in possibleParents)
        {
            if(is.na(scoreMat[i,t]))
            {
                resitot <- train_gam(X[,i],X[,currentTarget])$residuals
                scitot <- - log(var(X[,i])) - log(var(resitot))
                resttoi <- train_gam(X[,currentTarget],X[,i])$residuals
                scttoi <- - log(var(X[,currentTarget])) - log(var(resttoi)) # large score is good

                if(!silent)
                {
                    cat("\n\nmodel", i, " -> ", currentTarget, ":\nCAM score (large is good):",scitot)
                    cat("\n\nmodel", currentTarget, " -> ", i, ":\nCAM score (large is good):",scttoi)
                }
                
                scoreMat[i,t] <- scitot - scttoi
                if( (scitot - scttoi) > 0 )
                {
                    causalParents[i, t] <- TRUE    
                }

                if( i %in% parentsOf )
                {
                    iTarget <- which(parentsOf == i)
                    scoreMat[currentTarget,iTarget] <- scttoi - scitot
                    if( (scttoi - scitot) > 0 )
                    {
                        causalParents[currentTarget, iTarget] <- TRUE    
                    }                    
                }
            }
        }
    }
    
    
  
    return(list(causalParents = causalParents, scoreMat = scoreMat ))
}
