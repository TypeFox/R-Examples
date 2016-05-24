
bivariateANM <- function(X, parentsOf = dim(X)[2], variableSelMat = NULL, silent = TRUE)
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
                ititot <- indtestHsic(resitot,X[,i])$p.value
                resttoi <- train_gam(X[,currentTarget],X[,i])$residuals
                itttoi <- indtestHsic(resttoi,X[,currentTarget])$p.value
                if(!silent)
                {
                    cat("\n\nmodel", i, " -> ", currentTarget, ":\nHSIC score (p-val, i.e. large is good):",ititot)
                    cat("\n\nmodel", currentTarget, " -> ", i, ":\nHSIC score (p-val, i.e. large is good):",itttoi)
                }
                
                scoreMat[i,t] <- itttoi
                if( (ititot > 0.1) && (itttoi < (0.05/(p-1)) ) )
                {
                    causalParents[i, t] <- TRUE    
                }
                
                if( i %in% parentsOf )
                {
                    iTarget <- which(parentsOf == i)

                    scoreMat[currentTarget,iTarget] <- ititot
                    if( (itttoi > 0.1) && (ititot < (0.05/(p-1)) ) )
                    {
                        causalParents[currentTarget, iTarget] <- TRUE    
                    }
                }
            }
        }
    }
    
    
   
    return(list(causalParents = causalParents, scoreMat = scoreMat))
}

