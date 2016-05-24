


SPCALDA = function(X,Y,rho=exp(c((-2):6)),K=min(20,min(dim(X))), folds = NULL){   # input X=n\times p data matrix, Y label, rho = values of tuning parameters, K = number of output PC
     R = length(rho)
     ErrorCount = matrix(0, R, K)
     for (i in 1:length(folds)) {
        ind = folds[[i]]
        x = X[-ind,]; y = Y[-ind]; xt = X[ind,]; yt = Y[ind]
        RotationMatrix = RRRotation(X=x, Y=y, rho=rho, K=K)
        for (j in 1:R) {
           currentRotation = RotationMatrix[[j]]
           # rotate the data and apply nestedLDA
           Errors = nestedLDA(x%*%currentRotation,y,xt%*%currentRotation,yt)
           ErrorCount[j,] = ErrorCount[j,]+Errors
        }
     }
     minerror = min(ErrorCount); tuneposition = tiebreaker(ErrorCount)
     tunerho=rho[tuneposition[1]]; tuneK = tuneposition[2]
     tuneRotation = RRRotation(X, Y, tunerho, tuneK)[[1]]
     transfX = X%*%tuneRotation
     return(list(ob=lda(transfX,Y), tuneRotation=tuneRotation, ErrorCount=ErrorCount, minerror = minerror, rho = tunerho, K = tuneK ))

}

 