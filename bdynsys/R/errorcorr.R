# this scripts checks whether there are problematically high error correlations
# which will make it necessary to shift to the Bayesian only seemingly unrelated 
# regression procedure "baymodsur" (to be included in the package version 2.0)
# if the correlations are acceptable, the script re-calculates coefficients
# the user may compare how strong the corrected coefficients (under consideration of error correlations)
# deviate from the originally estimated and may decide to use the corrected 
# coefficiets in his/her result reports

# to call: errorcorr(pdata, 2, datap$logGDP, datap$EmanzV, f <- function(Y=c()) rbind(0.0012/Y[1]^2, + 0.0071*Y[1]^3), c(11), c(14), 2)
# with dx  = + 0.0012 /x^2 and dy = + 0.0071 x^3 
# call for indnr = 3: errorcorr(pdata, 3, datap$logGDP, datap$EmanzV, f <- function(Y=c()) rbind(0.0012/Y[1]^2, + 0.0071*Y[1]^3, 0.0025*Y[1]^2), c(29), c(34), 3, datap$DemocrH, c(28))
# call for indnr = 4: errorcorr(pdata, 4, datap$logGDP, datap$EmanzV, f <- function(Y=c()) rbind(0.0012/Y[1]^2, + 0.0071*Y[1]^3, 0.0025*Y[1]^2, + 0.123 + 0.076/Y[4]), c(83), c(90), 5, datap$DemocrH, c(82), datap$SeculTradV, c(1, 5)) 
errorcorr <- function(dataset, indnr, x, y, f, xterms, yterms, nrterms, z, zterms, v, vterms)
{
  if (indnr == 2)
  {
    procdata <- preprocess_data(indnr, x, y)
    xv <- procdata$allX
    yv <- procdata$allY
    chx <- procdata$tiChX
    chy <- procdata$tiChY
    
    # Specify the models    
    for (i in 1:length(xv))
    {
      Yprime <- f(rbind(xv[i], yv[i]))
      IncPredX <- Yprime[1]
      IncPredY <- Yprime[2]
    }
    
    # Calculating the actual errors
    errorX <- chx - IncPredX
    errorY <- chy - IncPredY
    
    errorXtmp <- errorX
    errorYtmp <- errorY
    
    errorX <- (errorXtmp - mean(errorXtmp, na.rm=TRUE))/sd(errorXtmp, na.rm=TRUE) 
    errorY <- (errorYtmp - mean(errorYtmp, na.rm=TRUE))/sd(errorYtmp, na.rm=TRUE)
    
    # Covariance of the errors
    covmat <- matrix(c(mean(errorX*errorX, na.rm=TRUE), mean(errorX*errorY, na.rm=TRUE), 
                      mean(errorY*errorX, na.rm=TRUE), mean(errorY*errorY, na.rm=TRUE)), 
                      nrow=2, byrow=TRUE)
    covmattmp <- covmat
    covmat <- solve(covmattmp)
    print(covmat)
    
    # computing the inverse of x and y
    invx <- xv^(-1)
    invy <- yv^(-1)
    
    # handling "Inf" as a result of inverse computation (/0), is this a problem?
    idx1 = which(is.infinite(invx))
    idx2 = which(is.infinite(invy)) 
    invx[idx1] <- NA
    invy[idx2] <- NA
    
    input <- cbind(rep(1, length(xv)), invx, invy, xv, yv, invx*invy, xv*invy, yv*invx,
                    xv*yv, xv^2, invx^2, yv^2, invy^2, xv^3, yv^3, invx^3, invy^3)
    Xterms <- input
    Yterms <- input
    
    # Selecting the terms from the models
    Xterms <- Xterms[, xterms]
    Yterms <- Yterms[, yterms]

    chmatr <- c(chx, chy)
    allincmatr <- chmatr
  
    # Reestimating the Beta coefficients
    xtkron <- matrix(c(covmat[1, 1]*t(Xterms), covmat[1, 2]*t(Xterms), 
                       covmat[2, 1]*t(Yterms), covmat[2, 2]*t(Yterms)), nrow=nrterms, byrow=TRUE)

    xtkronx <- matrix(c(covmat[1, 1]*t(Xterms)%*%Xterms, covmat[1, 2]*t(Xterms)%*%Yterms, 
                        covmat[2, 1]*t(Yterms)%*%Xterms, covmat[2, 2]*t(Yterms)%*%Yterms), 
                        nrow=nrterms, ncol=nrterms, byrow=TRUE)
    
    betapred <- ginv(xtkronx)%*%xtkron%*%allincmatr 
    print(betapred)
    #save(betapred, file = "Betapred.Rdata")
  }
  
  ########################## for indnr = 3 ####################################################################
  
  if (indnr == 3)
  {
    procdata <- preprocess_data(indnr, x, y, z)
    xv <- procdata$allX
    yv <- procdata$allY
    zv <- procdata$allZ
    chx <- procdata$tiChX
    chy <- procdata$tiChY
    chz <- procdata$tiChZ
    
    # Specify the models    
    for (i in 1:length(xv))
    {
      Yprime <- f(rbind(xv[i], yv[i], zv[i]))
      IncPredX <- Yprime[1]
      IncPredY <- Yprime[2]
      IncPredZ <- Yprime[3]
    }
    
    # Calculating the actual errors
    errorX <- chx - IncPredX
    errorY <- chy - IncPredY
    errorZ <- chz - IncPredZ
    
    errorXtmp <- errorX
    errorYtmp <- errorY
    errorZtmp <- errorZ
    
    errorX <- (errorXtmp - mean(errorXtmp, na.rm=TRUE))/sd(errorXtmp, na.rm=TRUE) 
    errorY <- (errorYtmp - mean(errorYtmp, na.rm=TRUE))/sd(errorYtmp, na.rm=TRUE)
    errorZ <- (errorZtmp - mean(errorZtmp, na.rm=TRUE))/sd(errorZtmp, na.rm=TRUE)
    
    # creating matrix with error covariances  
    covmat <- matrix(c(mean(errorX*errorX, na.rm=TRUE), mean(errorX*errorY, na.rm=TRUE), 
                      mean(errorX*errorZ, na.rm=TRUE), mean(errorY*errorX, na.rm=TRUE), 
                      mean(errorY*errorY, na.rm=TRUE), mean(errorY*errorZ, na.rm=TRUE), 
                      mean(errorZ*errorX, na.rm=TRUE), mean(errorZ*errorY, na.rm=TRUE), 
                      mean(errorZ*errorZ, na.rm=TRUE)), nrow=3, byrow=TRUE)
    covmattmp <- covmat
    covmat <- solve(covmat)
    print(covmat)
    
    # computing the inverse of x and y
    invx <- xv^(-1)
    invy <- yv^(-1)
    invz <- yv^(-1)
    
    # handling "Inf" as a result of inverse computation (/0), is this a problem?
    idx1 = which(is.infinite(invx))
    idx2 = which(is.infinite(invy)) 
    idx3 = which(is.infinite(invz)) 
    invx[idx1] <- NA
    invy[idx2] <- NA
    invz[idx3] <- NA
    
    input <- cbind(rep(1, length(xv)), invx, invy, invz, xv, yv, zv, invx*invy, invy*invz, invx*invz, 
                   xv*yv, yv*zv, xv*zv, xv*invy, yv*invx, xv*invz, zv*invx, yv*invz, zv*invy, xv*invy*invz,
                   yv*invx*invz, zv*invx*invy, xv*yv*invz, yv*zv*invx, xv*zv*invy, xv*yv*zv, invx*invy*invz, 
                   xv^2, invx^2, yv^2, invy^2, zv^2, invz^2, xv^3, yv^3, zv^3, invx^3, invy^3, invz^3)
    
    Xterms <- input
    Yterms <- input   
    Zterms <- input
    
    # Selecting the terms from the models
    Xterms <- Xterms[, xterms]
    Yterms <- Yterms[, yterms]
    Zterms <- Zterms[, zterms]
    
    chmatr <- c(chx, chy, chz)
    allincmatr <- chmatr
    
    # Reestimating the Beta coefficients    
    xtkron <- matrix(c(covmat[1, 1]*t(Xterms), covmat[1, 2]*t(Xterms), 
                      covmat[1, 3]*t(Xterms), covmat[2, 1]*t(Yterms), 
                      covmat[2, 2]*t(Yterms), covmat[2, 3]*t(Yterms),
                      covmat[3, 1]*t(Zterms), covmat[3, 2]*t(Zterms),
                      covmat[3, 3]*t(Zterms)), nrow=nrterms, byrow=TRUE)
    xtkronx <- matrix(c(covmat[1, 1]*t(Xterms)%*%Xterms, covmat[1, 2]*t(Xterms)%*%Yterms, 
                       covmat[1, 3]*t(Xterms)%*%Zterms, covmat[2, 1]*t(Yterms)%*%Xterms, 
                       covmat[2, 2]*t(Yterms)%*%Yterms, covmat[2, 3]*t(Yterms)%*%Zterms,
                       covmat[3, 1]*t(Zterms)%*%Xterms, covmat[3, 2]*t(Zterms)%*%Yterms,
                       covmat[3, 3]*t(Zterms)%*%Zterms), nrow=nrterms, ncol=nrterms, byrow=TRUE)

    betapred <- ginv(xtkronx)%*%xtkron%*%allincmatr
    print(betapred)
    #save(betapred, file = "Betapred.Rdata")   
  }

####################################### for indnr = 4 #####################################  
  
  if (indnr == 4)
  {
    procdata <- preprocess_data(indnr, x, y, z, v)
    xv <- procdata$allX
    yv <- procdata$allY
    zv <- procdata$allZ
    vv <- procdata$allV
    chx <- procdata$tiChX
    chy <- procdata$tiChY
    chz <- procdata$tiChZ
    chv <- procdata$tiChV
    
    # Specify the models    
    for (i in 1:length(xv))
    {
      Yprime <- f(rbind(xv[i], yv[i], zv[i], vv[i]))
      IncPredX <- Yprime[1]
      IncPredY <- Yprime[2]
      IncPredZ <- Yprime[3]
      IncPredV <- Yprime[4]
    }
    
    # Calculating the actual errors
    errorX <- chx - IncPredX
    errorY <- chy - IncPredY
    errorZ <- chz - IncPredZ
    errorV <- chv - IncPredV
    
    errorXtmp <- errorX
    errorYtmp <- errorY
    errorZtmp <- errorZ
    errorVtmp <- errorV
    
    errorX <- (errorXtmp - mean(errorXtmp, na.rm=TRUE))/sd(errorXtmp, na.rm=TRUE) 
    errorY <- (errorYtmp - mean(errorYtmp, na.rm=TRUE))/sd(errorYtmp, na.rm=TRUE)
    errorZ <- (errorZtmp - mean(errorZtmp, na.rm=TRUE))/sd(errorZtmp, na.rm=TRUE)
    errorV <- (errorVtmp - mean(errorVtmp, na.rm=TRUE))/sd(errorVtmp, na.rm=TRUE)
    
    # creating matrix with error covariance 
    covmat <- matrix(c(mean(errorX*errorX, na.rm=TRUE), mean(errorX*errorY, na.rm=TRUE), 
                       mean(errorX*errorZ, na.rm=TRUE), mean(errorX*errorV, na.rm=TRUE),
                       mean(errorY*errorX, na.rm=TRUE), mean(errorY*errorY, na.rm=TRUE), 
                       mean(errorY*errorZ, na.rm=TRUE), mean(errorY*errorV, na.rm=TRUE),
                       mean(errorZ*errorX, na.rm=TRUE), mean(errorZ*errorY, na.rm=TRUE), 
                       mean(errorZ*errorZ, na.rm=TRUE), mean(errorZ*errorV, na.rm=TRUE),
                       mean(errorV*errorX, na.rm=TRUE), mean(errorV*errorY, na.rm=TRUE),
                       mean(errorV*errorZ, na.rm=TRUE), mean(errorV*errorV, na.rm=TRUE)),
                     nrow=4, byrow=TRUE)
    covmattmp <- covmat
    covmat <- solve(covmat)
    print(covmat)
    
    # computing the inverse of x and y
    invx <- xv^(-1)
    invy <- yv^(-1)
    invz <- yv^(-1)
    invv <- vv^(-1)
    
    # handling "Inf" as a result of inverse computation (/0), is this a problem?
    idx1 = which(is.infinite(invx))
    idx2 = which(is.infinite(invy)) 
    idx3 = which(is.infinite(invz)) 
    idx4 = which(is.infinite(invv)) 
    invx[idx1] <- NA
    invy[idx2] <- NA
    invz[idx3] <- NA
    invv[idx4] <- NA
    
    input <- cbind(rep(1, length(xv)), invx, invy, invz, invv, xv, yv, zv, vv, invx*invy, invy*invz, invx*invz, 
                    invx*invv, invy*invv, invz*invv, xv*yv, yv*zv, xv*zv, xv*vv, yv*vv, zv*vv, xv*invy, yv*invx, 
                    xv*invz, zv*invx, yv*invz, zv*invy, xv*invv, vv*invx, yv*invv, vv*invy, zv*invv, vv*invz,
                    xv*invy*invz, yv*invx*invz, zv*invx*invy, vv*invx*invy, vv*invx*invz, vv*invy*invz, xv*invy*invv, 
                    xv*invz*invv, yv*invx*invv, yv*invz*invv, zv*invx*invv, zv*invy*invv, xv*yv*invz, yv*zv*invx, 
                    zv*xv*invy, xv*yv*invv, yv*zv*invv, zv*xv*invv, xv*vv*invz, yv*vv*invz, yv*vv*invx, zv*vv*invx, 
                    vv*xv*invy, vv*zv*invy, xv*yv*zv, xv*yv*vv, xv*vv*zv, vv*yv*zv, invx*invy*invz, invx*invy*invv, 
                    invx*invv*invz, invv*invy*invz, xv*invv*invy*invz, yv*invx*invv*invz, zv*invx*invy*invv, 
                    vv*invx*invy*invz, xv*yv*zv*invv, xv*yv*vv*invz, xv*vv*zv*invy, vv*yv*zv*invx, xv*yv*invv*invz,
                    xv*zv*invv*invy, xv*vv*invy*invz, yv*zv*invv*invx, yv*vv*invz*invx, zv*vv*invx*invy, 
                    invx*invy*invz*invv, xv*yv*zv*vv, xv^2, invx^2, yv^2, invy^2, zv^2, invz^2, vv^2, invv^2,
                    xv^3, yv^3, zv^3, vv^3, invx^3, invy^3, invz^3, invv^3)
    
    Xterms <- input    
    Yterms <- input   
    Zterms <- input    
    Vterms <- input
    
    # Selecting the terms from the models
    Xterms <- Xterms[, xterms]
    Yterms <- Yterms[, yterms]
    Zterms <- Zterms[, zterms]
    Vterms <- Vterms[, vterms]
    
    chmatr <-c(chx, chy, chz, chv)
    allincmatr <- chmatr
    
    # Reestimating the Beta coefficients
    xtkron <- matrix(c(covmat[1, 1]*t(Xterms), covmat[1, 2]*t(Xterms), 
                      covmat[1, 3]*t(Xterms), covmat[1, 4]*t(Xterms),
                      covmat[2, 1]*t(Yterms), covmat[2, 2]*t(Yterms), 
                      covmat[2, 3]*t(Yterms), covmat[2, 4]*t(Yterms),
                      covmat[3, 1]*t(Zterms), covmat[3, 2]*t(Zterms),
                      covmat[3, 3]*t(Zterms), covmat[3, 4]*t(Zterms),
                      covmat[4, 1]*t(Vterms), covmat[4, 2]*t(Vterms),
                      covmat[4, 3]*t(Vterms), covmat[4, 4]*t(Vterms)),
                     nrow=nrterms, byrow=TRUE)
    xtkronx <- matrix(c(covmat[1, 1]*t(Xterms)%*%Xterms, covmat[1, 2]*t(Xterms)%*%Yterms, 
                       covmat[1, 3]*t(Xterms)%*%Zterms, covmat[1, 4]*t(Xterms)%*%Vterms,
                       covmat[2, 1]*t(Yterms)%*%Xterms, covmat[2, 2]*t(Yterms)%*%Yterms, 
                       covmat[2, 3]*t(Yterms)%*%Zterms, covmat[2, 4]*t(Yterms)%*%Vterms,
                       covmat[3, 1]*t(Zterms)%*%Xterms, covmat[3, 2]*t(Zterms)%*%Yterms,
                       covmat[3, 3]*t(Zterms)%*%Zterms, covmat[3, 4]*t(Zterms)%*%Vterms,
                       covmat[4, 1]*t(Vterms)%*%Xterms, covmat[4, 2]*t(Vterms)%*%Yterms,
                       covmat[4, 3]*t(Vterms)%*%Zterms, covmat[4, 4]*t(Vterms)%*%Vterms),
                     nrow=nrterms, ncol=nrterms, byrow=TRUE)
    
    betapred <- ginv(xtkronx)%*%xtkron%*%allincmatr 
    print(betapred)
    #save(betapred, file = "Betapred.Rdata")    
  }
}
  