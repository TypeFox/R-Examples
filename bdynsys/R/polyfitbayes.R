# polyfitbayes corresponds to polyfitreg, the fitting is however not based
# on regression and maximum log likelihood but on Bayes Factor

polyfitbayes <- function(indnr, xv, yv, ch, sel, zv, vv)
{
  if (indnr == 2)
  {
    # Max number of iterations 
    maxiter = 200000; 
    
    # computing the inverse of x and y
    invx <- xv^(-1)
    invy <- yv^(-1)
    
    # handling "Inf" as a result of inverse computation (/0), is this a problem?
    idx1 = which(is.infinite(invx))
    idx2 = which(is.infinite(invy)) 
    invx[idx1] <- NA
    invy[idx2] <- NA
    
    # defining in- and output for the regression, computing the Ordinary Least Square Regression
    input <- cbind(rep(1, length(xv)), invx, invy, xv, yv, invx*invy, xv*invy, yv*invx,
                   xv*yv, xv^2, invx^2, yv^2, invy^2, xv^3, yv^3, invx^3, invy^3)
    inputs <- input[, sel]
    output <- ch
    
    # Monte Carlo Simulation, integration 
    paramspc <- seq(-0.1, 0.1, 0.001)
    randobsnum <- min(maxiter, (200^length(sel))) 
    indexchoose <- matrix(runif(randobsnum*length(sel)), randobsnum, length(sel))
    Cparam <- -0.1 + 200*0.001*indexchoose
    
    # computing log BayesFactor
    Plikeli <- c()
    logBF <- c()
    if (length(sel)==1)
    {
      for (tt in 1:nrow(Cparam))
      { 
        Plikeli[tt] <- exp(sum(-(output - inputs*(Cparam[tt,])^2/2/var(output))))
        logBF[tt] <- sum(-(output - inputs*(Cparam[tt,]))^2/(2*var(output))) 
      }
    } 
    if (length(sel)>1)
    {
      for (tt in 1:nrow(Cparam))
      {
        Plikeli[tt] <- exp(sum(-(output - inputs%*%(Cparam[tt,])^2/2/var(output))))
        logBF[tt] <- sum(-(output - inputs%*%(Cparam[tt,]))^2/(2*var(output))) 
      }
    }  
    
    # selecting best model based on BayesFactor
    bestm <- max(logBF)
    indexbest <- which.max(logBF)
    B <- Cparam[indexbest,]
    bestm <- max(logBF) + log(sum(exp(logBF - max(logBF)))/randobsnum)
    
    return(list(bestm, indexbest))
  }
  
  ############################# indnr == 2 ends, indnr == 3 begins #################################
  
  if (indnr == 3)
  {
    # Max number of iterations 
    maxiter = 200000; 
    
    # computing the inverse of x and y
    invx <- xv^(-1)
    invy <- yv^(-1)
    invz <- zv^(-1)
    
    # handling "Inf" as a result of inverse computation (/0), is this a problem?
    idx1 = which(is.infinite(invx))
    idx2 = which(is.infinite(invy)) 
    idx3 = which(is.infinite(invz)) 
    invx[idx1] <- NA
    invy[idx2] <- NA
    invz[idx3] <- NA
    
    # defining in- and output for the regression, computing the Ordinary Least Square Regression
    input <- cbind(rep(1, length(xv)), invx, invy, invz, xv, yv, zv, invx*invy, invy*invz, invx*invz, 
                   xv*yv, yv*zv, xv*zv, xv*invy, yv*invx, xv*invz, zv*invx, yv*invz, zv*invy, xv*invy*invz,
                   yv*invx*invz, zv*invx*invy, xv*yv*invz, yv*zv*invx, xv*zv*invy, xv*yv*zv, invx*invy*invz, 
                   xv^2, invx^2, yv^2, invy^2, zv^2, invz^2, xv^3, yv^3, zv^3, invx^3, invy^3, invz^3)
    inputs <- input[, sel]
    output <- ch
    
    # Monte Carlo Simulation, integration 
    paramspc <- seq(-0.1, 0.1, 0.001)
    randobsnum <- min(maxiter, (200^length(sel))) 
    indexchoose <- matrix(runif(randobsnum*length(sel)), randobsnum, length(sel))
    Cparam <- -0.1 + 200*0.001*indexchoose
    
    # computing log BayesFactor
    Plikeli <- c()
    logBF <- c()
    if (length(sel)==1)
    {
      for (tt in 1:nrow(Cparam))
      {
        Plikeli[tt] <- exp(sum(-(output - inputs*(Cparam[tt,])^2/2/var(output))))
        logBF[tt] <- sum(-(output - inputs*(Cparam[tt,]))^2/(2*var(output))) 
      }
    } 
    if (length(sel)>1)
    {
      for (tt in 1:nrow(Cparam))
      {
        Plikeli[tt] <- exp(sum(-(output - inputs%*%(Cparam[tt,])^2/2/var(output))))
        logBF[tt] <- sum(-(output - inputs%*%(Cparam[tt,]))^2/(2*var(output))) 
      }
    } 
    
    # selecting best model based on BayesFactor
    bestm <- max(logBF)
    indexbest <- which.max(logBF)
    B <- Cparam[indexbest,]
    bestm <- max(logBF) + log(sum(exp(logBF - max(logBF)))/randobsnum)
    
    return(list(bestm, indexbest))
  }
  
  ############################# indnr == 3 ends, indnr == 4 begins #################################
  
  if (indnr == 4)
  {
    # Max number of iterations 
    maxiter = 200000; 
    
    # computing the inverse of x and y
    invx <- xv^(-1)
    invy <- yv^(-1)
    invz <- zv^(-1)
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
    
    # defining in- and output for the regression
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
    inputs <- input[, sel]
    output <- ch
    
    # Monte Carlo Simulation, integration 
    paramspc <- seq(-0.1, 0.1, 0.001)
    randobsnum <- min(maxiter, (200^length(sel))) 
    indexchoose <- matrix(runif(randobsnum*length(sel)), randobsnum, length(sel))
    Cparam <- -0.1 + 200*0.001*indexchoose
    
    # computing log BayesFactor
    Plikeli <- c()
    logBF <- c()
    if (length(sel)==1)
    {
      for (tt in 1:nrow(Cparam))
      { 
        Plikeli[tt] <- exp(sum(-(output - inputs*(Cparam[tt,])^2/2/var(output))))
        logBF[tt] <- sum(-(output - inputs*(Cparam[tt,]))^2/(2*var(output))) 
      }
    } 
    if (length(sel)>1)
    {
      for (tt in 1:nrow(Cparam))
      {
        Plikeli[tt] <- exp(sum(-(output - inputs%*%(Cparam[tt,])^2/2/var(output))))
        logBF[tt] <- sum(-(output - inputs%*%(Cparam[tt,]))^2/(2*var(output))) 
      }
    } 
    
    # selecting best model based on BayesFactor
    bestm <- max(logBF)
    indexbest <- which.max(logBF)
    B <- Cparam[indexbest,]
    bestm <- max(logBF) + log(sum(exp(logBF - max(logBF)))/randobsnum)
    
    return(list(bestm, indexbest))
  }
}