HiddenF <-
function(ymtx){
  a <- nrow(ymtx)
  b <- ncol(ymtx)
  if(a > 20){cat("Due to computation time, \n HiddenF is unavailable for a > 20 rows\n");stop}
  else{
  tall <- maketall.fcn(ymtx)
  cc <- 2^(a-1)-1
  ra <- t(t(ymtx - apply(ymtx,1,mean) + mean(ymtx))-apply(ymtx,2,mean))
  PSS <- (1/(1-1/a))*apply(ra^2,1,sum)
  ind <- which.max(PSS)
  combmax <- ind
  PSSmax <- PSS[ind]
  if(a > 3){
    if(a %% 2 ==1){
      niter = round((a-1)/2)
      for(j in 2:niter){
        combs <- combn(a,j)
        PSS <- (1/(j*(1-j/a)))*apply(combs, 2, computePSS, ra)
        ind <- which.max(PSS)
        if(PSSmax < PSS[ind]){
          combmax <- combs[,ind]
          PSSmax <- PSS[ind]
        }
      }
    }
    else{
      niter = round(a/2 - 1)
      if(niter > 1){
        for(j in 2:niter){
          combs <- combn(a,j)
          PSS <- (1/(j*(1-j/a)))*apply(combs, 2, computePSS, ra)
          ind <- which.max(PSS)
          if(PSSmax < PSS[ind]){
            combmax <- combs[,ind]
            PSSmax <- PSS[ind]
          }
        }
      }
      niter <- niter + 1
      combs <- combn(a,niter)
      ncomb <- round(ncol(combs)/2)
      PSS <- (1/(niter*(1-niter/a)))*apply(combs[,1:ncomb], 2, computePSS, ra)
      ind <- which.max(PSS)
      if(PSSmax < PSS[ind]){
        combmax <- combs[,ind]
        PSSmax <- PSS[ind]
      }
    }
  }
  varY = (sum(ra^2) - PSSmax)/((b-1)*(a-2))
  config.vec <- 1*is.element(tall$rows,combmax)
  y.tmpout <- lm(tall$y~config.vec*tall$cols + tall$rows/config.vec)
  pvalue <- anova(y.tmpout)$P[4]
  adjpvalue <- cc*pvalue
  adjpvalue <- min(1,adjpvalue)
  #return(list(maxcomb = combmax, maxPSS = PSSmax, varY = varY))
  hfout <- list(adjpvalue=adjpvalue,config.vector=config.vec,tall=tall,cc=cc)
  class(hfout) <- "HiddenF"
  return(hfout)
}
}
