########## R function: SignifFeatureRegion ##########

# For determining the region of significant
# gradient for a particular bandwidth and
# significance level.

# Last changed: 18 JAN 2006

SignifFeatureRegion <- function(n, d, gcounts, gridsize, dest, bandwidth, signifLevel, range.x, grad=TRUE, curv=TRUE, neg.curv.only=TRUE)
{
  h <- bandwidth
  
  ESS <- n*dest$est*prod(h)*(sqrt(2*pi)^d)
  SigESS <- ESS >= 5

  Sig.scalar <- array(NA, dim=gridsize)
  Sig2.scalar <- array(NA, dim=gridsize)

  dest$est[dest$est<0] <- 0  
  ## constant for variance of gradient estimate
  Sig.scalar <- 1/2*(2*sqrt(pi))^(-d)*n^(-1)*prod(h)^(-1)*dest$est

  ##  constants for variance of curvature estimate  
  if (d==1)
    Sig2.scalar <- (8*sqrt(pi)*n*prod(h))^(-1)*dest$est
  else if (d==2)
    Sig2.scalar <- (16*pi*n*prod(h))^(-1)*dest$est
  else if (d==3)
    Sig2.scalar <- (32*pi^(3/2)*n*prod(h))^(-1)*dest$est
  else if (d==4)
    Sig2.scalar <- (64*pi^2*n*prod(h))^(-1)*dest$est

  ## Matrix square root - taken from Stephen Lake
  ## http://www5.biostat.wustl.edu/s-news/s-news-archive/200109/msg00067.html
  
  matrix.sqrt <- function(A)
  {
    sva <- svd(A)
    if (min(sva$d)>=0)
      Asqrt <- sva$u %*% diag(sqrt(sva$d)) %*% t(sva$v)
    else
      stop("Matrix square root is not defined")
    return(Asqrt)
  }

  
  if (d>1)
  {
    WaldGrad <- array(NA, dim=gridsize)
    WaldCurv <- array(NA, dim=gridsize)
    local.mode <- array(FALSE, dim=gridsize)
  }
  
  if (d==1)
  {
    if (grad)
    {
      obj1 <- drvkde(gcounts, drv=1, bandwidth=h, binned=TRUE, range.x=range.x, se=FALSE)
      fhat1 <- obj1$est

      Sig.inv12 <- 1/sqrt(Sig.scalar * h^(-2))
      WaldGrad <- (Sig.inv12 * fhat1)^2
    }

    if (curv)
    {
      obj2 <- drvkde(gcounts,drv=2,bandwidth=h,binned=TRUE,range.x=range.x, se=FALSE)
      fhat2 <- obj2$est
      
      Sig2.inv12 <- 1/sqrt(Sig2.scalar * 3*h^(-4))
      lambda1 <- Sig2.inv12 * fhat2
      WaldCurv <- lambda1^2
      local.mode <- (lambda1 < 0)
    }
  }
  
  if (d==2)
  {
    if (grad)
    { 
      obj10 <- drvkde(gcounts,drv=c(1,0),bandwidth=h,binned=TRUE,
                      range.x=range.x,se=FALSE)
      obj01 <- drvkde(gcounts,drv=c(0,1),bandwidth=h,binned=TRUE,
                      range.x=range.x,se=FALSE)
      fhat10 <- obj10$est
      fhat01 <- obj01$est

      for (i1 in 1:gridsize[1])
        for (i2 in 1:gridsize[2])
          if (SigESS[i1,i2])
          {  
            Sig.inv12 <- 1/sqrt(Sig.scalar[i1,i2] * h^(-2))
            WaldGrad[i1,i2] <- sum((Sig.inv12 * c(fhat10[i1,i2], fhat01[i1,i2]))^2) 
          }
    }


    if (curv)
    {         
      Sig2.mat <-
        matrix(c(3/h[1]^4, 0, 1/(h[1]^2*h[2]^2),
                 0, 1/(h[1]^2*h[2]^2), 0,
                 1/(h[1]^2*h[2]^2), 0, 3/h[2]^4),
               nrow=3, ncol=3)
      
      Sig2.mat.inv <- chol2inv(chol(Sig2.mat))

      obj20 <- drvkde(gcounts,drv=c(2,0),bandwidth=h,
                      binned=TRUE,range.x=range.x, se=FALSE)
      obj11 <- drvkde(gcounts,drv=c(1,1),bandwidth=h,
                      binned=TRUE,range.x=range.x, se=FALSE)
      obj02 <- drvkde(gcounts,drv=c(0,2),bandwidth=h,
                      binned=TRUE,range.x=range.x, se=FALSE)
      fhat20 <- obj20$est
      fhat11 <- obj11$est
      fhat02 <- obj02$est

      for (i1 in 1:gridsize[1])
        for (i2 in 1:gridsize[2])
          if (SigESS[i1,i2])
          {  
            Sig2.inv12 <- sqrt(1/Sig2.scalar[i1,i2])*matrix.sqrt(Sig2.mat.inv) 
            fhat.temp <- Sig2.inv12 %*%
              c(fhat20[i1,i2], fhat11[i1,i2], fhat02[i1,i2])

            WaldCurv[i1,i2] <- sum(fhat.temp^2)
          }
      
      lambda1 <- ((fhat20 + fhat02) - sqrt((fhat20-fhat02)^2 + 4*fhat11^2))/2
      lambda2 <- ((fhat20 + fhat02) + sqrt((fhat20-fhat02)^2 + 4*fhat11^2))/2
      local.mode <- (lambda1 < 0) & (lambda2 < 0)
    }
  }
    
 
  if (d==3)
  {
    if (grad)
    {    
      obj100 <- drvkde(gcounts,drv=c(1,0,0),bandwidth=h,
                       binned=TRUE,range.x=range.x, se=FALSE)
      obj010 <- drvkde(gcounts,drv=c(0,1,0),bandwidth=h,
                       binned=TRUE,range.x=range.x, se=FALSE)
      obj001 <- drvkde(gcounts,drv=c(0,0,1),bandwidth=h,
                       binned=TRUE,range.x=range.x, se=FALSE)
      
      fhat100 <- obj100$est
      fhat010 <- obj010$est
      fhat001 <- obj001$est
    
      for (i1 in 1:gridsize[1])
        for (i2 in 1:gridsize[2])
          for (i3 in 1:gridsize[3])
            if (SigESS[i1,i2,i3])
            {    
              Sig.inv12 <- 1/sqrt(Sig.scalar[i1,i2,i3] * h^(-2))
              WaldGrad[i1,i2,i3] <-
                sum((Sig.inv12 * c(fhat100[i1,i2,i3], fhat010[i1,i2,i3],
                                   fhat001[i1,i2,i3]))^2)
            }
    }

  
    if (curv)
    {
      obj200 <- drvkde(gcounts,drv=c(2,0,0),bandwidth=h,
                       binned=TRUE,range.x=range.x, se=FALSE)
      obj110 <- drvkde(gcounts,drv=c(1,1,0),bandwidth=h,
                       binned=TRUE,range.x=range.x, se=FALSE)
      obj101 <- drvkde(gcounts,drv=c(1,0,1),bandwidth=h,
                       binned=TRUE,range.x=range.x, se=FALSE)
      obj020 <- drvkde(gcounts,drv=c(0,2,0),bandwidth=h,
                       binned=TRUE,range.x=range.x, se=FALSE)
      obj011 <- drvkde(gcounts,drv=c(0,1,1),bandwidth=h,
                       binned=TRUE,range.x=range.x, se=FALSE)
      obj002 <- drvkde(gcounts,drv=c(0,0,2),bandwidth=h,
                       binned=TRUE,range.x=range.x, se=FALSE)
      fhat200 <- obj200$est
      fhat110 <- obj110$est
      fhat101 <- obj101$est
      fhat020 <- obj020$est
      fhat011 <- obj011$est
      fhat002 <- obj002$est
    
      Sig2.mat <-
        matrix(c(3/h[1]^4, 0, 0, 1/(h[1]*h[2])^2, 0, 1/(h[1]*h[3])^2,
                 0, 1/(h[1]*h[2])^2, 0, 0, 0, 0,
                 0, 0, 1/(h[1]*h[3])^2, 0, 0, 0,
                 1/(h[1]*h[2])^2, 0, 0, 3/h[2]^4, 0, 1/(h[2]*h[3])^2,
                 0, 0, 0, 0, 1/(h[2]*h[3])^2, 0,
                 1/(h[1]*h[3])^2, 0, 0, 1/(h[2]*h[3])^2, 0, 3/h[3]^4),
               nrow=6, ncol=6)
      
      Sig2.mat.inv <- chol2inv(chol(Sig2.mat))  

      ## at each grid point, find eigenvalues of vech'ed curvature matrix 
      for (i1 in 1:gridsize[1])
        for (i2 in 1:gridsize[2])
          for (i3 in 1:gridsize[3])
            if (SigESS[i1,i2,i3])
            {
              Sig2.inv12 <- sqrt(1/Sig2.scalar[i1,i2,i3]) *
                matrix.sqrt(Sig2.mat.inv) 
              fhat.temp <- Sig2.inv12 %*%
                c(fhat200[i1,i2,i3], fhat110[i1,i2,i3], fhat101[i1,i2,i3],
                  fhat020[i1,i2,i3], fhat011[i1,i2,i3], fhat002[i1,i2,i3])
              
              D2.mat <-
                matrix(c(fhat200[i1,i2,i3], fhat110[i1,i2,i3], fhat101[i1,i2,i3],
                         fhat110[i1,i2,i3], fhat020[i1,i2,i3], fhat011[i1,i2,i3],
                         fhat101[i1,i2,i3], fhat011[i1,i2,i3], fhat002[i1,i2,i3]),
                       nrow=3)
              lambda <- eigen(D2.mat, symmetric=TRUE, only.values=TRUE)$values
              
              WaldCurv[i1,i2,i3] <- sum(fhat.temp^2)
              
              local.mode[i1,i2,i3] <- all(lambda < 0)
          }
    }
    
  }
  
  if (d==4)
  {
    if (grad)
    {
      obj1000 <- drvkde(gcounts,drv=c(1,0,0,0),bandwidth=h,
                        binned=TRUE,range.x=range.x, se=FALSE)
      obj0100 <- drvkde(gcounts,drv=c(0,1,0,0),bandwidth=h,
                        binned=TRUE,range.x=range.x, se=FALSE)
      obj0010 <- drvkde(gcounts,drv=c(0,0,1,0),bandwidth=h,
                        binned=TRUE,range.x=range.x, se=FALSE)
      obj0001 <- drvkde(gcounts,drv=c(0,0,0,1),bandwidth=h,
                        binned=TRUE,range.x=range.x, se=FALSE)
      
      fhat1000 <- obj1000$est
      fhat0100 <- obj0100$est
      fhat0010 <- obj0010$est
      fhat0001 <- obj0001$est

      for (i1 in 1:gridsize[1])
        for (i2 in 1:gridsize[2])
          for (i3 in 1:gridsize[3])
            for (i4 in 1:gridsize[4])
              if (SigESS[i1,i2,i3,i4])
              {
                Sig.inv12 <- 1/sqrt(Sig.scalar[i1,i2,i3,i4] * h^(-2))
                WaldGrad[i1,i2,i3,i4] <-
                  sum((Sig.inv12*c(fhat1000[i1,i2,i3,i4],fhat0100[i1,i2,i3,i4],
                                   fhat0010[i1,i2,i3,i4],fhat0001[i1,i2,i3,i4]))^2)
              }
    }
    
    if (curv)
    {
      obj2000 <- drvkde(gcounts,drv=c(2,0,0,0),bandwidth=h,
                        binned=TRUE,range.x=range.x, se=FALSE)
      obj1100 <- drvkde(gcounts,drv=c(1,1,0,0),bandwidth=h,
                        binned=TRUE,range.x=range.x, se=FALSE)
      obj1010 <- drvkde(gcounts,drv=c(1,0,1,0),bandwidth=h,
                        binned=TRUE,range.x=range.x, se=FALSE)
      obj1001 <- drvkde(gcounts,drv=c(1,0,0,1),bandwidth=h,
                        binned=TRUE,range.x=range.x, se=FALSE)
      obj0200 <- drvkde(gcounts,drv=c(0,2,0,0),bandwidth=h,
                        binned=TRUE,range.x=range.x, se=FALSE)
      obj0110 <- drvkde(gcounts,drv=c(0,1,1,0),bandwidth=h,
                        binned=TRUE,range.x=range.x, se=FALSE)
      obj0101 <- drvkde(gcounts,drv=c(0,1,0,1),bandwidth=h,
                        binned=TRUE,range.x=range.x, se=FALSE)
      obj0020 <- drvkde(gcounts,drv=c(0,0,2,0),bandwidth=h,
                        binned=TRUE,range.x=range.x, se=FALSE)
      obj0011 <- drvkde(gcounts,drv=c(0,0,1,1),bandwidth=h,
                        binned=TRUE,range.x=range.x, se=FALSE)
      obj0002 <- drvkde(gcounts,drv=c(0,0,0,2),bandwidth=h,
                        binned=TRUE,range.x=range.x, se=FALSE)

      fhat2000 <- obj2000$est
      fhat1100 <- obj1100$est
      fhat1010 <- obj1010$est
      fhat1001 <- obj1001$est
      fhat0200 <- obj0200$est
      fhat0110 <- obj0110$est
      fhat0101 <- obj0101$est
      fhat0020 <- obj0020$est
      fhat0011 <- obj0011$est
      fhat0002 <- obj0002$est

      Sig2.mat <-
        matrix(c(3/h[1]^4, 0, 0, 0, 1/(h[1]*h[2])^2, 0, 0, 1/(h[1]*h[3])^2, 0, 1/(h[1]*h[4])^2,
                 0, 1/(h[1]*h[2])^2, 0, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 1/(h[1]*h[3])^2, 0, 0, 0, 0, 0, 0, 0,
                 0, 0, 0, 1/(h[1]*h[4])^2, 0, 0, 0, 0, 0, 0,
                 1/(h[1]*h[2])^2, 0, 0, 0, 3/h[2]^4, 0, 0, 1/(h[2]*h[3])^2, 0, 1/(h[2]*h[4])^2,
                 0, 0, 0, 0, 0, 1/(h[2]*h[3])^2, 0, 0, 0, 0,
                 0, 0, 0, 0 ,0, 0, 1/(h[2]*h[4])^2, 0, 0, 0,
                 1/(h[1]*h[3])^2, 0, 0, 0, 1/(h[2]*h[3])^2, 0, 0, 3/h[3]^4, 0, 1/(h[3]*h[4])^2,
                 0, 0, 0, 0, 0, 0, 0, 0, 1/(h[3]*h[4])^2, 0,
                 1/(h[1]*h[4])^2, 0, 0, 0, 1/(h[2]*h[4])^2, 0, 0, 1/(h[3]*h[4])^2, 0, 3/h[4]^4),
               nrow=10, ncol=10)
      
      Sig2.mat.inv <- chol2inv(chol(Sig2.mat))
     
      for (i1 in 1:gridsize[1])
        for (i2 in 1:gridsize[2])
          for (i3 in 1:gridsize[3])
            for (i4 in 1:gridsize[4])        
              if (SigESS[i1,i2,i3,i4])
              {
                Sig2.inv12 <- sqrt(1/Sig2.scalar[i1,i2,i3,i4]) *
                  matrix.sqrt(Sig2.mat.inv) 
                fhat.temp <- Sig2.inv12 %*%
                  c(fhat2000[i1,i2,i3,i4], fhat1100[i1,i2,i3,i4],
                    fhat1010[i1,i2,i3,i4], fhat1001[i1,i2,i3,i4],
                    fhat0200[i1,i2,i3,i4], fhat0110[i1,i2,i3,i4],
                    fhat0101[i1,i2,i3,i4], fhat0020[i1,i2,i3,i4],
                    fhat0011[i1,i2,i3,i4], fhat0002[i1,i2,i3,i4])
                
                D2.mat <-
                  matrix(c(fhat2000[i1,i2,i3,i4], fhat1100[i1,i2,i3,i4], fhat1010[i1,i2,i3,i4], fhat1001[i1,i2,i3,i4],
                           fhat1100[i1,i2,i3,i4], fhat0200[i1,i2,i3,i4], fhat0110[i1,i2,i3,i4], fhat0101[i1,i2,i3,i4],
                           fhat1010[i1,i2,i3,i4], fhat0110[i1,i2,i3,i4], fhat0020[i1,i2,i3,i4], fhat0011[i1,i2,i3,i4],
                           fhat1001[i1,i2,i3,i4], fhat0101[i1,i2,i3,i4], fhat0011[i1,i2,i3,i4], fhat0002[i1,i2,i3,i4]),
                         nrow=4)
                WaldCurv[i1,i2,i3,i4] <- sum(fhat.temp^2)
                lambda <- eigen(D2.mat, symmetric=TRUE, only.values=TRUE)$values
                
                local.mode[i1,i2,i3,i4] <- all(lambda < 0)
              }     
    }
  }


  ## multiple hypothesis testing - based on Hochberg's method
  ## - modified Bonferroni method using ordered p-values
  
  ## test statistic for gradient
  if (grad)
  {
    pval.Grad <- 1 - pchisq(WaldGrad, d)
    pval.Grad.ord <- pval.Grad[order(pval.Grad)]
    num.test <- sum(!is.na(pval.Grad.ord))

    if (num.test>=1)
      num.test.seq <- c(1:num.test, rep(NA, prod(gridsize) - num.test))
    else
      num.test.seq <- rep(NA, prod(gridsize))
    
    reject.nonzero <- ((pval.Grad.ord <= signifLevel/(num.test + 1 - num.test.seq)) &
                       (pval.Grad.ord > 0))  
    reject.nonzero.ind <- which(reject.nonzero)

    ## p-value == 0 => reject null hypotheses automatically
    SignifGrad <- array(FALSE, dim=gridsize)
    SignifGrad[which(pval.Grad==0, arr.ind=TRUE)] <- TRUE
    
    ## p-value > 0 then reject null hypotheses indicated in reject.nonzero.ind 
    for (i in reject.nonzero.ind)
      SignifGrad[which(pval.Grad==pval.Grad.ord[i], arr.ind=TRUE)] <- TRUE 
  }

  ## test statistic for curvature
  if (curv)
  {
    pval.Curv <- 1 - pchisq(WaldCurv, d*(d+1)/2)
    pval.Curv.ord <- pval.Curv[order(pval.Curv)]
    num.test <- sum(!is.na(pval.Curv.ord))

    if (num.test>=1)
      num.test.seq <- c(1:num.test, rep(NA, prod(gridsize) - num.test))
    else
      num.test.seq <- rep(NA, prod(gridsize))
    reject.nonzero <- ((pval.Curv.ord <= signifLevel/(num.test + 1 - num.test.seq)) &(pval.Curv.ord > 0))  
    reject.nonzero.ind <- which(reject.nonzero)

    SignifCurv <- array(FALSE, dim=gridsize)

    ## p-value == 0 => reject null hypotheses automatically
    SignifCurv[which(pval.Curv==0, arr.ind=TRUE)] <- TRUE

    ## p-value > 0 then reject null hypotheses indicated in reject.nonzero.ind
    for (i in reject.nonzero.ind)
      SignifCurv[which(pval.Curv==pval.Curv.ord[i], arr.ind=TRUE)] <- TRUE 

    if (neg.curv.only) SignifCurv <- SignifCurv & local.mode
  }
  
  if (grad & !curv)
    return(list(grad=SignifGrad))
  else if (!grad & curv)
    return(list(curv=SignifCurv))
  else if (grad & curv)
    return(list(grad=SignifGrad, curv=SignifCurv))
}


########## End of SignifFeatureRegion ##########


