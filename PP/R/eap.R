eap_4pl <- function(respm, thres, slopes, lowerA=NULL, upperA=NULL,
                    mu = NULL, sigma2 = NULL)
{
  
  nodmat <- matrix(0,nrow(respm),60)  
  weimat <- matrix(0,nrow(respm),60) 
  # #persons x #nodes
  
  for(qrun in 1:nrow(respm))
  {
    qres <- quadIT(nodes = 60,mu = mu[qrun],sigma = sqrt(sigma2[qrun]),absrange = 6)
    nodmat[qrun,]  <- qres[[1]]$nodes
    weimat[qrun,]  <- qres[[1]]$weights
  }
  
  
  FPL_eap <- sapply(1:nrow(respm),function(ep)
    {
      
      oben <- exp(t(outer(nodmat[ep,],thres[-1,],"-"))*slopes)
      oben2 <- oben * (upperA-lowerA)
      P <- lowerA + oben2/(1+oben)
      # item x nodes
      Q <- 1 - P
      
      xi <- respm[ep,]
      ximi <- 1-xi
      pxi <- P * xi
      qxi <- Q * ximi
      li  <- pxi + qxi
      
      Li <- apply(li,2,function(x)prod(x,na.rm=TRUE)) 
      EAP <- sum(Li * weimat[ep,] * nodmat[ep,])/sum((Li*weimat[ep,]))
      
      EAP_v <- sum(Li * weimat[ep,] * (nodmat[ep,]-EAP)^2)/sum((Li*weimat[ep,]))
      
      c(EAP,sqrt(EAP_v))
    })
  
  t(FPL_eap)   
}





eap_gpcm <- function(respm, thres, slopes, mu = NULL, sigma2 = NULL)
  {
    
    nodmat <- matrix(0,nrow(respm),60)  
    weimat <- matrix(0,nrow(respm),60) 
    # #persons x #nodes
    
    for(qrun in 1:nrow(respm))
    {
      qres <- quadIT(nodes = 60,mu = mu[qrun],sigma = sqrt(sigma2[qrun]),absrange = 6)
      nodmat[qrun,]  <- qres[[1]]$nodes
      weimat[qrun,]  <- qres[[1]]$weights
    }
    
    
  FPL_eap <- sapply(1:nrow(respm),function(ep)
      {
        
        Li  <- Likgpcm(respm[ep,],thres,slopes,nodmat[ep,])
  
        EAP <- sum(Li * weimat[ep,] * nodmat[ep,])/sum((Li*weimat[ep,]))
        
        EAP_v <- sum(Li * weimat[ep,] * (nodmat[ep,]-EAP)^2)/sum((Li*weimat[ep,]))
        
        c(EAP,sqrt(EAP_v))
      })
    
  t(FPL_eap)  
  
  }



## ------------------ mixed 4pl und gpcm

eap_mixed <- function(respm, thres, slopes, lowerA, upperA, mu = NULL, sigma2 = NULL,model2est)
  {
    
    nodmat <- matrix(0,nrow(respm),60)  
    weimat <- matrix(0,nrow(respm),60) 
    # #persons x #nodes
    
    wheregpcm <- model2est == "GPCM"
    where4pl  <- model2est == "4PL"
    
    respm_gpcm <- respm[,wheregpcm, drop=FALSE]
    respm_4pl  <- respm[,where4pl , drop=FALSE]
    
    
    for(qrun in 1:nrow(respm))
    {
      qres <- quadIT(nodes = 60,mu = mu[qrun],sigma = sqrt(sigma2[qrun]),absrange = 6)
      nodmat[qrun,]  <- qres[[1]]$nodes
      weimat[qrun,]  <- qres[[1]]$weights
    }
    
      # 4pl part
      thres4pl  <- thres[,where4pl,drop=FALSE]
      slopes4pl  <- slopes[where4pl]
      lowerA4pl <- lowerA[where4pl]
      upperA4pl <- upperA[where4pl]
      # gpcm part
      thresgpcm <- thres[,wheregpcm,drop=FALSE]
      slopegpcm <- slopes[wheregpcm]
      
    
    FPL_eap <- sapply(1:nrow(respm),function(ep)
    {
      
      oben <- exp(t(outer(nodmat[ep,],thres4pl[2,],"-"))*slopes4pl)
      oben2 <- oben * (upperA4pl - lowerA4pl)
      P <- lowerA4pl + oben2/(1+oben)
      # item x nodes
      Q <- 1 - P
      
      xi <- respm_4pl[ep,]
      ximi <- 1-xi
      pxi <- P * xi
      qxi <- Q * ximi
      li  <- pxi + qxi
      
      Li4pl <- apply(li,2,function(x)prod(x,na.rm=TRUE))
      
      Ligpcm    <- Likgpcm(respm_gpcm[ep,],thresgpcm,slopegpcm,nodmat[ep,])
      
      # all together
      Li <- Li4pl * Ligpcm
      
      EAP <- sum(Li * weimat[ep,] * nodmat[ep,])/sum((Li*weimat[ep,]))
      
      EAP_v <- sum(Li * weimat[ep,] * (nodmat[ep,]-EAP)^2)/sum((Li*weimat[ep,]))
      
      c(EAP,sqrt(EAP_v))
    })
    
    t(FPL_eap)  
    
  }









