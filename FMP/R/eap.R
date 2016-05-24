
eap<-function(data, bParams, NQuad = 21, priorVar =2, mintheta=-4, maxtheta=4){
  
  ##  Vers January 25, 2016
  ## eapMixed can be used to calculate eap theta estimates for item sets that fit 
  ##          monotonic polynomial IRT models of possibly different degree
  ##
  ## bParams an NItems by 9 matrix. Cols 1 : 8 include  coefficients of
  ##        the polynomial
  ##
  ## bParams col 9 should contain k (polynomial of degree 2k+1)

  if(is.data.frame(bParams)) bParams <- as.matrix(bParams)
  
  if(ncol(bParams)!=9) stop("\n\nbParams should have 9 columhns")
  
  breakPoints<-seq(mintheta,maxtheta,by=(maxtheta-mintheta)/NQuad)
  centers<-breakPoints+(maxtheta-mintheta)/(2*NQuad)
  #remove last center to create:
  #Quadrature points
  QuadPt<-centers[-length(centers)]
  
  
  NItem <- ncol(data)
  NSubj <- nrow(data)
  
  # inter QuadPt distance
  Quadinterval<-QuadPt[2] - QuadPt[1]
  # quadrature weights
  QuadWght <- dnorm(QuadPt, mean=0, sd=sqrt(priorVar))
  QuadWght <- QuadWght/sum(QuadWght)
  
  theta <- QuadPt
  thetaMat <- matrix(c(rep(1,NQuad), theta, theta^2, theta^3, theta^4, theta^5, theta^6, theta^7), nrow=NQuad, ncol=8)
  
  
  P <- function(m){
    1/(1+exp(-m))
  }
  
  eap <- vector(mode="numeric",length = NSubj)
  
  #initialize at 1
  EAPVec0 <- rep(1, length=NQuad)
  
  for(iSubj in 1:NSubj){
 
    EAPVec <- EAPVec0
   
    for(item in 1:NItem){
       mVec <- thetaMat %*% bParams[item,1:8]
       PVec <-P(mVec)   
       EAPVec <- EAPVec * ( PVec^data[iSubj,item]  *  (1-PVec)^(1 - data[iSubj,item ]))
    }
    eap[iSubj]<-sum(EAPVec *QuadWght * QuadPt)/sum(EAPVec * QuadWght)
  }
  
  #return eap estimates
  eap
  
}

