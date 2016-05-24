## FUP in R
# programmed by Niels Waller 
# October 14, 2014

###################################################################
FUP<-function(data, thetaInit, item, startvals, k = 0){
  
  
  # Prob of a keyed response 2PL ------------------------------##
  P <- function(m){
    1/(1+exp(-m))
  }
  
  ##-----------------------------------------------------------##
  fncFUP<- function(bParam, k=1, data, item, thetaInit){
    
       if(k>3) stop("\n\n*** FATAL ERROR:k must be <=3 ***\n\n")
    
       if(k==0){
      
         b0 = bParam[1]
         b1 = bParam[2]
      
         y = data[,item]
      
         m = b0 + b1*thetaInit 
      
         ##  avoid log(0)
         Pm <- P(m);
         Pm[Pm==1]   <-1-1e-12
         Pm[Pm<1e-12]<-1e-12
      
         # this is FHAT
         return(-mean( y*log(Pm) + (1 - y) * log(1 - Pm)  ))
      }#END k=0 
    
      if(k==1){
      
        b0 = bParam[1]
        b1 = bParam[2]
        b2 = bParam[3]
        b3 = bParam[4]
      
        y = data[,item]
      
        m = b0 + b1*thetaInit + b2*thetaInit^2 + b3*thetaInit^3
      
        ##  avoid log(0)
        Pm <- P(m);
      
        Pm[Pm=="NaN"]<-.5
        Pm[Pm==1]   <-1-1e-12
        Pm[Pm<1e-12]<-1e-12  
      
        # this is FHAT
        return(-mean( y*log(Pm) + (1 - y) * log(1 - Pm)  ))
      }#END k=1 
    
      if(k==2){
      
        b0 = bParam[1]
        b1 = bParam[2]
        b2 = bParam[3]
        b3 = bParam[4]
        b4 = bParam[5]
        b5 = bParam[6]
      
        y = data[,item]
      
        # create logit polynomial 
        m = b0 + b1*thetaInit + b2*thetaInit^2 + b3*thetaInit^3 + b4*thetaInit^4 + b5*thetaInit^5
      
        ##  avoid log(0)
        Pm <- P(m);
      
        Pm[Pm=="NaN"]<-.5
        Pm[Pm==1]   <-1-1e-12
        Pm[Pm<1e-12]<-1e-12  
      
        # this is FHAT
        return(-mean( y*log(Pm) + (1 - y) * log(1 - Pm)  ))
      }#END k=2 
    
    
      if(k==3){
      
        b0 = bParam[1]
        b1 = bParam[2]
        b2 = bParam[3]
        b3 = bParam[4]
        b4 = bParam[5]
        b5 = bParam[6]
        b6 = bParam[7]
        b7 = bParam[8]
      
        y = data[,item]
      
        # create logit polynomial 
        m = b0 + b1*thetaInit + b2*thetaInit^2 + b3*thetaInit^3 + 
            b4*thetaInit^4 + b5*thetaInit^5    +
            b6*thetaInit^6 + b7*thetaInit^7
      
        ##  avoid log(0)
        Pm <- P(m);
      
        Pm[Pm=="NaN"]<-.5
        Pm[Pm==1]   <-1-1e-12
        Pm[Pm<1e-12]<-1e-12  
      
        # this is FHAT
        return(-mean( y*log(Pm) + (1 - y) * log(1 - Pm)  ))
      }#END k=3 
    
  } #END fncFUP --------------------------------------------##
  
 
  
    out <- optim(par=startvals, 
                 fn=fncFUP, 
                 method="BFGS",
                 k=k,
                 data=data,
                 item=item,
                 thetaInit=thetaInit,
                 control=list(factr=1e-12,
                              maxit=1000))
    
  ## compute scaled (pseudo) AIC and BIC
     NSubj <- nrow(data)
     q <- 2 * k + 2
     AIC <- 2 * out$value + (2/NSubj) * q
     BIC <- 2 * out$value + (log(NSubj)/NSubj) * q
  
 
   list( b = out$par,
        FHAT=out$value,
        counts=out$counts,
        AIC = AIC,
        BIC = BIC,
        convergence = out$convergence)
  
} #END FMP




###################################################################
#                 END OF FUNCTION DEFINITIONS
###################################################################



