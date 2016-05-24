## FMP 
# programmed by Niels Waller 
# January 27, 2016
# bugs in T2 T3 fixed on February 8, 2016



#########################################################################
FMP<-function(data, thetaInit, item, startvals, k=0, eps=1e-6){
  
	if(k > 3) stop("\nFatal Error: k must equal 0, 1, 2, or 3.\n")
  
    T1 <- matrix(0,3,1)
    T2 <- matrix(0,5,3); 
    T2[1,1]<-T2[2,2]<-T2[3,3]<-1
    T3 <- matrix(0,7,5)
    T3[1,1]<-T3[2,2]<-T3[3,3]<-T3[4,4]<-T3[5,5]<-1
   
  
    # Prob of a keyed response 2PL
    P <- function(m){
        1/(1+exp(-m))
    }  
  
    #------------------------------------------------------------------##
    # gammaTob:  generate final b coefficients from final gamma estimates
    gammaTob<-function(gammaParam, k=1){
    
           if(k>3) stop("\n\n*** FATAL ERROR:k must be <=3 ***\n\n")
    
           if(k==0){
             Ksi = gammaParam[1]
             lambda = gammaParam[2]
             b0 = Ksi
             b1 = lambda
             return(c(b0,b1))
           } 
    
           if(k==1){
             Ksi = gammaParam[1]
             lambda = gammaParam[2]
             alpha1<-gammaParam[3]
             beta1 <-gammaParam[4]
             phi1<-alpha1^2 + beta1
             T1 <-as.vector(c(1, -2*alpha1, phi1))
             asup0 = lambda
             asup1 = T1 * asup0
             b0 = Ksi
             b1 = asup1[1]/1
             b2 = asup1[2]/2
             b3 = asup1[3]/3
            return( c(b0,b1,b2,b3))
          }  
    
          if(k==2){
            Ksi = gammaParam[1]
            lambda = gammaParam[2]
            alpha1 <- gammaParam[3]
            beta1  <- gammaParam[4]
            alpha2 <- gammaParam[5]
            beta2  <- gammaParam[6]
            phi1   <- alpha1^2 + beta1
            phi2   <- alpha2^2 + beta2
            T1 <-as.vector(c(1, -2*alpha1, phi1))
            T2[2,1]<-T2[3,2]<-T2[4,3] <- -2*alpha2
            T2[3,1] <- T2[4,2] <- T2[5,3]<- phi2
            # see (A6) Liang & Browne
            asup1 = (T2 %*% T1) * lambda
            b0 = Ksi
            b1 = asup1[1]/1
            b2 = asup1[2]/2
            b3 = asup1[3]/3
            b4 = asup1[4]/3
            b5 = asup1[5]/3
            return( c(b0,b1,b2,b3,b4,b5))
          }  
    
          if(k==3){
            Ksi = gammaParam[1]
            lambda = gammaParam[2]
      
            alpha1 <- gammaParam[3]
            beta1  <- gammaParam[4]
            alpha2 <- gammaParam[5]
            beta2  <- gammaParam[6]
            alpha3 <- gammaParam[7]
            beta3  <- gammaParam[8]
      
            phi1 <- alpha1^2 + beta1
            phi2 <- alpha2^2 + beta2
            phi3 <- alpha3^2 + beta3
      
            T1 <-as.vector(c(1, -2*alpha1, phi1))
            T2[2,1]<-T2[3,2]<-T2[4,3] <- -2*alpha2
            T2[3,1] <- T2[4,2] <- T2[5,3]<- phi2
            T3[2,1]<-T3[3,2]<-T3[4,3]<-T3[5,4]<-T3[6,5]<- -2*alpha3
            T3[3,1]<-T3[4,2]<-T3[5,3]<-T3[6,4]<-T3[7,5]<- phi3
      
            # see (A6) Liang & Browne
            asup3 = (T3 %*% T2 %*% T1) * lambda
      
            b0 = Ksi
            b1 = asup3[1]/1
            b2 = asup3[2]/2
            b3 = asup3[3]/3
            b4 = asup3[4]/4
            b5 = asup3[5]/5
            b6 = asup3[6]/6
            b7 = asup3[7]/7
            return( c(b0,b1,b2,b3,b4,b5,b6,b7))
         }  
  } # END gammaTob -----------------------------------------------##
  
  
  
  
  
  ###################################################################
  ## fncFMP
    fncFMP<- function(gammaParam, k=1, data, item, thetaInit){
    
          if(k>3) stop("\n\n*** FATAL ERROR:k must be <=3 ***\n\n")
    
          if(k==0){
                Ksi = gammaParam[1]
                lambda = gammaParam[2]
         
                b0 = Ksi
                b1 = lambda
      
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
                Ksi = gammaParam[1]
                lambda = gammaParam[2]
      
                alpha1<-gammaParam[3]
                beta1 <-gammaParam[4]
      
                phi1<-alpha1^2 + beta1
      
               T1 <-as.vector(c(1, -2*alpha1, phi1))
      
               asup0 = lambda
               asup1 = T1 * asup0
               b0 = Ksi
               b1 = asup1[1]/1
               b2 = asup1[2]/2
               b3 = asup1[3]/3
      
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
                Ksi = gammaParam[1]
                lambda = gammaParam[2]
      
                alpha1 <- gammaParam[3]
                beta1  <- gammaParam[4]
                alpha2 <- gammaParam[5]
                beta2  <- gammaParam[6]
      
                phi1 <- alpha1^2 + beta1
                phi2 <- alpha2^2 + beta2
      
                T1 <-as.vector(c(1, -2*alpha1, phi1))
                T2[2,1]<-T2[3,2]<-T2[4,3] <- -2*alpha2
                T2[3,1] <- T2[4,2] <- T2[5,3]<- phi2
      
                # see (A6) Liang & Browne
                asup2 = (T2 %*% T1) * lambda
      
                b0 = Ksi
                b1 = asup2[1]/1
                b2 = asup2[2]/2
                b3 = asup2[3]/3
                b4 = asup2[4]/4
                b5 = asup2[5]/5
      
                y = data[,item]
      
                # create logit polynomial 
                m = b0 + b1*thetaInit + b2*thetaInit^2 + b3*thetaInit^3 +    b4*thetaInit^4 + b5*thetaInit^5
      
                ##  avoid log(0)
                Pm <- P(m);
      
                Pm[Pm=="NaN"]<-.5
                Pm[Pm==1]   <-1-1e-12
                Pm[Pm<1e-12]<-1e-12  
      
                # this is FHAT
                return(-mean( y*log(Pm) + (1 - y) * log(1 - Pm)  ))
          }#END k=2 
    
    
          if(k==3){
                Ksi = gammaParam[1]
                lambda = gammaParam[2]
      
                alpha1 <- gammaParam[3]
                beta1  <- gammaParam[4]
                alpha2 <- gammaParam[5]
                beta2  <- gammaParam[6]
                alpha3 <- gammaParam[7]
                beta3  <- gammaParam[8]
      
                phi1 <- alpha1^2 + beta1
                phi2 <- alpha2^2 + beta2
                phi3 <- alpha3^2 + beta3
      
                T1 <-as.vector(c(1, -2*alpha1, phi1))
                T2[2,1]<-T2[3,2]<-T2[4,3] <- -2*alpha2
                T2[3,1] <- T2[4,2] <- T2[5,3]<- phi2
                T3[2,1]<-T3[3,2]<-T3[4,3]<-T3[5,4]<-T3[6,5]<- -2*alpha3
                T3[3,1]<-T3[4,2]<-T3[5,3]<-T3[6,4]<-T3[7,5]<- phi3
                # see (A6) Liang & Browne
                asup3 = (T3 %*% T2 %*% T1) * lambda
      
                b0 = Ksi
                b1 = asup3[1]/1
                b2 = asup3[2]/2
                b3 = asup3[3]/3
                b4 = asup3[4]/4
                b5 = asup3[5]/5
                b6 = asup3[6]/6
                b7 = asup3[7]/7
      
      
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
    
  } #END fncFMP ##################################################

  
  
  
  
  if(k>3) stop("\n\n*** FATAL ERROR:k must be <= 3 ***\n\n")
  
    if(k==0) Constraints = c(-Inf,0)
    if(k==1) Constraints = c(-Inf,0,-Inf,0)
    if(k==2) Constraints = c(-Inf,0,-Inf,0,-Inf,0)
    if(k==3) Constraints = c(-Inf,0,-Inf,0,-Inf,0,-Inf,0)
  
    Nparam<-length(startvals)
    
    out <- optim(par=startvals, 
                 fn=fncFMP, 
                 method="L-BFGS-B",
                 lower= Constraints,
                 k=k,
                 data=data,
                 item=item,
                 thetaInit=thetaInit,
                 control=list(ndeps=rep(eps,Nparam),
                              maxit=2000))
    
  ## compute scaled (pseudo) AIC and BIC
     NSubj <- nrow(data)
     q <- 2 * k + 2
     AIC <- 2 * out$value + (2/NSubj) * q
     BIC <- 2 * out$value + (log(NSubj)/NSubj) * q
  
     
   
     list( b = gammaTob(out$par,k),
         gamma=out$par, 
         FHAT=out$value,
         counts=out$counts,
         AIC = AIC,
         BIC = BIC,
         convergence = out$convergence)
  
} #END FMP




###################################################################
#                 END OF FUNCTION DEFINITIONS
###################################################################



