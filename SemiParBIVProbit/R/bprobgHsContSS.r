bprobgHsContSS <- function(params, respvec, VC, sp = NULL, qu.mag = NULL, AT = FALSE){

  eta1 <- VC$X1%*%params[1:VC$X1.d2]
  eta2 <- VC$X2%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]
  etad <- etas <- l.ln <- etan <- NULL 

  epsilon <- 0.0000001 # 0.9999999 0.0001 # sqrt(.Machine$double.eps)
  max.p   <- 0.9999999
  
  
if(is.null(VC$X3)){  
  sigma2.st <- etas <- params[(VC$X1.d2 + VC$X2.d2 + 1)]
  teta.st   <- etad <- params[(VC$X1.d2 + VC$X2.d2 + 2)]
} 

if(!is.null(VC$X3)){  
  sigma2.st <- etas <- VC$X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]
  teta.st   <- etad <- VC$X4%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2)] #remove and change teta.st
}  
  
  
    sstr1 <- esp.tr(sigma2.st, VC$margins[2])  
    sigma2.st <- sstr1$vrb.st 
    sigma2    <- sstr1$vrb 


    eta2 <- eta.tr(eta2, VC$margins[2])
    

 dHs  <- distrHs(respvec$y2, eta2, sigma2, sigma2.st, margin2=VC$margins[2], naive = FALSE)
  
 pdf2                         <- dHs$pdf2
 p2                           <- dHs$p2 
 derpdf2.dereta2              <- dHs$derpdf2.dereta2 
 derpdf2.dersigma2.st         <- dHs$derpdf2.dersigma2.st 
 derp2.dersigma.st            <- dHs$derp2.dersigma.st
 derp2.dereta2                <- dHs$derp2.dereta2
 der2p2.dereta2eta2           <- dHs$der2p2.dereta2eta2 
 der2pdf2.dereta2             <- dHs$der2pdf2.dereta2
 der2p2.dersigma2.st2         <- dHs$der2p2.dersigma2.st2
 der2pdf2.dersigma2.st2       <- dHs$der2pdf2.dersigma2.st2
 
 der2p2.dereta2dersigma2.st   <- dHs$der2p2.dereta2dersigma2.st            
 der2pdf2.dereta2dersigma2.st <- dHs$der2pdf2.dereta2dersigma2.st  
 
  
  pd1 <- probm(eta1, VC$margins[1], bc = TRUE) 
  
  p1 <- 1 - pd1$pr

########################################################################################################  
  
resT    <- teta.tr(VC, teta.st)
teta.st <- resT$teta.st
teta    <- resT$teta
    

########################################################################################################

dH <- copgHs(p1[VC$inde], p2, eta1 = NULL, eta2 = NULL, teta, teta.st, VC$BivD)

h <- dH$c.copula.be2  

l.par1          <- log(p1)
l.par1[VC$inde] <- log(pdf2 * (1-h)) 

l.par <- VC$weights*l.par1   

  
########################################################################################################
 
  c.copula2.be2    <- dH$c.copula2.be2 
  c.copula2.be1be2 <- dH$c.copula2.be1be2   
  c.copula2.be2th  <- dH$c.copula2.be2th 
  
  derp1.dereta1    <- pd1$derp1.dereta1 
  
  derh.dereta1     <- c.copula2.be1be2 * derp1.dereta1[VC$inde]
  derh.dereta2     <- c.copula2.be2 * derp2.dereta2
    
  dl.dbe11          <- 1/p1* derp1.dereta1 
  dl.dbe11[VC$inde] <- -c.copula2.be1be2 / (1-h) * derp1.dereta1[VC$inde]  
  dl.dbe1           <- VC$weights*dl.dbe11  
  
  dl.dbe2      <- VC$weights[VC$inde]*( 1/(pdf2* (1-h)) * (derpdf2.dereta2      * (1-h) - pdf2 * derh.dereta2 )  )
  
  dl.dsigma.st <- VC$weights[VC$inde]*( 1/(pdf2* (1-h)) * (derpdf2.dersigma2.st * (1-h) - pdf2 * c.copula2.be2 * derp2.dersigma.st )  )
    
  dl.dteta.st  <- VC$weights[VC$inde]*( 1/(1-h) * (-c.copula2.be2th) )                       
 
 
######################################################################################################## 
 
  BITS <- copgHsCont(p1[VC$inde], p2, teta, teta.st, VC)
  
  der2h.derp2p2              <- BITS$der2h.derp2p2 
  der2h.derteta.teta.st      <- BITS$der2h.derteta.teta.st  
  derteta.derteta.st         <- BITS$derteta.derteta.st 
  der2teta.derteta.stteta.st <- BITS$der2teta.derteta.stteta.st  
  der2h.derp1p2              <- BITS$der2h.derp1p2  
  der2h.derp1teta            <- BITS$der2h.derp1teta                                     
  der2h.derp2teta            <- BITS$der2h.derp2teta  
  der2h.derp1p1              <- BITS$der2h.derp1p1
  
 
  der2p1.dereta1eta1 <- pd1$der2p1.dereta1eta1     
                    
  der2h.dereta2.dereta2         <- der2h.derp2p2*derp2.dereta2^2 + c.copula2.be2*der2p2.dereta2eta2                                        
  der2h.derteta.st2             <- der2h.derteta.teta.st*derteta.derteta.st^2 + c.copula2.be2th/derteta.derteta.st* der2teta.derteta.stteta.st                                                                                        
  der2h.derp2dersigma2.st       <- der2h.derp2p2*derp2.dersigma.st 
  
  #der2h.derp2dernu.st           <- der2h.derp2p2*derp2.dernu.st 

  der2h.dersigma2.st2           <- der2h.derp2dersigma2.st*derp2.dersigma.st + c.copula2.be2*der2p2.dersigma2.st2
  
  derh.dersigma2.st             <- c.copula2.be2 * derp2.dersigma.st
  
  der2h.dereta1.dereta2         <- der2h.derp1p2*derp1.dereta1[VC$inde]*derp2.dereta2                                                                   
  der2h.dereta1.derteta.st      <- der2h.derp1teta*derp1.dereta1[VC$inde]*derteta.derteta.st  
  der2h.dereta1.dersigma2.st    <- der2h.derp1p2 * derp2.dersigma.st*derp1.dereta1[VC$inde]  
  
  
  der2h.dereta2.derteta.st      <- der2h.derp2teta*derp2.dereta2*derteta.derteta.st 
  
  der2h.derteta.st.dersigma2.st <- der2h.derp2teta* derteta.derteta.st*derp2.dersigma.st 
  
  
  der2h.dereta2.dersigma2.st    <- der2h.derp2dersigma2.st*derp2.dereta2 + c.copula2.be2*der2p2.dereta2dersigma2.st  
  
  
  
  
d2l.be1.be11          <- -1/p1^2*derp1.dereta1*derp1.dereta1 + 1/p1*der2p1.dereta1eta1  
d2l.be1.be11[VC$inde] <- -(derh.dereta1/(1-h)^2 * c.copula2.be1be2 + 1/(1-h)* der2h.derp1p1 * derp1.dereta1[VC$inde])*derp1.dereta1[VC$inde] - c.copula2.be1be2/(1-h)*der2p1.dereta1eta1[VC$inde]   
d2l.be1.be1           <- -VC$weights*d2l.be1.be11   
   


  d2l.be2.be2      <- -VC$weights[VC$inde]*(( -(derpdf2.dereta2 * (1-h) - pdf2 * derh.dereta2 ) / (pdf2*(1-h))^2  ) * ( (derpdf2.dereta2  * (1-h) - pdf2 * derh.dereta2 ))
                                    + 1/(pdf2*(1-h)) * (der2pdf2.dereta2 * (1-h) - derpdf2.dereta2 * derh.dereta2 - (derpdf2.dereta2 * derh.dereta2 + pdf2*der2h.dereta2.dereta2) ) )
                                    
  d2l.rho.rho      <- -VC$weights[VC$inde]*((-c.copula2.be2th/(h-1)^2) * c.copula2.be2th + (1/(h-1) * der2h.derteta.st2))
  
  d2l.sigma.sigma  <- -VC$weights[VC$inde]*(( -(derpdf2.dersigma2.st * (1-h) - pdf2 * derh.dersigma2.st) / (pdf2*(1-h))^2  ) * ( (derpdf2.dersigma2.st  * (1-h) - pdf2 * derh.dersigma2.st ))
                                    + 1/(pdf2*(1-h)) * ( der2pdf2.dersigma2.st2 * (1-h) - derpdf2.dersigma2.st * derh.dersigma2.st - (derpdf2.dersigma2.st * derh.dersigma2.st + pdf2*der2h.dersigma2.st2) ) )
                                                                 
                                                                                                                                                     
  d2l.be1.be2      <- -VC$weights[VC$inde]*(-( derh.dereta2/(1-h)^2 * c.copula2.be1be2 + 1/(1-h)*der2h.derp1p2*derp2.dereta2)*derp1.dereta1[VC$inde])
  
  d2l.be1.rho      <- -VC$weights[VC$inde]*(-( c.copula2.be2th/(1-h)^2 * c.copula2.be1be2 + 1/(1-h)*derteta.derteta.st * der2h.derp1teta)*derp1.dereta1[VC$inde])
  
  d2l.be1.sigma    <- -VC$weights[VC$inde]*(-( derh.dersigma2.st/(1-h)^2 * c.copula2.be1be2 + 1/(1-h)*der2h.derp1p2*derp2.dersigma.st)*derp1.dereta1[VC$inde])
  
  
  d2l.be2.rho      <- -VC$weights[VC$inde]*( -c.copula2.be2*derp2.dereta2/((h-1)^2)*c.copula2.be2th + 1/(h-1)*der2h.dereta2.derteta.st)
  
  d2l.be2.sigma    <- -VC$weights[VC$inde]*(( -(derpdf2.dersigma2.st * (1-h) - pdf2 * derh.dersigma2.st ) / (pdf2*(1-h))^2  ) * ( (derpdf2.dereta2  * (1-h) - pdf2 * derh.dereta2 ))
                                    + 1/(pdf2*(1-h)) * (der2pdf2.dereta2dersigma2.st * (1-h) - derpdf2.dereta2 * derh.dersigma2.st - (derpdf2.dersigma2.st * derh.dereta2 + pdf2*der2h.dereta2.dersigma2.st) ) )
                                    
                                    
  d2l.rho.sigma    <- -VC$weights[VC$inde]*( -c.copula2.be2*derp2.dersigma.st/((h-1)^2)*c.copula2.be2th + 1/(h-1)*der2h.derteta.st.dersigma2.st)
  
  
                         
                                
if( is.null(VC$X3) ){


  be1.be1   <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
  be2.be2   <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be1.be2   <- crossprod(VC$X1[VC$inde,]*c(d2l.be1.be2),VC$X2)  
  be1.rho   <- t(t(rowSums(t(VC$X1[VC$inde,]*c(d2l.be1.rho))))) 
  be1.sigma <- t(t(rowSums(t(VC$X1[VC$inde,]*c(d2l.be1.sigma)))))  
  be2.rho   <- t(t(rowSums(t(VC$X2*c(d2l.be2.rho)))))
  be2.sigma <- t(t(rowSums(t(VC$X2*c(d2l.be2.sigma))))) 
     
  


         
 
   H <- rbind( cbind( be1.be1     ,  be1.be2     ,  be1.sigma   ,    be1.rho   ), 
               cbind( t(be1.be2)  ,  be2.be2     ,  be2.sigma   ,    be2.rho   ), 
               cbind( t(be1.sigma),  t(be2.sigma),  sum(d2l.sigma.sigma) ,   sum(d2l.rho.sigma) ),
               cbind( t(be1.rho)  ,  t(be2.rho)  ,  sum(d2l.rho.sigma), sum(d2l.rho.rho)   ) 
               
               
             )  
         
         
         
         
         
  G   <- -c( colSums( c(dl.dbe1)*VC$X1 ) ,
             colSums( c(dl.dbe2)*VC$X2 ) ,
             sum( dl.dsigma.st ),
             sum( dl.dteta.st )    )
    
}




if( !is.null(VC$X3) ){

  be1.be1   <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
  be2.be2   <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be1.be2   <- crossprod(VC$X1[VC$inde,]*c(d2l.be1.be2),VC$X2)
  
  be1.rho   <- crossprod(VC$X1[VC$inde,]*c(d2l.be1.rho),  VC$X4)                                     
  be1.sigma <- crossprod(VC$X1[VC$inde,]*c(d2l.be1.sigma),VC$X3)  
                                

  be2.rho   <- crossprod(VC$X2*c(d2l.be2.rho),  VC$X4)                                     
  be2.sigma <- crossprod(VC$X2*c(d2l.be2.sigma),VC$X3) 

  
  
  sigma.sigma <- crossprod(VC$X3*c(d2l.sigma.sigma),VC$X3)
  
  
  sigma.rho   <- crossprod(VC$X3*c(d2l.rho.sigma),VC$X4) 
 
  rho.rho     <- crossprod(VC$X4*c(d2l.rho.rho),  VC$X4) 

  
  

  
  

  H <- rbind( cbind( be1.be1     ,  be1.be2     ,  be1.sigma       , be1.rho   ), 
              cbind( t(be1.be2)  ,  be2.be2     ,  be2.sigma       , be2.rho   ), 
              cbind( t(be1.sigma),  t(be2.sigma),  sigma.sigma     , sigma.rho ),
              cbind( t(be1.rho)  ,  t(be2.rho)  ,  t(sigma.rho)    , rho.rho   ) 
              
              
             )  
            
   
  G   <- -c( colSums(      c(dl.dbe1)*VC$X1 ) ,
             colSums(      c(dl.dbe2)*VC$X2 ) ,
             colSums( c(dl.dsigma.st)*VC$X3 ) ,
             colSums(  c(dl.dteta.st)*VC$X4 ) )   
   
    
}




      res <- -sum(l.par)




if( ( VC$l.sp1==0 && VC$l.sp2==0 && VC$l.sp3==0 && VC$l.sp4==0 && VC$l.sp5==0 ) || VC$fp==TRUE) ps <- list(S.h = 0, S.h1 = 0, S.h2 = 0) else ps <- pen(params, qu.mag, sp, VC)

 
if(VC$extra.regI == "pC") H <- regH(H, type = 1)
  
         S.res <- res
         res   <- S.res + ps$S.h1
         G     <- G + ps$S.h2
         H     <- H + ps$S.h  
        
if(VC$extra.regI == "sED") H <- regH(H, type = 2)   
  
  


  
if(!is.null(VC$X3)){ 

teta.stt <- etast <- etadt <- etant <- rep(0, length(eta1))
etast[VC$inde]    <- etas    ; etas     <- etast
etadt[VC$inde]    <- etad    ; etad     <- etadt
teta.stt[VC$inde] <- teta.st ; teta.st	<- teta.stt

}

p2t <- dl.dbe2t <- dl.dsigma.stt <- dl.dteta.stt <- rep(0, length(eta1))

p2t[VC$inde]           <- p2            ;p2           <- p2t          
dl.dbe2t[VC$inde]      <- dl.dbe2	;dl.dbe2      <- dl.dbe2t     
dl.dsigma.stt[VC$inde] <- dl.dsigma.st 	;dl.dsigma.st <- dl.dsigma.stt
dl.dteta.stt[VC$inde]  <- dl.dteta.st	;dl.dteta.st  <- dl.dteta.stt


 if(VC$margins[2] == "LN"){
    
 dHs  <- distrHsAT(exp(respvec$y2), eta2, sigma2, 1, margin2 = VC$margins[2])
 pdf2 <- dHs$pdf2
 p2   <- dHs$p2 
 dH   <- copgHsAT(p1[VC$inde], p2, teta, VC$BivD)
 h    <- dH$c.copula.be2  
 
 l.par1 <- log(p1)
 l.par1[VC$inde] <- log(pdf2 * (1-h)) 
  
 l.ln <- -sum(VC$weights*l.par1)  
 
 }



         list(value=res, gradient=G, hessian=H, S.h=ps$S.h, l=S.res, l.par=l.par, ps = ps, etas = etas, l.ln=l.ln, 
              eta1 = eta1, eta2 = eta2, etad = etad, etan = etan,
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.dsigma.st = dl.dsigma.st, 
              dl.dteta.st = dl.dteta.st,
              BivD=VC$BivD, p1=1-p1, p2=p2, theta.star = teta.st)      

}






     























