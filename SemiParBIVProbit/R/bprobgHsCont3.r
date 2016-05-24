bprobgHsCont3 <- function(params, respvec, VC, sp = NULL, qu.mag = NULL, AT = FALSE){


  eta1 <- VC$X1%*%params[1:VC$X1.d2]
  eta2 <- VC$X2%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]
  etad <- etas <- etan <- NULL 

  epsilon <- 0.0000001 # 0.9999999 0.0001 # sqrt(.Machine$double.eps)
  max.p   <- 0.9999999
  
  
if(is.null(VC$X3)){  
  sigma2.st <- etas <- params[(VC$X1.d2 + VC$X2.d2 + 1)]
  nu.st     <- etan <- params[(VC$X1.d2 + VC$X2.d2 + 2)]
  teta.st   <- etad <- params[(VC$X1.d2 + VC$X2.d2 + 3)]
} 

if(!is.null(VC$X3)){  
  sigma2.st <- etas <- VC$X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]
  nu.st     <- etan <- VC$X4%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2)]
  teta.st   <- etad <- VC$X5%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2+VC$X5.d2)]
}  
  
  
    sstr1 <- esp.tr(sigma2.st, VC$margins[2])  
    sigma2.st <- sstr1$vrb.st 
    sigma2    <- sstr1$vrb 

    
sstr1 <- esp.tr(nu.st, VC$margins[2])  
nu.st <- sstr1$vrb.st 
nu    <- sstr1$vrb 

    eta2 <- eta.tr(eta2, VC$margins[2])
    

 dHs  <- distrHs(respvec$y2, eta2, sigma2, sigma2.st, nu, nu.st, margin2=VC$margins[2], naive = FALSE)
  
 pdf2                         <- dHs$pdf2
 p2                           <- dHs$p2 
 derpdf2.dereta2              <- dHs$derpdf2.dereta2 
 derpdf2.dersigma2.st         <- dHs$derpdf2.dersigma2.st 
 derp2.dersigma.st            <- dHs$derp2.dersigma.st
 derpdf2.dernu.st             <- dHs$derpdf2.dernu.st 
 derp2.dernu.st               <- dHs$derp2.nu.st
 derp2.dereta2                <- dHs$derp2.dereta2
 der2p2.dereta2eta2           <- dHs$der2p2.dereta2eta2 
 der2pdf2.dereta2             <- dHs$der2pdf2.dereta2
 der2p2.dersigma2.st2         <- dHs$der2p2.dersigma2.st2
 der2pdf2.dersigma2.st2       <- dHs$der2pdf2.dersigma2.st2
 
 
 
 der2p2.dernu.st2         <- dHs$der2p2.dernu.st2
 der2pdf2.dernu.st2       <- dHs$der2pdf2.dernu.st2
 
 der2p2.dereta2dersigma2.st   <- dHs$der2p2.dereta2dersigma2.st            
 der2pdf2.dereta2dersigma2.st <- dHs$der2pdf2.dereta2dersigma2.st  
 
 
 der2p2.dereta2dernu.st   <- dHs$der2p2.dereta2dernu.st            
 der2pdf2.dereta2dernu.st <- dHs$der2pdf2.dereta2dernu.st 
  
 der2p2.dersigma2.stdernu.st   <- dHs$der2p2.dersigma2.stdernu.st            
 der2pdf2.dersigma2.stdernu.st <- dHs$der2pdf2.sigma2.st2dernu.st 
  
  pd1 <- probm(eta1, VC$margins[1], bc = TRUE) 
  
  p1 <- 1 - pd1$pr #   pnorm(-eta1)

########################################################################################################  
  
resT    <- teta.tr(VC, teta.st)
teta.st <- resT$teta.st
teta    <- resT$teta
    

########################################################################################################

dH <- copgHs(p1,p2,eta1=NULL,eta2=NULL,teta,teta.st,VC$BivD)

h <- dH$c.copula.be2  

l.par <- VC$weights*( respvec$cy*log(h) + respvec$y1*log(1 - h) + log(pdf2) )
  
########################################################################################################
 
  c.copula2.be2    <- dH$c.copula2.be2 
  c.copula2.be1be2 <- dH$c.copula2.be1be2   
  c.copula2.be2th  <- dH$c.copula2.be2th 
  
  derp1.dereta1    <- pd1$derp1.dereta1 
  
  derh.dereta1     <- c.copula2.be1be2 * derp1.dereta1
  derh.dereta2     <- c.copula2.be2 * derp2.dereta2

  dl.dbe1      <- VC$weights*( derh.dereta1 *(respvec$cy/h - respvec$y1/(1-h)) ) 
  dl.dbe2      <- VC$weights*( derh.dereta2 * (respvec$cy/h - respvec$y1/(1-h)) + derpdf2.dereta2/pdf2 )
  dl.dsigma.st <- VC$weights*( c.copula2.be2 * derp2.dersigma.st *(respvec$cy/h - respvec$y1/(1-h)) + derpdf2.dersigma2.st/pdf2 )
  dl.dnu.st    <- VC$weights*( c.copula2.be2 * derp2.dernu.st *(respvec$cy/h - respvec$y1/(1-h)) + derpdf2.dernu.st/pdf2 )
  dl.dteta.st  <- VC$weights*( c.copula2.be2th*(respvec$cy/h - respvec$y1/(1-h)) )                     
 
 
######################################################################################################## 
 
  BITS <- copgHsCont(p1, p2, teta, teta.st, VC)
  
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
  
  der2h.derp2dernu.st           <- der2h.derp2p2*derp2.dernu.st 

  der2h.dersigma2.st2           <- der2h.derp2dersigma2.st*derp2.dersigma.st + c.copula2.be2*der2p2.dersigma2.st2
  der2h.dernu.st2               <- der2h.derp2dernu.st*derp2.dernu.st        + c.copula2.be2*der2p2.dernu.st2
  
  derh.dersigma2.st             <- c.copula2.be2 * derp2.dersigma.st
  derh.dernu.st                 <- c.copula2.be2 * derp2.dernu.st
  
  der2h.dereta1.dereta2         <- der2h.derp1p2*derp1.dereta1*derp2.dereta2                                                                   
  der2h.dereta1.derteta.st      <- der2h.derp1teta*derp1.dereta1*derteta.derteta.st  
  der2h.dereta1.dersigma2.st    <- der2h.derp1p2 * derp2.dersigma.st*derp1.dereta1  
  
  der2h.dereta1.dernu.st        <- der2h.derp1p2 * derp2.dernu.st*derp1.dereta1 
  
  der2h.dereta2.derteta.st      <- der2h.derp2teta*derp2.dereta2*derteta.derteta.st 
  
  der2h.derteta.st.dersigma2.st <- der2h.derp2teta* derteta.derteta.st*derp2.dersigma.st 
  der2h.derteta.st.dernu.st     <- der2h.derp2teta* derteta.derteta.st*derp2.dernu.st 
  
  der2h.dersigma2.st.dernu.st   <- der2h.derp2dernu.st*derp2.dersigma.st + c.copula2.be2*der2p2.dersigma2.stdernu.st
  
  der2h.dereta2.dersigma2.st    <- der2h.derp2dersigma2.st*derp2.dereta2 + c.copula2.be2*der2p2.dereta2dersigma2.st  
  der2h.dereta2.dernu.st        <- der2h.derp2dernu.st*derp2.dereta2 + c.copula2.be2*der2p2.dereta2dernu.st 
  
  der2h.dereta1.dereta1         <- der2h.derp1p1*derp1.dereta1^2 + c.copula2.be1be2*der2p1.dereta1eta1      
  
  
  d2l.be1.be1      <- -VC$weights*(der2h.dereta1.dereta1 *(respvec$cy/h - respvec$y1/(1-h)) - derh.dereta1^2 * (respvec$cy/h^2 + respvec$y1/(1-h)^2) )
  d2l.be2.be2      <- -VC$weights*(der2h.dereta2.dereta2 *(respvec$cy/h - respvec$y1/(1-h)) - derh.dereta2^2 * (respvec$cy/h^2 + respvec$y1/(1-h)^2) + (der2pdf2.dereta2*pdf2-(derpdf2.dereta2)^2)/(pdf2)^2 )
  d2l.rho.rho      <- -VC$weights*(der2h.derteta.st2*(respvec$cy/h - respvec$y1/(1-h)) - c.copula2.be2th^2 * (respvec$cy/h^2 + respvec$y1/(1-h)^2) )
  d2l.sigma.sigma  <- -VC$weights*(der2h.dersigma2.st2 *(respvec$cy/h - respvec$y1/(1-h)) - derh.dersigma2.st^2 * (respvec$cy/h^2 + respvec$y1/(1-h)^2) + (der2pdf2.dersigma2.st2*pdf2-(derpdf2.dersigma2.st)^2)/(pdf2)^2 )
  d2l.nu.nu        <- -VC$weights*(der2h.dernu.st2 *(respvec$cy/h - respvec$y1/(1-h)) - derh.dernu.st^2 * (respvec$cy/h^2 + respvec$y1/(1-h)^2) + (der2pdf2.dernu.st2*pdf2-(derpdf2.dernu.st)^2)/(pdf2)^2 )
                                                                                                                                                   
  d2l.be1.be2      <- -VC$weights*(der2h.dereta1.dereta2 *(respvec$cy/h - respvec$y1/(1-h)) - derh.dereta1*derh.dereta2 * (respvec$cy/h^2 + respvec$y1/(1-h)^2)  )
  d2l.be1.rho      <- -VC$weights*(der2h.dereta1.derteta.st *(respvec$cy/h - respvec$y1/(1-h)) - derh.dereta1 * c.copula2.be2th* (respvec$cy/h^2 + respvec$y1/(1-h)^2) )
  d2l.be1.sigma    <- -VC$weights*(der2h.dereta1.dersigma2.st *(respvec$cy/h - respvec$y1/(1-h)) - derh.dereta1 * derh.dersigma2.st* (respvec$cy/h^2 + respvec$y1/(1-h)^2) )
  d2l.be1.nu       <- -VC$weights*(der2h.dereta1.dernu.st *(respvec$cy/h - respvec$y1/(1-h)) - derh.dereta1 * derh.dernu.st* (respvec$cy/h^2 + respvec$y1/(1-h)^2) )

  d2l.be2.rho      <- -VC$weights*(der2h.dereta2.derteta.st *(respvec$cy/h - respvec$y1/(1-h)) - derh.dereta2*c.copula2.be2th * (respvec$cy/h^2 + respvec$y1/(1-h)^2)  )
  d2l.be2.sigma    <- -VC$weights*(der2h.dereta2.dersigma2.st *(respvec$cy/h - respvec$y1/(1-h)) - derh.dereta2*derh.dersigma2.st * (respvec$cy/h^2 + respvec$y1/(1-h)^2) + (der2pdf2.dereta2dersigma2.st*pdf2-(derpdf2.dereta2*derpdf2.dersigma2.st))/(pdf2)^2 ) 
  d2l.be2.nu       <- -VC$weights*(der2h.dereta2.dernu.st *(respvec$cy/h - respvec$y1/(1-h)) - derh.dereta2*derh.dernu.st * (respvec$cy/h^2 + respvec$y1/(1-h)^2) + (der2pdf2.dereta2dernu.st*pdf2-(derpdf2.dereta2*derpdf2.dernu.st))/(pdf2)^2 ) 

  d2l.rho.sigma    <- -VC$weights*(der2h.derteta.st.dersigma2.st *(respvec$cy/h - respvec$y1/(1-h)) - derh.dersigma2.st * c.copula2.be2th* (respvec$cy/h^2 + respvec$y1/(1-h)^2) )
  d2l.rho.nu       <- -VC$weights*(der2h.derteta.st.dernu.st *(respvec$cy/h - respvec$y1/(1-h)) - derh.dernu.st * c.copula2.be2th* (respvec$cy/h^2 + respvec$y1/(1-h)^2) )

  d2l.sigma.nu     <- -VC$weights*(der2h.dersigma2.st.dernu.st *(respvec$cy/h - respvec$y1/(1-h)) - derh.dersigma2.st * derh.dernu.st * (respvec$cy/h^2 + respvec$y1/(1-h)^2) + (der2pdf2.dersigma2.stdernu.st*pdf2-(derpdf2.dersigma2.st*derpdf2.dernu.st))/(pdf2)^2 )
                                   
if( is.null(VC$X3) ){


  be1.be1   <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
  be2.be2   <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be1.be2   <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)
  be1.rho   <- t(t(rowSums(t(VC$X1*c(d2l.be1.rho)))))
  be1.sigma <- t(t(rowSums(t(VC$X1*c(d2l.be1.sigma))))) 
  be2.rho   <- t(t(rowSums(t(VC$X2*c(d2l.be2.rho)))))
  be2.sigma <- t(t(rowSums(t(VC$X2*c(d2l.be2.sigma))))) 
  be2.nu <- t(t(rowSums(t(VC$X2*c(d2l.be2.nu)))))
  be1.nu <- t(t(rowSums(t(VC$X1*c(d2l.be1.nu)))))   
  


         
 
   H <- rbind( cbind( be1.be1     ,  be1.be2     ,  be1.sigma   , be1.nu,    be1.rho   ), 
               cbind( t(be1.be2)  ,  be2.be2     ,  be2.sigma   , be2.nu,    be2.rho   ), 
               cbind( t(be1.sigma),  t(be2.sigma),  sum(d2l.sigma.sigma) , sum(d2l.sigma.nu),  sum(d2l.rho.sigma) ),
               cbind( t(be1.nu),     t(be2.nu),     sum(d2l.sigma.nu) ,    sum(d2l.nu.nu), sum(d2l.rho.nu)  ),
               cbind( t(be1.rho)  ,  t(be2.rho)  ,  sum(d2l.rho.sigma), sum(d2l.rho.nu), sum(d2l.rho.rho)   ) 
               
               
             )  
         
         
         
         
         
  G   <- -c( colSums( c(dl.dbe1)*VC$X1 ) ,
             colSums( c(dl.dbe2)*VC$X2 ) ,
             sum( dl.dsigma.st ),
             sum( dl.dnu.st ),
             sum( dl.dteta.st )    )
    
}




if( !is.null(VC$X3) ){

  be1.be1   <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
  be2.be2   <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be1.be2   <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)
  
  be1.rho   <- crossprod(VC$X1*c(d2l.be1.rho),  VC$X5)                                     
  be1.sigma <- crossprod(VC$X1*c(d2l.be1.sigma),VC$X3)  
  be1.nu    <- crossprod(VC$X1*c(d2l.be1.nu),VC$X4)                                   

  be2.rho   <- crossprod(VC$X2*c(d2l.be2.rho),  VC$X5)                                     
  be2.sigma <- crossprod(VC$X2*c(d2l.be2.sigma),VC$X3) 
  be2.nu    <- crossprod(VC$X2*c(d2l.be2.nu),VC$X4)  
  
  
  sigma.sigma <- crossprod(VC$X3*c(d2l.sigma.sigma),VC$X3)
  sigma.nu    <- crossprod(VC$X3*c(d2l.sigma.nu),VC$X4)   
  
  sigma.rho   <- crossprod(VC$X3*c(d2l.rho.sigma),VC$X5) 
 
  rho.rho     <- crossprod(VC$X5*c(d2l.rho.rho),    VC$X5) 
  rho.nu      <- crossprod(VC$X4*c(d2l.rho.nu),    VC$X5) 
  nu.nu      <- crossprod(VC$X4*c(d2l.nu.nu),    VC$X4)    
  

  
  

  H <- rbind( cbind( be1.be1     ,  be1.be2     ,  be1.sigma   , be1.nu      , be1.rho   ), 
              cbind( t(be1.be2)  ,  be2.be2     ,  be2.sigma   , be2.nu      , be2.rho   ), 
              cbind( t(be1.sigma),  t(be2.sigma),  sigma.sigma , sigma.nu    , sigma.rho ),
              cbind( t(be1.nu)   ,  t(be2.nu)   ,  t(sigma.nu) , nu.nu       , rho.nu    ),
              cbind( t(be1.rho)  ,  t(be2.rho)  ,  t(sigma.rho), t(rho.nu)   , rho.rho   ) 
              
              
             )  
            
   
  G   <- -c( colSums(      c(dl.dbe1)*VC$X1 ) ,
             colSums(      c(dl.dbe2)*VC$X2 ) ,
             colSums( c(dl.dsigma.st)*VC$X3 ) ,
             colSums(    c(dl.dnu.st)*VC$X4 ) ,
             colSums(  c(dl.dteta.st)*VC$X5 ) )   
   
    
}




      res <- -sum(l.par)




if( ( VC$l.sp1==0 && VC$l.sp2==0 && VC$l.sp3==0 && VC$l.sp4==0 && VC$l.sp5==0 ) || VC$fp==TRUE) ps <- list(S.h = 0, S.h1 = 0, S.h2 = 0) else ps <- pen(params, qu.mag, sp, VC)

 
if(VC$extra.regI == "pC") H <- regH(H, type = 1)
  
         S.res <- res
         res   <- S.res + ps$S.h1
         G     <- G + ps$S.h2
         H     <- H + ps$S.h  
        
if(VC$extra.regI == "sED") H <- regH(H, type = 2)   
  
  


  
  
  

         list(value=res, gradient=G, hessian=H, S.h=ps$S.h, l=S.res, l.par=l.par, ps = ps, etas = etas,
              eta1 = eta1, eta2 = eta2, etad = etad, etan = etan,
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.dsigma.st = dl.dsigma.st, 
              dl.dnu.st = dl.dnu.st, dl.dteta.st = dl.dteta.st,
#d2l.be1.be1     = d2l.be1.be1    , 
#d2l.be2.be2     = d2l.be2.be2    ,
#d2l.rho.rho     = d2l.rho.rho    ,
#d2l.sigma.sigma = d2l.sigma.sigma,
#d2l.nu.nu       = d2l.nu.nu      ,          
#d2l.be1.be2     = d2l.be1.be2    ,
#d2l.be1.rho     = d2l.be1.rho    ,
#d2l.be1.sigma   = d2l.be1.sigma  ,
#d2l.be1.nu      = d2l.be1.nu     ,
#d2l.be2.rho     = d2l.be2.rho    ,
#d2l.be2.sigma   = d2l.be2.sigma  ,
#d2l.be2.nu      = d2l.be2.nu     ,
#d2l.rho.sigma   = d2l.rho.sigma  ,
#d2l.rho.nu      = d2l.rho.nu     ,
#d2l.sigma.nu    = d2l.sigma.nu   ,            
BivD=VC$BivD, p1=1-p1, p2=p2, theta.star = teta.st)      

}




     























