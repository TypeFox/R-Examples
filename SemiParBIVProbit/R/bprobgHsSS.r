bprobgHsSS <- function(params, respvec, VC, sp = NULL, qu.mag = NULL){

  epsilon <- 0.0000001
  max.p   <- 0.9999999
  
  eta1 <- VC$X1%*%params[1:VC$X1.d2]
  eta2 <- VC$X2%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]
  etad <- NULL ## New bit

  pd1 <- probm(eta1, VC$margins[1], only.pr = FALSE)
  pd2 <- probm(eta2, VC$margins[2], only.pr = FALSE)
  
  p1 <- pd1$pr; d.n1 <- pd1$d.n 
  p2 <- pd2$pr; d.n2 <- pd2$d.n  

  der2p1.dereta12 <- pd1$der2p.dereta
  der2p2.dereta22 <- pd2$der2p.dereta


  if(is.null(VC$X3))  teta.st <- etad <- params[(VC$X1.d2+VC$X2.d2+1)]
  if(!is.null(VC$X3)) teta.st <- etad <- VC$X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]
  
  
########################################################################################################  
  
resT    <- teta.tr(VC, teta.st)
teta.st <- resT$teta.st
teta    <- resT$teta


  p11 <- BiCDF(p1[VC$inde], p2, VC$nC, teta)

########################################################################################################

  p10 <- pmax(p1[VC$inde] - p11, epsilon)
  p0  <- pmax(1 - p1, epsilon)


  l.par1          <- respvec$cy1*log(p0) 
  l.par1[VC$inde] <- respvec$y1.y2*log(p11) + respvec$y1.cy2*log(p10) 
  l.par           <- VC$weights*l.par1 


dH <- copgHs(p1[VC$inde], p2, eta1=NULL, eta2=NULL, teta, teta.st, VC$BivD)

c.copula.be1   <- dH$c.copula.be1
c.copula.be2   <- dH$c.copula.be2
c.copula.theta <- dH$c.copula.theta 
         
                                
c.copula2.be1 <- dH$c.copula2.be1  
c.copula2.be2 <- dH$c.copula2.be2 


bit1.b1b1 <- c.copula2.be1*d.n1[VC$inde]^2 + c.copula.be1*der2p1.dereta12[VC$inde]                                                                                                                                
bit2.b1b1 <- -c.copula2.be1*d.n1[VC$inde]^2  + (1-c.copula.be1)*der2p1.dereta12[VC$inde]
bit3.b1b1 <- (-der2p1.dereta12*p0-d.n1^2)/p0^2
bit1.b2b2 <- c.copula2.be2*d.n2^2 + c.copula.be2*der2p2.dereta22


bit2.b2b2 <- -bit1.b2b2


c.copula2.be1be2 <- dH$c.copula2.be1be2
bit1.b1b2 <- c.copula2.be1be2*d.n1[VC$inde]*d.n2
bit2.b1b2 <- -bit1.b1b2

c.copula2.be1th <- dH$c.copula2.be1th 
bit1.b1th <- c.copula2.be1th*d.n1[VC$inde]
bit2.b1th <- -bit1.b1th 

c.copula2.be2th <- dH$c.copula2.be2th
bit1.b2th <- c.copula2.be2th*d.n2
bit2.b2th <- -bit1.b2th 

bit1.th2 <- dH$bit1.th2
bit2.th2 <- -bit1.th2



  dl.dbe11          <- d.n1*respvec$cy1/-p0 
  dl.dbe11[VC$inde] <- d.n1[VC$inde]*( (respvec$y1.y2*c.copula.be1/p11) + (respvec$y1.cy2*(1-c.copula.be1)/p10)  ) 

  dl.dbe1 <-  VC$weights*dl.dbe11
                                        
  dl.dbe2 <-  VC$weights[VC$inde]*d.n2*( (respvec$y1.y2*c.copula.be2/p11)  +
                                (respvec$y1.cy2*(c.copula.be2)/(-p10)) )

  dl.drho <-  VC$weights[VC$inde]*( respvec$y1.y2*c.copula.theta/p11 + respvec$y1.cy2*(-c.copula.theta)/p10  )
  
  
if(VC$hess==TRUE){

  d2l.be1.be11 <- bit3.b1b1 
  d2l.be1.be11[VC$inde] <- (respvec$y1.y2*(bit1.b1b1*p11-(c.copula.be1*d.n1[VC$inde])^2)/p11^2+respvec$y1.cy2*(bit2.b1b1*p10-((1-c.copula.be1)*d.n1[VC$inde])^2)/p10^2 ) 


  d2l.be1.be1  <- -VC$weights*d2l.be1.be11

  d2l.be2.be2  <- -VC$weights[VC$inde]*(respvec$y1.y2*(bit1.b2b2*p11 - (c.copula.be2*d.n2)^2 )/p11^2+
                              respvec$y1.cy2*(bit2.b2b2*p10-(-c.copula.be2*d.n2)^2)/p10^2 )

  d2l.be1.be2  <- -VC$weights[VC$inde]*(respvec$y1.y2*(bit1.b1b2*p11-(c.copula.be1*d.n1[VC$inde]*c.copula.be2*d.n2))/p11^2+
                              respvec$y1.cy2*(bit2.b1b2*p10-((1-c.copula.be1)*d.n1[VC$inde])*(-c.copula.be2*d.n2))/p10^2)

  d2l.be1.rho  <- -VC$weights[VC$inde]*(respvec$y1.y2*(bit1.b1th*p11-(c.copula.be1*d.n1[VC$inde]*c.copula.theta))/p11^2+
                              respvec$y1.cy2*(bit2.b1th*p10-((1-c.copula.be1)*d.n1[VC$inde])*(-c.copula.theta))/p10^2 )

  d2l.be2.rho  <- -VC$weights[VC$inde]*(respvec$y1.y2*(bit1.b2th*p11-(c.copula.be2*d.n2*c.copula.theta))/p11^2+
                              respvec$y1.cy2*(bit2.b2th*p10-(-c.copula.be2*d.n2)*(-c.copula.theta))/p10^2 )

  d2l.rho.rho  <- -VC$weights[VC$inde]*(respvec$y1.y2*(bit1.th2*p11-c.copula.theta^2)/p11^2+
                              respvec$y1.cy2*(bit2.th2*p10-(-c.copula.theta)^2)/p10^2 )

}

if(VC$hess==FALSE){

fi <- p11/p1[VC$inde]
se <- p10/p1[VC$inde]


eta2a <- rep(0,length(eta1))
eta2a[VC$inde] <- eta2

p2r <- probm(eta2a, VC$margins[2], only.pr = FALSE)$pr

p11a <- BiCDF(p1, p2r, VC$nC, teta)
p10a <- pmax(p1 - p11a, epsilon)

c.copula.be1a  <- copgHs(p1, p2r, eta1=NULL, eta2=NULL, teta, teta.st, VC$BivD)$c.copula.be1
c.copula2.be1a <- copgHs(p1, p2r, eta1=NULL, eta2=NULL, teta, teta.st, VC$BivD)$c.copula2.be1

bit1.b1b1a     <- c.copula2.be1a*d.n1^2 + c.copula.be1a*der2p1.dereta12 
bit2.b1b1a     <- -c.copula2.be1a*d.n1^2  + (1-c.copula.be1a)*der2p1.dereta12


  d2l.be1.be1  <- -VC$weights*( p0*bit3.b1b1 + (bit1.b1b1a*p11a-(c.copula.be1a*d.n1)^2)/p11a + (bit2.b1b1a*p10a-((1-c.copula.be1a)*d.n1)^2)/p10a )
  
  d2l.be2.be2  <- -VC$weights[VC$inde]*(fi*(bit1.b2b2*p11 - (c.copula.be2*d.n2)^2 )/p11^2+
                            se*(bit2.b2b2*p10-(-c.copula.be2*d.n2)^2)/p10^2 )#*(1-respvec$cy1)

  d2l.be1.be2  <- -VC$weights[VC$inde]*(fi*(bit1.b1b2*p11-(c.copula.be1*d.n1[VC$inde]*c.copula.be2*d.n2))/p11^2+
                            se*(bit2.b1b2*p10-((1-c.copula.be1)*d.n1[VC$inde])*(-c.copula.be2*d.n2))/p10^2)#*(1-respvec$cy1)

  d2l.be1.rho  <- -VC$weights[VC$inde]*(fi*(bit1.b1th*p11-(c.copula.be1*d.n1[VC$inde]*c.copula.theta))/p11^2+
                            se*(bit2.b1th*p10-((1-c.copula.be1)*d.n1[VC$inde])*(-c.copula.theta))/p10^2 )#*(1-respvec$cy1)

  d2l.be2.rho  <- -VC$weights[VC$inde]*(fi*(bit1.b2th*p11-(c.copula.be2*d.n2*c.copula.theta))/p11^2+
                            se*(bit2.b2th*p10-(-c.copula.be2*d.n2)*(-c.copula.theta))/p10^2 )#*(1-respvec$cy1)

  d2l.rho.rho  <- -VC$weights[VC$inde]*(fi*(bit1.th2*p11-c.copula.theta^2)/p11^2+
                            se*(bit2.th2*p10-(-c.copula.theta)^2)/p10^2 )#*(1-respvec$cy1)

}

      
      
      if( is.null(VC$X3) ){
      
        be1.be1 <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
        be2.be2 <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
        be1.be2 <- crossprod(VC$X1[VC$inde,]*c(d2l.be1.be2),VC$X2)
        be1.rho <- t(t(rowSums(t(VC$X1[VC$inde,]*c(d2l.be1.rho)))))
        be2.rho <- t(t(rowSums(t(VC$X2*c(d2l.be2.rho)))))
        
        H <- rbind( cbind( be1.be1    , be1.be2    , be1.rho ), 
                    cbind( t(be1.be2) , be2.be2    , be2.rho ), 
                    cbind( t(be1.rho) , t(be2.rho) , sum(d2l.rho.rho) ) 
                  ) 
                  
                  
               
               G   <- -c( colSums( c(dl.dbe1)*VC$X1 ),
                          colSums( c(dl.dbe2)*VC$X2 ),
                          sum( dl.drho )  )
          
      }
      
      if( !is.null(VC$X3) ){
      
        be1.be1 <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
        be2.be2 <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
        be1.be2 <- crossprod(VC$X1[VC$inde,]*c(d2l.be1.be2),VC$X2)
        be1.rho <- crossprod(VC$X1[VC$inde,]*c(d2l.be1.rho),VC$X3)
        be2.rho <- crossprod(VC$X2*c(d2l.be2.rho),VC$X3)
        rho.rho <- crossprod(VC$X3*c(d2l.rho.rho),VC$X3)
        
        H <- rbind( cbind( be1.be1    , be1.be2    , be1.rho ), 
                    cbind( t(be1.be2) , be2.be2    , be2.rho ), 
                    cbind( t(be1.rho) , t(be2.rho) , rho.rho ) 
                  ) 
                  
                  
               
               G   <- -c( colSums( c(dl.dbe1)*VC$X1 ),
                          colSums( c(dl.dbe2)*VC$X2 ),
                          colSums( c(dl.drho)*VC$X3 )  )
          
      }

    
    
 res <- -sum(l.par)   
    
if( ( VC$l.sp1==0 && VC$l.sp2==0 && VC$l.sp3==0) || VC$fp==TRUE) ps <- list(S.h = 0, S.h1 = 0, S.h2 = 0) else ps <- pen(params, qu.mag, sp, VC)


if(VC$extra.regI == "pC" && VC$hess==FALSE) H <- regH(H, type = 1)
  
         S.res <- res
         res   <- S.res + ps$S.h1
         G     <- G + ps$S.h2
         H     <- H + ps$S.h  
        
if(VC$extra.regI == "sED") H <- regH(H, type = 2)   
  

if(!is.null(VC$X3)){ 

etadt <- rep(0, length(eta1))
etadt[VC$inde]    <- etad    ; etad     <- etadt

}

dl.dbe2t <- dl.drhot <- rep(0, length(eta1))

dl.dbe2t[VC$inde]  <- dl.dbe2	;dl.dbe2  <- dl.dbe2t     
dl.drhot[VC$inde]  <- dl.drho	;dl.drho  <- dl.drhot


         list(value=res, gradient=G, hessian=H, S.h=ps$S.h, l=S.res, l.par=l.par, ps = ps,
              p11=p11, p10=p10, p0=p0, eta1=eta1, eta2=eta2, etad=etad,
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.drho=dl.drho,
              BivD=VC$BivD, p1=p1)     

}




     























