bprobgHsPO <- function(params, respvec, VC, sp = NULL, qu.mag = NULL){



  epsilon <- 0.0000001 # 0.9999999 0.0001 # sqrt(.Machine$double.eps)
  max.p   <- 0.9999999

  eta1 <- VC$X1%*%params[1:VC$X1.d2]
  eta2 <- VC$X2%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]
  etad <- NULL ## New bit


  pd1 <- probm(eta1, VC$margins[1], only.pr = FALSE)
  pd2 <- probm(eta2, VC$margins[2], only.pr = FALSE)
  
  p1 <- pd1$pr; d.n1 <- pd1$d.n 
  p2 <- pd2$pr; d.n2 <- pd2$d.n  

  if(is.null(VC$X3))  teta.st <- etad <- params[(VC$X1.d2+VC$X2.d2+1)]
  if(!is.null(VC$X3)) teta.st <- etad <- VC$X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]


########################################################################################################  
  
resT    <- teta.tr(VC, teta.st)
teta.st <- resT$teta.st
teta    <- resT$teta


p11 <-  BiCDF(p1, p2, VC$nC, teta)

########################################################################################################

  cp11 <- pmax(1 - p11, epsilon)
  
  l.par <- VC$weights*( respvec$y1*log(p11) + respvec$cy*log(cp11) )


dH <- copgHs(p1,p2,eta1=NULL,eta2=NULL,teta,teta.st,VC$BivD)

c.copula.be1   <- dH$c.copula.be1
c.copula.be2   <- dH$c.copula.be2
c.copula.theta <- dH$c.copula.theta 
        
c.copula2.be1 <- dH$c.copula2.be1   
c.copula2.be2 <- dH$c.copula2.be2 

bit1.b1b1 <- c.copula2.be1*(d.n1)^2-c.copula.be1*d.n1*eta1
bit1.b2b2 <- c.copula2.be2*(d.n2)^2-c.copula.be2*d.n2*eta2

c.copula2.be1be2 <- dH$c.copula2.be1be2
bit1.b1b2 <- c.copula2.be1be2 * d.n1 *d.n2

c.copula2.be1th <- dH$c.copula2.be1th 
bit1.b1th <- c.copula2.be1th*d.n1

c.copula2.be2th <- dH$c.copula2.be2th
bit1.b2th <- c.copula2.be2th*d.n2

bit1.th2 <- dH$bit1.th2


  dl.dbe1 <-  VC$weights*d.n1*( respvec$y1*c.copula.be1/p11   - respvec$cy*c.copula.be1/cp11   )  
                                          
  dl.dbe2 <-  VC$weights*d.n2*( respvec$y1*c.copula.be2/p11   - respvec$cy*c.copula.be2/cp11   )
                                                  
  dl.drho <-       VC$weights*( respvec$y1*c.copula.theta/p11 - respvec$cy*c.copula.theta/cp11 ) 

if(VC$hess==TRUE){

  d2l.be1.be1  <- -VC$weights*(respvec$y1*(bit1.b1b1*p11-(c.copula.be1*d.n1)^2)/p11^2 -
                           respvec$cy*(bit1.b1b1*cp11+(c.copula.be1*d.n1)^2)/cp11^2 )

  d2l.be2.be2  <- -VC$weights*(respvec$y1*(bit1.b2b2*p11 - (c.copula.be2*d.n2)^2 )/p11^2 -
                           respvec$cy*(bit1.b2b2*cp11 + (c.copula.be2*d.n2)^2 )/cp11^2 )

  d2l.be1.be2  <- -VC$weights*(respvec$y1*(bit1.b1b2*p11-(c.copula.be1*d.n1*c.copula.be2*d.n2))/p11^2 -
                           respvec$cy*(bit1.b1b2*cp11+(c.copula.be1*d.n1*c.copula.be2*d.n2))/cp11^2)

  d2l.be1.rho  <- -VC$weights*(respvec$y1*(bit1.b1th*p11-(c.copula.be1*d.n1*c.copula.theta))/p11^2 -
                           respvec$cy*(bit1.b1th*cp11+(c.copula.be1*d.n1*c.copula.theta))/cp11^2)

  d2l.be2.rho  <- -VC$weights*(respvec$y1*(bit1.b2th*p11-(c.copula.be2*d.n2*c.copula.theta))/p11^2 -
                           respvec$cy*(bit1.b2th*cp11+(c.copula.be2*d.n2*c.copula.theta))/cp11^2)

  d2l.rho.rho  <- -VC$weights*(respvec$y1*(bit1.th2*p11-c.copula.theta^2)/p11^2 -
                           respvec$cy*(bit1.th2*cp11+c.copula.theta^2)/cp11^2)

}  
  
if(VC$hess==FALSE){
                          
  d2l.be1.be1  <- -VC$weights*((bit1.b1b1*p11-(c.copula.be1*d.n1)^2)/p11 -
                            (bit1.b1b1*cp11+(c.copula.be1*d.n1)^2)/cp11 )

  d2l.be2.be2  <- -VC$weights*((bit1.b2b2*p11 - (c.copula.be2*d.n2)^2 )/p11 -
                            (bit1.b2b2*cp11 + (c.copula.be2*d.n2)^2 )/cp11 )

  d2l.be1.be2  <- -VC$weights*((bit1.b1b2*p11-(c.copula.be1*d.n1*c.copula.be2*d.n2))/p11 -
                            (bit1.b1b2*cp11+(c.copula.be1*d.n1*c.copula.be2*d.n2))/cp11)

  d2l.be1.rho  <- -VC$weights*((bit1.b1th*p11-(c.copula.be1*d.n1*c.copula.theta))/p11 -
                            (bit1.b1th*cp11+(c.copula.be1*d.n1*c.copula.theta))/cp11)

  d2l.be2.rho  <- -VC$weights*((bit1.b2th*p11-(c.copula.be2*d.n2*c.copula.theta))/p11 -
                            (bit1.b2th*cp11+(c.copula.be2*d.n2*c.copula.theta))/cp11)

  d2l.rho.rho  <- -VC$weights*((bit1.th2*p11-c.copula.theta^2)/p11 -
                            (bit1.th2*cp11+c.copula.theta^2)/cp11)                          
}

       
if( is.null(VC$X3) ){

  be1.be1 <- crossprod(VC$X1*c(d2l.be1.be1),VC$X1)
  be2.be2 <- crossprod(VC$X2*c(d2l.be2.be2),VC$X2)
  be1.be2 <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)
  be1.rho <- t(t(rowSums(t(VC$X1*c(d2l.be1.rho)))))
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
  be1.be2 <- crossprod(VC$X1*c(d2l.be1.be2),VC$X2)
  be1.rho <- crossprod(VC$X1*c(d2l.be1.rho),VC$X3)
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

if( ( VC$l.sp1==0 && VC$l.sp2==0 && VC$l.sp3==0 ) || VC$fp==TRUE) ps <- list(S.h = 0, S.h1 = 0, S.h2 = 0) else ps <- pen(params, qu.mag, sp, VC)


if(VC$extra.regI == "pC" && VC$hess==FALSE) H <- regH(H, type = 1)
  
         S.res <- res
         res   <- S.res + ps$S.h1
         G     <- G + ps$S.h2
         H     <- H + ps$S.h  
        
if(VC$extra.regI == "sED") H <- regH(H, type = 2)  
   

           

         list(value=res, gradient=G, hessian=H, S.h=ps$S.h, l=S.res, l.par=l.par, ps = ps,
              p11=p11, cp11=cp11, eta1=eta1, eta2=eta2, etad=etad,
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.drho=dl.drho,
              #d2l.be1.be1=d2l.be1.be1, d2l.be2.be2=d2l.be2.be2, 
              #d2l.be1.be2=d2l.be1.be2, d2l.be1.rho=d2l.be1.rho,
              #d2l.be2.rho=d2l.be2.rho, d2l.rho.rho=d2l.rho.rho, 
              BivD=VC$BivD, p1=p1, p2=p2)      

}




     























