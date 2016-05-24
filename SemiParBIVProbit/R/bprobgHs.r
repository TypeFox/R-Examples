bprobgHs <- function(params, respvec, VC, sp = NULL, qu.mag = NULL, AT = FALSE){

 
  epsilon <- 0.0000001 
  max.p   <- 0.9999999

  eta1 <- VC$X1%*%params[1:VC$X1.d2]
  eta2 <- VC$X2%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]
  etad <- NULL 
  
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
    
p11 <- BiCDF(p1, p2, VC$nC, teta)

########################################################################################################

  p10 <- pmax(p1 - p11, epsilon)
  p01 <- pmax(p2 - p11, epsilon)
  p00 <- pmax(1 - p11 - p10 - p01, epsilon)

  l.par <- VC$weights*( respvec$y1.y2*log(p11) + respvec$y1.cy2*log(p10) + respvec$cy1.y2*log(p01) + respvec$cy1.cy2*log(p00) )

dH <- copgHs(p1,p2,eta1=NULL,eta2=NULL,teta,teta.st,VC$BivD)

c.copula.be1   <- dH$c.copula.be1
c.copula.be2   <- dH$c.copula.be2
c.copula.theta <- dH$c.copula.theta 
  
c.copula2.be1 <- dH$c.copula2.be1   
c.copula2.be2 <- dH$c.copula2.be2 


bit1.b1b1 <- c.copula2.be1*d.n1^2 + c.copula.be1*der2p1.dereta12                                                                                                                                
bit2.b1b1 <- -c.copula2.be1*d.n1^2  + (1-c.copula.be1)*der2p1.dereta12
bit3.b1b1 <- -bit1.b1b1
bit4.b1b1 <- -bit2.b1b1


bit1.b2b2 <- c.copula2.be2*d.n2^2 + c.copula.be2*der2p2.dereta22                                                                                                                             
bit2.b2b2 <- -bit1.b2b2
bit3.b2b2 <- -c.copula2.be2*d.n2^2 + (1-c.copula.be2)*der2p2.dereta22
bit4.b2b2 <- -bit3.b2b2

c.copula2.be1be2 <- dH$c.copula2.be1be2
bit1.b1b2 <- c.copula2.be1be2*d.n1*d.n2
bit2.b1b2 <- -bit1.b1b2
bit3.b1b2 <- -bit1.b1b2
bit4.b1b2 <- bit1.b1b2

c.copula2.be1th <- dH$c.copula2.be1th 
bit1.b1th <- c.copula2.be1th*d.n1
bit2.b1th <- -bit1.b1th 
bit3.b1th <- -bit1.b1th 
bit4.b1th <- bit1.b1th 

c.copula2.be2th <- dH$c.copula2.be2th
bit1.b2th <- c.copula2.be2th*d.n2
bit2.b2th <- -bit1.b2th 
bit3.b2th <- -bit1.b2th 
bit4.b2th <- bit1.b2th 

if(AT==TRUE) bit1.th2 <- dH$bit1.th2ATE else bit1.th2 <- dH$bit1.th2

bit2.th2 <- -bit1.th2 
bit3.th2 <- -bit1.th2 
bit4.th2 <- bit1.th2 

  dl.dbe1 <-  VC$weights*d.n1*( (respvec$y1.y2*c.copula.be1/p11)  +
                      (respvec$y1.cy2*(1-c.copula.be1)/p10) +
                      (respvec$cy1.y2*c.copula.be1/(-p01)) +
                      (respvec$cy1.cy2*(c.copula.be1-1)/p00) )
                                
  dl.dbe2 <-  VC$weights*d.n2*( (respvec$y1.y2*c.copula.be2/p11)  +
                            (respvec$y1.cy2*c.copula.be2/(-p10)) +
                                (respvec$cy1.y2*(1-c.copula.be2)/(p01)) +
                                (respvec$cy1.cy2*(c.copula.be2-1)/p00) )

  dl.drho <- VC$weights*( respvec$y1.y2*c.copula.theta/p11+respvec$y1.cy2*(-c.copula.theta)/p10 + 
                       respvec$cy1.y2*(-c.copula.theta)/p01+respvec$cy1.cy2*c.copula.theta/p00 ) 


add.b  <- 1
if(AT==TRUE){
    if(VC$BivD %in% c("N") )                                    add.b <- 1/cosh(teta.st)^2
    if(VC$BivD %in% c("C0", "C180","J0", "J180","G0", "G180") ) add.b <-  exp(teta.st)     
    if(VC$BivD %in% c("C90","C270","J90","J270","G90","G270") ) add.b <- -exp(teta.st)        
}



if(VC$hess==TRUE){

  d2l.be1.be1  <- -VC$weights*(respvec$y1.y2*(bit1.b1b1*p11-(c.copula.be1*d.n1)^2)/p11^2+
                              respvec$y1.cy2*(bit2.b1b1*p10-((1-c.copula.be1)*d.n1)^2)/p10^2+
                              respvec$cy1.y2*(bit3.b1b1*p01-(-c.copula.be1*d.n1)^2)/p01^2+
                              respvec$cy1.cy2*(bit4.b1b1*p00-((c.copula.be1-1)*d.n1)^2)/p00^2 )

  d2l.be2.be2  <- -VC$weights*(respvec$y1.y2*(bit1.b2b2*p11 - (c.copula.be2*d.n2)^2 )/p11^2+
                              respvec$y1.cy2*(bit2.b2b2*p10-(-c.copula.be2*d.n2)^2)/p10^2+
                              respvec$cy1.y2*(bit3.b2b2*p01- ((1-c.copula.be2)*d.n2)^2 )/p01^2+
                              respvec$cy1.cy2*(bit4.b2b2*p00-((c.copula.be2-1)*d.n2)^2)/p00^2 )

  d2l.be1.be2  <- -VC$weights*(respvec$y1.y2*(bit1.b1b2*p11-(c.copula.be1*d.n1*c.copula.be2*d.n2))/p11^2+
                              respvec$y1.cy2*(bit2.b1b2*p10-((1-c.copula.be1)*d.n1)*(-c.copula.be2*d.n2))/p10^2+
                              respvec$cy1.y2*(bit3.b1b2*p01-(-c.copula.be1*d.n1*(1-c.copula.be2)*d.n2))/p01^2+
                              respvec$cy1.cy2*(bit4.b1b2*p00-((c.copula.be1-1)*d.n1*((c.copula.be2-1)*d.n2)))/p00^2 )

  d2l.be1.rho  <- -VC$weights*(respvec$y1.y2*(bit1.b1th*p11-(c.copula.be1*d.n1*c.copula.theta))/p11^2+
                              respvec$y1.cy2*(bit2.b1th*p10-((1-c.copula.be1)*d.n1)*(-c.copula.theta))/p10^2+
                             respvec$cy1.y2*(bit3.b1th*p01-(-c.copula.be1*d.n1)*(-c.copula.theta))/p01^2+
                              respvec$cy1.cy2*(bit4.b1th*p00-((c.copula.be1-1)*d.n1)*c.copula.theta)/p00^2 )/add.b

  d2l.be2.rho  <- -VC$weights*(respvec$y1.y2*(bit1.b2th*p11-(c.copula.be2*d.n2*c.copula.theta))/p11^2+
                              respvec$y1.cy2*(bit2.b2th*p10-(-c.copula.be2*d.n2)*(-c.copula.theta))/p10^2+
                             respvec$cy1.y2*(bit3.b2th*p01-((1-c.copula.be2)*d.n2)*(-c.copula.theta))/p01^2+
                              respvec$cy1.cy2*(bit4.b2th*p00-((c.copula.be2-1)*d.n2)*c.copula.theta)/p00^2 )/add.b
                              
  d2l.rho.rho  <- (-VC$weights*(   respvec$y1.y2*(bit1.th2*p11-( c.copula.theta/add.b)^2)/p11^2+
                               respvec$y1.cy2*(bit2.th2*p10-(-c.copula.theta/add.b)^2)/p10^2+
                               respvec$cy1.y2*(bit3.th2*p01-(-c.copula.theta/add.b)^2)/p01^2+
                              respvec$cy1.cy2*(bit4.th2*p00-( c.copula.theta/add.b)^2)/p00^2 ) )                            
}     


if(VC$hess==FALSE && VC$end==0){

  d2l.be1.be1  <- -VC$weights*(  (bit1.b1b1*p11-(c.copula.be1*d.n1)^2)/p11+
                              (bit2.b1b1*p10-((1-c.copula.be1)*d.n1)^2)/p10+
                              (bit3.b1b1*p01-(-c.copula.be1*d.n1)^2)/p01+
                              (bit4.b1b1*p00-((c.copula.be1-1)*d.n1)^2)/p00 )

  d2l.be2.be2  <- -VC$weights*((bit1.b2b2*p11 - (c.copula.be2*d.n2)^2 )/p11+
                              (bit2.b2b2*p10-(-c.copula.be2*d.n2)^2)/p10+
                              (bit3.b2b2*p01- ((1-c.copula.be2)*d.n2)^2 )/p01+
                              (bit4.b2b2*p00-((c.copula.be2-1)*d.n2)^2)/p00 )

  d2l.be1.be2  <- -VC$weights*((bit1.b1b2*p11-(c.copula.be1*d.n1*c.copula.be2*d.n2))/p11+
                              (bit2.b1b2*p10-((1-c.copula.be1)*d.n1)*(-c.copula.be2*d.n2))/p10+
                              (bit3.b1b2*p01-(-c.copula.be1*d.n1*(1-c.copula.be2)*d.n2))/p01+
                              (bit4.b1b2*p00-((c.copula.be1-1)*d.n1*((c.copula.be2-1)*d.n2)))/p00 )

  d2l.be1.rho  <- -VC$weights*((bit1.b1th*p11-(c.copula.be1*d.n1*c.copula.theta))/p11+
                              (bit2.b1th*p10-((1-c.copula.be1)*d.n1)*(-c.copula.theta))/p10+
                              (bit3.b1th*p01-(-c.copula.be1*d.n1)*(-c.copula.theta))/p01+
                              (bit4.b1th*p00-((c.copula.be1-1)*d.n1)*c.copula.theta)/p00 )/add.b

  d2l.be2.rho  <- -VC$weights*(  (bit1.b2th*p11-(c.copula.be2*d.n2*c.copula.theta))/p11+
                              (bit2.b2th*p10-(-c.copula.be2*d.n2)*(-c.copula.theta))/p10+
                              (bit3.b2th*p01-((1-c.copula.be2)*d.n2)*(-c.copula.theta))/p01+
                              (bit4.b2th*p00-((c.copula.be2-1)*d.n2)*c.copula.theta)/p00 )/add.b

  d2l.rho.rho  <- -VC$weights*(  (bit1.th2*p11-( c.copula.theta/add.b)^2)/p11+
                              (bit2.th2*p10-(-c.copula.theta/add.b)^2)/p10+
                              (bit3.th2*p01-(-c.copula.theta/add.b)^2)/p01+
                              (bit4.th2*p00-( c.copula.theta/add.b)^2)/p00 )
}  

if(VC$hess==FALSE && (VC$end==1 || VC$end==2)){

if(VC$end==1){

fi <- p11/p1
se <- p10/p1
th <- p01/(1-p1)
fo <- p00/(1-p1)

resp1 <- respvec$y1
resp2 <- respvec$y1
resp3 <- 1 - respvec$y1 
resp4 <- 1 - respvec$y1 


}

if(VC$end==2){

fi <- p11/p2
se <- p10/(1-p2)
th <- p01/p2
fo <- p00/(1-p2)

resp1 <- respvec$y2
resp2 <- 1 - respvec$y2
resp3 <- respvec$y2 
resp4 <- 1 - respvec$y2

}

  d2l.be1.be1  <- -VC$weights*(  fi*(bit1.b1b1*p11-(c.copula.be1*d.n1)^2)/p11^2*resp1+
                              se*(bit2.b1b1*p10-((1-c.copula.be1)*d.n1)^2)/p10^2*resp2+
                              th*(bit3.b1b1*p01-(-c.copula.be1*d.n1)^2)/p01^2*resp3+
                              fo*(bit4.b1b1*p00-((c.copula.be1-1)*d.n1)^2)/p00^2*resp4 )

  d2l.be2.be2  <- -VC$weights*(  fi*(bit1.b2b2*p11 - (c.copula.be2*d.n2)^2 )/p11^2*resp1+
                              se*(bit2.b2b2*p10-(-c.copula.be2*d.n2)^2)/p10^2*resp2+
                              th*(bit3.b2b2*p01- ((1-c.copula.be2)*d.n2)^2 )/p01^2*resp3+
                              fo*(bit4.b2b2*p00-((c.copula.be2-1)*d.n2)^2)/p00^2*resp4 )

  d2l.be1.be2  <- -VC$weights*(  fi*(bit1.b1b2*p11-(c.copula.be1*d.n1*c.copula.be2*d.n2))/p11^2*resp1+
                              se*(bit2.b1b2*p10-((1-c.copula.be1)*d.n1)*(-c.copula.be2*d.n2))/p10^2*resp2+
                              th*(bit3.b1b2*p01-(-c.copula.be1*d.n1*(1-c.copula.be2)*d.n2))/p01^2*resp3+
                              fo*(bit4.b1b2*p00-((c.copula.be1-1)*d.n1*((c.copula.be2-1)*d.n2)))/p00^2*resp4 )

  d2l.be1.rho  <- -VC$weights*(  fi*(bit1.b1th*p11-(c.copula.be1*d.n1*c.copula.theta))/p11^2*resp1+
                              se*(bit2.b1th*p10-((1-c.copula.be1)*d.n1)*(-c.copula.theta))/p10^2*resp2+
                              th*(bit3.b1th*p01-(-c.copula.be1*d.n1)*(-c.copula.theta))/p01^2*resp3+
                              fo*(bit4.b1th*p00-((c.copula.be1-1)*d.n1)*c.copula.theta)/p00^2*resp4 )/add.b

  d2l.be2.rho  <- -VC$weights*(  fi*(bit1.b2th*p11-(c.copula.be2*d.n2*c.copula.theta))/p11^2*resp1+
                              se*(bit2.b2th*p10-(-c.copula.be2*d.n2)*(-c.copula.theta))/p10^2*resp2+
                              th*(bit3.b2th*p01-((1-c.copula.be2)*d.n2)*(-c.copula.theta))/p01^2*resp3+
                              fo*(bit4.b2th*p00-((c.copula.be2-1)*d.n2)*c.copula.theta)/p00^2*resp4 )/add.b

  d2l.rho.rho  <- -VC$weights*(  fi*(bit1.th2*p11-( c.copula.theta/add.b)^2)/p11^2*resp1+
                              se*(bit2.th2*p10-(-c.copula.theta/add.b)^2)/p10^2*resp2+
                              th*(bit3.th2*p01-(-c.copula.theta/add.b)^2)/p01^2*resp3+
                              fo*(bit4.th2*p00-( c.copula.theta/add.b)^2)/p00^2*resp4 )
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


#nr <- 20000
#nc <- 1000
#X <- matrix(runif(nr*nc), nr, nc)
#system.time(crossprod(X))
#system.time(crossprod(X[1:(nr/2),]) + crossprod(X[(nr/2 + 1):nr,]))


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
              p11=p11, p10=p10, p01=p01, p00=p00, eta1=eta1, eta2=eta2, etad=etad, 
              dl.dbe1=dl.dbe1, dl.dbe2=dl.dbe2, dl.drho=dl.drho,
              #d2l.be1.be1=d2l.be1.be1, d2l.be2.be2=d2l.be2.be2, 
              #d2l.be1.be2=d2l.be1.be2, d2l.be1.rho=d2l.be1.rho,
              #d2l.be2.rho=d2l.be2.rho, d2l.rho.rho=d2l.rho.rho, 
              BivD=VC$BivD, p1=p1, p2=p2, theta.star = teta.st)      

}




     























