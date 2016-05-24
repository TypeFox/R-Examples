bprobgHsContUniv3 <- function(params, respvec, VC, sp = NULL, qu.mag = NULL, AT = FALSE){

epsilon <- 0.0000001 # 0.9999999 0.0001 # sqrt(.Machine$double.eps)

weights <- VC$weights

if(VC$Cont == "NO"){

X2 <- VC$X2
X3 <- VC$X3
X4 <- VC$X4
y2 <- respvec$y2

if(VC$ccss == "yes"){weights <- weights[VC$inde]} # ; y2 <- y2[VC$inde]; X2 <- X2[VC$inde,]; X3 <- X3[VC$inde,]; X4 <- X4[VC$inde,]} 

eta2 <- X2%*%params[1:VC$X2.d2]
eta2 <- eta.tr(eta2, VC$margins[2])
    
if(is.null(VC$X3))  sigma2.st <- params[(VC$X2.d2 + 1)] 
if(!is.null(VC$X3)) sigma2.st <- X3%*%params[(VC$X2.d2+1):(VC$X2.d2+VC$X3.d2)]

if(is.null(VC$X4))  nu.st <- params[(VC$X2.d2 + 2)] 
if(!is.null(VC$X4)) nu.st <- X4%*%params[(VC$X2.d2+VC$X3.d2+1):(VC$X2.d2+VC$X3.d2+VC$X4.d2)]

sstr1 <- esp.tr(sigma2.st, VC$margins[2])  
sigma2.st <- sstr1$vrb.st 
sigma2    <- sstr1$vrb 
    
sstr1 <- esp.tr(nu.st, VC$margins[2])  
nu.st <- sstr1$vrb.st 
nu    <- sstr1$vrb 

dHs  <- distrHs(y2, eta2, sigma2, sigma2.st, nu, nu.st, margin2=VC$margins[2], naive = TRUE)

}




if(VC$Cont == "YES" && respvec$univ == 2){ # interested in first equation

eta2 <- VC$X1%*%params[1:VC$X1.d2] # this is eta1
eta2 <- eta.tr(eta2, VC$margins[1])
    
if(is.null(VC$X3))  sigma2.st <- params[(VC$X1.d2 + 1)] 
if(!is.null(VC$X3)) sigma2.st <- VC$X3%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X3.d2)]

if(is.null(VC$X4))  nu.st <- params[(VC$X1.d2 + 2)] 
if(!is.null(VC$X4)) nu.st <- VC$X5%*%params[(VC$X1.d2+VC$X3.d2+1):(VC$X1.d2+VC$X3.d2+VC$X5.d2)]

sstr1 <- esp.tr(sigma2.st, VC$margins[1])  
sigma2.st <- sstr1$vrb.st 
sigma2    <- sstr1$vrb 
    
sstr1 <- esp.tr(nu.st, VC$margins[1])  
nu.st <- sstr1$vrb.st 
nu    <- sstr1$vrb           
 
dHs  <- distrHs(respvec$y1, eta2, sigma2, sigma2.st, nu, nu.st, margin2=VC$margins[1], naive = TRUE)

}


if(VC$Cont == "YES" && respvec$univ == 3){ # interested in second equation

eta2 <- VC$X2%*%params[1:VC$X2.d2]
eta2 <- eta.tr(eta2, VC$margins[2])
    

if(is.null(VC$X3))  sigma2.st <- params[(VC$X2.d2 + 1)] 
if(!is.null(VC$X3)) sigma2.st <- VC$X4%*%params[(VC$X2.d2+1):(VC$X2.d2+VC$X4.d2)]

if(is.null(VC$X4))  nu.st <- params[(VC$X2.d2 + 2)] 


if(!is.null(VC$X4)){

if( VC$margins[1] %in% c("DAGUM","SM")  ) nu.st <- VC$X6%*%params[(VC$X2.d2+VC$X4.d2+1):(VC$X2.d2+VC$X4.d2+VC$X6.d2)] 
if(!(VC$margins[1] %in% c("DAGUM","SM"))) nu.st <- VC$X5%*%params[(VC$X2.d2+VC$X4.d2+1):(VC$X2.d2+VC$X4.d2+VC$X5.d2)] 

}

sstr1 <- esp.tr(sigma2.st, VC$margins[2])  
sigma2.st <- sstr1$vrb.st 
sigma2    <- sstr1$vrb 
    
sstr1 <- esp.tr(nu.st, VC$margins[2])  
nu.st <- sstr1$vrb.st 
nu    <- sstr1$vrb          
 
dHs  <- distrHs(respvec$y2, eta2, sigma2, sigma2.st, nu, nu.st, margin2=VC$margins[2], naive = TRUE)

}




########################################################################################################

pdf2                         <- dHs$pdf2
derpdf2.dereta2              <- dHs$derpdf2.dereta2 
derpdf2.dersigma2.st         <- dHs$derpdf2.dersigma2.st  
der2pdf2.dereta2             <- dHs$der2pdf2.dereta2
der2pdf2.dersigma2.st2       <- dHs$der2pdf2.dersigma2.st2         
der2pdf2.dereta2dersigma2.st <- dHs$der2pdf2.dereta2dersigma2.st  


der2pdf2.dereta2dernu.st    <- dHs$der2pdf2.dereta2dernu.st   
der2pdf2.dersigma2.stdernu.st <- dHs$der2pdf2.sigma2.st2dernu.st
derpdf2.dernu.st            <- dHs$derpdf2.dernu.st           
der2pdf2.dernu.st2          <- dHs$der2pdf2.dernu.st2         
    
########################################################################################################

l.par <- weights*log(pdf2)
res   <- -sum(l.par)
  
########################################################################################################
 
dl.dbe       <- -weights*( derpdf2.dereta2/pdf2 )
dl.dsigma.st <- -weights*( derpdf2.dersigma2.st/pdf2 )
dl.dnu.st    <- -weights*( derpdf2.dernu.st/pdf2 )   
          
          
d2l.be.be        <- -weights*( (der2pdf2.dereta2*pdf2 - (derpdf2.dereta2)^2)/pdf2^2 )
d2l.sigma.sigma  <- -weights*( (der2pdf2.dersigma2.st2*pdf2-(derpdf2.dersigma2.st)^2)/pdf2^2 )
d2l.be.sigma     <- -weights*( (der2pdf2.dereta2dersigma2.st*pdf2 - derpdf2.dereta2*derpdf2.dersigma2.st)/pdf2^2 )
d2l.be.nu        <- -weights*( (der2pdf2.dereta2dernu.st*pdf2 - derpdf2.dereta2*derpdf2.dernu.st)/pdf2^2 )
     
d2l.sigma.nu     <- -weights*( (der2pdf2.dersigma2.stdernu.st*pdf2-(derpdf2.dersigma2.st*derpdf2.dernu.st))/pdf2^2 )
  d2l.nu.nu      <- -weights*( (der2pdf2.dernu.st2*pdf2-(derpdf2.dernu.st)^2)/pdf2^2 )
     
     
     
########################################################################################################
     

if(VC$Cont == "NO"){

if( is.null(VC$X3) ){

  G   <- c( colSums( c(dl.dbe)*X2 ) ,
            sum( dl.dsigma.st ), sum( dl.dnu.st ) )
                
  be.be    <- crossprod(X2*c(d2l.be.be),X2)
  be.sigma <- t(t(rowSums(t(X2*c(d2l.be.sigma))))) 
  be.nu    <- t(t(rowSums(t(X2*c(d2l.be.nu))))) 

  H <- rbind( cbind( be.be      , be.sigma,             be.nu             ), 
              cbind( t(be.sigma), sum(d2l.sigma.sigma), sum(d2l.sigma.nu) ),
              cbind( t(be.nu),    sum(d2l.sigma.nu),    sum(d2l.nu.nu)    )  )     }


if( !is.null(VC$X3) ){

  G   <- c( colSums( c(dl.dbe)*X2      ),
            colSums( c(dl.dsigma.st)*X3),
            colSums( c(dl.dnu.st)*X4   )  )
                
  be.be    <- crossprod(X2*c(d2l.be.be),X2)
  be.sigma <- crossprod(X2*c(d2l.be.sigma),X3)
  be.nu    <- crossprod(X2*c(d2l.be.nu),X4)
  
  si.sigma <- crossprod(X3*c(d2l.sigma.sigma),X3)  
  si.nu    <- crossprod(X4*c(d2l.nu.nu),X4)    
  
  sa.nu    <- crossprod(X3*c(d2l.sigma.nu),X4)    
  
  H <- rbind( cbind( be.be      , be.sigma , be.nu ), 
              cbind( t(be.sigma), si.sigma,  sa.nu ),
              cbind( t(be.nu),    t(sa.nu),  si.nu ) )       }


if( (VC$l.sp2==0 && VC$l.sp3==0 && VC$l.sp4==0) || VC$fp==TRUE) ps <- list(S.h = 0, S.h1 = 0, S.h2 = 0) else ps <- pen(params, qu.mag, sp, VC, respvec$univ)

}




if(VC$Cont == "YES" && respvec$univ == 2){

if( is.null(VC$X3) ){

  G   <- c( colSums( c(dl.dbe)*VC$X1 ) ,
            sum( dl.dsigma.st ), sum( dl.dnu.st ) )
                
  be.be    <- crossprod(VC$X1*c(d2l.be.be),VC$X1)
  be.sigma <- t(t(rowSums(t(VC$X1*c(d2l.be.sigma))))) 
  be.nu    <- t(t(rowSums(t(VC$X1*c(d2l.be.nu))))) 

  H <- rbind( cbind( be.be      , be.sigma,             be.nu             ), 
              cbind( t(be.sigma), sum(d2l.sigma.sigma), sum(d2l.sigma.nu) ),
              cbind( t(be.nu),    sum(d2l.sigma.nu),    sum(d2l.nu.nu)    )  )     }


if( !is.null(VC$X3) ){

  G   <- c( colSums( c(dl.dbe)*VC$X1      ),
            colSums( c(dl.dsigma.st)*VC$X3),
            colSums( c(dl.dnu.st)*VC$X5   )  )
                
  be.be    <- crossprod(VC$X1*c(d2l.be.be),VC$X1)
  be.sigma <- crossprod(VC$X1*c(d2l.be.sigma),VC$X3)
  be.nu    <- crossprod(VC$X1*c(d2l.be.nu),VC$X5)
  
  si.sigma <- crossprod(VC$X3*c(d2l.sigma.sigma),VC$X3)  
  si.nu    <- crossprod(VC$X5*c(d2l.nu.nu),VC$X5)    
  
  sa.nu    <- crossprod(VC$X3*c(d2l.sigma.nu),VC$X5)    
  
  H <- rbind( cbind( be.be      , be.sigma , be.nu ), 
              cbind( t(be.sigma), si.sigma,  sa.nu ),
              cbind( t(be.nu),    t(sa.nu),  si.nu ) )       }



if( (VC$l.sp1==0 && VC$l.sp3==0 && VC$l.sp5==0) || VC$fp==TRUE) ps <- list(S.h = 0, S.h1 = 0, S.h2 = 0) else ps <- pen(params, qu.mag, sp, VC, respvec$univ)

}






if(VC$Cont == "YES" && respvec$univ == 3){



if( is.null(VC$X3) ){

  G   <- c( colSums( c(dl.dbe)*VC$X2 ) ,
            sum( dl.dsigma.st ), sum( dl.dnu.st ) )
                
  be.be    <- crossprod(VC$X2*c(d2l.be.be),VC$X2)
  be.sigma <- t(t(rowSums(t(VC$X2*c(d2l.be.sigma))))) 
  be.nu    <- t(t(rowSums(t(VC$X2*c(d2l.be.nu))))) 

  H <- rbind( cbind( be.be      , be.sigma,             be.nu             ), 
              cbind( t(be.sigma), sum(d2l.sigma.sigma), sum(d2l.sigma.nu) ),
              cbind( t(be.nu),    sum(d2l.sigma.nu),    sum(d2l.nu.nu)    )  )  
              
              
}



if( !is.null(VC$X3) ){



if( VC$margins[1] %in% c("DAGUM","SM")  ){

  G   <- c( colSums( c(dl.dbe)*VC$X2      ),
            colSums( c(dl.dsigma.st)*VC$X4),
            colSums( c(dl.dnu.st)*VC$X6   )  )
                
  be.be    <- crossprod(VC$X2*c(d2l.be.be),VC$X2)
  be.sigma <- crossprod(VC$X2*c(d2l.be.sigma),VC$X4)
  be.nu    <- crossprod(VC$X2*c(d2l.be.nu),VC$X6)
  
  si.sigma <- crossprod(VC$X4*c(d2l.sigma.sigma),VC$X4)  
  si.nu    <- crossprod(VC$X6*c(d2l.nu.nu),VC$X6)    
  
  sa.nu    <- crossprod(VC$X4*c(d2l.sigma.nu),VC$X6)    
  
  H <- rbind( cbind( be.be      , be.sigma , be.nu ), 
              cbind( t(be.sigma), si.sigma,  sa.nu ),
              cbind( t(be.nu),    t(sa.nu),  si.nu ) )   }



if( !(VC$margins[1] %in% c("DAGUM","SM"))  ){

  G   <- c( colSums( c(dl.dbe)*VC$X2      ),
            colSums( c(dl.dsigma.st)*VC$X4),
            colSums( c(dl.dnu.st)*VC$X5   )  )
                
  be.be    <- crossprod(VC$X2*c(d2l.be.be),VC$X2)
  be.sigma <- crossprod(VC$X2*c(d2l.be.sigma),VC$X4)
  be.nu    <- crossprod(VC$X2*c(d2l.be.nu),VC$X5)
  
  si.sigma <- crossprod(VC$X4*c(d2l.sigma.sigma),VC$X4)  
  si.nu    <- crossprod(VC$X5*c(d2l.nu.nu),VC$X5)    
  
  sa.nu    <- crossprod(VC$X4*c(d2l.sigma.nu),VC$X5)    
  
  H <- rbind( cbind( be.be      , be.sigma , be.nu ), 
              cbind( t(be.sigma), si.sigma,  sa.nu ),
              cbind( t(be.nu),    t(sa.nu),  si.nu ) )   }


}


if( VC$margins[1] %in% c("DAGUM","SM")  ) { if( (VC$l.sp2==0 && VC$l.sp4==0 && VC$l.sp6==0) || VC$fp==TRUE) ps <- list(S.h = 0, S.h1 = 0, S.h2 = 0) else ps <- pen(params, qu.mag, sp, VC, respvec$univ) }
if( !(VC$margins[1] %in% c("DAGUM","SM"))  ) { if( (VC$l.sp2==0 && VC$l.sp4==0 && VC$l.sp5==0) || VC$fp==TRUE) ps <- list(S.h = 0, S.h1 = 0, S.h2 = 0) else ps <- pen(params, qu.mag, sp, VC, respvec$univ) }


}



if(VC$extra.regI == "pC") H <- regH(H, type = 1)
  
         S.res <- res
         res   <- S.res + ps$S.h1
         G     <- G + ps$S.h2
         H     <- H + ps$S.h  
        
if(VC$extra.regI == "sED") H <- regH(H, type = 2) 
  

list(value=res, gradient=G, hessian=H, S.h=ps$S.h, l=S.res, l.par=l.par, 
     ps = ps, sigma2.st = sigma2.st, nu.st = nu.st,
     BivD=VC$BivD, eta2 = eta2, sigma2 = sigma2, nu = nu)      

}
