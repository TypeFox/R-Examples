bprobgHsContUniv <- function(params, respvec, VC, sp = NULL, qu.mag = NULL, AT = FALSE){

epsilon <- 0.0000001 # 0.9999999 0.0001 # sqrt(.Machine$double.eps)

weights <- VC$weights

#########################################################################

if(VC$Cont == "NO"){

X2 <- VC$X2
X3 <- VC$X3
y2 <- respvec$y2

if(VC$ccss == "yes"){weights <- weights[VC$inde]} #; y2 <- y2[VC$inde]; X2 <- X2[VC$inde,]; X3 <- X3[VC$inde,] } 


eta2 <- X2%*%params[1:VC$X2.d2]
eta2 <- eta.tr(eta2, VC$margins[2])

if(is.null(VC$X3))  sigma2.st <- params[(VC$X2.d2 + 1)] 
if(!is.null(VC$X3)) sigma2.st <- X3%*%params[(VC$X2.d2+1):(VC$X2.d2+VC$X3.d2)] 
                    
    sstr1 <- esp.tr(sigma2.st, VC$margins[2])  
    sigma2.st <- sstr1$vrb.st 
    sigma2    <- sstr1$vrb 


dHs  <- distrHs(y2, eta2, sigma2, sigma2.st, nu = 1, nu.st = 1, margin2=VC$margins[2], naive = TRUE)

}

#########################################################################
#########################################################################


if(VC$Cont == "YES" && respvec$univ == 2){ # interested in first equation


eta2 <- VC$X1%*%params[1:VC$X1.d2] # this is eta1 but changed to eta2 for convenience
eta2 <- eta.tr(eta2, VC$margins[1])

if(is.null(VC$X3))  sigma2.st <- params[(VC$X1.d2 + 1)] 
if(!is.null(VC$X3)) sigma2.st <- VC$X3%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X3.d2)]

    sstr1 <- esp.tr(sigma2.st, VC$margins[1])  
    sigma2.st <- sstr1$vrb.st 
    sigma2    <- sstr1$vrb 


dHs  <- distrHs(respvec$y1, eta2, sigma2, sigma2.st, nu = 1, nu.st = 1, margin2=VC$margins[1], naive = TRUE)

}


if(VC$Cont == "YES" && respvec$univ == 3){ # interested in second equation

eta2 <- VC$X2%*%params[1:VC$X2.d2]
eta2 <- eta.tr(eta2, VC$margins[2])

if(is.null(VC$X3))  sigma2.st <- params[(VC$X2.d2 + 1)] 
if(!is.null(VC$X3)) sigma2.st <- VC$X4%*%params[(VC$X2.d2+1):(VC$X2.d2+VC$X4.d2)]


    sstr1 <- esp.tr(sigma2.st, VC$margins[2])  
    sigma2.st <- sstr1$vrb.st 
    sigma2    <- sstr1$vrb 

dHs  <- distrHs(respvec$y2, eta2, sigma2, sigma2.st, nu = 1, nu.st = 1, margin2=VC$margins[2], naive = TRUE)


}

#########################################################################
#########################################################################


pdf2                         <- dHs$pdf2
derpdf2.dereta2              <- dHs$derpdf2.dereta2 
derpdf2.dersigma2.st         <- dHs$derpdf2.dersigma2.st  
der2pdf2.dereta2             <- dHs$der2pdf2.dereta2
der2pdf2.dersigma2.st2       <- dHs$der2pdf2.dersigma2.st2         
der2pdf2.dereta2dersigma2.st <- dHs$der2pdf2.dereta2dersigma2.st  
  
    
########################################################################################################

l.par <- weights*log(pdf2)
res   <- -sum(l.par)
  
########################################################################################################
 
dl.dbe       <- -weights*( derpdf2.dereta2/pdf2 )
dl.dsigma.st <- -weights*( derpdf2.dersigma2.st/pdf2 )
                     
d2l.be.be        <- -weights*( (der2pdf2.dereta2*pdf2 - (derpdf2.dereta2)^2)/pdf2^2 )
d2l.sigma.sigma  <- -weights*( (der2pdf2.dersigma2.st2*pdf2-(derpdf2.dersigma2.st)^2)/pdf2^2 )
d2l.be.sigma     <- -weights*( (der2pdf2.dereta2dersigma2.st*pdf2 - derpdf2.dereta2*derpdf2.dersigma2.st)/pdf2^2 )
     
########################################################################################################
     

if(VC$Cont == "NO"){


if( is.null(VC$X3) ){


#if(VC$ccss == "no") X2 <- VC$X2

  G   <- c( colSums( c(dl.dbe)*X2 ) ,
            sum( dl.dsigma.st ) )
                
  be.be    <- crossprod(X2*c(d2l.be.be),X2)
  be.sigma <- t(t(rowSums(t(X2*c(d2l.be.sigma))))) 

  H <- rbind( cbind( be.be      , be.sigma             ), 
              cbind( t(be.sigma), sum(d2l.sigma.sigma) )  ) 
}


if( !is.null(VC$X3) ){

#if(VC$ccss == "no"){ X2 <- VC$X2; X3 <- VC$X3} 


  G   <- c( colSums( c(dl.dbe)*X2        ) ,
            colSums( c(dl.dsigma.st)*X3  ) )
                
  be.be    <- crossprod(X2*c(d2l.be.be),X2)
  be.sigma <- crossprod(X2*c(d2l.be.sigma),X3)
  si.sigma <- crossprod(X3*c(d2l.sigma.sigma),X3)  

  H <- rbind( cbind( be.be      , be.sigma  ), 
              cbind( t(be.sigma), si.sigma  )  ) 
              
}

if( (VC$l.sp2==0 && VC$l.sp3==0 ) || VC$fp==TRUE) ps <- list(S.h = 0, S.h1 = 0, S.h2 = 0) else ps <- pen(params, qu.mag, sp, VC, respvec$univ)

}



if(VC$Cont == "YES" && respvec$univ == 2){


if( is.null(VC$X3) ){

  G   <- c( colSums( c(dl.dbe)*VC$X1 ) ,
            sum( dl.dsigma.st )
          )
                
  be.be    <- crossprod(VC$X1*c(d2l.be.be),VC$X1)
  be.sigma <- t(t(rowSums(t(VC$X1*c(d2l.be.sigma))))) 

  H <- rbind( cbind( be.be      , be.sigma             ), 
              cbind( t(be.sigma), sum(d2l.sigma.sigma) )  
            ) 
}


if( !is.null(VC$X3) ){

  G   <- c( colSums( c(dl.dbe)*VC$X1        ) ,
            colSums( c(dl.dsigma.st)*VC$X3  )
          )
                
  be.be    <- crossprod(VC$X1*c(d2l.be.be),VC$X1)
  be.sigma <- crossprod(VC$X1*c(d2l.be.sigma),VC$X3)
  si.sigma <- crossprod(VC$X3*c(d2l.sigma.sigma),VC$X3)  

  H <- rbind( cbind( be.be      , be.sigma  ), 
              cbind( t(be.sigma), si.sigma  )  
            )  
}



if( (VC$l.sp1==0 && VC$l.sp3==0 ) || VC$fp==TRUE) ps <- list(S.h = 0, S.h1 = 0, S.h2 = 0) else ps <- pen(params, qu.mag, sp, VC, respvec$univ)

}





if(VC$Cont == "YES" && respvec$univ == 3){


if( is.null(VC$X3) ){

  G   <- c( colSums( c(dl.dbe)*VC$X2 ) ,
            sum( dl.dsigma.st )
          )
                
  be.be    <- crossprod(VC$X2*c(d2l.be.be),VC$X2)
  be.sigma <- t(t(rowSums(t(VC$X2*c(d2l.be.sigma))))) 

  H <- rbind( cbind( be.be      , be.sigma             ), 
              cbind( t(be.sigma), sum(d2l.sigma.sigma) )  
            ) 
}


if( !is.null(VC$X3) ){

  G   <- c( colSums( c(dl.dbe)*VC$X2        ) ,
            colSums( c(dl.dsigma.st)*VC$X4  )
          )
                
  be.be    <- crossprod(VC$X2*c(d2l.be.be),VC$X2)
  be.sigma <- crossprod(VC$X2*c(d2l.be.sigma),VC$X4)
  si.sigma <- crossprod(VC$X4*c(d2l.sigma.sigma),VC$X4)  

  H <- rbind( cbind( be.be      , be.sigma  ), 
              cbind( t(be.sigma), si.sigma  )  
            )  
}



if( (VC$l.sp2==0 && VC$l.sp4==0 ) || VC$fp==TRUE) ps <- list(S.h = 0, S.h1 = 0, S.h2 = 0) else ps <- pen(params, qu.mag, sp, VC, respvec$univ)

}




if(VC$extra.regI == "pC") H <- regH(H, type = 1)
  
         S.res <- res
         res   <- S.res + ps$S.h1
         G     <- G + ps$S.h2
         H     <- H + ps$S.h  
        
if(VC$extra.regI == "sED") H <- regH(H, type = 2) 
  

list(value=res, gradient=G, hessian=H, S.h=ps$S.h, l=S.res, l.par=l.par, ps = ps, sigma2.st = sigma2.st,
     BivD=VC$BivD, eta2 = eta2, sigma2 = sigma2, nu = NULL)      


}


