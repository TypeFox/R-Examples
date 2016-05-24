ghss <- function(params, dat, VC, qu.mag=NULL, sp=NULL){

  etatheta <- etasqv <- etak <- NULL

  X1 <- as.matrix(VC$X1)
  X2 <- as.matrix(VC$X2)
  
  if(!is.null(VC$X3)){
  X3 <- as.matrix(VC$X3)
  X4 <- as.matrix(VC$X4)
  }


  eta1 <- X1%*%params[1:VC$X1.d2]  
  eta2 <- X2%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]   
  
  if(is.null(VC$X4))  teta.st <- params[(VC$X1.d2+VC$X2.d2+2)]
  if(!is.null(VC$X4)) teta.st <- etatheta <- X4%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2)]

  eps <- sqrt(.Machine$double.eps)

  if(VC$BivD %in% c("N","FGM","AMH") ){ 
                         teta <- tanh(teta.st)
		         teta <- ifelse(teta < -0.9999999, -0.9999999, teta)
                         teta <- ifelse(teta > 0.9999999 ,  0.9999999, teta)
  }
 
  
  if( !(VC$BivD %in% c("F","N","FGM","AMH")) ) {
  
    teta.st <- ifelse( teta.st > 709, 709, teta.st )  
    teta.st <- ifelse( teta.st < -20, -20, teta.st )  
  
  }
  
  if(VC$BivD=="F")                                      teta <- teta.st + eps
  if(VC$BivD %in% c("C0", "C90", "C180", "C270") )      teta <- exp(teta.st) + eps
  if(VC$BivD %in% c("J0", "J90", "J180", "J270") )      teta <- exp(teta.st) + 1 + eps 
  if(VC$BivD %in% c("G0", "G90", "G180", "G270") )      teta <- exp(teta.st) + 1 



  i1 <- dat[,1]
  i0 <- 1-i1
  ind <- i1==0 

  F1 <- pnorm(-eta1)				
  F1 <- ifelse(F1>0.00000001,F1,0.00000001)
  F1 <- ifelse(F1<0.99999999,F1,0.99999999)
  ph <- dnorm(-eta1)


if(VC$margins[2]=="N"){


    if(is.null(VC$X3))  sqv.st <- params[(VC$X1.d2+VC$X2.d2+1)]
    if(!is.null(VC$X3)) sqv.st <- etasqv <- X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]

    sqv.st <- ifelse( sqv.st > 709, 709, sqv.st )  
    sqv.st <- ifelse( sqv.st < -20, -20, sqv.st )  

    sqv <- exp(sqv.st)


    if(VC$BivD=="N"){ 
 
      a  <- sqrt(1-teta^2)
      e2 <- dat[,2]-eta2  
      A  <- (eta1+(teta/sqv)*e2)/a
      l1 <- ifelse(eta1<37.5,dnorm(-eta1)/pnorm(-eta1),eta1)
      l2 <- ifelse(A>-37.5,dnorm(A)/pnorm(A),-A)
      ec <- exp(2*teta.st) 

      PA    <- -(pnorm(A)*dnorm(A)*A+dnorm(A)^2)/pnorm(A)^2 
      Peta1 <- -(pnorm(-eta1)*dnorm(-eta1)*(-eta1)+dnorm(-eta1)^2)/pnorm(-eta1)^2  


      l.par <- VC$weights*(i0*log(pnorm(-eta1)) + i1*(log(pnorm(A))-1/2*log(2*pi)-log(sqv)-1/2*(e2/sqv)^2)) 
 

      dl.1 <- -VC$weights*( i0*l1 - i1*(l2/a) ) 
      dl.2 <- VC$weights*i1*( e2/sqv^2 - l2*( teta/(sqv*a) )   )
      dl.3 <- VC$weights*i1*( -1 + (e2/sqv)^2 - l2*((teta*e2)/(a*sqv)) ) 
      dl.4 <- VC$weights*i1*( l2*( ((2*ec*e2)/((ec+1)*sqv) - (2*teta*e2*ec)/((ec+1)*sqv))/a  - 0.5*( (eta1+(teta*e2)/sqv)*( -4*teta*ec/(ec+1) + 4*teta^2*ec/(ec+1) ) )/a^3 )  )   


      d2l.11 <- -VC$weights*( i0*Peta1 + i1*(PA/a^2) )  
      d2l.12 <- -VC$weights*( -i1*PA*(teta/(sqv*a^2)) )   
      d2l.13 <- -VC$weights*( -i1*PA*((teta*e2)/(sqv*a^2)) ) 
      d2l.14 <- -VC$weights*( i1*(       (PA/a)*(  ( (2*ec*e2)/((ec+1)*sqv) - (2*(ec-1)*e2*ec)/((ec+1)^2*sqv) )/a      
               				- 0.5*(eta1+(teta*e2)/sqv)*((-4*teta*ec)/(ec+1)+(4*teta^2*ec)/(ec+1) )/a^3 
         				  )  
   				 -0.5*(l2/a^3)*( (-4*teta*ec)/(ec+1) + (4*teta^2*ec)/(ec+1) )           
                     )    )
      d2l.22 <- -VC$weights*( i1*( PA*(teta/(sqv*a))^2 - 1/sqv^2 ) ) 
      d2l.23 <- -VC$weights*( i1*( PA*(teta^2*e2)/(sqv^2*a^2) + (teta*l2)/(sqv*a) - (2*e2)/sqv^2 ) )
      d2l.24 <- -VC$weights*( i1*(     -PA*(    
    				( (  (2*e2*ec)/((ec+1)*sqv) - (2*teta*e2*ec)/((ec+1)*sqv)  )/a 
      				   -0.5/a^3*( ( eta1+(teta*e2)/sqv )*( -4*teta*ec/(ec+1) + 4*teta^2*ec/(ec+1) ))   )*(ec-1)  
     				)/((ec+1)*sqv*a)
   				- (2*l2*ec)/(sqv*a*(ec+1))      
  				 + (2*l2*teta*ec)/(sqv*a*(ec+1))  
  				 + 0.5*l2*teta/(sqv*a^3)*( -4*teta*ec/(ec+1) + 4*teta^2*ec/(ec+1)   ) 
                     )    )
      d2l.33 <- -VC$weights*( i1*(  -2*(e2/sqv)^2 + l2*((e2*teta)/(a*sqv)) + PA*((teta*e2)/(sqv*a))^2  )   )
      d2l.34 <- -VC$weights*( i1*(   -PA*(    (( (2*e2*ec)/((ec+1)*sqv) -( 2*teta*e2*ec)/((ec+1)*sqv) )/a 
         				 -0.5/a^3*( ( eta1+(teta*e2)/sqv )*( -4*teta*ec/(ec+1) + 4*teta^2*ec/(ec+1) )))*(ec-1)*e2  
      				)/((ec+1)*sqv*a)
   				- (2*l2*ec*e2)/(sqv*a*(ec+1))      
  				 + (2*l2*teta*e2*ec)/(sqv*a*(ec+1))  
  				 + 0.5*l2*teta*e2/(sqv*a^3)*( -4*teta*ec/(ec+1) + 4*teta^2*ec/(ec+1)   ) 
                     )    )
      d2l.44 <- -VC$weights*( i1*(   PA*(       (
        				 ( (2*e2*ec)/((ec+1)*sqv) -(2*teta*e2*ec)/((ec+1)*sqv) )/a 
       				 -0.5/a^3*( ( eta1+(teta*e2)/sqv )*( -4*teta*ec/(ec+1) + 4*teta^2*ec/(ec+1) )) 
          				   )^2
     				 )
 				+ l2*(
 		 1/a*( 4*ec*e2/((ec+1)*sqv) - 8*ec^2*e2/((ec+1)^2*sqv) + 8*teta*ec^2*e2/((ec+1)^2*sqv)  - 4*teta*ec*e2/((ec+1)*sqv)  ) 
  		-1/a^3*( (2*ec*e2/((ec+1)*sqv) - 2*ec*e2*teta/((ec+1)*sqv) )*( -4*teta*ec/(ec+1) + 4*teta^2*ec/(ec+1) )  ) 
  		+3/(4*a^5)*( ( eta1+(teta*e2)/sqv )*( -4*teta*ec/(ec+1) + 4*teta^2*ec/(ec+1) )^2 ) 
  		-0.5/a^3*( (eta1+(teta*e2)/sqv )*(  -8*ec^2/(ec+1)^2 + 32*teta*ec^2/(ec+1)^2 - 8*teta*ec/(ec+1) - 24*teta^2*ec^2/(ec+1)^2  + 8*teta^2*ec/(ec+1)    ) 
       		    )
     		 )
                 )    )

    }
    
    if(VC$BivD!="N"){


        e2 <- (dat[,2]-eta2)/sqv
        F2 <- pnorm(e2)
        F2 <- ifelse(F2>0.00000001,F2,0.00000001)
        F2 <- ifelse(F2<0.99999999,F2,0.99999999)
        f2 <- dnorm(e2)/sqv;
        f2 <- ifelse(f2>0.00000001,f2,0.00000001)


        bits <- bitsgHs(cop=VC$BivD,margin=VC$margins[2],i1=i1,F1=F1,F2=F2,f2=f2,eta1=eta1,ph=ph,teta=teta,e2=e2,sqv=sqv,ver=0)

        z <- bits$z; p <- bits$p; h <- bits$h; b <- bits$b; P <- bits$P; A <- bits$A; h14 <- bits$h14; E <- bits$E; B <- bits$B; h44 <- bits$h44 

        z[ind] <- p[ind] <- h[ind] <- b[ind] <- P[ind] <- A[ind] <- h14[ind]<- E[ind] <- B[ind] <- h44[ind] <- 0


        l.par <- VC$weights*(i0*log(F1) + i1*(log(1-z)+log(f2)))

        # Gradient:

        dl.1  <- VC$weights*( -i0*F1^(-1) + i1*p )*ph                  # dl.dbe1
        dl.2  <- VC$weights*i1*( h + e2/sqv )                          # dl.dbe2
        dl.3  <- VC$weights*i1*( h*e2*sqv + e2^2 - 1 )                 # dl.dsqv.st
        dl.4  <- VC$weights*i1*b                                       # dl.dteta.st


        # (Minus) Hessian:

        d2l.11 <- -VC$weights*ph*( -i0*(ph/F1-eta1)/F1 + i1*P )        # d2l.be1.be1
        d2l.12 <- -VC$weights*A                                        # d2l.be1.be2
        d2l.13 <- -VC$weights*A*e2*sqv                                 # d2l.be1.sqv.st
        d2l.14 <- -VC$weights*h14                                      # d2l.be1.teta.st

        d2l.22 <- -VC$weights*i1*( h*E - sqv^(-2) )                    # d2l.be2.be2
        d2l.23 <- -VC$weights*i1*( h*(E*e2*sqv-1) - 2*e2/sqv )         # d2l.be2.sqv.st
        d2l.24 <- -VC$weights*B                                        # d2l.be2.teta.st 

        d2l.33 <- -VC$weights*i1*e2*(h*(E*e2*sqv^2-sqv)-2*e2)          # d2l.sqv.st.sqv.st
        d2l.34 <- -VC$weights*e2*B*sqv                                 # d2l.sqv.st.teta.st

        d2l.44 <- -VC$weights*h44                                      # d2l.teta.st.teta.st

    }

}



if(VC$margins[2]=="G"){

    if(is.null(VC$X3))  k.st <- params[(VC$X1.d2+VC$X2.d2+1)]
    if(!is.null(VC$X3)) k.st <- etak <- X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]

    k.st <- ifelse( k.st > 709, 709, k.st )  
    k.st <- ifelse( k.st < -20, -20, k.st )  

    k <- exp(k.st)
    i2 <- dat[,2]
    F2 <- pgamma(i2,shape=k,rate=k*exp(-eta2));		
    F2 <- ifelse(F2>0.00000001,F2,0.00000001)
    F2 <- ifelse(F2<0.99999999,F2,0.99999999)
    f2 <- dgamma(i2,shape=k,rate=k*exp(-eta2));		
    f2[i2==0] <- 0
    f2 <- ifelse(f2>0.00000001,f2,0.00000001)

    ye <- i2*exp(-eta2)

#    G <- gamma(k)
    psi <- psigamma(k)
    psiprim <- psigamma(k,deriv=1)
#    fun1 <- function(t,k) { t^(k-1)*exp(-t)*log(t) }
#    Gkf  <- function(k, y) { sapply(y, function(y) integrate(fun1,y,Inf,k)$value ) } # fingers crossed here!
#    Gk   <- Gkf(k,k*ye)
#    fun2 <- function(t,k){t^(k-1)*exp(-t)*(log(t))^2}
#    G2kf <- function(k, y) { sapply(y, function(y) integrate(fun2,y,Inf,k)$value ) }
#    G2k  <- G2kf(k,k*ye)
#    dF2k <- i2*f2*k^(-1)+psi*(1-F2)-Gk/G
    df2kbyf2 <- -psi+log(k*i2)+1-eta2-ye;  
    df2kbyf2[i1==0]=0;
#    d2F2k <-i2*k^(-1)*f2*(2*df2kbyf2+ye-1-1/k) + (1-F2)*(psiprim-psi^2) + 2*psi*Gk*G^(-1)-G2k/G 
    fp <- function(s) pgamma(i2, shape=s, rate=s*exp(-eta2))
    eps <- 1e-06
    dF2k <- (fp(k+eps/2) - fp(k-eps/2))/eps
    d2F2k <- (fp(k+2*eps) - fp(k+eps) - fp(k+eps) + fp(k))/(eps*eps)


    bits <- bitsgHs(cop=VC$BivD,margin=VC$margins[2],i1=i1,F1=F1,F2=F2,f2=f2,eta1=eta1,ph=ph,teta=teta,ver=0)

    z <- bits$z; p <- bits$p; h <- bits$h; b <- bits$b; P <- bits$P; A <- bits$A; h14 <- bits$h14; E <- bits$E; B <- bits$B; h44 <- bits$h44 

    z[ind] <- p[ind] <- h[ind] <- b[ind] <- P[ind] <- A[ind] <- h14[ind] <- E[ind] <- B[ind] <- h44[ind] <- 0


    l.par <- VC$weights*(i0*log(F1) + i1*(log(1-z)+log(f2)))


    # Gradient:

    dl.1  <- VC$weights*( -i0*F1^(-1) + i1*p )*ph                   # dl.dbe1
    dl.2  <- VC$weights*i1*( i2*f2*h - k*(1-ye) )                   # dl.dbe2
    dl.3  <- VC$weights*i1*( -dF2k*h + df2kbyf2 )*k                 # dl.dk.st 
    dl.4  <- VC$weights*i1*b                                        # dl.dteta.st

 
    # (Minus) Hessian:

    d2l.11  <- -VC$weights*ph*( -i0*(ph/F1-eta1)/F1 + i1*P )        # d2l.be1.be1
    d2l.12  <- -VC$weights*i2*f2*A                                  # d2l.be1.be2
    d2l.13  <-  VC$weights*A*dF2k*k                                  # d2l.be1.k.st
    d2l.14  <- -VC$weights*h14                                      # be1.teta.st

    d2l.22  <- -VC$weights*i1*i2*( f2*h*(i2*f2*E+k*ye-k) - k*exp(-eta2) )    # d2l.be2.be2
    d2l.23  <- -VC$weights*i1*( -i2*f2*h*(dF2k*E-df2kbyf2) - 1 + ye )*k      # d2l.be2.k.st
    d2l.24  <- -VC$weights*i2*B*f2                                           # d2l.be2.teta.st 

    d2l.33  <- -VC$weights*i1*( -dF2k*h + df2kbyf2 - h*k*(d2F2k-dF2k^2*E) + 1 - k*psiprim )*k    # d2l.k.st.k.st
    d2l.34  <-  VC$weights*dF2k*B*k                                                                 # d2l.k.st.teta.st

    d2l.44  <- -VC$weights*h44                                       # d2l.teta.st.teta.st

} 




if( is.null(VC$X3) && is.null(VC$X4)  ){

  H11 <- crossprod(X1*c(d2l.11),X1)
  H12 <- crossprod(X1*c(d2l.12),X2) 
  H13 <- t(t(rowSums(t(X1*c(d2l.13)))))
  H14 <- t(t(rowSums(t(X1*c(d2l.14)))))

  H22 <- crossprod(X2*c(d2l.22),X2) 
  H23 <- t(t(rowSums(t(X2*c(d2l.23)))))
  H24 <- t(t(rowSums(t(X2*c(d2l.24)))))

  H <- rbind( cbind( H11    , H12    , H13  ,  H14 ), 
              cbind( t(H12) , H22    , H23  ,  H24 ),
              cbind( t(H13) , t(H23) , sum(d2l.33), sum(d2l.34) ) ,
              cbind( t(H14) , t(H24) , sum(d2l.34), sum(d2l.44) )
            ) 

 

  G   <- -c(colSums( c(dl.1)*X1 ) ,
            colSums( c(dl.2)*X2 )    ,
            sum( dl.3 ) ,  
            sum( dl.4 )   )
            
            
} else{


  H11 <- crossprod(X1*c(d2l.11),X1)
  H12 <- crossprod(X1*c(d2l.12),X2) 
  H13 <- crossprod(X1*c(d2l.13),X3) 
  H14 <- crossprod(X1*c(d2l.14),X4) 

  H22 <- crossprod(X2*c(d2l.22),X2) 
  H23 <- crossprod(X2*c(d2l.23),X3) 
  H24 <- crossprod(X2*c(d2l.24),X4) 
  
  H33 <- crossprod(X3*c(d2l.33),X3) 
  H34 <- crossprod(X3*c(d2l.34),X4) 
  H44 <- crossprod(X4*c(d2l.44),X4)   
  

  H <- rbind( cbind( H11    , H12    , H13  ,  H14 ), 
              cbind( t(H12) , H22    , H23  ,  H24 ),
              cbind( t(H13) , t(H23) , H33  ,  H34 ) ,
              cbind( t(H14) , t(H24) , t(H34), H44 )
            ) 

  G   <- -c(colSums( c(dl.1)*X1 ) ,
            colSums( c(dl.2)*X2 ) ,
            colSums( c(dl.3)*X3 ) ,  
            colSums( c(dl.4)*X4 )  )


}



 res <- -sum(l.par)

if( ( VC$l.sp1==0 && VC$l.sp2==0 && VC$l.sp3==0 && VC$l.sp4==0 && VC$l.sp5==0 ) || VC$fp==TRUE) ps <- list(S.h = 0, S.h1 = 0, S.h2 = 0) else ps <- pen(params, qu.mag, sp, VC)
   
   
  S.res <- res
  res <- S.res + ps$S.h1
  G   <- G + ps$S.h2
  H   <- H + ps$S.h  

  list(value=res, gradient=G, hessian=H, S.h=ps$S.h, S.h2 = ps$S.h2, l=S.res, eta1=eta1, eta2=eta2,   
       etatheta = etatheta, etasqv = etasqv, etak = etak,
       dl.dbe1=dl.1, dl.dbe2=dl.2, l.par=l.par,
       dl.dsqv.st=dl.3,
       dl.dcor.st=dl.4, 
       d2l.be1.be1=d2l.11, d2l.be1.be2=d2l.12, d2l.be2.be2=d2l.22,
       d2l.be1.sqv.st=d2l.13,
       d2l.be1.cor.st=d2l.14,
       d2l.be2.sqv.st=d2l.23, 
       d2l.be2.cor.st=d2l.24,
       d2l.sqv.st.sqv.st=d2l.33,
       d2l.sqv.st.cor.st=d2l.34,    
       d2l.cor.st.cor.st=d2l.44 )

}






