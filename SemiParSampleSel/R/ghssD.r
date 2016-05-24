ghssD <- function(params, dat, VC, qu.mag=NULL, sp=NULL) {

    etatheta <- etasqv <- etanu <- NULL
  
    bd <- VC$bd
    
    X1 <- as.matrix(VC$X1)
    X2 <- as.matrix(VC$X2)
    
    if(!is.null(VC$X3)){
    
    X3 <- as.matrix(VC$X3)
    
    if(!is.null(VC$X4)) X4 <- as.matrix(VC$X4)
    if(!is.null(VC$X5)) X5 <- as.matrix(VC$X5)    
    
    }


  eta1 <- X1%*%params[1:VC$X1.d2]  
  eta2 <- X2%*%params[(VC$X1.d2+1):(VC$X1.d2+VC$X2.d2)]   
  
  
if(is.null(VC$X3)){  
  
  if(VC$margins[2] %in% c("P", "BI", "GEOM", "LG", "YULE"))                                                     theta.star <- params[(VC$X1.d2+VC$X2.d2+1)]
  if(VC$margins[2] %in% c("NB", "PIG", "NBII", "BB", "WARING", "ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2") )  theta.star <- params[(VC$X1.d2+VC$X2.d2+2)]
  if(VC$margins[2] %in% c("D", "S", "ZABB", "ZIBB", "ZANBI", "ZINBI", "ZIPIG") )                                theta.star <- params[(VC$X1.d2+VC$X2.d2+3)]
    
}

if(!is.null(VC$X3)){
    
  if(VC$margins[2] %in% c("P", "BI", "GEOM", "LG", "YULE"))                                                     theta.star <- etatheta <- X3%*%params[(VC$X1.d2+VC$X2.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2)]
  if(VC$margins[2] %in% c("NB", "PIG", "NBII", "BB", "WARING", "ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2") )  theta.star <- etatheta <- X4%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2)]
  if(VC$margins[2] %in% c("D", "S", "ZABB", "ZIBB", "ZANBI", "ZINBI", "ZIPIG") )                                theta.star <- etatheta <- X5%*%params[(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2+1):(VC$X1.d2+VC$X2.d2+VC$X3.d2+VC$X4.d2+VC$X5.d2)]  
       
}    
    
  precision <- 10^(-8)
    
  marginComp <- marginBitsD(params=params, eta1=eta1, eta2=eta2, X3=X3, X4=X4, dat=dat, VC=VC, precision=precision, univariate=FALSE, bd=bd)
  
  
  i0 <- marginComp$i0 
  i1 <- marginComp$i1 
  i2 <- marginComp$i2 
  ind <- marginComp$ind 
  F1 <- marginComp$F1 
  dF1 <- marginComp$dF1 
  d2F1ddelta1delta1 <- marginComp$d2F1ddelta1delta1 
  mu <- marginComp$mu 
  sigma <- marginComp$sigma 
  nu <- marginComp$nu 
  F2 <- marginComp$F2 
  f2 <- marginComp$f2 
  F22 <- marginComp$F22 
  df2 <- marginComp$df2 
  df2.sigma <- marginComp$df2.sigma 
  df2.nu <- marginComp$df2.nu 
  dF2 <- marginComp$dF2 
  dF22 <- marginComp$dF22 
  dF2.sigma <- marginComp$dF2.sigma 
  dF22.sigma <- marginComp$dF22.sigma 
  dF2.nu <- marginComp$dF2.nu 
  dF22.nu <- marginComp$dF22.nu 
  d2f2delta22 <- marginComp$d2f2delta22 
  d2f2sigma2 <- marginComp$d2f2sigma2 
  d2f2delta2sigma <- marginComp$d2f2delta2sigma 
  d2f2nu2 <- marginComp$d2f2nu2 
  d2f2delta2nu <- marginComp$d2f2delta2nu 
  d2f2nusigma <- marginComp$d2f2nusigma
  d2f2sigmanu <- marginComp$d2f2sigmanu
  d2F2ddelta22 <- marginComp$d2F2ddelta22 
  d2F22ddelta22 <- marginComp$d2F22ddelta22 
  d2F2dsigma2 <- marginComp$d2F2dsigma2 
  d2F22dsigma2 <- marginComp$d2F22dsigma2 
  d2F2ddelta2dsigma <- marginComp$d2F2ddelta2dsigma 
  d2F22ddelta2dsigma <- marginComp$d2F22ddelta2dsigma 
  d2F2dnu2 <- marginComp$d2F2dnu2 
  d2F22dnu2 <- marginComp$d2F22dnu2 
  d2F2ddelta2dnu <- marginComp$d2F2ddelta2dnu 
  d2F22ddelta2dnu <- marginComp$d2F22ddelta2dnu
  d2F2dnudsigma <- marginComp$d2F2dnudsigma
  d2F22dnudsigma <- marginComp$d2F22dnudsigma
  d2F2dsigmadnu <- marginComp$d2F2dsigmadnu 
  d2F22dsigmadnu <- marginComp$d2F22dsigmadnu
  etasqv <- marginComp$etasqv
  etanu <- marginComp$etanu
  
  copulaComp <- copulaBitsD(VC=VC, theta.star=theta.star, F1=F1, F2=F2, F22=F22, f2=f2, precision=precision)
  
  
  theta <- copulaComp$theta 
  lx <- copulaComp$lx
  dcop1 <- copulaComp$dcop1
  dcop11 <- copulaComp$dcop11 
  dcop2 <- copulaComp$dcop2 
  dcop22 <- copulaComp$dcop22 
  dcop.theta1 <- copulaComp$dcop.theta1 
  dcop.theta2 <- copulaComp$dcop.theta2 
  d2Cdcop12 <- copulaComp$d2Cdcop12 
  d2Cdcop112 <- copulaComp$d2Cdcop112 
  d2Cdcop22 <- copulaComp$d2Cdcop22 
  d2Cdcop222 <- copulaComp$d2Cdcop222
  d2Cdcop1dcop2 <- copulaComp$d2Cdcop1dcop2 
  d2Cdcop11dcop22 <- copulaComp$d2Cdcop11dcop22 
  d2Cdcop.theta12 <- copulaComp$d2Cdcop.theta12 
  d2Cdcop.theta22 <- copulaComp$d2Cdcop.theta22 
  d2Cdcop1dcop.theta1 <- copulaComp$d2Cdcop1dcop.theta1 
  d2Cdcop11dcop.theta2 <- copulaComp$d2Cdcop11dcop.theta2
  d2Cdcop2dcop.theta1 <- copulaComp$d2Cdcop2dcop.theta1 
  d2Cdcop22dcop.theta2 <- copulaComp$d2Cdcop22dcop.theta2
  theta.append <- copulaComp$theta.append
  theta.append.der <- copulaComp$theta.append.der
    
  


  
#-------------------------------------------------------

    # Likelihood

    l.par <- VC$weights*(i0*log(F1) + i1*log(lx))


    # Gradient components:

    dl.1  <- VC$weights*(  i0*1/as.vector(F1)*dF1 + i1*as.vector(1/lx)*as.vector(-dcop1+dcop11)*dF1 )                  # dl.dbe1
    dl.2  <- VC$weights*( i1*as.vector(1/lx)*(df2-as.vector(dcop2)*dF2+as.vector(dcop22)*dF22) )                 # dl.dbe2
    dl.5  <- VC$weights*( i1*(1/lx)*(-dcop.theta1+dcop.theta2)*theta.append )                          # dl.dteta.st
     
    
  if (VC$margins[2] %in% c("P", "BI", "GEOM", "LG", "YULE")) {
     dl.3  <- NULL                  # dl.sigma.st 
     dl.4  <- NULL                    # dl.nu.st
  } else if (VC$margins[2] %in% c("NB", "NBII", "BB", "WARING")) {
      dl.3  <- VC$weights*( i1*as.vector(1/lx)*(df2.sigma-as.vector(dcop2)*dF2.sigma+as.vector(dcop22)*dF22.sigma)*sigma )         # dl.sigma.st 
      dl.4  <- NULL                   # dl.nu.st
  } else if (VC$margins[2]=="D") {
      dl.3  <- VC$weights*( i1*as.vector(1/lx)*(df2.sigma-as.vector(dcop2)*dF2.sigma+as.vector(dcop22)*dF22.sigma)*sigma )                    # dl.sigma.st 
      dl.4  <- VC$weights*( i1*as.vector(1/lx)*(df2.nu-as.vector(dcop2)*dF2.nu+as.vector(dcop22)*dF22.nu)*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2) )      # dl.nu.st
  } else if (VC$margins[2]=="PIG") {
     dl.3  <- VC$weights*( (i1*as.vector(1/lx)*(df2.sigma-as.vector(dcop2)*dF2.sigma+as.vector(dcop22)*dF22.sigma)*sigma) )      # dl.sigma.st 
     dl.4  <- NULL                   # dl.nu.st
  } else if (VC$margins[2]=="S") {
     dl.3  <- VC$weights*( i1*as.vector(1/lx)*(df2.sigma-as.vector(dcop2)*dF2.sigma+as.vector(dcop22)*dF22.sigma)*sigma )       # dl.sigma.st 
     dl.4  <- VC$weights*( i1*as.vector(1/lx)*(df2.nu-as.vector(dcop2)*dF2.nu+as.vector(dcop22)*dF22.nu)  )                     # dl.nu.st  
  } else if (VC$margins[2] %in% c("ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2")) {
      dl.3  <- VC$weights*( (i1*as.vector(1/lx)*(df2.sigma-as.vector(dcop2)*dF2.sigma+as.vector(dcop22)*dF22.sigma)*((sigma/(1-sigma)*(1+sigma/(1-sigma))-(sigma/(1-sigma))^2)/(1+sigma/(1-sigma))^2)) )          
      dl.4  <- NULL
  } else if (VC$margins[2] %in% c("ZABB", "ZIBB", "ZANBI", "ZINBI", "ZIPIG")) {
      dl.3  <- VC$weights*( i1*as.vector(1/lx)*(df2.sigma-as.vector(dcop2)*dF2.sigma+as.vector(dcop22)*dF22.sigma)*sigma )                  
      dl.4  <- VC$weights*( i1*as.vector(1/lx)*(df2.nu-as.vector(dcop2)*dF2.nu+as.vector(dcop22)*dF22.nu)*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2))      
  }
    
    
    # (Minus) Hessian components:

  
  
    d2l.11  <- -VC$weights*( i0*as.vector(-1/(F1^2)*(dF1)^2) + i0*as.vector(1/F1*d2F1ddelta1delta1) + i1*as.vector(-1/(lx^2)*(-dcop1+dcop11)^2*(dF1)^2) + i1*as.vector((1/lx)*((-d2Cdcop12+d2Cdcop112)*dF1^2 + (-dcop1+dcop11)*d2F1ddelta1delta1))  )        # d2l.be1.be1
    d2l.12  <- -VC$weights*i1*( as.vector((-1/(lx^2)*(-dcop1*dF1+dcop11*dF1)*(df2-dcop2*dF2+dcop22*dF22))+1/lx*(-d2Cdcop1dcop2*dF2*dF1+d2Cdcop11dcop22*dF22*dF1))  )                                  # d2l.be1.be2
    d2l.15  <- -VC$weights*i1*( as.vector(-1/(lx^2)*(-dcop.theta1+dcop.theta2)*(-dcop1*dF1+dcop11*dF1)+1/lx*(-d2Cdcop.theta12*dF1+d2Cdcop.theta22*dF1))*(theta.append)  )                                     # be1.teta.st
  
    d2l.22  <- -VC$weights*i1*( as.vector(-1/(lx^2)*(df2-dcop2*dF2+dcop22*dF22)^2) + as.vector(1/lx*(d2f2delta22-(d2Cdcop22*dF2^2+dcop2*d2F2ddelta22)+(d2Cdcop222*dF22^2+dcop22*d2F22ddelta22)))  )  # d2l.be2.be2
    d2l.25  <- -VC$weights*i1*( as.vector(-1/(lx^2)*(-dcop.theta1+dcop.theta2)*(df2-dcop2*dF2+dcop22*dF22)+1/lx*(-d2Cdcop1dcop.theta1*dF2+d2Cdcop11dcop.theta2*dF22))*(theta.append)  )                                           # d2l.be2.teta.st

    d2l.55  <- -VC$weights*i1*((-1/(lx^2)*(-dcop.theta1+dcop.theta2)^2+1/lx*(-d2Cdcop2dcop.theta1+d2Cdcop22dcop.theta2))*theta.append^2 + (1/lx)*(-dcop.theta1+dcop.theta2)*theta.append.der  )                                      # d2l.teta.st.teta.st

  
  
  
  if (VC$margins[2] %in% c("P", "BI", "GEOM", "LG", "YULE")) {
    
    
    d2l.13  <-  NULL                                 # d2l.be1.sigma.st
    d2l.14  <-  NULL                                     # d2l.be1.nu.st
    
    d2l.23  <- NULL                                           # d2l.be2.k.st
    d2l.24  <- NULL                                           # d2l.be2.nu.st
    
    d2l.33  <- NULL    # d2l.sigma.st.k.st
    d2l.34  <- NULL                                                                  # d2l.sigma.st.nu.st
    d2l.35  <- NULL                                                                 # d2l.sigma.st.teta.st
    
    d2l.44  <- NULL                                      # d2l.nu.st.nu.st
    d2l.45  <- NULL                                       # d2l.nu.st.teta.st
      
      
    
  } else if (VC$margins[2] %in% c("NB", "NBII", "BB", "WARING")) {
    
    
    d2l.13  <- -VC$weights*i1*(as.vector((-1/(lx^2)*(-dcop1*dF1+dcop11*dF1)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma))+1/lx*(-d2Cdcop1dcop2*dF2.sigma*dF1+d2Cdcop11dcop22*dF22.sigma*dF1))*sigma)                                   # d2l.be1.sigma.st
    d2l.14  <- NULL                                      # d2l.be1.nu.st
    
    d2l.23  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2-dcop2*dF2+dcop22*dF22)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)) + as.vector(1/lx*(d2f2delta2sigma-(d2Cdcop22*dF2*dF2.sigma+dcop2*d2F2ddelta2dsigma)+(d2Cdcop222*dF22*dF22.sigma+dcop22*d2F22ddelta2dsigma))))*sigma)                                            # d2l.be2.k.st
    d2l.24  <-  NULL                                          # d2l.be2.nu.st
    
    d2l.33  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)^2)+as.vector(1/lx*(d2f2sigma2-(d2Cdcop22*dF2.sigma^2+dcop2*d2F2dsigma2)+(d2Cdcop222*dF22.sigma^2+dcop22*d2F22dsigma2))))*sigma^2+(as.vector(1/lx)*(df2.sigma-as.vector(dcop2)*dF2.sigma+as.vector(dcop22)*dF22.sigma))*sigma )    # d2l.sigma.st.k.st
    d2l.34  <-  NULL                                                                 # d2l.sigma.st.nu.st
    d2l.35  <- -VC$weights*i1*(as.vector(-1/(lx^2)*(-dcop.theta1+dcop.theta2)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)+1/lx*(-d2Cdcop1dcop.theta1*dF2.sigma+d2Cdcop11dcop.theta2*dF22.sigma))*(theta.append)*sigma)                                                                # d2l.sigma.st.teta.st
    
    d2l.44  <-  NULL                                       # d2l.nu.st.nu.st
    d2l.45  <-  NULL                                     # d2l.nu.st.teta.st
    
    
  } else if (VC$margins[2]=="D") {
    
    d2l.13  <- -VC$weights*i1*(as.vector((-1/(lx^2)*(-dcop1*dF1+dcop11*dF1)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma))+1/lx*(-d2Cdcop1dcop2*dF2.sigma*dF1+d2Cdcop11dcop22*dF22.sigma*dF1))*sigma)                                   # d2l.be1.sigma.st
    d2l.14  <- -VC$weights*i1*(as.vector((-1/(lx^2)*(-dcop1*dF1+dcop11*dF1)*(df2.nu-dcop2*dF2.nu+dcop22*dF22.nu))+1/lx*(-d2Cdcop1dcop2*dF2.nu*dF1+d2Cdcop11dcop22*dF22.nu*dF1))*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2))                                      # d2l.be1.nu.st
    
    d2l.23  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2-dcop2*dF2+dcop22*dF22)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)) + as.vector(1/lx*(d2f2delta2sigma-(d2Cdcop22*dF2*dF2.sigma+dcop2*d2F2ddelta2dsigma)+(d2Cdcop222*dF22*dF22.sigma+dcop22*d2F22ddelta2dsigma))))*sigma)                                            # d2l.be2.k.st
    d2l.24  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2-dcop2*dF2+dcop22*dF22)*(df2.nu-dcop2*dF2.nu+dcop22*dF22.nu))+as.vector(1/lx*(d2f2delta2nu-(d2Cdcop22*dF2*dF2.nu+dcop2*d2F2ddelta2dnu)+(d2Cdcop222*dF22*dF22.nu+dcop22*d2F22ddelta2dnu))))*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2) )                                           # d2l.be2.nu.st
    
    d2l.33  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)^2)+as.vector(1/lx*(d2f2sigma2-(d2Cdcop22*dF2.sigma^2+dcop2*d2F2dsigma2)+(d2Cdcop222*dF22.sigma^2+dcop22*d2F22dsigma2))))*sigma^2+(as.vector(1/lx)*(df2.sigma-as.vector(dcop2)*dF2.sigma+as.vector(dcop22)*dF22.sigma))*sigma )    # d2l.sigma.st.k.st
    d2l.34  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)*(df2.nu-dcop2*dF2.nu+dcop22*dF22.nu))+as.vector(1/lx*(d2f2nusigma-(d2Cdcop22*dF2.sigma*dF2.nu+dcop2*d2F2dnudsigma)+(d2Cdcop222*dF22.sigma*dF22.nu+dcop22*d2F22dnudsigma))))*sigma*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2) )                                                                 # d2l.sigma.st.nu.st
    d2l.35  <- -VC$weights*i1*(as.vector(-1/(lx^2)*(-dcop.theta1+dcop.theta2)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)+1/lx*(-d2Cdcop1dcop.theta1*dF2.sigma+d2Cdcop11dcop.theta2*dF22.sigma))*(theta.append)*sigma)                                                                # d2l.sigma.st.teta.st
    
    d2l.44  <-  -VC$weights*i1*((as.vector(-1/(lx^2)*(df2.nu-dcop2*dF2.nu+dcop22*dF22.nu)^2)+as.vector(1/lx*(d2f2nu2-(d2Cdcop22*dF2.nu^2+dcop2*d2F2dnu2)+(d2Cdcop222*dF22.nu^2+dcop22*d2F22dnu2))))*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2)^2+(as.vector(1/lx)*(df2.nu-as.vector(dcop2)*dF2.nu+as.vector(dcop22)*dF22.nu))*((nu/(1-nu))*(1+(nu/(1-nu)))^2-2*(1+(nu/(1-nu)))*(nu/(1-nu))^2)/(1+nu/(1-nu))^4  )                                       # d2l.nu.st.nu.st
    d2l.45  <-  -VC$weights*i1*(as.vector(-1/(lx^2)*(-dcop.theta1+dcop.theta2)*(df2.nu-dcop2*dF2.nu+dcop22*dF22.nu)+1/lx*(-d2Cdcop1dcop.theta1*dF2.nu+d2Cdcop11dcop.theta2*dF22.nu))*theta.append*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2)  )                                     # d2l.nu.st.teta.st
      
    
  } else if (VC$margins[2]=="PIG") {
    
    d2l.13  <- -VC$weights*i1*(as.vector((-1/(lx^2)*(-dcop1*dF1+dcop11*dF1)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma))+1/lx*(-d2Cdcop1dcop2*dF2.sigma*dF1+d2Cdcop11dcop22*dF22.sigma*dF1))*sigma)                                   # d2l.be1.sigma.st
    d2l.14  <- NULL                                      # d2l.be1.nu.st
    
    d2l.23  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2-dcop2*dF2+dcop22*dF22)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)) + as.vector(1/lx*(d2f2delta2sigma-(d2Cdcop22*dF2*dF2.sigma+dcop2*d2F2ddelta2dsigma)+(d2Cdcop222*dF22*dF22.sigma+dcop22*d2F22ddelta2dsigma))))*sigma)                                            # d2l.be2.k.st
    d2l.24  <-  NULL                                          # d2l.be2.nu.st
    
    d2l.33  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)^2)+as.vector(1/lx*(d2f2sigma2-(d2Cdcop22*dF2.sigma^2+dcop2*d2F2dsigma2)+(d2Cdcop222*dF22.sigma^2+dcop22*d2F22dsigma2))))*sigma^2+(as.vector(1/lx)*(df2.sigma-as.vector(dcop2)*dF2.sigma+as.vector(dcop22)*dF22.sigma))*sigma )    # d2l.sigma.st.k.st
    d2l.34  <-  NULL                                                                 # d2l.sigma.st.nu.st
    d2l.35  <- -VC$weights*i1*(as.vector(-1/(lx^2)*(-dcop.theta1+dcop.theta2)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)+1/lx*(-d2Cdcop1dcop.theta1*dF2.sigma+d2Cdcop11dcop.theta2*dF22.sigma))*(theta.append)*sigma)                                                                # d2l.sigma.st.teta.st
    
    d2l.44  <-  NULL                                       # d2l.nu.st.nu.st
    d2l.45  <-  NULL                                     # d2l.nu.st.teta.st
                                             
    
  } else if (VC$margins[2]=="S") {
    
    d2l.13  <- -VC$weights*i1*(as.vector((-1/(lx^2)*(-dcop1*dF1+dcop11*dF1)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma))+1/lx*(-d2Cdcop1dcop2*dF2.sigma*dF1+d2Cdcop11dcop22*dF22.sigma*dF1))*sigma)                                   # d2l.be1.sigma.st
    d2l.14  <- -VC$weights*i1*( as.vector((-1/(lx^2)*(-dcop1*dF1+dcop11*dF1)*(df2.nu-dcop2*dF2.nu+dcop22*dF22.nu))+1/lx*(-d2Cdcop1dcop2*dF2.nu*dF1+d2Cdcop11dcop22*dF22.nu*dF1)) )                                      # d2l.be1.nu.st
    
    d2l.23  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2-dcop2*dF2+dcop22*dF22)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)) + as.vector(1/lx*(d2f2delta2sigma-(d2Cdcop22*dF2*dF2.sigma+dcop2*d2F2ddelta2dsigma)+(d2Cdcop222*dF22*dF22.sigma+dcop22*d2F22ddelta2dsigma))))*sigma)                                            # d2l.be2.k.st
    d2l.24  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2-dcop2*dF2+dcop22*dF22)*(df2.nu-dcop2*dF2.nu+dcop22*dF22.nu))+as.vector(1/lx*(d2f2delta2nu-(d2Cdcop22*dF2*dF2.nu+dcop2*d2F2ddelta2dnu)+(d2Cdcop222*dF22*dF22.nu+dcop22*d2F22ddelta2dnu))))  )                                           # d2l.be2.nu.st
    
    d2l.33  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)^2)+as.vector(1/lx*(d2f2sigma2-(d2Cdcop22*dF2.sigma^2+dcop2*d2F2dsigma2)+(d2Cdcop222*dF22.sigma^2+dcop22*d2F22dsigma2))))*sigma^2+(as.vector(1/lx)*(df2.sigma-as.vector(dcop2)*dF2.sigma+as.vector(dcop22)*dF22.sigma))*sigma )    # d2l.sigma.st.k.st
    d2l.34  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)*(df2.nu-dcop2*dF2.nu+dcop22*dF22.nu))+as.vector(1/lx*(d2f2sigmanu-(d2Cdcop22*dF2.sigma*dF2.nu+dcop2*d2F2dsigmadnu)+(d2Cdcop222*dF22.sigma*dF22.nu+dcop22*d2F22dsigmadnu))))*sigma )                                                                 # d2l.sigma.st.nu.st
    d2l.35  <- -VC$weights*i1*(as.vector(-1/(lx^2)*(-dcop.theta1+dcop.theta2)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)+1/lx*(-d2Cdcop1dcop.theta1*dF2.sigma+d2Cdcop11dcop.theta2*dF22.sigma))*(theta.append)*sigma)                                                                # d2l.sigma.st.teta.st
    
    d2l.44  <-  -VC$weights*i1*(as.vector(-1/(lx^2)*(df2.nu-dcop2*dF2.nu+dcop22*dF22.nu)^2)+as.vector(1/lx*(d2f2nu2-(d2Cdcop22*dF2.nu^2+dcop2*d2F2dnu2)+(d2Cdcop222*dF22.nu^2+dcop22*d2F22dnu2)))  )                                       # d2l.nu.st.nu.st
    d2l.45  <-  -VC$weights*i1*(as.vector(-1/(lx^2)*(-dcop.theta1+dcop.theta2)*(df2.nu-dcop2*dF2.nu+dcop22*dF22.nu)+1/lx*(-d2Cdcop1dcop.theta1*dF2.nu+d2Cdcop11dcop.theta2*dF22.nu))*theta.append )                                     # d2l.nu.st.teta.st
      
    
    
  } else if (VC$margins[2] %in% c("ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2")) {
    
    d2l.13  <- -VC$weights*i1*(as.vector((-1/(lx^2)*(-dcop1*dF1+dcop11*dF1)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma))+1/lx*(-d2Cdcop1dcop2*dF2.sigma*dF1+d2Cdcop11dcop22*dF22.sigma*dF1))*((sigma/(1-sigma)*(1+sigma/(1-sigma))-(sigma/(1-sigma))^2)/(1+sigma/(1-sigma))^2))                                   # d2l.be1.sigma.st
    d2l.14  <- NULL                                      # d2l.be1.nu.st
    
    d2l.23  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2-dcop2*dF2+dcop22*dF22)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)) + as.vector(1/lx*(d2f2delta2sigma-(d2Cdcop22*dF2*dF2.sigma+dcop2*d2F2ddelta2dsigma)+(d2Cdcop222*dF22*dF22.sigma+dcop22*d2F22ddelta2dsigma))))*((sigma/(1-sigma)*(1+sigma/(1-sigma))-(sigma/(1-sigma))^2)/(1+sigma/(1-sigma))^2))                                            # d2l.be2.k.st
    d2l.24  <-  NULL                                          # d2l.be2.nu.st
    
    d2l.33  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)^2)+as.vector(1/lx*(d2f2sigma2-(d2Cdcop22*dF2.sigma^2+dcop2*d2F2dsigma2)+(d2Cdcop222*dF22.sigma^2+dcop22*d2F22dsigma2))))*((sigma/(1-sigma)*(1+sigma/(1-sigma))-(sigma/(1-sigma))^2)/(1+sigma/(1-sigma))^2)^2 + (as.vector(1/lx)*(df2.sigma-as.vector(dcop2)*dF2.sigma+as.vector(dcop22)*dF22.sigma))*((sigma/(1-sigma))*(1+(sigma/(1-sigma)))^2-2*(1+(sigma/(1-sigma)))*(sigma/(1-sigma))^2)/(1+sigma/(1-sigma))^4 )    # d2l.sigma.st.k.st
    d2l.34  <-  NULL                                                                 # d2l.sigma.st.nu.st
    d2l.35  <- -VC$weights*i1*(as.vector(-1/(lx^2)*(-dcop.theta1+dcop.theta2)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)+1/lx*(-d2Cdcop1dcop.theta1*dF2.sigma+d2Cdcop11dcop.theta2*dF22.sigma))*(theta.append)*((sigma/(1-sigma)*(1+sigma/(1-sigma))-(sigma/(1-sigma))^2)/(1+sigma/(1-sigma))^2))                                                                # d2l.sigma.st.teta.st
    
    d2l.44  <-  NULL                                       # d2l.nu.st.nu.st
    d2l.45  <-  NULL                                     # d2l.nu.st.teta.st
                                             
    
  } else if (VC$margins[2] %in% c("ZABB", "ZIBB", "ZANBI", "ZINBI", "ZIPIG")) {
      
    
    d2l.13  <- -VC$weights*i1*(as.vector((-1/(lx^2)*(-dcop1*dF1+dcop11*dF1)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma))+1/lx*(-d2Cdcop1dcop2*dF2.sigma*dF1+d2Cdcop11dcop22*dF22.sigma*dF1))*sigma)                                   # d2l.be1.sigma.st
    d2l.14  <- -VC$weights*i1*(as.vector((-1/(lx^2)*(-dcop1*dF1+dcop11*dF1)*(df2.nu-dcop2*dF2.nu+dcop22*dF22.nu))+1/lx*(-d2Cdcop1dcop2*dF2.nu*dF1+d2Cdcop11dcop22*dF22.nu*dF1))*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2))                                      # d2l.be1.nu.st
    
    d2l.23  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2-dcop2*dF2+dcop22*dF22)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)) + as.vector(1/lx*(d2f2delta2sigma-(d2Cdcop22*dF2*dF2.sigma+dcop2*d2F2ddelta2dsigma)+(d2Cdcop222*dF22*dF22.sigma+dcop22*d2F22ddelta2dsigma))))*sigma)                                            # d2l.be2.k.st
    d2l.24  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2-dcop2*dF2+dcop22*dF22)*(df2.nu-dcop2*dF2.nu+dcop22*dF22.nu))+as.vector(1/lx*(d2f2delta2nu-(d2Cdcop22*dF2*dF2.nu+dcop2*d2F2ddelta2dnu)+(d2Cdcop222*dF22*dF22.nu+dcop22*d2F22ddelta2dnu))))*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2) )                                           # d2l.be2.nu.st
    
    d2l.33  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)^2)+as.vector(1/lx*(d2f2sigma2-(d2Cdcop22*dF2.sigma^2+dcop2*d2F2dsigma2)+(d2Cdcop222*dF22.sigma^2+dcop22*d2F22dsigma2))))*sigma^2+(as.vector(1/lx)*(df2.sigma-as.vector(dcop2)*dF2.sigma+as.vector(dcop22)*dF22.sigma))*sigma )    # d2l.sigma.st.k.st
    d2l.34  <- -VC$weights*i1*((as.vector(-1/(lx^2)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)*(df2.nu-dcop2*dF2.nu+dcop22*dF22.nu))+as.vector(1/lx*(d2f2sigmanu-(d2Cdcop22*dF2.sigma*dF2.nu+dcop2*d2F2dsigmadnu)+(d2Cdcop222*dF22.sigma*dF22.nu+dcop22*d2F22dsigmadnu))))*sigma*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2) )                                                                 # d2l.sigma.st.nu.st
    d2l.35  <- -VC$weights*i1*(as.vector(-1/(lx^2)*(-dcop.theta1+dcop.theta2)*(df2.sigma-dcop2*dF2.sigma+dcop22*dF22.sigma)+1/lx*(-d2Cdcop1dcop.theta1*dF2.sigma+d2Cdcop11dcop.theta2*dF22.sigma))*(theta.append)*sigma)                                                                # d2l.sigma.st.teta.st
    
    d2l.44  <-  -VC$weights*i1*((as.vector(-1/(lx^2)*(df2.nu-dcop2*dF2.nu+dcop22*dF22.nu)^2)+as.vector(1/lx*(d2f2nu2-(d2Cdcop22*dF2.nu^2+dcop2*d2F2dnu2)+(d2Cdcop222*dF22.nu^2+dcop22*d2F22dnu2))))*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2)^2+(as.vector(1/lx)*(df2.nu-as.vector(dcop2)*dF2.nu+as.vector(dcop22)*dF22.nu))*((nu/(1-nu))*(1+(nu/(1-nu)))^2-2*(1+(nu/(1-nu)))*(nu/(1-nu))^2)/(1+nu/(1-nu))^4  )                                       # d2l.nu.st.nu.st
    d2l.45  <-  -VC$weights*i1*(as.vector(-1/(lx^2)*(-dcop.theta1+dcop.theta2)*(df2.nu-dcop2*dF2.nu+dcop22*dF22.nu)+1/lx*(-d2Cdcop1dcop.theta1*dF2.nu+d2Cdcop11dcop.theta2*dF22.nu))*theta.append*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2)  )                                     # d2l.nu.st.teta.st
      
    
  }


  
  
  
  
  
  
  
  
  
#---------------------------------------------------------------------

  
  res <- -sum(l.par)
  
  ######################################################################
  
  # Creating gradient and Hessian
  
  if (VC$margins[2] %in% c("P", "BI", "GEOM", "LG", "YULE")) {
     
     
   if(is.null(VC$X3)){  
     
    H11 <- crossprod(X1*c(d2l.11),X1)
    H12 <- crossprod(X1*c(d2l.12),X2) 
    H15 <- t(t(rowSums(t(X1*c(d2l.15)))))

    H22 <- crossprod(X2*c(d2l.22),X2) 
    H25 <- t(t(rowSums(t(X2*c(d2l.25)))))

    H <- rbind( cbind( H11    , H12    ,  H15  ), 
              cbind( t(H12) , H22    ,  H25  ),
              cbind( t(H15) , t(H25) ,  sum(d2l.55) )
            ) 


    G   <- -c( colSums( c(dl.1)*X1 ) ,
               colSums( c(dl.2)*X2 )    ,   
               sum( dl.5 )                 )
    
    }
    
   if(!is.null(VC$X3)){  
     
    H11 <- crossprod(X1*c(d2l.11),X1)
    H12 <- crossprod(X1*c(d2l.12),X2)      
    H15 <- crossprod(X1*c(d2l.15),X3)    

    H22 <- crossprod(X2*c(d2l.22),X2) 
    H25 <- crossprod(X2*c(d2l.25),X3)
    H55 <- crossprod(X3*c(d2l.55),X3)  
    

    H <- rbind( cbind( H11    , H12    ,  H15  ), 
                cbind( t(H12) , H22    ,  H25  ),
                cbind( t(H15) , t(H25) ,  H55  )
              ) 


    G   <- -c( colSums( c(dl.1)*X1 ) ,
               colSums( c(dl.2)*X2 ) ,   
               colSums( c(dl.5)*X3 )   )
    
    }    
    
    
    
    
  } 
  
  
  ######################################################################
  
  
  if (VC$margins[2] %in% c("NB", "PIG", "NBII", "BB", "WARING", "ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2") ) {
      
      

   if(is.null(VC$X3)){        
      
    H11 <- crossprod(X1*c(d2l.11),X1)
    H12 <- crossprod(X1*c(d2l.12),X2) 
    H13 <- t(t(rowSums(t(X1*c(d2l.13)))))
    H15 <- t(t(rowSums(t(X1*c(d2l.15)))))

    H22 <- crossprod(X2*c(d2l.22),X2) 
    H23 <- t(t(rowSums(t(X2*c(d2l.23)))))
    H25 <- t(t(rowSums(t(X2*c(d2l.25)))))

    H <- rbind( cbind( H11    , H12    , H13  , H15  ), 
              cbind( t(H12) , H22    , H23  ,  H25  ),
              cbind( t(H13) , t(H23) , sum(d2l.33),  sum(d2l.35) ) ,
              cbind( t(H15) , t(H25) , sum(d2l.35),  sum(d2l.55) )
            ) 


    G  <- -c( colSums( c(dl.1)*X1 ) ,
              colSums( c(dl.2)*X2 )    ,
              sum( dl.3 ) ,    
              sum( dl.5 )                 )


}
      
      
      
   if(!is.null(VC$X3)){        
      
    H11 <- crossprod(X1*c(d2l.11),X1)
    H12 <- crossprod(X1*c(d2l.12),X2) 
    H13 <- crossprod(X1*c(d2l.13),X3) 
    H15 <- crossprod(X1*c(d2l.15),X4) 
  
    H22 <- crossprod(X2*c(d2l.22),X2) 
    H23 <- crossprod(X2*c(d2l.23),X3) 
    H25 <- crossprod(X2*c(d2l.25),X4) 
    
    H33 <- crossprod(X3*c(d2l.33),X3) 
    H55 <- crossprod(X4*c(d2l.55),X4) 
    H35 <- crossprod(X3*c(d2l.35),X4) 
    
   
    H <- rbind( cbind( H11    , H12    , H13  , H15  ), 
                cbind( t(H12) , H22    , H23  , H25  ),
                cbind( t(H13) , t(H23) , H33  , H35  ),
                cbind( t(H15) , t(H25) , t(H35)  , H55  )
            ) 


   G   <- -c( colSums( c(dl.1)*X1 ) ,
              colSums( c(dl.2)*X2 ) ,
              colSums( c(dl.3)*X3 ) ,    
              colSums( c(dl.5)*X4 )  )

}




  } 
  
  ######################################################################
  
  if (VC$margins[2] %in% c("D", "S", "ZABB", "ZIBB", "ZANBI", "ZINBI", "ZIPIG") ) {
      
      
    if(is.null(VC$X3)){       
      
    H11 <- crossprod(X1*c(d2l.11),X1)
    H12 <- crossprod(X1*c(d2l.12),X2) 
    H13 <- t(t(rowSums(t(X1*c(d2l.13)))))
    H14 <- t(t(rowSums(t(X1*c(d2l.14)))))
    H15 <- t(t(rowSums(t(X1*c(d2l.15)))))

    H22 <- crossprod(X2*c(d2l.22),X2) 
    H23 <- t(t(rowSums(t(X2*c(d2l.23)))))
    H24 <- t(t(rowSums(t(X2*c(d2l.24)))))
    H25 <- t(t(rowSums(t(X2*c(d2l.25)))))

    H <- rbind( cbind( H11    , H12    , H13  ,  H14 , H15  ), 
              cbind( t(H12) , H22    , H23  ,  H24 , H25  ),
              cbind( t(H13) , t(H23) , sum(d2l.33), sum(d2l.34), sum(d2l.35) ) ,
              cbind( t(H14) , t(H24) , sum(d2l.34), sum(d2l.44), sum(d2l.45) ) ,
              cbind( t(H15) , t(H25) , sum(d2l.35), sum(d2l.45), sum(d2l.55) )
            ) 


   G   <- -c( colSums( c(dl.1)*X1 ) ,
              colSums( c(dl.2)*X2 )    ,
              sum( dl.3 ) ,  
              sum( dl.4 ) ,   
              sum( dl.5 )                 )
  
 
 }
 
 
 
       
     if(!is.null(VC$X3)){       
       
     H11 <- crossprod(X1*c(d2l.11),X1)
     H12 <- crossprod(X1*c(d2l.12),X2) 
     H13 <- crossprod(X1*c(d2l.13),X3) 
     H14 <- crossprod(X1*c(d2l.14),X4) 
     H15 <- crossprod(X1*c(d2l.15),X5) 
 
     H22 <- crossprod(X2*c(d2l.22),X2) 
     H23 <- crossprod(X2*c(d2l.23),X3) 
     H24 <- crossprod(X2*c(d2l.24),X4) 
     H25 <- crossprod(X2*c(d2l.25),X5) 
 
     H33 <- crossprod(X3*c(d2l.33),X3) 
     H44 <- crossprod(X4*c(d2l.44),X4) 
     H55 <- crossprod(X5*c(d2l.55),X5) 
     
     H34 <- crossprod(X3*c(d2l.34),X4)  
     H35 <- crossprod(X3*c(d2l.35),X5)   
     H45 <- crossprod(X4*c(d2l.45),X5)        
 
     H <- rbind( cbind( H11    , H12    , H13  ,  H14 ,   H15 ) , 
                 cbind( t(H12) , H22    , H23  ,  H24 ,   H25 ) ,
                 cbind( t(H13) , t(H23) , H33,    H34 ,   H35 ) ,
                 cbind( t(H14) , t(H24) , t(H34), H44 ,   H45 ) ,
                 cbind( t(H15) , t(H25) , t(H35), t(H45), H55 )
             ) 
 
 
     G  <- -c( colSums( c(dl.1)*VC$X1 ) ,
               colSums( c(dl.2)*VC$X2 ) ,
               colSums( c(dl.3)*VC$X3 ) ,  
               colSums( c(dl.4)*VC$X4 ) ,   
               colSums( c(dl.5)*VC$X5 )  )
   
  
  }
 
 
 
 
 
 
 
 
 }
  
######################################################################
 
  
  if( ( VC$l.sp1==0 && VC$l.sp2==0 && VC$l.sp3==0 && VC$l.sp4==0 && VC$l.sp5==0 ) || VC$fp==TRUE) ps <- list(S.h = 0, S.h1 = 0, S.h2 = 0) else ps <- pen(params, qu.mag, sp, VC)
     
     
    S.res <- res
    res <- S.res + ps$S.h1
    G   <- G + ps$S.h2
    H   <- H + ps$S.h  

  list(value=res, gradient=G, hessian=H, S.h=ps$S.h, S.h2 = ps$S.h2, l=S.res, eta1=eta1, eta2=eta2,
       dl.dbe1=dl.1, dl.dbe2=dl.2, l.par=l.par, etatheta = etatheta, etasqv = etasqv , etanu = etanu, 
       dl.dsqv.st=dl.3,
       dl.dsh.st=dl.4,
       dl.dcor.st=dl.5, 
       d2l.be1.be1=d2l.11, d2l.be1.be2=d2l.12, d2l.be2.be2=d2l.22,
       d2l.be1.sqv.st=d2l.13,
       d2l.be1.sh.st=d2l.14,
       d2l.be1.cor.st=d2l.15,
       d2l.be2.sqv.st=d2l.23, 
       d2l.be2.sh.st=d2l.24,
       d2l.be2.cor.st=d2l.25,
       d2l.sqv.st.sqv.st=d2l.33,
       d2l.sqv.st.sh.st=d2l.34,
       d2l.sqv.st.cor.st=d2l.35,
       d2l.sh.st.sh.st=d2l.44, 
       d2l.sh.st.cor.st=d2l.45,
       d2l.cor.st.cor.st=d2l.55  )

}



