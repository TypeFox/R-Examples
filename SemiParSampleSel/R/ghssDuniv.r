ghssDuniv <- function (params, dat, VC, qu.mag = NULL, sp = NULL) {
    etasqv <- etanu <- NULL
    bd <- VC$bd
    X2 <- as.matrix(VC$X2)
    if (!is.null(VC$X3)) {
        X3 <- as.matrix(VC$X3)
        if (!is.null(VC$X4)) 
            X4 <- as.matrix(VC$X4)
    }
    eta1 <- NULL
    eta2 <- X2 %*% params[1:VC$X2.d2]
    precision <- 10^(-8)
    
    marginComp <- marginBitsD(params=params, eta1=eta1, eta2=eta2, X3=X3, X4=X4, dat=dat, VC=VC, precision=precision, univariate=TRUE, bd=bd)
  
    i0 <- marginComp$i0 
    i1 <- marginComp$i1 
    i2 <- marginComp$i2 
    ind <- marginComp$ind  
    mu <- marginComp$mu 
    sigma <- marginComp$sigma 
    nu <- marginComp$nu 
    F2 <- marginComp$F2 
    f2 <- marginComp$f2  
    df2 <- marginComp$df2 
    df2.sigma <- marginComp$df2.sigma 
    df2.nu <- marginComp$df2.nu  
    d2f2delta22 <- marginComp$d2f2delta22 
    d2f2sigma2 <- marginComp$d2f2sigma2 
    d2f2delta2sigma <- marginComp$d2f2delta2sigma 
    d2f2nu2 <- marginComp$d2f2nu2 
    d2f2delta2nu <- marginComp$d2f2delta2nu
    d2f2nusigma <- marginComp$d2f2nusigma
    d2f2sigmanu <- marginComp$d2f2sigmanu
    etasqv <- marginComp$etasqv
    etanu <- marginComp$etanu
  
    l.par <- VC$weights*i1*(log(f2))
    dl.2  <- VC$weights*i1*(1/f2)*(df2)               
    if (VC$margins[2] %in% c("P", "BI", "GEOM", "LG", "YULE")) {
     dl.3  <- NULL                  
     dl.4  <- NULL                 
  } else if (VC$margins[2] %in% c("NB", "NBII", "BB", "WARING")) {
      dl.3  <- VC$weights*i1*(1/f2)*(df2.sigma*sigma)         
      dl.4  <- NULL                  
  } else if (VC$margins[2]=="D") {
      dl.3  <- VC$weights*i1*(1/f2)*(df2.sigma*sigma)                    
      dl.4  <- VC$weights*i1*(1/f2)*(df2.nu*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2))     
  } else if (VC$margins[2]=="PIG") {
     dl.3  <- VC$weights*i1*(1/f2)*(df2.sigma*sigma)       
     dl.4  <- NULL                   
  } else if (VC$margins[2]=="S") {
     dl.3  <- VC$weights*i1*(1/f2)*(df2.sigma*sigma)       
     dl.4  <- VC$weights*i1*(1/f2)*(df2.nu)
  } else if (VC$margins[2] %in% c("ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2")) {
      dl.3  <- VC$weights*i1*(1/f2)*(df2.sigma*((sigma/(1-sigma)*(1+sigma/(1-sigma))-(sigma/(1-sigma))^2)/(1+sigma/(1-sigma))^2))         
      dl.4  <- NULL
  } else if (VC$margins[2] %in% c("ZABB", "ZIBB", "ZANBI", "ZINBI", "ZIPIG")) {
      dl.3  <- VC$weights*i1*(1/f2)*(df2.sigma*sigma)                    
      dl.4  <- VC$weights*i1*(1/f2)*(df2.nu*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2))     
  }
  
  
  
    d2l.22  <- -VC$weights*i1*(-(1/f2^2)*df2^2+(1/f2)*d2f2delta22)  
  if (VC$margins[2] %in% c("P", "BI", "GEOM", "LG", "YULE")) {
    d2l.23  <- NULL                                       
    d2l.24  <- NULL                                         
    d2l.33  <- NULL   
    d2l.34  <- NULL           
    d2l.44  <- NULL                                      
  } else if (VC$margins[2] %in% c("NB", "NBII", "BB", "WARING")) {
    d2l.23  <- -VC$weights*i1*( -(1/f2^2)*df2.sigma*df2+(1/f2)*d2f2delta2sigma )*sigma             
    d2l.24  <-  NULL                  
    d2l.33  <- -VC$weights*i1*((-(1/f2^2)*df2.sigma^2+(1/f2)*d2f2sigma2)*sigma^2 + (1/f2)*df2.sigma*sigma)   
    d2l.34  <-  NULL                                                              
    d2l.44  <-  NULL                            
  } else if (VC$margins[2]=="D") {
    d2l.23  <- -VC$weights*i1*( -(1/f2^2)*df2.sigma*df2+(1/f2)*d2f2delta2sigma )*sigma                                     
    d2l.24  <- -VC$weights*i1*( (-(1/f2^2)*df2.nu*df2+(1/f2)*d2f2delta2nu)*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2) )         
    d2l.33  <- -VC$weights*i1*((-(1/f2^2)*df2.sigma^2+(1/f2)*d2f2sigma2)*sigma^2 + (1/f2)*df2.sigma*sigma)    # d2l.sigma.st.k.st
    d2l.34  <- -VC$weights*i1*( (-(1/f2^2)*df2.nu*df2.sigma+(1/f2)*d2f2nusigma)*sigma*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2) )
    d2l.44  <-  -VC$weights*i1*((-(1/f2^2)*df2.nu^2+(1/f2)*d2f2nu2)*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2)^2+(1/f2)*(df2.nu)*((nu/(1-nu))*(1+(nu/(1-nu)))^2-2*(1+(nu/(1-nu)))*(nu/(1-nu))^2)/(1+nu/(1-nu))^4  )
  } else if (VC$margins[2]=="PIG") {
    d2l.23  <- -VC$weights*i1*( -(1/f2^2)*df2.sigma*df2+(1/f2)*d2f2delta2sigma  )*sigma                  
    d2l.24  <-  NULL                                        
    d2l.33  <- -VC$weights*i1*( (-(1/f2^2)*df2.sigma^2+(1/f2)*d2f2sigma2)*sigma^2 + (1/f2)*df2.sigma*sigma ) 
    d2l.34  <-  NULL    
    d2l.44  <-  NULL
  } else if (VC$margins[2]=="S") {
    d2l.23  <- -VC$weights*i1*( -(1/f2^2)*df2.sigma*df2+(1/f2)*d2f2delta2sigma  )*sigma                                           
    d2l.24  <- -VC$weights*i1*(   (-(1/f2^2)*df2.nu*df2+(1/f2)*d2f2delta2nu)   )                                         
    d2l.33  <- -VC$weights*i1*( (-(1/f2^2)*df2.sigma^2+(1/f2)*d2f2sigma2)*sigma^2 + (1/f2)*df2.sigma*sigma )    
    d2l.34  <- -VC$weights*i1*( (-(1/f2^2)*df2.nu*df2.sigma+(1/f2)*d2f2sigmanu)*sigma )                            
    d2l.44  <-  -VC$weights*i1*( (-(1/f2^2)*df2.nu^2+(1/f2)*d2f2nu2)  )
  } else if (VC$margins[2] %in% c("ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2")) {
    d2l.23  <- -VC$weights*i1*( -(1/f2^2)*df2.sigma*df2+(1/f2)*d2f2delta2sigma  )*((sigma/(1-sigma)*(1+sigma/(1-sigma))-(sigma/(1-sigma))^2)/(1+sigma/(1-sigma))^2)                  
    d2l.24  <-  NULL                                        
    d2l.33  <- -VC$weights*i1*( (-(1/f2^2)*df2.sigma^2+(1/f2)*d2f2sigma2)*((sigma/(1-sigma)*(1+sigma/(1-sigma))-(sigma/(1-sigma))^2)/(1+sigma/(1-sigma))^2)^2 + (1/f2)*df2.sigma*((sigma/(1-sigma))*(1+(sigma/(1-sigma)))^2-2*(1+(sigma/(1-sigma)))*(sigma/(1-sigma))^2)/(1+sigma/(1-sigma))^4 ) 
    d2l.34  <-  NULL    
    d2l.44  <-  NULL  
  } else if (VC$margins[2] %in% c("ZABB", "ZIBB", "ZANBI", "ZINBI", "ZIPIG")) {
    d2l.23  <- -VC$weights*i1*( -(1/f2^2)*df2.sigma*df2+(1/f2)*d2f2delta2sigma )*sigma                                     
    d2l.24  <- -VC$weights*i1*( (-(1/f2^2)*df2.nu*df2+(1/f2)*d2f2delta2nu)*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2) )         
    d2l.33  <- -VC$weights*i1*((-(1/f2^2)*df2.sigma^2+(1/f2)*d2f2sigma2)*sigma^2 + (1/f2)*df2.sigma*sigma)    # d2l.sigma.st.k.st
    d2l.34  <- -VC$weights*i1*( (-(1/f2^2)*df2.nu*df2.sigma+(1/f2)*d2f2sigmanu)*sigma*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2) )
    d2l.44  <-  -VC$weights*i1*((-(1/f2^2)*df2.nu^2+(1/f2)*d2f2nu2)*((nu/(1-nu)*(1+nu/(1-nu))-(nu/(1-nu))^2)/(1+nu/(1-nu))^2)^2+(1/f2)*(df2.nu)*((nu/(1-nu))*(1+(nu/(1-nu)))^2-2*(1+(nu/(1-nu)))*(nu/(1-nu))^2)/(1+nu/(1-nu))^4  )
  }
  
  
    res <- -sum(l.par)
    if (VC$margins[2] %in% c("P", "BI", "GEOM", "LG", "YULE")) {
        H22 <- crossprod(X2 * c(d2l.22), X2)
        H <- H22
        G <- -colSums(c(dl.2) * X2)
    }
    if (VC$margins[2] %in% c("NB", "PIG", "NBII", "BB", "WARING", "ZABI", "ZIBI", "ZALG", "ZAP", "ZIP", "ZIP2") ) {
        if (is.null(VC$X3)) {
            H22 <- crossprod(X2 * c(d2l.22), X2)
            H23 <- t(t(rowSums(t(X2 * c(d2l.23)))))
            H <- rbind(cbind(H22, H23), cbind(t(H23), sum(d2l.33)))
            G <- -c(colSums(c(dl.2) * X2), sum(dl.3))
        }
        if (!is.null(VC$X3)) {
            H22 <- crossprod(X2 * c(d2l.22), X2)
            H23 <- crossprod(X2 * c(d2l.23), X3)
            H33 <- crossprod(X3 * c(d2l.33), X3)
            H <- rbind(cbind(H22, H23), cbind(t(H23), H33))
            G <- -c(colSums(c(dl.2) * X2), colSums(c(dl.3) * 
                X3))
        }
    }
    if (VC$margins[2] %in% c("D", "S", "ZABB", "ZIBB", "ZANBI", "ZINBI", "ZIPIG")) {
        if (is.null(VC$X3)) {
            H22 <- crossprod(X2 * c(d2l.22), X2)
            H23 <- t(t(rowSums(t(X2 * c(d2l.23)))))
            H24 <- t(t(rowSums(t(X2 * c(d2l.24)))))
            H <- rbind(cbind(H22, H23, H24), cbind(t(H23), sum(d2l.33), 
                sum(d2l.34)), cbind(t(H24), sum(d2l.34), sum(d2l.44)))
            G <- -c(colSums(c(dl.2) * X2), sum(dl.3), sum(dl.4))
        }
        if (!is.null(VC$X3)) {
            H22 <- crossprod(X2 * c(d2l.22), X2)
            H23 <- crossprod(X2 * c(d2l.23), X3)
            H24 <- crossprod(X2 * c(d2l.24), X4)
            H33 <- crossprod(X3 * c(d2l.33), X3)
            H44 <- crossprod(X4 * c(d2l.44), X4)
            H34 <- crossprod(X3 * c(d2l.34), X4)
            H <- rbind(cbind(H22, H23, H24), cbind(t(H23), H33, 
                H34), cbind(t(H24), t(H34), H44))
            G <- -c(colSums(c(dl.2) * VC$X2), colSums(c(dl.3) * 
                VC$X3), colSums(c(dl.4) * VC$X4))
        }
    }
    if ((VC$l.sp2 == 0 && VC$l.sp3 == 0 && VC$l.sp4 == 0) || 
        VC$fp == TRUE) 
        ps <- list(S.h = 0, S.h1 = 0, S.h2 = 0)
    else ps <- pen(params, qu.mag, sp, VC, eq1 = "no")
    S.res <- res
    res <- S.res + ps$S.h1
    G <- G + ps$S.h2
    H <- H + ps$S.h
    list(value = res, gradient = G, hessian = H, S.h = ps$S.h, 
        S.h2 = ps$S.h2, l = S.res, eta2 = eta2, dl.dbe2 = dl.2, 
        l.par = l.par, etasqv = etasqv, etanu = etanu, dl.dsqv.st = dl.3, 
        dl.dsh.st = dl.4, d2l.be2.be2 = d2l.22, d2l.be2.sqv.st = d2l.23, 
        d2l.be2.sh.st = d2l.24, d2l.sqv.st.sqv.st = d2l.33, d2l.sqv.st.sh.st = d2l.34, 
        d2l.sh.st.sh.st = d2l.44)
}