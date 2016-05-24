triprobgHs <- function (params, respvec, VC, sp = NULL, qu.mag = NULL, AT = FALSE){

  eta1 <- VC$X1%*%params[1:VC$X1.d2]
  eta2 <- VC$X2%*%params[(VC$X1.d2 + 1):(VC$X1.d2 + VC$X2.d2)]
  eta3 <- VC$X3%*%params[(VC$X1.d2 + VC$X2.d2 + 1):(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)]
  
  etad <- p111 <- A <- NULL
  
  theta12.st <- params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2+1)]    
  theta12    <- tanh(theta12.st)    
  theta13.st <- params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2+2)]    
  theta13    <- tanh(theta13.st) 
  theta23.st <- params[(VC$X1.d2 + VC$X2.d2 + VC$X3.d2+3)]    
  theta23    <- tanh(theta23.st)
  
  p1 <- 10
  p2 <- 10
  p3 <- 10
  
  p11 <- 10
  p13 <- 10
  p23 <- 10
  
  Sigma <-  matrix( c( 1,        theta12, theta13,
                       theta12,        1, theta23,
                       theta13,  theta23,        1), 3 , 3) 
                   
                   
  eS <- eigen(Sigma)                 
  check.eigen <- any(eS$values < 0)
  
  if(check.eigen == TRUE){
  
    D.dash <- diag(abs(eS$values), 3, 3)
    P      <- eS$vectors
    R.dash <- P %*% D.dash %*% t(P)
    D1 <- diag(1/sqrt( diag(R.dash) ), 3, 3)
    Sigma <- D1 %*% R.dash %*% D1
    
    theta12 <- Sigma[1,2]
    theta13 <- Sigma[1,3]
    theta23 <- Sigma[2,3]   
    
    theta12.st <- atanh(theta12) 
    theta13.st <- atanh(theta13)
    theta23.st <- atanh(theta23)
    
    params <- c(params[1:(VC$X1.d2 + VC$X2.d2 + VC$X3.d2)],theta12.st,theta13.st,theta23.st)
  
  } else Sigma <- Sigma 
  

  for(i in 1:VC$n) p111[i] <- 10
  
  p011 <- mm(p23 - p111)
  p101 <- mm(p13 - p111)
  p110 <- mm(p11 - p111)
  p100 <- mm(p1 - p11 - p101)
  p010 <- mm(p2 - p11 - p011)
  p001 <- mm(p3 - p23 - p101)
  p000 <- mm(1 - p111 - p011 - p101 - p110 - p001 - p010 - p100)
  
  ##########################################################################################
  
  l.par <- VC$weights * (respvec$y1.y2.y3    * log(p111) + respvec$y1.y2.cy3  * log(p110) +
                         respvec$cy1.y2.y3   * log(p011) + respvec$cy1.y2.cy3 * log(p010) +
                         respvec$cy1.cy2.cy3 * log(p000) + respvec$cy1.cy2.y3 * log(p001) + 
                         respvec$y1.cy2.cy3  * log(p100) + respvec$y1.cy2.y3  * log(p101)   ) 
   
  res <- -sum(l.par) 
   
  ##########################################################################################
                                         
  TIn <- list(eta1 = eta1, eta2 = eta2, eta3 = eta3, 
              theta12 = theta12, theta13 = theta13, theta23 = theta23, 
              theta12.st = theta12.st, theta13.st = theta13.st, theta23.st = theta23.st, 
              p111 = p111, p110 = p110, p011 = p011, p010 = p010, 
              p000 = p000, p001 = p001, p100 = p100, p101 = p101)
              
  gTRI <- g.tri(respvec = respvec, VC = VC, TIn = TIn)

  G <- -c( colSums(c(gTRI$dl.de1) * VC$X1), 
           colSums(c(gTRI$dl.de2) * VC$X2), 
           colSums(c(gTRI$dl.de3) * VC$X3), 
           sum(gTRI$dl.dtheta12.st), 
           sum(gTRI$dl.dtheta13.st), 
           sum(gTRI$dl.dtheta23.st) )
           
  ##########################################################################################
               
  LgTRI <- list(p12.g = gTRI$p12.g, p13.g = gTRI$p13.g, p23.g = gTRI$p23.g, 
                p12.g.c = gTRI$p12.g.c, p13.g.c = gTRI$p13.g.c, p23.g.c = gTRI$p23.g.c, 
                d11.12 = gTRI$d11.12, d11.13 = gTRI$d11.13, d11.23 = gTRI$d11.23,
                p.1.11 = gTRI$p.1.11, p.1.10 = gTRI$p.1.10, p.1.00 = gTRI$p.1.00, 
                p.1.01 = gTRI$p.1.01, p.2.11 = gTRI$p.2.11, p.2.10 = gTRI$p.2.10, 
                p.2.00 = gTRI$p.2.00, p.2.01 = gTRI$p.2.01, p.3.11 = gTRI$p.3.11, 
                p.3.10 = gTRI$p.3.10, p.3.00 = gTRI$p.3.00, p.3.01 = gTRI$p.3.01,
                mean.12 = gTRI$mean.12,
                mean.13 = gTRI$mean.13,
                mean.23 = gTRI$mean.23, sd.12 = gTRI$sd.12,
                sd.13 = gTRI$sd.13,
                sd.23 = gTRI$sd.23,
                dl.dtheta12 = gTRI$dl.dtheta12, dl.dtheta13 = gTRI$dl.dtheta13, dl.dtheta23 = gTRI$dl.dtheta23) 
                
  HTRI <- H.tri(respvec = respvec, VC = VC, TIn = TIn, LgTRI = LgTRI)
  H    <- - rbind(HTRI$h1, HTRI$h2, HTRI$h3, HTRI$h4, HTRI$h5, HTRI$h6)
  
##########################################################################################
##########################################################################################

# && VC$l.sp4==0 && VC$l.sp5==0 && VC$l.sp6==0
# at the moment the setup below excludes the case of varying corrs
# and most probably we will not do it

if( (VC$l.sp1==0 && VC$l.sp2==0 && VC$l.sp3==0 && VC$l.sp4==0 ) || VC$fp==TRUE) ps <- list(S.h = 0, S.h1 = 0, S.h2 = 0) else ps <- pen(params, qu.mag, sp, VC)

if (VC$extra.regI == "pC") H <- regH(H, type = 1)

  S.res <- res
  res   <- S.res + ps$S.h1
  G     <- G + ps$S.h2
  H     <- H + ps$S.h
  
  
if (VC$extra.regI == "sED") H <- regH(H, type = 2)


  list(value = res, gradient = G, hessian = H, S.h = ps$S.h, qu.mag = qu.mag, qu.mag.a = ps$qu.mag.a,
       l = S.res, l.par = l.par, ps = ps, 
       eta1 = eta1, eta2 = eta2, eta3 = eta3, 
       p111 = p111,
       p011 = p011,
       p101 = p101,
       p110 = p110,
       p100 = p100,
       p010 = p010,
       p001 = p001,
       p000 = p000,
       theta12 = theta12,
       theta13 = theta13,
       theta23 = theta23,
       dl.de1 = gTRI$dl.de1,
       dl.de2 = gTRI$dl.de2,
       dl.de3 = gTRI$dl.de3,
       dl.dtheta12.st = gTRI$dl.dtheta12.st, 
       dl.dtheta13.st = gTRI$dl.dtheta13.st, 
       dl.dtheta23.st = gTRI$dl.dtheta23.st)
              
}
