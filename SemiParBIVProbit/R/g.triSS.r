g.triSS <- function(params, respvec, VC, TIn){
    
  mean1 <- TIn$theta12 * TIn$eta1
  mean2 <- TIn$theta13 * TIn$eta1
  mean3 <- TIn$theta12 * TIn$eta2
  mean4 <- TIn$theta23 * TIn$eta2
  mean5 <- TIn$theta13 * TIn$eta3
  mean6 <- TIn$theta23 * TIn$eta3
  
  var1 <- 1 - TIn$theta12^2   
  var2 <- 1 - TIn$theta13^2
  var3 <- 1 - TIn$theta23^2
  
  cov1 <- TIn$theta23 - TIn$theta12 * TIn$theta13
  cov2 <- TIn$theta13 - TIn$theta12 * TIn$theta23
  cov3 <- TIn$theta12 - TIn$theta13 * TIn$theta23
  
  d.1 <- dnorm(TIn$eta1)
  d.2 <- dnorm(TIn$eta2)
  d.3 <- dnorm(TIn$eta3)
  
  p.1.11  <- mm(pbinorm( TIn$eta2,   TIn$eta3, mean1 =   mean1, mean2 =   mean2, var1 = var1, var2 = var2, cov12 =   cov1) )            
  p.1.10  <- mm(pbinorm( TIn$eta2,  -TIn$eta3, mean1 =   mean1, mean2 =  -mean2, var1 = var1, var2 = var2, cov12 =  -cov1) ) 

  p.2.11  <- mm(pbinorm( TIn$eta1,   TIn$eta3, mean1 =   mean3, mean2 =   mean4, var1 = var1, var2 = var3, cov12 =   cov2) )
  p.2.10  <- mm(pbinorm( TIn$eta1,  -TIn$eta3, mean1 =   mean3, mean2 =  -mean4, var1 = var1, var2 = var3, cov12 =  -cov2) )

  p.3.11  <- mm(pbinorm( TIn$eta1,   TIn$eta2, mean1 =   mean5, mean2 =   mean6, var1 = var2, var2 = var3, cov12 =   cov3) )
  p.3.10  <- mm(pbinorm( TIn$eta1,  -TIn$eta2, mean1 =   mean5, mean2 =  -mean6, var1 = var2, var2 = var3, cov12 =  -cov3) ) 

  upst.1   <- mm( pnorm( (TIn$eta2  - TIn$theta12 * TIn$eta1)/sqrt(1 - TIn$theta12^2)) )
  upst.2   <- mm( pnorm( (TIn$eta1  - TIn$theta12 * TIn$eta2)/sqrt(1 - TIn$theta12^2)) )
  
  dl.de1 <- - respvec$cy1/TIn$p0 * d.1 +
    respvec$y1.cy2/TIn$p10 * (d.1 - d.1 * upst.1) +
    respvec$y1.y2.cy3/TIn$p110 * d.1 * p.1.10 +
    respvec$y1.y2.y3/TIn$p111  * d.1 * p.1.11 
  
  dl.de2 <- - respvec$y1.cy2/TIn$p10 * d.2 * upst.2 +
    respvec$y1.y2.cy3/TIn$p110 * d.2 * p.2.10 +
    respvec$y1.y2.y3/TIn$p111 * d.2 * p.2.11 
    
    dl.de3 <- - respvec$y1.y2.cy3/TIn$p110 * d.3 * p.3.11 +
    respvec$y1.y2.y3/TIn$p111 * d.3 * p.3.11 
    

  mean.12 <- ( TIn$eta1 * (TIn$theta13 - TIn$theta12 * TIn$theta23) + TIn$eta2 * (TIn$theta23 - TIn$theta12 * TIn$theta13) )/( 1 - TIn$theta12^2 )
  mean.13 <- ( TIn$eta1 * (TIn$theta12 - TIn$theta13 * TIn$theta23) + TIn$eta3 * (TIn$theta23 - TIn$theta12 * TIn$theta13) )/( 1 - TIn$theta13^2 )
  mean.23 <- ( TIn$eta2 * (TIn$theta12 - TIn$theta13 * TIn$theta23) + TIn$eta3 * (TIn$theta13 - TIn$theta12 * TIn$theta23) )/( 1 - TIn$theta23^2 )
  
  deno <- 1 - TIn$theta12^2 - TIn$theta13^2 - TIn$theta23^2 + 2 * TIn$theta12 * TIn$theta13 * TIn$theta23
  sd.12   <- sqrt( deno / ( 1 - TIn$theta12^2 ) )
  sd.13   <- sqrt( deno / ( 1 - TIn$theta13^2 ) )
  sd.23   <- sqrt( deno / ( 1 - TIn$theta23^2 ) )
  
  p12.g <- probm( (TIn$eta3 - mean.12)/sd.12, "probit")$pr
  p13.g <- probm( (TIn$eta2 - mean.13)/sd.13, "probit")$pr
  p23.g <- probm( (TIn$eta1 - mean.23)/sd.23, "probit")$pr
  
  p12.g.c <- mm(1 - p12.g)
  p13.g.c <- mm(1 - p13.g)
  p23.g.c <- mm(1 - p23.g)
  
  dtheta12.dtheta12.st <- 4 * exp( 2 * TIn$theta12.st )/( exp(2 * TIn$theta12.st) + 1 )^2
  dtheta13.dtheta13.st <- 4 * exp( 2 * TIn$theta13.st )/( exp(2 * TIn$theta13.st) + 1 )^2
  dtheta23.dtheta23.st <- 4 * exp( 2 * TIn$theta23.st )/( exp(2 * TIn$theta23.st) + 1 )^2
  
  d11.12 <- dbinorm( TIn$eta1,  TIn$eta2, cov12 =  TIn$theta12)
  d11.13 <- dbinorm( TIn$eta1,  TIn$eta3, cov12 =  TIn$theta13)
  d11.23 <- dbinorm( TIn$eta2,  TIn$eta3, cov12 =  TIn$theta23)
  
  
  dl.dtheta12 <- - respvec$y1.cy2/TIn$p10   * d11.12 +
    respvec$y1.y2.cy3/TIn$p110 * d11.12 * p12.g.c +
    respvec$y1.y2.y3/TIn$p111   * d11.12 * p12.g 
  
  dl.dtheta13 <- - respvec$y1.y2.cy3/TIn$p110 * d11.13 * p13.g +
    respvec$y1.y2.y3/TIn$p111 * d11.13 * p13.g
  
  dl.dtheta23 <- - respvec$y1.y2.cy3/TIn$p110 * d11.23 * p23.g +
    respvec$y1.y2.y3/TIn$p111 * d11.23 * p23.g
  
  dl.dtheta12.st <- dl.dtheta12 * dtheta12.dtheta12.st
  dl.dtheta13.st <- dl.dtheta13 * dtheta13.dtheta13.st
  dl.dtheta23.st <- dl.dtheta23 * dtheta23.dtheta23.st
  
  GTRIVec <- list(p12.g = p12.g, p13.g = p13.g, p23.g = p23.g, 
                  p12.g.c = p12.g.c, p13.g.c = p13.g.c, p23.g.c = p23.g.c, 
                  d11.12 = d11.12, d11.13 = d11.13, d11.23 = d11.23, 
                  p.1.11 = p.1.11, p.1.10 = p.1.10,  
                  p.2.11 = p.2.11, p.2.10 = p.2.10,  
                  p.3.11 = p.3.11, p.3.10 = p.3.10, 
                  dl.de1 = VC$weights*dl.de1, dl.de2 = VC$weights*dl.de2, dl.de3 = VC$weights*dl.de3, 
                  dl.dtheta12.st = VC$weights*dl.dtheta12.st, dl.dtheta13.st = VC$weights*dl.dtheta13.st,
                  dl.dtheta23.st = VC$weights*dl.dtheta23.st,
                  mean.12 = mean.12,
                  mean.13 = mean.13,
                  mean.23 = mean.23, 
                  sd.12 = sd.12,
                  sd.13 = sd.13,
                  sd.23 = sd.23,
                  upst.1 = upst.1,
                  upst.2 = upst.2,
                  dl.dtheta12 =dl.dtheta12, dl.dtheta13 = dl.dtheta13, dl.dtheta23 = dl.dtheta23)
  
  GTRIVec
  
}