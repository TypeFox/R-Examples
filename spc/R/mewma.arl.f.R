# Computation of MEWMA ARLs (multivariate mean monitoring), returns function
mewma.arl.f <- function(l, cE, p, delta=0, r=20, ntype=NULL, qm0=20, qm1=qm0) {
  if ( l<=0 | l>1 )		stop("l has to be between 0 and 1")
  if ( cE<=0 )			stop("threshold c has to be positive")
  if ( p<1 )			stop("wrong dimension parameter")
  if ( delta<0 )		stop("wrong magnitude value")
  
  if ( r<4 )			stop("resolution too small")
  if ( qm0<5 )			stop("more quadrature nodes needed")
  if ( qm1<5 )			stop("more quadrature nodes needed")
  
  if ( is.null(ntype) ) {
    if ( delta <1e-10 ) {
      ntype <- "gl2"
    } else {
      #if ( p %in% c(2,4) ) {
      if ( p==2 ) {
        ntype <- "gl3"
      } else {
        ntype <- "gl5"
      }
    }
  }
  
  # collocation basis of Chebshev polynomials
  Tn <- Vectorize(function(z, n) {
    if ( n==0 ) result <- 1
    if ( n==1 ) result <- z
    if ( n==2 ) result <- 2*z^2 - 1
    if ( n==3 ) result <- 4*z^3 - 3*z
    if ( n==4 ) result <- 8*z^4 - 8*z^2 + 1
    if ( n==5 ) result <- 16*z^5 - 20*z^3 + 5*z
    if ( n>5 )  result <- cos( n*acos(z) ) 
    result
  })
  
  qtyp <- pmatch(tolower(ntype), c("gl", "co", "ra", "cc", "mc", "sr", "co2", "gl2", "gl3", "gl4", "gl5", "co3", "co4")) - 1
  if ( is.na(qtyp) )		stop("invalid type of numerical algorithm")

  if ( abs(delta) < 1e-10 ) { # in-control    
    LENGTH <- 3*r
    zeug <- .C("mewma_arl_f", as.double(l), as.double(cE), as.integer(p), as.double(delta),
                              as.integer(r), as.integer(qtyp), as.integer(qm0), as.integer(qm1),
                              ans=double(length=LENGTH), PACKAGE="spc")$ans
                              
    g <- zeug[1:r]                          
    w <- zeug[1:r + r]
    z <- zeug[1:r + 2*r]
    
    # helper functions
    cE <- cE * l/(2-l)
    l2 <- ( (1-l)/l )^2 
    fchi <- function(a, u) dchisq( u/l^2, p, ncp=l2*a ) / l^2
    FCHI <- function(a, u) pchisq( u/l^2, p, ncp=l2*a )
    
    if ( qtyp %in% c(0,2,5) ) arl <- Vectorize(function(a) 1 + sum( w * fchi(a, z) * g ), "a") # ordinary GL or Radau or Simpson rule Nystroem 
    if ( qtyp==7 ) arl <- Vectorize(function(a) 1 + sum( w * fchi(a, z^2) * g * 2*z ), "a") # GL Nystroem with ()^2 substitution
    if ( qtyp==1 ) arl <- Vectorize(function(a) sum( Tn( (2*a-cE)/cE, 0:(r-1) ) * g ), "a") # collocation
    if ( qtyp==3 ) arl <- Vectorize(function(a) 1 + sum( w * fchi(a, z^2) * g ) * cE/2 , "a") # Clenshaw-Curtis
    if ( qtyp==4 ) arl <- Vectorize(function(a) 1 + sum( ( FCHI(a, z^2) - FCHI(a, c(0, z[-length(z)]^2)) ) * g), "a") # Markov chain (Runger/Prabhu)
    
  } else { # out-of-control
    if ( qtyp==4 ) {
      cE_ <- sqrt( cE * l/(2-l) )
      w  <- 2*cE_/( 2*r + 1 )
      ii <- function(ix, iy) (ix-r)^2 + iy^2 < cE_^2/w^2
      CIRC <- as.vector( t(outer( 0:(2*r), 0:r, ii)) )
      dQ <- sum(CIRC)
      LENGTH <- dQ
    } else {
      r2 <- r^2
      LENGTH <- r2 + 4*r
    }
    zeug <- .C("mewma_arl_f", as.double(l), as.double(cE), as.integer(p), as.double(delta),
                              as.integer(r), as.integer(qtyp), as.integer(qm0), as.integer(qm1),
                              ans=double(length=LENGTH), PACKAGE="spc")$ans                              
    if ( qtyp!=4 ) {
      g  <- zeug[1:r2]
      w0 <- zeug[1:r + r2]
      z0 <- zeug[1:r + r2 + r]
      w1 <- zeug[1:r + r2 + 2*r]
      z1 <- zeug[1:r + r2 + 3*r]
    } else {
      g <- zeug[1:dQ]
    }
    
    # helpers
    l2 <- ( (1-l)/l )^2
    h <- cE * l/(2-l)
    rdc <- l * sqrt( delta/h )
    sig <- l / sqrt( h )
    
    if ( qtyp %in% c(0, 2, 3, 5) ) arl <- Vectorize(function(a, b) { # ordinary GL or Radau or CC or Simpson rule Nystroem 
      if ( abs(h-a) < 1e-10 ) a_ <- 1 else a_ <- ( a - b^2/delta ) / ( h - b^2/delta )
      b_ <- b / sqrt( delta * h )
      m <- rdc + (1-l)*b_
      eta <- l2 * h * (1-b_^2) * a_
      if ( eta < 1e-10 ) eta <- 0
      result <- 1
      for ( i in 1:r ) {
        korr <- h * (1-z1[i]^2) / l^2 
        outer <- korr * w1[i] * dnorm( z1[i], mean=m, sd=sig)
        inner <- sum( w0 * dchisq( korr*z0, p-1, eta) * g[ (i-1)*r + 1:r ] )
        result <- result + inner * outer
      }
      result
    })
    
    if ( qtyp==7 ) arl <- Vectorize(function(a, b) { # GL Nystroem with ()^2 substitution
      if ( abs(h-a) < 1e-10 ) a_ <- 1 else a_ <- ( a - b^2/delta ) / ( h - b^2/delta )
      b_ <- b / sqrt( delta * h )
      m <- rdc + (1-l)*b_
      eta <- l2 * h * (1-b_^2) * a_
      if ( eta < 1e-10 ) eta <- 0
      result <- 1
      for ( i in 1:r ) {
        korr <- h * (1-z1[i]^2) / l^2 
        outer <- korr * w1[i] * dnorm( z1[i], mean=m, sd=sig)
        inner <- sum( w0 * dchisq( korr*z0^2, p-1, eta) * 2*z0 * g[ (i-1)*r + 1:r ] )
        result <- result + inner * outer
      }
      result
    })
        
    if ( qtyp==8 ) arl <- Vectorize(function(a, b) { # GL Nystroem with ()^2 plus sin() substitution
      if ( abs(h-a) < 1e-10 ) a_ <- 1 else a_ <- ( a - b^2/delta ) / ( h - b^2/delta )
      b_ <- b / sqrt( delta * h )
      m <- rdc + (1-l)*b_
      eta <- l2 * h * (1-b_^2) * a_
      if ( eta < 1e-10 ) eta <- 0
      result <- 1
      for ( i in 1:r ) {
        korr <- h * ( 1 - sin(z1[i])^2 ) / l^2 
        outer <- korr * w1[i] * dnorm( sin(z1[i]), mean=m, sd=sig) * cos(z1[i])
        inner <- sum( w0 * dchisq( korr*z0^2, p-1, eta) * 2*z0 * g[ (i-1)*r + 1:r ] )
        result <- result + inner * outer
      }
      result
    })
    
    if ( qtyp==9 ) arl <- Vectorize(function(a, b) { # GL Nystroem with ()^2 plus tan() substitution
      if ( abs(h-a) < 1e-10 ) a_ <- 1 else a_ <- ( a - b^2/delta ) / ( h - b^2/delta )
      b_ <- b / sqrt( delta * h )
      m <- rdc + (1-l)*b_
      eta <- l2 * h * (1-b_^2) * a_
      if ( eta < 1e-10 ) eta <- 0
      result <- 1
      for ( i in 1:r ) {
        korr <- h * ( 1 - tan(z1[i])^2 ) / l^2 
        outer <- korr * w1[i] * dnorm( tan(z1[i]), mean=m, sd=sig) / cos(z1[i])^2
        inner <- sum( w0 * dchisq( korr*z0^2, p-1, eta) * 2*z0 * g[ (i-1)*r + 1:r ] )
        result <- result + inner * outer
      }
      result
    })   
 
    if ( qtyp==10 ) arl <- Vectorize(function(a, b) { # GL Nystroem with ()^2 plus sinh() substitution
      norm <- sinh(1)
      if ( abs(h-a) < 1e-10 ) a_ <- 1 else a_ <- ( a - b^2/delta ) / ( h - b^2/delta )
      b_ <- b / sqrt( delta * h )
      m <- rdc + (1-l)*b_
      eta <- l2 * h * (1-b_^2) * a_
      if ( eta < 1e-10 ) eta <- 0
      result <- 1
      for ( i in 1:r ) {
        korr <- h * ( 1 - (sinh(z1[i])/norm)^2 ) / l^2         
        outer <- korr * w1[i] * dnorm( sinh(z1[i])/norm, mean=m, sd=sig) * cosh(z1[i])/norm
        inner <- sum( w0 * dchisq( korr*z0^2, p-1, eta) * 2*z0 * g[ (i-1)*r + 1:r ] )
        result <- result + inner * outer
      }
      result
    })
    
    if ( qtyp %in% c(1, 6, 11, 12) ) arl <- Vectorize(function(a, b) { # collocation 
      if ( abs(h-a) < 1e-10 ) a_ <- 1 else a_ <- ( a - b^2/delta ) / ( h - b^2/delta )
      b_ <- b / sqrt( delta * h )
      result <- 0
      for ( i in 1:r ) {
        outer <- Tn( 2*a_-1, i-1 ) 
        inner <- sum( Tn( b_, 0:(r-1)) * g[ (i-1)*r + 1:r ] )
        result <- result + inner * outer
      }
      result
    })
    
    if ( qtyp==4 ) arl <- Vectorize(function(a, b) { # Markov chain (Runger/Prabhu)
      #a   <- sqrt(a)
      cE_ <- sqrt( cE * l/(2-l) )
      w  <- 2*cE_/( 2*r + 1 )
      wl <- w^2 / l^2
      ii <- function(ix, iy) (ix-r)^2 + iy^2 < cE_^2/w^2
      Vf <- function(iy,jy) pchisq( (jy+.5)^2*wl, p-1, ncp=(iy*w)^2*l2) - as.numeric(jy>0)*pchisq( (jy-.5)^2*wl, p-1, ncp=(iy*w)^2*l2)
      Hf <- function(ix,jx) pnorm( (-cE_+(jx+1)*w-(1-l)*(-cE_+(ix+.5)*w) )/l, mean=delta) - pnorm( (-cE_+jx*w-(1-l)*(-cE_+(ix+.5)*w) )/l, mean=delta)
      CIRC <- as.vector( t(outer( 0:(2*r), 0:r, ii)) )
      Vv <- Vf( a/w, 0:r )
      Hv <- Hf( (b+cE_)/w-.5, 0:(2*r) )
      dQ <- sum(CIRC)
      Qv <- as.vector( Vv %o% Hv )
      Qv_ <- Qv[ CIRC ]
      result <- 1 + sum( Qv_ * g )
    })
    
  } 
  
  arl
}
