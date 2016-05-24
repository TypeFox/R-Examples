## All the functions here convert things to
##  penAA, penAB, penBB, returned as a vector.

MOI_ADDITIVE <- 1;
MOI_MULTIPLICATIVE <- 2;
MOI_RECESSIVE <- 3;
MOI_DOMINANT <- 4;

## helper functions
h_K <- function( paa, pab, pbb, p, K )
  return( paa*p^2 + 2*pab*p*(1-p) + pbb*(1-p)^2 - K )

h_model <- function( paa, pab, pbb, model ) {
  if( model==MOI_ADDITIVE )
    return( paa + pbb - 2*pab )
    #return( paa-pab - (pab-pbb) )
  if( model==MOI_DOMINANT )
    return( paa-pab )
  if( model==MOI_RECESSIVE )
    return( pab-pbb )
  if( model==MOI_MULTIPLICATIVE )
    return( paa/pab - pab/pbb )

  stop( "h_model MOI not understood." )
}

h_orAllelic <- function( paa, pab, pbb, p,  orAllelic ) {
  a <- paa*p^2 + pab*p*(1-p)
  b <- pbb*(1-p)^2 + pab*p*(1-p)
  c <- (1-paa)*p^2 + (1-pab)*p*(1-p)
  d <- (1-pbb)*(1-p)^2 + (1-pab)*p*(1-p)

  return( orAllelic - a*d/b/c )
}

solve_allelicOR <- function( pen, afreq, popPrev, orAllelic, moi ) {
  paa <- pen[1]
  pab <- pen[2]
  pbb <- pen[3]
  return( h_K( paa, pab, pbb, afreq, popPrev )^2
          + h_model( paa, pab, pbb, moi )^2
          + h_orAllelic( paa, pab, pbb, afreq, orAllelic )^2 )
}

## export - solving for pen** given allelic OR + info
pen_allelicOR <- function( moi, afreq, popPrev, orAllelic ) {
  res <- optim( rep(0.3,3), fn=solve_allelicOR,
                moi=moi, afreq=afreq, popPrev=popPrev, orAllelic=orAllelic )
  #print( str( res ) )
  return( res$par )  ## value is value at
}
## export - reverse solving for allelic OR given pen** and afreq
allelic_OR <- function( paa, pab, pbb, p ) {
  return( -h_orAllelic( paa, pab, pbb, p, 0 ) )
}

##################
## OR functions ##

## export - calculating OR1, OR2
calc_or <- function( paa, pab, pbb, p ) {
  case0 <- pbb*(1-p)^2
  case1 <- 2*pab*p*(1-p)
  case2 <- paa*p^2
  cont0 <- (1-pbb)*(1-p)^2
  cont1 <- 2*(1-pab)*p*(1-p)
  cont2 <- (1-paa)*p^2

  #or1 <- 1/( case0 * cont1 / case1 / cont0 )  ## OR from options, OR_heterozygous, OR1 output
  #or2 <- 1/( case0 * cont2 / case2 / cont0 )  ## OR_homozygous, OR2 output
  or1 <- (case1 * cont0) / (case0 * cont1)
  or2 <- (case2 * cont0) / (case0 * cont2)

  return( c(or1=or1, or2=or2) )
}

##h_or <- function( paa, pab, pbb, afreq, or1 ) {
##  #return( or1 - calc_or(paa,pab,pbb,afreq)[1] )
##  return( or1 - (pab*(1-pbb))/(pbb*(1-pab)) )
##}

#h_or2 <- function( paa, pab, pbb, afreq, or1 ) {
#  return( or1 - (paa*(1-pab))/(pab*(1-paa)) )
#}

h_or <- function( paa, pab, pbb, afreq, or1, moi ) {
  return(  ( or1-(pab*(1-pbb))/(pbb*(1-pab)) ) )

  ## STUPID, It was correct all along!
  
#  if( moi==MOI_ADDITIVE )
#    return( abs( 2*or1 - (paa*(1-pbb))/(pbb*(1-paa)) ) + abs( or1-(pab*(1-pbb))/(pbb*(1-pab)) )  )
    #return( ( 2*or1 - (paa*(1-pbb))/(pbb*(1-paa)) )^2 + ( or1-(pab*(1-pbb))/(pbb*(1-pab)) )^2  )
#    return( ( or1-(paa*(1-pab))/(pab*(1-paa)) )^2 + ( or1-(pab*(1-pbb))/(pbb*(1-pab)) )^2 )
  #  return( ( 2*or1 - (paa*(1-pbb))/(pbb*(1-paa)) )^2 + ( or1-(paa*(1-pab))/(pab*(1-paa)) )^2 +  ( or1-(pab*(1-pbb))/(pbb*(1-pab)) )^2  )
  #  return( abs( or1-(paa*(1-pab))/(pab*(1-paa)) ) + abs( or1-(pab*(1-pbb))/(pbb*(1-pab)) ) )
}

solve_or1 <- function( pen, afreq, popPrev, or1, moi ) {
  paa <- pen[1]
  pab <- pen[2]
  pbb <- pen[3]
  return( h_K( paa, pab, pbb, afreq, popPrev )^2
          + h_model( paa, pab, pbb, moi )^2
          + h_or( paa, pab, pbb, afreq, or1, moi )^2 )
  #return( abs( h_K( paa, pab, pbb, afreq, popPrev ) )
  #        + abs( h_model( paa, pab, pbb, moi ) )
  #        + abs( h_or( paa, pab, pbb, afreq, or1 ) ) )
}

## export - solving for pen** given OR
#pen_or1 <- function( moi, afreq, popPrev, or1 ) {
#  res <- optim( rep(0.3,3), fn=solve_or1,
#                moi=moi, afreq=afreq, popPrev=popPrev, or1=or1 )
#  return( res$par )
#}

pen_or1_a <- function( moi, afreq, popPrev, or1 ) {
  bestPen <- NULL
  bestVal <- -1
  
  NRUNS <- 10
  #if( popPrev < 0.01 ) NRUNS <- 100
  RUNIF_MAX <- 0.5
  if( popPrev < 0.01 ) RUNIF_MAX <- 0.05

  for( i in 1:NRUNS ) {
    #cat( i, "" )
    res <- optim( runif(n=3)*RUNIF_MAX, fn=solve_or1, moi=moi, afreq=afreq, popPrev=popPrev, or1=or1 )
    pen <- res$par
    val <- solve_or1( pen=pen, afreq=afreq, popPrev=popPrev, or1=or1, moi=moi )
    if( bestVal==-1 || val<bestVal ) {
      if( all(pen>0 & pen<1) ) {
        bestPen <- pen
        bestVal <- val
        #cat( "bestPen" ); print( bestPen );
      }
    }
  }
  return( bestPen )
}

## -- redoing this function, now using the rootSolve routine
pen_or1 <- function( moi, afreq, popPrev, or1 ) {
  if( moi == MOI_RECESSIVE )
    return( rev( pen_or1( moi=MOI_DOMINANT, afreq=1-afreq, popPrev=popPrev, or1=1/or1 ) ) )

  require( rootSolve )

  solve_or_multiroot <- function( pen, afreq, popPrev, or1, moi ) {
    paa <- pen[1]
    pab <- pen[2]
    pbb <- pen[3]

    return( c( h_K( paa, pab, pbb, afreq, popPrev ),
               h_model( paa, pab, pbb, moi ),
               h_or( paa, pab, pbb, afreq, or1 ) ) )
  }

  NRUNS <- 1000
  if( popPrev < 0.01 ) NRUNS <- 100
  RUNIF_MAX <- 0.5
  if( popPrev < 0.01 ) RUNIF_MAX <- 0.05

  NSUCRUNS <- 1 #100
  SUCCESSES <- 0
  bestroot <- NULL
  

  for( i in 1:NRUNS ) {
    try( {
      res <- multiroot(f = solve_or_multiroot,
                      start = runif(n=3)*RUNIF_MAX,
                      maxiter=10000,
                      rtol=1e-10, atol=1e-10, ctol=1e-10,
#                      rtol=1e-6, atol=1e-8, ctol=1e-8,
                      useFortran=FALSE,
                      positive=FALSE,
                      afreq=afreq, popPrev=popPrev, or1=or1, moi=moi)
      root <- res$root
      if( all( root>=0 & root<=1 ) ) {
        #cat( i ); print( res )
        #return( res$root )

        if( SUCCESSES==0 || sum(abs(res$root))<sum(abs(bestroot)) ) {
          bestroot <- res$root
          SUCCESSES <- SUCCESSES + 1
          if( SUCCESSES > NSUCRUNS )
            return( bestroot )
        }
      }
      #print( str( res$root ) )
      #print( class( res$root ) )
      #stop()
      #return( res$root )
    }, silent=TRUE )
  }
  if( SUCCESSES > 0 )
    return( bestroot )
    
  ##stop( "ALL FAILED" )
  return( pen_or1_a( moi=moi, afreq=afreq, popPrev=popPrev, or1=or1 ) )
}

##################
## AR functions ##

## export - calculting AR from pens
calc_ar <- function( paa, pab, pbb, p ) {
  case0 <- pbb*(1-p)^2
  case1 <- 2*pab*p*(1-p)
  case2 <- paa*p^2
  cont0 <- (1-pbb)*(1-p)^2
  cont1 <- 2*(1-pab)*p*(1-p)
  #cont2 <- (1-paa)*p^2

  caseA <- case2 + 0.5*case1
  caseB <- case0 + 0.5*case1
  #contA <- cont2 + 0.5*cont1  ## NOT USED!
  contB <- cont0 + 0.5*cont1

  K <- caseA+caseB

  ar <- ( 1 - caseB/K/(contB+caseB) )  ## attributable risk
  return( ar )
}

h_ar <- function( paa, pab, pbb, afreq, ar ) {
  return( ar - calc_ar(paa,pab,pbb,afreq) )
}

solve_ar <- function( pen, afreq, popPrev, ar, moi ) {
  paa <- pen[1]
  pab <- pen[2]
  pbb <- pen[3]
  return( h_K( paa, pab, pbb, afreq, popPrev )^2
          + h_model( paa, pab, pbb, moi )^2
          + h_ar( paa, pab, pbb, afreq, ar )^2 )
}

## export - solving for pen** given AR
pen_ar <- function( moi, afreq, popPrev, ar ) {
  res <- optim( rep(0.3,3), fn=solve_ar,
                moi=moi, afreq=afreq, popPrev=popPrev, ar=ar )
  return( res$par )
}

## export -- population prevalence!
calc_popPrev <- function( paa, pab, pbb, p ) {
  return( paa*p^2 + 2*pab*p*(1-p) + pbb*(1-p)^2 )
}

###########
## DEBUG ##
###########
debug_powerConversions <- function() {
  pen_allelicOR( moi=MOI_ADDITIVE, afreq=0.2, popPrev=0.3, orAllelic=2.471 ) ## should return 0.6271, 0.4227, 0.2182
  allelic_OR( 0.6271, 0.4227, 0.2182, 0.2 )

  calc_or( 0.6271, 0.4227, 0.2182, 0.2 )
  pen_or1( moi=MOI_ADDITIVE, afreq=0.2, popPrev=0.3, or1=2.622 )

  calc_ar( 0.6271, 0.4227, 0.2182, 0.2 )
  pen_ar( moi=MOI_ADDITIVE, afreq=0.2, popPrev=0.3, ar=0.2726 )
}
