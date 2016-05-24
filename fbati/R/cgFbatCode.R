## Moved into where it was used now...
#require( rootSolve )


## THIS is now set in the functions, and would have to be manually changed in
##  copyFBAT
##  haplotypeDensity
##  condGeneP2
#FBAT <- "fbat_2.02c"
##FBAT <- "/tmp/th/fbat_2.02c"

## for simulations -- trying to copy it to the local hard drive if possible...
copyFBAT <- function(FBAT="/tmp/th/fbat_2.02c") {
  FBATSOURCE <- "~/bin/fbat_2.02c"
  if( !file.exists(FBAT) ) {
    #file.copy( from=FBATSOURCE, to=FBAT )
    system( paste( "cp", FBATSOURCE, FBAT ) )
  }
}

condGeneFBATControl_load <- function( filename ) {
  i = as.integer(0)
  res <- .C("condGeneFBATControl_load", filename, i, DUP=TRUE )
  return( res[[2]] )
}

condGeneFBATControl_free <- function( reference ) {
  .C( "condGeneFBATControl_free", as.integer(reference), DUP=TRUE ) #DUP=FALSE )
}

condGeneFBATControl_print <- function( reference ) {
  .C( "condGeneFBATControl_print", as.integer(reference), DUP=TRUE ) #DUP=FALSE )
}

condGeneFBATControl_linkTrait <- function( reference, pid, trait ) {
  #print( reference )
  #print( head( pid ) )
  #print( head( trait ) )
  .C( "condGeneFBATControl_linkTrait", as.integer(reference), as.integer(pid), as.double(trait), as.integer(length(pid)), DUP=TRUE, NAOK=TRUE) #DUP=FALSE, NAOK=TRUE )
}

condGeneFBATControl_numFam <- function( reference ) {
  ret <- as.integer( 0 );
  res <- .C( "condGeneFBATControl_numFam", as.integer(reference), ret, DUP=TRUE ) #DUP=FALSE )
  ##cat( "n =", res[[2]], "\n" ) ## DEBUG ONLY
  return( res[[2]] )
}

condGeneFBATControl_proportionInformative <- function( reference ) {
  informative <- as.double(-1)
  ##res <- .C( "condGeneFBATControl_proportionInformative", as.integer(reference), informative, DUP=FALSE )  ## codetools fun
  #.C( "condGeneFBATControl_proportionInformative", as.integer(reference), informative, DUP=FALSE )
  #return( informative )
  res = .C( "condGeneFBATControl_proportionInformative", as.integer(reference), informative, DUP=TRUE )
  return(res[[2]])
}

condGeneFBATControl_removeUnphased <- function( reference ) {
  ##cat( "# Families before removing unphased:", condGeneFBATControl_numFam( reference ), "\n" )
  .C( "condGeneFBATControl_removeUnphased", as.integer(reference), DUP=TRUE ) #DUP=FALSE )
  ##cat( "# Families after removing unphased:", condGeneFBATControl_numFam( reference ), "\n" )
}

## for qtl
condGeneFBATControl_centerTrait <- function( reference, center=0, mean=TRUE ) {
  #.C( "condGeneFBATControl_centerTrait", as.integer(reference), as.double(center), as.integer(mean), DUP=FALSE )
  #return( center ) ## should be the mean computed by the function
  res = .C( "condGeneFBATControl_centerTrait", as.integer(reference), as.double(center), as.integer(mean), DUP=TRUE )
  return(res[[2]])
}

## binary trait
condGeneFBATControl_uimc <-
  function( reference,
            bm, bc0, bc1,
            analyze_allele_index, conditional_allele_index,
            onlyComputeConditional=TRUE ) {
  ##
  n <- condGeneFBATControl_numFam( reference )
  if( n==0 ) stop( "No informative families." );
  na <- length(analyze_allele_index)
  nc <- length(conditional_allele_index)
  ret_analyze <- as.double( rep( 0, n*na ) )
  ret_conditional0 <- as.double( rep( 0, n*nc ) )
  ret_conditional1 <- as.double( rep( 0, n*nc ) )

  #cat( "uimc bm" ); print( bm )
  #cat( "uimc bc0" ); print( bc0 )
  #cat( "uimc bc1" ); print( bc1 )

  #cat( "condGeneFBATControl_uimc before .C\n" )

  res = .C( "condGeneFBATControl_uimc",
      as.integer(reference),
      as.double(bm), as.double(bc0), as.double(bc1),
      as.integer(analyze_allele_index), as.integer(length(analyze_allele_index)),
      as.integer(conditional_allele_index), as.integer(length(conditional_allele_index)),
      as.integer(onlyComputeConditional),
      ret_analyze,
      ret_conditional0, ret_conditional1,
      DUP=TRUE, NAOK=TRUE ) #DUP=FALSE, NAOK=TRUE )
  ##return( list( ret_analyze=ret_analyze, ret_conditional0=ret_conditional0, ret_conditional1=ret_conditional1 ) )

  #cat( "condGeneFBATControl_uimc after .C\n" )

  #return( list( analyze=matrix(ret_analyze,nrow=n),
  #              conditional0=matrix(ret_conditional0,nrow=n),
  #              conditional1=matrix(ret_conditional1,nrow=n) ) )
  return(list(
    analyze=matrix(res[[10]], nrow=n),
    conditional0=matrix(res[[1]], nrow=n),
    conditional1=matrix(res[[1]], nrow=n)))

}

condGeneFBATControl_uimc_nuis <-
  function( beta,
            reference,
            analyze_allele_index, conditional_allele_index,
            onlyComputeConditional=TRUE ) {
  #cat( "beta before" )
  #print( beta )

  ## Pad in beta when doing only conditional
  R <- length(analyze_allele_index) + length(conditional_allele_index)*2
  if( length(beta) != R )
    beta <- c( rep(0,length(analyze_allele_index)), beta )
  if( length(beta) != R )
    stop( "condGeneFBATControl_uimc_nuis malformed beta" )

  #cat( "beta after" )
  #print( beta )

  ##
  na <- length(analyze_allele_index)
  nc <- length(conditional_allele_index)
  bm <- beta[1:na]
  bc0 <- beta[1:nc + na]
  bc1 <- beta[1:nc + na + nc]

  #cat( "bm" ); print( bm )
  #cat( "bc0" ); print( bc0 )
  #cat( "bc1" ); print( bc1 )

  res <- condGeneFBATControl_uimc( reference=reference,
                                   bm=bm, bc0=bc0, bc1=bc1,
                                   analyze_allele_index=analyze_allele_index, conditional_allele_index=conditional_allele_index,
                                   onlyComputeConditional=onlyComputeConditional )

  ret <- rep( 0, length(beta) )
  if( na>0 )
    for( i in 1:na )
      ret[i] <- sum( res$analyze[,i], na.rm=TRUE )
  for( i in 1:nc ) {
    ret[i+na]    <- sum( res$conditional0[,i], na.rm=TRUE )
    ret[i+na+nc] <- sum( res$conditional1[,i], na.rm=TRUE )
  }

  if( onlyComputeConditional )
   ret <- ret[(na+1):length(ret)]

  #cat( "nuisance ret" )
  #print( ret )

  return( ret )
}

condGeneFBATControl_imc <-
  function( reference,
            bm, bc0, bc1,
            analyze_allele_index, conditional_allele_index ) {
  ##
  na <- length(analyze_allele_index)
  nc <- length(conditional_allele_index)
  R <- na + nc*2

  ##print( "imc analyze_allele_index" )
  ##print( analyze_allele_index )

  ret_I <- as.double( rep( 0, R^2 ) )

  res = .C( "condGeneFBATControl_imc",
      as.integer(reference),
      as.double(bm), as.double(bc0), as.double(bc1),
      as.integer(analyze_allele_index), as.integer(length(analyze_allele_index)),
      as.integer(conditional_allele_index), as.integer(length(conditional_allele_index)),
      ret_I,
      DUP=TRUE, NAOK=TRUE ) ##DUP=FALSE, NAOK=TRUE )

  #return( matrix( ret_I, nrow=R ) )
  return(matrix(res[[9]], nrow=R))
}

condGeneFBATControl_robustStat <-
  function( reference,
            analyze_allele_index,
            conditional_allele_index ) {
  ##
  n <- condGeneFBATControl_numFam( reference )
  if( n==0 ) stop( "No informative families." );
  na <- length(analyze_allele_index)
  ret_analyze <- as.double( rep( 0, n*na ) )

  res = .C( "condGeneFBATControl_robustStat",
      as.integer( reference ),
      as.integer(analyze_allele_index), as.integer(length(analyze_allele_index)),
      as.integer(conditional_allele_index), as.integer(length(conditional_allele_index)),
      ret_analyze,
      DUP=TRUE, NAOK=TRUE) #DUP=FALSE, NAOK=TRUE )

  #return( matrix( ret_analyze, nrow=n ) )
  return(matrix(res[[6]], nrow=n))
}

condGeneFBATControl_contsUimc <-
  function( beta,
            reference,
            alpha, sigma,
            analyze_allele_index, conditional_allele_index,
            onlyComputeConditional=TRUE,
            ignoreBtX=TRUE ) {
  if( is.null(alpha) ) alpha <- 0

  ##
  n <- condGeneFBATControl_numFam( reference )
  if( n==0 ) stop( "No informative families." );
  nn <- length(analyze_allele_index) + length(conditional_allele_index)*2
  ret_b <- as.double( rep( 0, n*nn ) )

  R <- length(analyze_allele_index) + length(conditional_allele_index)*2
  if( length(beta) != R )
    beta <- c( rep(0, R-length(beta)), beta )

  res = .C( "condGeneFBATControl_contsUimc",
      as.integer(reference),
      as.double(alpha), as.double(sigma), as.double(beta),
      as.integer(analyze_allele_index), as.integer(length(analyze_allele_index)),
      as.integer(conditional_allele_index),
      as.integer(length(conditional_allele_index)),
      as.integer(onlyComputeConditional),
      as.integer(ignoreBtX),
      ret_b,
      DUP=TRUE, NAOK=TRUE ) #DUP=FALSE, NAOK=TRUE )

  #return( matrix( ret_b, nrow=n ) )
  return(matrix(res[[11]], nrow=n))
}

condGeneFBATControl_contsUimc_nuis <-
  function( beta,
            reference,
            alpha, sigma,
            analyze_allele_index, conditional_allele_index,
            onlyComputeConditional=TRUE,
            ignoreBtX=TRUE ) {
  if( is.null(alpha) ) alpha <- 0
  ##cat( "beta before " )
  ##print( beta )

  nn <- length(analyze_allele_index) + length(conditional_allele_index)*2
  if( length(beta) != nn )
    beta <- c( rep(0, nn-length(beta)), beta )
  res <- condGeneFBATControl_contsUimc( beta=beta, reference=reference, alpha=alpha, sigma=sigma, analyze_allele_index=analyze_allele_index, conditional_allele_index=conditional_allele_index, onlyComputeConditional=TRUE, ignoreBtX=ignoreBtX )

  ret <- rep( 0, length(beta) )
  for( i in 1:length(beta) )   ## We could make this loop faster...
    ret[i] <- sum( res[,i], na.rm=TRUE )
  ##cat( "return sum " )
  ##print( ret )
  ##cat( "actual return sum " )
  ##print( ret[-1:analyze_allele_index] )

  #if( length(analyze_allele_index) == 0 )
  #  return( ret ) ## special case now...
  #
  #return( ret[-1:length(analyze_allele_index)] )

  return( ret[ (length(analyze_allele_index)+1) : length(ret) ] )
}

condGeneFBATControl_contsImc <-
  function( beta,
            reference,
            alpha, sigma,
            analyze_allele_index, conditional_allele_index,
            ignoreBtX=TRUE ) {
  if( is.null(alpha) ) alpha <- 0

  R <- length(analyze_allele_index) + length(conditional_allele_index)*2
  if( length(beta) != R )
    beta <- c( rep(0, R-length(beta)), beta )

  ret_I <- as.double( rep( 0.0, R*R ) )
  #print( "regular" )
  #print( (alpha) )
  #print( (sigma) )
  #print( (beta) )
  #print( "as.double" )
  #print( as.double(alpha) )
  #print( as.double(sigma) )
  #print( as.double(beta) )
  res = .C( "condGeneFBATControl_contsImc",
      as.integer(reference),
      as.double(alpha), as.double(sigma), as.double(beta),
      as.integer(analyze_allele_index), as.integer(length(analyze_allele_index)),
      as.integer(conditional_allele_index), as.integer(length(conditional_allele_index)),
      as.integer(ignoreBtX),
      ret_I,
      DUP=TRUE, NAOK=TRUE) #DUP=FALSE, NAOK=TRUE )

  #return( matrix( ret_I, nrow=R ) )
  return(matrix(res[[10]], nrow=R))
}

condGeneFBATControl_numInfFam <-
  function( reference ) {

  numInf <- as.integer(0)

  res = .C( "condGeneFBATControl_numInfFam",
      as.integer(reference), numInf,
      DUP=TRUE ) #DUP=FALSE )
  #return( numInf )
  return(res[[2]])
} ## condGeneFBATControl_numInfFam

condGeneFBATControl_pids <-
  function( reference ) {

  n <- condGeneFBATControl_numFam( reference )
  pid <- as.integer( rep( 0, n ) )
  res = .C( "condGeneFBATControl_pids",
      as.integer(reference), pid,
      DUP=TRUE ) #DUP=FALSE )
  #return( pid )
  return(res[[2]])
}

condGeneFBATControl_estEqNuis <-
  function( referenceCondition, offset=0 ) { ## if offset is 0, continuous, otherwise, it's the offset.

  nc <- length(referenceCondition)
  nc2 <- nc*2
  ret_lhs <- as.double( rep(0,nc2*nc2) )
  ret_rhs <- as.double( rep(0,nc2) )

  res = .C( "condGeneFBATControl_estEqNuis",
      as.integer(referenceCondition), as.integer(length(referenceCondition)),
      as.double(offset),
      ret_lhs, ret_rhs,
      DUP=TRUE, NAOK=TRUE ) #DUP=FALSE, NAOK=TRUE )
  ret_lhs = res[[4]]
  ret_rhs = res[[5]]

  ret_lhs <- matrix( ret_lhs, ncol=nc2, byrow=FALSE )
  #ret_lhs <- matrix( ret_lhs, ncol=nc2 )

  #cat( "ret_lhs\n" ); print( ret_lhs );
  #cat( "ret_rhs\n" ); print( ret_rhs );

  solve.svd <- getFromNamespace( "solve.svd", "fbati" )
  bc <- solve.svd( svd(ret_lhs), ret_rhs ) ## Will do generalized inverse if necessary...

  #print( "########" )
  #print( bc )

  bc <- as.vector( bc )

  #attr(bc,"rank") <- NULL
  #attr(bc,"rcond") <- NULL

  return( bc )
}

condGeneFBATControl_estEqNuisUpdate <-
  function( referenceCondition, bc ) {
  .C( "condGeneFBATControl_estEqNuisUpdate",
      as.integer(referenceCondition), as.integer(length(referenceCondition)),
      as.double(bc),
      DUP=TRUE)#DUP=FALSE )
  return( invisible() ) ## no return
}

condGeneFBATControl_estEqNuisUpdate2 <-
  function( referenceCondition, bc ) {
  .C( "condGeneFBATControl_estEqNuisUpdate2",
      as.integer(referenceCondition), as.integer(length(referenceCondition)),
      as.double(bc),
      DUP=TRUE) #DUP=FALSE )
  return( invisible() ) ## no return
}

condGeneFBATControl_estEq <-
  function( referenceAnalyze, referenceCondition, bc, offset=0.0, verbose=FALSE ) {

  ## Some useful constants
  na <- length(referenceAnalyze)
  nc <- length(referenceCondition)
  nc2 <- nc * 2
  np <- condGeneFBATControl_numFam( referenceCondition[1] )
  R <- na + nc2

  ## The return value, will be transformed into a matrix
  ret_uij <- as.double( rep( 0, (na+nc2)*np ) )
  ret_xmxc <- as.double( rep(0, na*nc2) )
  ret_xcxc <- as.double( rep(0, nc2*nc2) )

  ## Call the C Code
  res = .C( "condGeneFBATControl_estEq",
      as.integer(referenceAnalyze), as.integer(length(referenceAnalyze)),
      as.integer(referenceCondition), as.integer(length(referenceCondition)),
      as.double(bc),  as.double(offset),
      ret_uij, ret_xmxc, ret_xcxc,
      DUP=TRUE, NAOK=TRUE) #DUP=FALSE, NAOK=TRUE )
  ret_uij = res[[7]]
  ret_xmxc = res[[8]]
  ret_xcxc = res[[9]]

  ## Form the matrix
  ret_uij <- matrix( ret_uij, nrow=np )
  ret_xmxc <- matrix( ret_xmxc, nrow=na )
  ret_xcxc <- matrix( ret_xcxc, nrow=nc2 )

  #cat( "dims ret_uij=", dim(ret_uij), ", ret_xmxc", dim(ret_xmxc), ", ret_xcxc", dim(ret_xcxc), "\n" )

  ## Compute the test statistic
  #cat( "na", na, "\n")
  uim <- t(ret_uij[,1:na])
  ####if( !is.matrix(uim) ) uim <- t(as.matrix(uim))
  uic <- t(ret_uij[,(na+1):R])  ## Will always be a matrix (conditioning on at least dim 2)
  #cat( "dims, uim=", dim(uim), ", uic=", dim(uic), "\n" )

  if( verbose ) {
    print( "condGeneFBATControl_estEq::ret_xmxc" )
    print( ret_xmxc )
    print( "condGeneFBATControl_estEq::ret_xcxc" )
    print( ret_xcxc )
  }

  ##ti = uim - ret_xmxc %*% solve( ret_xcxc, uic )
  solve.svd <- getFromNamespace( "solve.svd", "fbati" )
  #print( "dims" )
  #print( dim( solve.svd( svd(ret_xcxc), uic )))
  #print( dim( ret_xmxc %*% solve.svd( svd(ret_xcxc), uic )))
  #print( dim(uim))
  ti = t(   uim - ret_xmxc %*% solve.svd( svd(ret_xcxc), uic )   )   ## Bloody hell, transpose!
  #print( "dim(ti)" )
  #print( dim(ti) )

  #print( "ti" )
  #print( ti )

  #print( "sum(ti)")
  #print( sum(ti) )
  #print( "sum(uim)" )
  #print( sum(uim) )
  #print( "sum(ti^2)" )
  #print( sum(ti^2) )
  #print( "sum(ti)^2 / sum(ti^2)" )
  #print( sum(ti)^2 / sum(ti^2) )

  #print( "uim" )
  #print( uim )
  #print( "uic" )
  #print( uic )
  #print( "ti" )
  #print( ti )
  #stop()

  u <- rep( 0, na )
  uu <- matrix( 0, nrow=na, ncol=na )
  for( k1 in 1:na ) {
    u[k1] <- sum( ti[,k1], na.rm=TRUE )
    for( k2 in 1:na ) {
      uu[k1,k2] <- sum( ti[,k1] * ti[,k2], na.rm=TRUE )
    }
  }

  if( verbose ) {
    print( "condGeneFBATControl_estEq::u" )
    print( u )
    print( "condGeneFBATControl_estEq::uu" )
    print( uu )
  }

  pvalue <- 1; rank <- 1;
  solve.svd <- getFromNamespace( "solve.svd", "fbati" )
  try( {
        ## Compute the test statistic
        post <- solve.svd( svd(uu), u )
        #print( post )
        rank <- attr( post, "rank" )
        #print( rank )
        stat <- u %*% post
        #cat( "stat", stat, "\n" )
        ## and compute the p-value
        pvalue <- pchisq( stat, df=rank, lower.tail=FALSE )
      }
  )

  #cat( "pvalue=", pvalue, "\n")

  #numInf <- rep( 0, na )
  #for( k1 in 1:na )
  #  numInf[k1] <- sum( abs(ret_uij[,k1]) >= sqrt(.Machine$double.eps), na.rm=TRUE )
  #
  #return( list( pvalue=pvalue, rank=rank, numInf=paste(numInf,collapse=",") ) )

  #numInf_na <- rep(0,na)
  #for( k1 in 1:na )
  #  numInf_na[k1] <- sum( abs(ret_uij[,k1]) >= sqrt(.Machine$double.eps), na.rm=TRUE )
  #numInf_nc <- rep(0,nc2)
  #for( k1 in 1:nc2 )
  #  numInf_nc[k1] <- sum( abs(ret_uij[,k1+na]) >= sqrt(.Machine$double.eps), na.rm=TRUE )
  #
  #return( list( pvalue=pvalue, rank=rank, numInf=paste( paste(numInf_na,collapse=","), paste(numInf_nc,collapse=","), sep="|" ) ) )

  numInf <- ""
  for( ki in 1:na ) {
    temp <- sum( abs(ret_uij[,k1]) >= sqrt(.Machine$double.eps), na.rm=TRUE )
    if( k1!=na ) {
      numInf <- paste(numInf, temp, sep=",")
    }else{
      numInf <- as.character(temp)
    }
  }
  for( k1 in 1:nc ) {
    temp1 <- sum( abs(ret_uij[,(k1-1)*2+1]) >= sqrt(.Machine$double.eps), na.rm=TRUE )
    temp2 <- sum( abs(ret_uij[,(k1-1)*2+2]) >= sqrt(.Machine$double.eps), na.rm=TRUE )
    temp <- paste( "(", temp1, ",", temp2, ")", sep="" )
    if( k1==1 ) {
      numInf <- paste( numInf, temp, sep="|" )
    }else{
      numInf <- paste( numInf, temp, sep="," )
    }
  }

  return( list( pvalue=pvalue, rank=rank, numInf=numInf ) )
}

#void condGeneFBATControl_dUdBc(
#    int *referenceAnalyze, int *referenceAnalyzeSize,
#    int *referenceCondition, int *referenceConditionSize,
#    int *analyzeAlleleIndex, int *analyzeAlleleIndexSize,
#    int *conditionAlleleIndex, int *conditionAlleleIndexSize, // for referenceAnalyze
#    int *conditionAlleleIndex2, int *conditionAlleleIndexSize2, // for referenceCondition
#    double *bc,
#    double *ret_m, double *ret_c ) {
condGeneFBATControl_dUdBc <-
  function( referenceAnalyze,
      referenceCondition,
      analyzeAlleleIndex, conditionAlleleIndex, # for referenceAnalyze
      conditionAlleleIndex2, # for reference condition
      bc ) {
  ## Create the return pieces
  #na <- length( analyzeAlleleIndex )
  na <- length( referenceAnalyze )
  #cat( "na", na, "\n" )
  nc <- length( conditionAlleleIndex )
  nc2 <- nc*2
  ret_m <- as.double( rep(0,na*nc2) )
  ret_c <- as.double( rep(0,nc2*nc2) )

  ## Call the function
  ret = .C( "condGeneFBATControl_dUdBc",
      as.integer(referenceAnalyze), as.integer(length(referenceAnalyze)),
      as.integer(referenceCondition), as.integer(length(referenceCondition)),
      as.integer(analyzeAlleleIndex), as.integer(length(analyzeAlleleIndex)),
      as.integer(conditionAlleleIndex), as.integer(length(conditionAlleleIndex)),
      as.integer(conditionAlleleIndex2), as.integer(length(conditionAlleleIndex2)),
      as.double(bc),
      ret_m, ret_c,
      DUP=TRUE) #DUP=FALSE )
  ret_m = ret[[12]]
  ret_c = ret[[13]]


  umc <- matrix( ret_m, nrow=na, ncol=nc2 )
  ucc <- matrix( ret_c, nrow=nc2, ncol=nc2 )
  return( list(umc=umc,ucc=ucc) )
}

condGeneFBATControl_varExplConts <-
    function( reference, betaEst ) {
  ret_varExpl <- as.double(0.0)
  res = .C( "condGeneFBATControl_varExplConts",
      as.integer(reference), as.integer(length(reference)),
      as.double(betaEst),
      ret_varExpl,
      DUP=TRUE ) #DUP=FALSE )
  #return( ret_varExpl )
  return(res[[4]])
}

condGeneFBATControl_varContsMean <-
    function( reference, betaEst ) {
  ret_var <- as.double(0.0)
  res = .C( "condGeneFBATControl_varContsMean",
      as.integer(reference), as.integer(length(reference)),
      as.double(betaEst),
      ret_var,
      DUP=TRUE ) #DUP=FALSE )
  #return( as.double(ret_var) )
  return(res[[4]])
}
condGeneFBATControl_varContsModel <-
    function( reference, betaEst ) {
  ret_var <- as.double(0.0)
  res = .C( "condGeneFBATControl_varContsModel",
      as.integer(reference), as.integer(length(reference)),
      as.double(betaEst),
      ret_var,
      DUP=TRUE) #DUP=FALSE )
  #return( ret_var )
  return(res[[4]])
}


condGeneFBATControl_backupTrait <- function( reference ) {
  .C( "condGeneFBATControl_backupTrait",
      as.integer(reference), as.integer(length(reference)),
      DUP=TRUE)#DUP=FALSE)
}
condGeneFBATControl_restoreTrait <- function( reference ) {
  .C( "condGeneFBATControl_restoreTrait",
      as.integer(reference), as.integer(length(reference)),
      DUP=TRUE)#DUP=FALSE)
}


##############################################################
## PBAT WARNING: *NIX _ONLY                                  #
##   I don't know what the hell happened to the windows PBAT #
##   version, but it doesn't work properly there (back when  #
##   trying to control the damn thing for power).            #
## Send a prespecified listing of commands to any program.   #
##############################################################
programControl <- function( program, commands, filename='systemControl.sh', intern=TRUE, output="", ... ){
  ## Creates a shell
  file <- file( filename, 'w' );

  cat( '\n', program, ' <<EOF\n', sep='', file=file );

  for( i in 1:length(commands) )
    cat( commands[i], '\n', sep='', file=file );

  cat( '\nEOF', file=file );
  close( file );

  sh <- 'sh';
  if( getFromNamespace("isWindows","pbatR")() ) {
    sh <- SH.EXE; ##sh <- "sh.exe";

    if( output=="" )
      return( system( paste( sh, filename ), intern=FALSE, ... ) ); ## confuses our output pipes otherwise
    return( system( paste( sh, filename, ">", output ), intern=FALSE, ... ) );
  }

  if( output=="" )
    return( system( paste( sh, filename ), intern=intern, ... ) );
  return( system( paste( sh, filename, ">", output ), intern=intern, ... ) );
}

## ROOT solving via a different routine
## -- This isn't really necessary -- for the most part, if multiroot fails, it's gone, and not going to get any better.
mmlroot <- function( f, start, maxiter, rtol, atol, ctol, useFortran, reference, analyze_allele_index, conditional_allele_index, onlyComputeConditional,
                     minRR=rep(-4,length(start)), maxRR=rep(4,length(start)), tryRR=5, ## bounds for random restart
                     verbose=FALSE, ... ) {
  for( t in 1:tryRR ) {
    res <- NULL
    try(  res <- multiroot( f, start=start, maxiter=maxiter, rtol=rtol, atol=atol, ctol=ctol, useFortran=useFortran, reference=reference, analyze_allele_index=analyze_allele_index, conditional_allele_index=conditional_allele_index, onlyComputeConditional=onlyComputeConditional, ... ),  silent=!verbose  )
    #print( res )
    if( !is.null(res) )
      return( res )

    if( verbose ) cat( t, " ", sep="" )

    start <- runif( n=length(start), min=minRR, max=maxRR )
  }

  return( NULL )
}


## ROOT solving via the optim routine...
## How about trying it with the optim routine...
mroot <- function( f, ## The function
                   initial,  ## initial value
                   method="CG",
                   ... ) { ## Other parameters
  of <- function( x, ... ){
    res <- f( x, ... )
    return( sum( abs( res ) ) )
    #return( sum( res^2 ) )
  }

  ##print( "mroot" )

  res <- NULL
  try(  res <- optim( initial, of, gr=NULL, ..., method=method ), silent=TRUE  )

  return( res )
}

mrroot <- function( f, ## the function
                    initial, ## the initial value,
                    minRR=rep(-4,length(initial)), maxRR=rep(4,length(initial)), tryRR=5, ## bounds for random restart
                    method="Nelder-Mead",  ##"BFGS",
                    MAXITER=1000, TOL=sqrt(.Machine$double.eps),
                    verbose=FALSE,
                    ... ) {
  ## New, try to give some result
  best <- NULL

  if( verbose ) cat( "mrroot entered: " )
  for( t in 1:tryRR ) {
    res <- NULL
    try(  res <- mroot( f, initial, method, control=list(maxit=MAXITER,reltol=TOL), ... ) , silent=!verbose  )
    ##print( res )
    if( !is.null(res) && res$convergence==0 && abs(res$value)<0.001 ) {
      if( verbose ) {
        cat( "converged:" )
        print( res )
      }
      return( res )
    }

    ## failed -- was it the best so far?
    if( is.null(best) ) best <- res
    if( !is.null(res) && !is.null(best$value) && !is.na(best$value) && abs(res$value)<abs(best$value) )
      best <- res

    ## failed, try a random restart
    if( verbose ) {
      if( is.null(res) ) {
        cat( t, "(BAD-START) ", sep="" )  ## NOT debug
      }else{
        CONVSTR <- "?"
        if( res$convergence==1 ) {
          CONVSTR <- "MAXITER"
        }else if( res$convergence==10) {
          CONVSTR <- "NM-DEGEN"
        }else{
          try( CONVSTR <- res$message, silent=!verbose )
        }
        cat( t, "(", res$value, ",", CONVSTR, ") ", sep="" )
      }
    }
    initial <- runif( n=length(initial), min=minRR, max=maxRR )
  }

  if( verbose )
    cat( "FAILED\n" )

  msg <- paste( "Convergence failed -- should converge to zero; the best found converged to ", best$value, ".", sep="" )
  warning( msg )
  cat( msg, "\n", sep="" )
  return(best)

  #return( NULL )
}

## Improved version always tries a certain number of times
##  to make sure it doesn't converge to a better value...
mrroot2 <- function( f, ## the function
                     initial, ## the initial value,
                     minRR=rep(-4,length(initial)), maxRR=rep(4,length(initial)),  tryRR=5, minRRsuccess=1, ## bounds for random restart
                     method="Nelder-Mead",  ##"BFGS",
                     MAXITER=1000, TOL=sqrt(.Machine$double.eps),
                     verbose=FALSE,
                    ... ) {
  ## for storing the best result
  best <- NULL

  ## No successes
  successes <- 0

  for( t in 1:tryRR ) {
    ## Try to run the algorithm
    res <- NULL
    try(  res <- mroot( f, initial, method, control=list(maxit=MAXITER,reltol=TOL), ... ) , silent=!verbose  )

    ## Was anything returned?
    if( !is.null(res) ) {
      ##try( cat( "res$value", res$value, "best$value", best$value, "\n" ) )
      ## Was it the best so far?
      if( is.null(best) ) {
        ##cat( "INSTALLING BEST\n" )
        #print( res )
        best <- res
      }else if( abs(res$value) < abs(best$value) ) {
        ##cat( "REPLACING BEST\n" )
        #print( res )
        best <- res
      }

      ## Was it a success?
      if( res$convergence==0 && abs(res$value)<0.001 ) {
        successes <- successes + 1

        ## Have we met the minimum number of successes?
        if( successes>=minRRsuccess )
          return( best )
      }else{
        ##cat( "FAILED TO CONVERGE..." )
      }
    }

    ## And a random start
    initial <- runif( n=length(initial), min=minRR, max=maxRR )
  } ## t

  ## it didn't quite work the way that we wanted it to...
  msg <- NULL
  if( successes > 0 ) {
    msg <- paste( "Failed to converge required number of times for nuisance parameter (ran", successes, " instead of ", minRRsuccess, " times), you may want to run the routine again to ensure results are stable.", sep="" )
  }else if( is.null(best) ) {
    msg <- paste( "Failed to converge, or even start solving for the nuisance parameters. Is there any chance you are trying this on a case with a very, very small number of families? Otherwise, this is a very bad sign." )
  }else{
    msg <- paste( "Failed to converge for nuisance parameter -- the best solution should be approximately zero, but yields ", best$value, " instead." )
  }
  print( msg )
  warning( msg )

  ## and still return the best
  return( best )
}

mroptim <- function( f,
                     initial,
                     minRR=rep(-4,length(initial)), maxRR=rep(4,length(initial)), tryRR=5, ## bounds and # of times to try for random restart
                     method="Nelder-Mead",
                     verbose=FALSE,
                     ... ) {
  for( t in 1:tryRR ) {
    res <- optim( fn=f, par=initial, gr=NULL, ..., method=method )
    if( res$convergence==0 ) {
      if( verbose ) {
        cat( "converged:" )
        print( res )
      }
      return( res )
    }

    ## failed, try a random restart
    if( verbose ) cat( t, "" )
    initial <- runif( n=length(initial), min=minRR, max=maxRR )
  }
  if( verbose )
    cat( "FAILED\n" )
  return( NULL ) ## failed!
}

## Remove the .a and .b from markers...
fixMarkerNames <- function( names )
{
  pop <- function(x) return( x[1:(length(x)-1)] )
  ## pops off the .a/.b
  for( i in 1:length(names) )
    names[i] <- paste( pop(unlist(strsplit(names[i],"\\."))), collapse="." )
  return(names)
}

haplotypeDensity <- function( data, markerCol, traitCol, tempPrefix, FBAT="/tmp/th/fbat_2.02c" ) {
  usingAffectionStatus <- ( traitCol <= 6 )

  ## Argh... inverse
  markerNames <- unique( fixMarkerNames( names(data)[markerCol] ) )
  traitName <- names(data)[traitCol]
  #cat( "markerNames, traitName\n" )
  #print( markerNames )
  #print( traitName )
  #stop()

  ## Are any of the children missing? Need complete data here...
  ## NO -- we really don't! How does that affect this?
  ## TRY moving the badChildren piece until later
  badChildren <- rep(FALSE, nrow(data))
  ##for( mi in 1:length(markerCol) )
  ##  badChildren <- badChildren | data[[markerCol[mi]]]==0  ## AHHH markerCol[mi], not mi!!!
  ##badChildren <- badChildren & data$idfath!=0 & data$idmoth!=0 ## Not killing parents
  ## Subset out markers of interest, and kill badChildren
  ped <- data[ !badChildren, c(1:6,markerCol) ]
  class(ped) <- c( "ped", "data.frame" )
  phe <- NULL
  if( !usingAffectionStatus ) {
    phe <- data[ !badChildren, c(1:2,traitCol) ]
    class(phe) <- c( "phe", "data.frame" )
  }

  ## write out data to disk
  pedFilename <- paste( tempPrefix, "_fbat.ped", sep="" )
  pheFilename <- paste( tempPrefix, "_fbat.phe", sep="" )
  write.ped( pedFilename, ped )
  if( !usingAffectionStatus )
    write.phe( pheFilename, phe )

  ## Create commands for controlling FBAT
  ##markerStr <- paste( "_m", paste( markerCol, collapse="m" ), sep="" )
  markerStr <- paste( "_m" )
  logfile <- paste( tempPrefix, markerStr, ".log", sep="" )
  #stop( logfile )
  fbatCommands <- NULL
  traitCommands <- NULL
  if( !usingAffectionStatus )
    traitCommands <- c( paste( "load", pheFilename ),
                        paste( "trait", traitName ) )
  fbatCommands <- c( paste( "load", pedFilename ),
                     traitCommands,
##                     paste( "maxcmh", "100000" ),
                     paste( "log", logfile ),
                     paste( "viewhap", paste( markerNames, collapse=" " ) ),
                     "log off",
                     "quit" )

  ## Get our FBAT on
  progContFilename <- paste( tempPrefix, markerStr, "_systemControl.sh", sep="")
  programControl( program=FBAT, commands=fbatCommands, filename=progContFilename, intern=FALSE, output="/dev/null" )

  ## Load in the output
  reference <- condGeneFBATControl_load( logfile )

  ## Pre-linking in the trait
  ##print( head( ped, n=20 ) )
  for( mi in 1:length(markerCol) )
    badChildren <- badChildren | data[[ markerCol[mi] ]]==0
  badChildren <- badChildren & data$idfath!=0 & data$idmoth!=0 ## Not killing parents
  ped <- ped[ !badChildren, ]
  if( !usingAffectionStatus ) {
    phe <- phe[ !badChildren, ]
  }

  ##print( head( ped, n=20 ) )
  ##print( head( phe, n=20 ) )


  ## Link in the trait
  children <- ped$idfath!=0 & ped$idmoth!=0
  traitData <- NULL
  if( usingAffectionStatus ) {
    traitData <- ped[[traitCol]]
    traitData[ traitData == 0 ] <- NA
    traitData <- traitData - 1
  }else{
    traitData <- phe[[3]]
  }

  ##print( cbind( ped[children,1], traitData[children] ) )

  condGeneFBATControl_linkTrait( reference, ped[children,1], traitData[children] )

  ## And finally return the reference!
  return( reference )
}

resolveTraitType <- function( phe=NULL, trait="AffectionStatus" ) {
  ## Least resource intensive to do this special case
  if( trait=="AffectionStatus" )
    return( "binary" )

  ## Handle trait info
  traitCol <- which( names(phe) == trait )
  if( length(unique(phe[[traitCol]]))<=2 )
    return( "binary" )
  ##traitType <- "continuous"
  return( "continuous" )
}

## This is exported, so we should make sure to be a little more robust
condGene <- function( ped=NULL, phe=NULL, data=mergePhePed( ped, phe ),
                      trait="AffectionStatus", traitType="auto",
                      markerAnalyze=NULL, markerCondition=NULL,
                      sigma=1, alpha=NULL,
                      ignoreBtX=TRUE,
                      tempPrefix="temp",
                      removeUnphased=FALSE,
                      verbose=FALSE ) { ## essentially debug mode
                      ##RETURN_BC_ESTIMATE=FALSE, ## DEBUG ONLY!!! KILL ME
  ######################
  ## Precompute setup ##
  ######################

  solve.svd <- getFromNamespace( "solve.svd", "fbati" )

  ## if there wasn't a phenotype, then data is a little messed up... should fix this... -- or is this an R issue?
  if( is.null(data) )
    data <- ped

  ## We want the markerCondition to be the first two markers in marker
  if( is.null(markerAnalyze) || is.null(markerCondition) )
    stop( "markerCondition/markerAnalyze must be non-null." )

  ## Find the marker locations
  allMarkerNames       <- fixMarkerNames( names( ped ) )
  markerAnalyzeIndexTemp   <- match( markerAnalyze,   allMarkerNames ) ## Represents first index only - second is this +1
  markerConditionIndexTemp <- match( markerCondition, allMarkerNames )
  n <- length(markerAnalyzeIndexTemp)*2
  markerAnalyzeIndex <- rep(0,n)
  markerAnalyzeIndex[seq(from=1,to=n,by=2)] <- markerAnalyzeIndexTemp
  markerAnalyzeIndex[seq(from=2,to=n,by=2)] <- markerAnalyzeIndexTemp+1
  nC <- length(markerConditionIndexTemp)*2
  markerConditionIndex <- rep(0,nC)
  markerConditionIndex[seq(from=1,to=nC,by=2)] <- markerConditionIndexTemp
  markerConditionIndex[seq(from=2,to=nC,by=2)] <- markerConditionIndexTemp+1

  if( length(markerAnalyzeIndex)/2 != length(markerAnalyze) )
    stop( "Some markers in markerAnalyze do not exist." )
  if( length(markerConditionIndex)/2 != length(markerCondition) )
    stop( "Some markers in markerCondition do not exist." )

  ## Recode the pedigree alleles (in data) to all be 1/2
  #marker <- c(markerAnalyze, markerCondition)
  marker <- c( markerAnalyzeIndex, markerConditionIndex )
  for( i in seq( from=1, to=length(marker), by=2 ) ) {
    m1 <- marker[i]
    m2 <- marker[i+1]
    alleles <- sort(  setdiff( unique( c(data[[m1]], data[[m2]]) ), 0 )  ) ## sort is really only for debugging purposes -- won't change our simulated data then.
    ##alleles <- setdiff( unique( c(data[[m1]], data[[m2]]) ), 0 )
    if( length(alleles)>2 ) {
      stop( "biallelic loci only" )
    }else if( length(alleles) == 1 ) {
      data[[m1]][ data[[m1]]==alleles ] <- 1
      data[[m2]][ data[[m2]]==alleles ] <- 1
    }else if( length(alleles) == 2 ) {
      ## The case we expect most of the time
      a1 <- data[[m1]]==alleles[1]
      a2 <- data[[m1]]==alleles[2]
      data[[m1]][a1] <- 1
      data[[m1]][a2] <- 2
      a1 <- data[[m2]]==alleles[1]
      a2 <- data[[m2]]==alleles[2]
      data[[m2]][a1] <- 1
      data[[m2]][a2] <- 2
    }
  }

  if( verbose ) cat( "condGene:: Alleles have been recoded.\n" );

  ## nuclify merge the data
  data <- nuclifyMerged( data )

  if( verbose ) cat( "condGene:: nuclifyMerged\n" );

  ## Handle trait info
  traitCol <- trait
  if( is.character(trait) )
    traitCol <- which( names(data) == trait )
  if( length(traitCol) != 1 )
    stop( paste( "Trait '", trait, "' was not found.", sep="" ) )

  data_trait <- data[[traitCol]]

  if( traitType=="auto" ) {
    if( trait=="AffectionStatus" || length(unique(data_trait))<=2 ) {
      traitType <- "binary"
    }else{
      traitType <- "continuous"
    }

    if( trait=="AffectionStatus" ) {
      data_trait[ data_trait==0 ] <- NA
      data_trait <- data_trait - 1
      #print( data_trait )
    }
  }

  if( verbose ) cat( "Trait type: ", traitType, "\n", sep="" )

  #######################
  ## Computation piece ##
  #######################

  ## Compute the haplotype density
  if( verbose ) cat( "About to call upon FBAT.\n" )
  reference <- haplotypeDensity( data, c(markerAnalyzeIndex,markerConditionIndex), traitCol, tempPrefix )
  if( verbose ) cat( "Haplotype density read in from FBAT.\n" )
  ##condGeneFBATControl_print( reference ) ## DEBUG
  numInfFam <- condGeneFBATControl_numInfFam( reference )

  ## Solve for the nuisance parameters
  TOL <- .Machine$double.eps
  analyze_allele_index <- 1:length(markerAnalyze) - 1
  conditional_allele_index <- 1:length(markerCondition) + length(markerAnalyze) - 1

  umcMat <- NULL
  uuModel <- NULL
  if( traitType == "binary" ) {
    if( verbose ) cat( "Entered binary specific code.\n" )
    #print( analyze_allele_index )
    #print( conditional_allele_index )
    #stop()
    ##initial <- rep( 0, length(markerAnalyze) + 2*length(markerCondition) )
    initial <- rep( 0, 2*length(markerCondition) )
    bcSolve <- NULL
    bc <- NULL
    if( verbose ) cat( "About to try multiroot.\n" )
    try({
      ##require( rootSolve )
      bcSolve <- multiroot( f=condGeneFBATControl_uimc_nuis, start=initial, maxiter=1000, rtol=TOL, atol=TOL, ctol=TOL, useFortran=FALSE, reference=reference, analyze_allele_index=analyze_allele_index, conditional_allele_index=conditional_allele_index, onlyComputeConditional=TRUE )
      bc <- bcSolve$root
    }, silent=TRUE )
    if( is.null(bcSolve) ) {
      if( verbose ) cat( "multiroot failed, falling back on slower convergence method.\n" )
      try({
        bcSolve <- mrroot( f=condGeneFBATControl_uimc_nuis, initial=initial, reference=reference, analyze_allele_index=analyze_allele_index, conditional_allele_index=conditional_allele_index, onlyComputeConditional=TRUE, ignoreBtX )
        bc <- bcSolve$par
      }, silent=TRUE)
    }
    ##print( bcSolve )
    ##print( bc )
    bc0 <- bc[ 1 : length(conditional_allele_index) ]
    bc1 <- bc[ (length(conditional_allele_index)+1) : length(bc) ]

    if( verbose ) {
      cat( "Nuisance parameters have been solved for:\n" )
      print( bc0 )
      print( bc1 )
      cat( "Number of inf fams used: ", condGeneFBATControl_numInfFam(reference), "\n" )
    }

    ## Model based, empirical variance
    umc <- condGeneFBATControl_uimc( reference,
                                     rep(0,length(analyze_allele_index)), bc0, bc1,
                                     analyze_allele_index, conditional_allele_index,
                                     FALSE )
    #print( umc )
    umcMat <- cbind( umc$analyze, umc$conditional0, umc$conditional1 )
    #print( head( umcMat ) )
    #print( umcMat )
    #stop()


    ## Model based, model variance
    uuModel <- condGeneFBATControl_imc( reference=reference,
                                        bm=rep(0,length(analyze_allele_index)), bc0=bc0, bc1=bc1,
                                        analyze_allele_index=analyze_allele_index, conditional_allele_index=conditional_allele_index )
  }else{ ## traitType=qtl
    if( verbose ) cat( "Entered qtl specific code.\n" )
    ## First center the qtl by the sample mean
    if( is.null( alpha ) ) {
      center <- condGeneFBATControl_centerTrait( reference )
      if( verbose ) cat( "QTL center", center, "\n" )
      alpha <- 0
    }

    initial <- rep( 0, length(conditional_allele_index)*2 )

    bc <- NULL
    bcSolve <- NULL

    try( {
      ##require( rootSolve )
      bcSolve <- multiroot( f=condGeneFBATControl_contsUimc_nuis, start=initial, maxiter=100, rtol=TOL, atol=TOL, ctol=TOL, useFortran=FALSE, reference=reference, alpha=alpha, sigma=sigma, analyze_allele_index=analyze_allele_index, conditional_allele_index=conditional_allele_index, onlyComputeConditional=TRUE )
      bc <- bcSolve$root
    }, silent=TRUE )
    ##cat( "multiroot bcSolve" ); print( bc ); bcSolve <- NULL;
    if( is.null( bcSolve ) ) {
      try( {
        bcSolve <- mrroot( f=condGeneFBATControl_contsUimc_nuis, initial=initial, reference=reference, alpha=alpha, sigma=sigma, analyze_allele_index=analyze_allele_index, conditional_allele_index=conditional_allele_index, onlyComputeConditional=TRUE )
        bc <- bcSolve$par
      }, silent=TRUE )
    }
    ##print( bcSolve )
    if( verbose ) {
      cat( "bc" )
      print( bc )
      cat( "Number of inf fams used: ", condGeneFBATControl_numInfFam(reference), "\n" )
    }
    ##stop()

    umcMat <- condGeneFBATControl_contsUimc( beta=bc,
                                             reference=reference,
                                             alpha=alpha, sigma=sigma,
                                             analyze_allele_index=analyze_allele_index, conditional_allele_index=conditional_allele_index,
                                             onlyComputeConditional=FALSE )

    ## Model based, model variance
    uuModel <- condGeneFBATControl_contsImc( beta=bc, reference=reference, alpha=alpha, sigma=sigma, analyze_allele_index=analyze_allele_index, conditional_allele_index=conditional_allele_index )
    ##cat( "Model variance\n" )
    ##print( uuModel )
  }

  #if( RETURN_BC_ESTIMATE )
  #  return( bc )
  #cat( "bc non-pairwise" ); print( bc );

  #print( "condGene umcMat" );
  #print( umcMat )

  ######################################
  ## Then compute the test statistics ##

  if( verbose ) cat( "About to compute the test statistic.\n" )

  nA <- length(analyze_allele_index)
  nC <- length(conditional_allele_index)
  n <- nA + 2*nC
  u <- rep( 0, n )
  uu <- matrix( 0, nrow=n, ncol=n )
  for( k1 in 1:n ) {
    u[k1] <- sum( umcMat[,k1], na.rm=TRUE )
    for( k2 in k1:n ) {
      uu[k1,k2] <- sum( umcMat[,k1] * umcMat[,k2], na.rm=TRUE )
      uu[k2,k1] <- uu[k1,k2]
    }
  }
  u[(nA+1):n] <- 0 ## ADDITION -- THEORETICALLY ZERO
  ##print( "u, uu" )
  ##print( u )
  ##print( uu )
  ##stop()
  pvalueEmp <- 1
  rankEmp <- NA
  try( {
    post <- solve.svd( svd(uu), u )
    stat <- u %*% post

    ## These pieces are potentially not as stable (may need to try using Monte-Carlo...)
    ##rankEmp <- attr( post, "rank" ) - nC*2
    uInv <- solve.svd( svd(uu) )
    uInvA <- uInv[1:length(analyze_allele_index), 1:length(analyze_allele_index)]
    #uInvA <- matrix( uInvA )
    uInvAInv <- solve.svd( svd(uInvA) )
    rankEmp <- attr(uInvAInv, "rank" )

    pvalueEmp <- pchisq( stat, df=rankEmp, lower.tail=FALSE )
  }, silent=TRUE )
  ##cat( "pvalueEmp", pvalueEmp, "\n" )

  ## Model based, model variance
  uu <- uuModel
  if( verbose ) {
    cat( "condGene u, uu\n" )
    print( u )
    print( uu )
  }
  pvalueModel <- 1
  rankModel <- NA
  try( {
    post <- solve.svd( svd(uu), u )
    stat <- u %*% post

    ##rankModel <- attr( post, "rank" ) - nC*2
    uInv <- solve.svd( svd(uu) )
    uInvA <- uInv[1:length(analyze_allele_index), 1:length(analyze_allele_index)]
    #uInvA <- matrix( uInvA )
    uInvAInv <- solve.svd( svd(uInvA) )
    rankModel <- attr(uInvAInv, "rank" )

    pvalueModel <- pchisq( stat, df=rankModel, lower.tail=FALSE )
  }, silent=TRUE )
  ##cat( "pvalueModel", pvalueModel, "\n" )
  ##cat( "rankModel", rankModel, "\n" )

  ## Model free, empirical variance
  umcMat <- condGeneFBATControl_robustStat( reference, analyze_allele_index, conditional_allele_index )
  #print( umcMat )
  n <- ncol(umcMat)
  u <- rep( 0, n )
  uu <- matrix( 0, nrow=n, ncol=n )
  for( k1 in 1:n ) {
    u[k1] <- sum( umcMat[,k1], na.rm=TRUE )
    for( k2 in k1:n ) {
      uu[k1,k2] <- sum( umcMat[,k1] * umcMat[,k2], na.rm=TRUE )
      uu[k2,k1] <- uu[k1,k2]
    }
  }
  #print( u )
  #print( uu )
  pvalueRobust <- 1
  rankRobust <- NA
  try( {
    post <- solve.svd( svd(uu), u )
    rankRobust <- attr( post, "rank" )
    stat <- u %*% post
    pvalueRobust <- pchisq( stat, df=rankRobust, lower.tail=FALSE )
  }, silent=TRUE )
  ##cat( "pvalueRobust", pvalueRobust, "\n" )

  ## Free the data!
  condGeneFBATControl_free( reference )

  return( data.frame( pvalueEmp=pvalueEmp, rankEmp=rankEmp,
                      pvalueModel=pvalueModel, rankModel=rankModel,
                      pvalueRobust=pvalueRobust, rankRobust=rankRobust,
                      numInfFam=numInfFam ) )
}



## Conditional gene test, but do it PAIRWISE!!!
## -- I thought previously that it would be less powerful, but now I'm convinced
##     that actually it should _always_ be more powerful? Very confusing,
##     more things will contribute to this one...
## -- So the idea then is to go ahead and get all of these coded in,
##     and then compare the two to see what happens
## -- Even if they are the same, this one is greatly preferred,
##     as otherwise reconstructing the haplotype density fails
##     miserably!
condGeneP <- function( ped=NULL, phe=NULL, data=mergePhePed( ped, phe ),
                      trait="AffectionStatus", traitType="auto",
                      markerAnalyze=NULL, markerCondition=NULL,
                      sigma=1, alpha=NULL,
                      ignoreBtX=TRUE,
                      tempPrefix="temp",
                      removeUnphased=FALSE,
                      verbose=FALSE ) { ## essentially debug mode
  ######################
  ## Precompute setup ##
  ######################
  #markerNames <- unique( fixMarkerNames( #names(data)[markerCol] ) )

  solve.svd <- getFromNamespace( "solve.svd", "fbati" )

  ## if there wasn't a phenotype, then data is a little messed up... should fix this... -- or is this an R issue?
  if( is.null(data) )
    data <- ped

  ## We want the markerCondition to be the first two markers in marker
  if( is.null(markerAnalyze) || is.null(markerCondition) )
    stop( "markerCondition/markerAnalyze must be non-null." )

  ## Find the marker locations
  allMarkerNames       <- fixMarkerNames( names( ped ) )
  markerAnalyzeIndexTemp   <- match( markerAnalyze,   allMarkerNames ) ## Represents first index only - second is this +1
  markerConditionIndexTemp <- match( markerCondition, allMarkerNames )
  n <- length(markerAnalyzeIndexTemp)*2
  markerAnalyzeIndex <- rep(0,n)
  markerAnalyzeIndex[seq(from=1,to=n,by=2)] <- markerAnalyzeIndexTemp
  markerAnalyzeIndex[seq(from=2,to=n,by=2)] <- markerAnalyzeIndexTemp+1
  nC <- length(markerConditionIndexTemp)*2
  markerConditionIndex <- rep(0,nC)
  markerConditionIndex[seq(from=1,to=nC,by=2)] <- markerConditionIndexTemp
  markerConditionIndex[seq(from=2,to=nC,by=2)] <- markerConditionIndexTemp+1

  if( length(markerAnalyzeIndex)/2 != length(markerAnalyze) )
    stop( "Some markers in markerAnalyze do not exist." )
  if( length(markerConditionIndex)/2 != length(markerCondition) )
    stop( "Some markers in markerCondition do not exist." )

  ## Recode the pedigree alleles (in data) to all be 1/2
  #marker <- c(markerAnalyze, markerCondition)
  marker <- c( markerAnalyzeIndex, markerConditionIndex )
  for( i in seq( from=1, to=length(marker), by=2 ) ) {
    m1 <- marker[i]
    m2 <- marker[i+1]
    alleles <- sort(  setdiff( unique( c(data[[m1]], data[[m2]]) ), 0 )  ) ## sort is really only for debugging purposes -- won't change our simulated data then.
    ##alleles <- setdiff( unique( c(data[[m1]], data[[m2]]) ), 0 )
    if( length(alleles)>2 ) {
      stop( "biallelic loci only" )
    }else if( length(alleles) == 1 ) {
      data[[m1]][ data[[m1]]==alleles ] <- 1
      data[[m2]][ data[[m2]]==alleles ] <- 1
    }else if( length(alleles) == 2 ) {
      ## The case we expect most of the time
      a1 <- data[[m1]]==alleles[1]
      a2 <- data[[m1]]==alleles[2]
      data[[m1]][a1] <- 1
      data[[m1]][a2] <- 2
      a1 <- data[[m2]]==alleles[1]
      a2 <- data[[m2]]==alleles[2]
      data[[m2]][a1] <- 1
      data[[m2]][a2] <- 2
    }
  }

  if( verbose ) cat( "condGeneP:: Alleles have been recoded.\n" );

  ## nuclify merge the data
  data <- nuclifyMerged( data )

  if( verbose ) cat( "condGeneP:: nuclifyMerged\n" );

  ## Handle trait info
  traitCol <- trait
  if( is.character(trait) )
    traitCol <- which( names(data) == trait )
  if( length(traitCol) != 1 )
    stop( paste( "Trait '", trait, "' was not found.", sep="" ) )

  data_trait <- data[[traitCol]]

  ##print( trait )

  if( traitType=="auto" ) {
    if( trait=="AffectionStatus" || length(unique(data_trait))<=2 ) {
      traitType <- "binary"
    }else{
      traitType <- "continuous"
    }

    if( trait=="AffectionStatus" ) {
      data_trait[ data_trait==0 ] <- NA
      data_trait <- data_trait - 1
      #print( data_trait )
    }
  }

  if( verbose ) cat( "Trait type: ", traitType, "\n", sep="" )

  #######################
  ## Computation piece ##
  #######################

  ################################################
  ## First solve for the nuisance parameters... ##

  ## Compute the haplotype density
  conditional_allele_index <- markerConditionIndex
  reference <- haplotypeDensity( data, markerConditionIndex, traitCol, tempPrefix )
  ##condGeneFBATControl_print( reference ) ## DEBUG
  ##print( condGeneFBATControl_pids( reference ) ) ## DEBUG
  n <- condGeneFBATControl_numFam( reference )

  ## Now solve for the nuisance parameters
  bc <- NULL
  analyze_allele_index <- c()
  conditional_allele_index <- 1:length(markerCondition) - 1
  TOL <- .Machine$double.eps
  umcNuis <- NULL
  CENTER_QTL <- FALSE
  if( traitType == "binary" ) {
    initial <- rep( 0, 2*length(markerCondition) )
    bcSolve <- NULL
    try({
      ##require( rootSolve )
      bcSolve <- multiroot( f=condGeneFBATControl_uimc_nuis, start=initial, maxiter=100, rtol=TOL, atol=TOL, ctol=TOL, useFortran=FALSE, reference=reference, analyze_allele_index=analyze_allele_index, conditional_allele_index=conditional_allele_index, onlyComputeConditional=TRUE ) ## onlyComputeConditional irrelevant here...
      bc <- bcSolve$root
    }, silent=TRUE )
    if( is.null(bcSolve) ) {
      if( verbose ) cat( "multiroot failed, falling back on slower convergence method.\n" )
      try( {
        bcSolve <- mrroot( f=condGeneFBATControl_uimc_nuis, initial=initial, reference=reference, analyze_allele_index=analyze_allele_index, conditional_allele_index=conditional_allele_index, onlyComputeConditional=TRUE )
        bc <- bcSolve$par
      }, silent=TRUE );
      if( is.null( bcSolve ) )
        return( data.frame( pvalue=1, rank=0,
                            pvalueR=1, rankR=0 ) )
    }
    #cat( "bc before" ); print( bc );
    #bc <- condGene( ped=ped, phe=phe,
    #                trait=trait, traitType=traitType,
    #                markerAnalyze=markerAnalyze, markerCondition=markerCondition,
    #                sigma=sigma, alpha=alpha,
    #                ignoreBtX=ignoreBtX,
    #                tempPrefix=tempPrefix,
    #                removeUnphased=removeUnphased,
    #                RETURN_BC_ESTIMATE=TRUE,
    #                verbose=verbose ) ## DEBUG ONLY -- KILL ME
    #cat( "bc after" ); print( bc );
    umcNuis <- condGeneFBATControl_uimc( reference,
                                         c(), bc[1:length(conditional_allele_index)], bc[(length(conditional_allele_index)+1):length(bc)],
                                         analyze_allele_index, conditional_allele_index,
                                         FALSE )   ## $conditional0, $conditional1
    umcNuis <- cbind( umcNuis$conditional0, umcNuis$conditional1 )
  }else{ ## traitType == qtl
    ## First center the qtl by the sample mean
    if( is.null( alpha ) ) {
      center <- condGeneFBATControl_centerTrait( reference )
      if( verbose ) cat( "QTL center", center, "\n" )
      alpha <- 0
      CENTER_QTL <- TRUE
    }

    ## Then solve for the nuisance parameter bc
    initial <- rep( 0, 2*length(conditional_allele_index) )
    bcSolve <- NULL
    try( {
      ##require( rootSolve )
      bcSolve <- multiroot( f=condGeneFBATControl_contsUimc_nuis, start=initial, maxiter=100, rtol=TOL, atol=TOL, ctol=TOL, useFortran=FALSE, reference=reference, alpha=alpha, sigma=sigma, analyze_allele_index=analyze_allele_index, conditional_allele_index=conditional_allele_index, onlyComputeConditional=TRUE, ignoreBtX=ignoreBtX )
      bc <- bcSolve$root
    }, silent=TRUE )
    ##cat( "multiroot bcSolve" ); print( bc ); bcSolve <- NULL;
    if( is.null( bcSolve ) ) {
      if( verbose ) cat( "multiroot failed, falling back on slower convergence method.\n" )
      try( {
        bcSolve <- mrroot( f=condGeneFBATControl_contsUimc_nuis, initial=initial, reference=reference, alpha=alpha, sigma=sigma, analyze_allele_index=analyze_allele_index, conditional_allele_index=conditional_allele_index, onlyComputeConditional=TRUE, ignoreBtX=ignoreBtX )
        bc <- bcSolve$par
      }, silent=TRUE )
      if( is.null( bcSolve ) )
        return( data.frame( pvalue=1, rank=0,
                            pvalueR=1, rankR=0 ) )
    }

    ## only needed for the empirical variance
    umcNuis <- condGeneFBATControl_contsUimc( beta=bc,
                                              reference=reference,
                                              alpha=alpha, sigma=sigma,
                                              analyze_allele_index=analyze_allele_index, conditional_allele_index=conditional_allele_index,
                                              onlyComputeConditional=FALSE )
  }
  if( verbose ) {
    cat( "bc" )
    print( bc )
    cat( "Number of inf fams used: ", condGeneFBATControl_numInfFam(reference), "\n" )
  }
  condGeneFBATControl_free( reference ) ## Free the data here

  #############################################
  ## Now handle each of the pairwise fits... ##
  #analyze_allele_index <- 1:length(markerAnalyze) - 1
  #conditional_allele_index <- 1:length(markerCondition) + length(markerAnalyze) - 1
  analyze_allele_index <- 0 ## pairwise, so always moving across this marker...
  conditional_allele_index <- 1:length(markerCondition) + length(analyze_allele_index) - 1

  ##R <- length(analyze_allele_index) + 2*length(conditional_allele_index)
  R <- length(markerAnalyze) + 2*length(markerCondition)
  A <- length(markerAnalyze)
  umcMat <- matrix( 0, nrow=n, ncol=R )
  umcMat[, (R-ncol(umcNuis)+1):R ] <- umcNuis

  umcRMat <- matrix( 0, nrow=n, ncol=length(markerAnalyze) )

  for( i in 1:length(markerAnalyze) ) {
    ## Pull out the pairwise pedigree
    #pairwisePed <- data[ , c(  1:6, markerAnalyzeIndex[i*2-c(1,0)], markerConditionIndex  ) ]
    #print( head( pairwisePed ) )

    ## Get the haplotype distribution
    reference <- haplotypeDensity( data, c(markerAnalyzeIndex[i*2-c(1,0)], markerConditionIndex), traitCol, tempPrefix )

    ##condGeneFBATControl_print( reference ) ## DEBUG
    ##print( condGeneFBATControl_pids( reference ) ) ## DEBUG

    ## Compute the pieces of the statistic
    if( traitType=="binary" ) {
      bc0 <- bc[1:length(conditional_allele_index)]
      bc1 <- bc[(length(conditional_allele_index)+1):length(bc)]
      umcAddi <- condGeneFBATControl_uimc( reference,
                                           rep(0,length(analyze_allele_index)), bc0, bc1,
                                           analyze_allele_index, conditional_allele_index,
                                           FALSE )
      #print( str( umcAddi ) )
      #print( head( umcMat ) )
      #stop()
      umcMat[ , i ] <- umcAddi$analyze
    }else{ ## qtl
      ## First center the qtl by the sample mean
      if( CENTER_QTL ) {
        center <- condGeneFBATControl_centerTrait( reference )
        if( verbose ) cat( "QTL center (analyze parameter", i, ")", center, "\n" )
      }

      umcAddi <- condGeneFBATControl_contsUimc( beta=bc,
                                                reference=reference,
                                                alpha=alpha, sigma=sigma,
                                                analyze_allele_index=analyze_allele_index, conditional_allele_index=conditional_allele_index,
                                                onlyComputeConditional=FALSE )
      umcMat[ , i ] <- umcAddi[ , 1 ]
    }

    ##print( umcMat )
    ##stop()

    ## And the robust statistic piece
    umcRMat[ , i ] <- condGeneFBATControl_robustStat( reference, analyze_allele_index, conditional_allele_index )

    ## And finally free the memory
    condGeneFBATControl_free( reference )
  } ## i

  #print( "condGeneP umcMat" )
  #print( umcMat )

  ################################
  ## Compute the test statistic ##
  if( verbose ) cat( "About to compute the test statistic.\n" )

  ## Model-based test using empirical variance
  u <- rep( 0, R )
  uu <- matrix( 0, nrow=R, ncol=R )
  for( k1 in 1:R ) {
    u[k1] <- sum( umcMat[,k1], na.rm=TRUE )
    for( k2 in k1:R ) {
      uu[k1,k2] <- sum( umcMat[,k1] * umcMat[,k2], na.rm=TRUE )
      uu[k2,k1] <- uu[k1,k2]
    }
  }
  u[(A+1):R] <- 0 ## ADDITION -- THEORETICALLY ZERO

  if( verbose ) {
    print( "condGeneP Pairwise u, uu" )
    print( u )
    print( uu )
  }

  pvalue <- 1 ## only computes empirical
  rank <- 1
  try( {
    ## Compute the test statistic
    post <- solve.svd( svd( uu ), u )
    stat <- u %*% post

    ## Semi-unstable computing the rank...
    uInv <- solve.svd( svd( uu ) )
    uInvA <- uInv[ 1:length(markerAnalyze), 1:length(markerAnalyze) ]
    uInvAInv <- solve.svd( svd( uInvA ) )
    rank <- attr( uInvAInv, "rank" )

    ## and compute the pvalue
    pvalue <- pchisq( stat, df=rank, lower.tail=FALSE )
  }, silent=TRUE )

  ## Robust test statistic (empirical variance)
  na <- length(markerAnalyze)
  u <- rep( 0, na )
  uu <- matrix( 0, nrow=na, ncol=na )
  for( k1 in 1:na ) {
    u[k1] <- sum( umcRMat[,k1], na.rm=TRUE )
    for( k2 in k1:na ) {
      uu[k1,k2] <- sum( umcRMat[,k1] * umcRMat[,k2], na.rm=TRUE )
      uu[k2,k1] <- uu[k1,k2]
    }
  }

  pvalueR <- 1
  rankR <- 1
  try( {
    ## Compute the test statistic
    post <- solve.svd( svd(uu), u )
    stat <- u %*% post
    ##cat( "stat", stat, "\n" )

    ## Get the rank (should be stable)
    rankR <- attr( post, "rank" )

    ## and compute the pvalue
    pvalueR <- pchisq( stat, df=rankR, lower.tail=FALSE )
  }, silent=TRUE )

  ##stop()

  #############
  ## Return! ##
  return( data.frame( pvalue=pvalue, rank=rank,
                      pvalueR=pvalueR, rankR=rankR ) )
}






## NEW: alpha is the offset for a dichotomous trait
##      for the robust test _only_
## Redo the pairwise, so that it truly is pairwise; reestimate the nuisance parameters for everything
condGeneP2 <- function( ped=NULL, phe=NULL, data=mergePhePed( ped, phe ),
                        trait="AffectionStatus", traitType="auto",
                        markerAnalyze=NULL, markerCondition=NULL,
                        sigma=1, alpha=NULL,
                        ignoreBtX=TRUE,
                        tempPrefix="temp",
                        removeUnphased=FALSE,
                        MAXITER=1000, TOL=sqrt(.Machine$double.eps),
                        verbose=FALSE,
                        FBAT="/tmp/th/fbat_2.02c" ) { ## essentially debug mode
  ######################
  ## Precompute setup ##
  ######################

  solve.svd <- getFromNamespace( "solve.svd", "fbati" )

  ## if there wasn't a phenotype, then data is a little messed up... should fix this... -- or is this an R issue?
  if( is.null(data) )
    data <- ped

  ## We want the markerCondition to be the first two markers in marker
  if( is.null(markerAnalyze) || is.null(markerCondition) )
    stop( "markerCondition/markerAnalyze must be non-null." )

  ## Find the marker locations
  allMarkerNames       <- fixMarkerNames( names( ped ) )
  markerAnalyzeIndexTemp   <- match( markerAnalyze,   allMarkerNames ) ## Represents first index only - second is this +1
  markerConditionIndexTemp <- match( markerCondition, allMarkerNames )
  n <- length(markerAnalyzeIndexTemp)*2
  markerAnalyzeIndex <- rep(0,n)
  markerAnalyzeIndex[seq(from=1,to=n,by=2)] <- markerAnalyzeIndexTemp
  markerAnalyzeIndex[seq(from=2,to=n,by=2)] <- markerAnalyzeIndexTemp+1
  nC <- length(markerConditionIndexTemp)*2
  markerConditionIndex <- rep(0,nC)
  markerConditionIndex[seq(from=1,to=nC,by=2)] <- markerConditionIndexTemp
  markerConditionIndex[seq(from=2,to=nC,by=2)] <- markerConditionIndexTemp+1

  if( length(markerAnalyzeIndex)/2 != length(markerAnalyze) )
    stop( "Some markers in markerAnalyze do not exist." )
  if( length(markerConditionIndex)/2 != length(markerCondition) )
    stop( "Some markers in markerCondition do not exist." )

  ## Recode the pedigree alleles (in data) to all be 1/2
  #marker <- c(markerAnalyze, markerCondition)
  marker <- c( markerAnalyzeIndex, markerConditionIndex )
  for( i in seq( from=1, to=length(marker), by=2 ) ) {
    m1 <- marker[i]
    m2 <- marker[i+1]
    alleles <- sort(  setdiff( unique( c(data[[m1]], data[[m2]]) ), 0 )  ) ## sort is really only for debugging purposes -- won't change our simulated data then.
    ##alleles <- setdiff( unique( c(data[[m1]], data[[m2]]) ), 0 )
    if( length(alleles)>2 ) {
      stop( "biallelic loci only" )
    }else if( length(alleles) == 1 ) {
      data[[m1]][ data[[m1]]==alleles ] <- 1
      data[[m2]][ data[[m2]]==alleles ] <- 1
    }else if( length(alleles) == 2 ) {
      ## The case we expect most of the time
      a1 <- data[[m1]]==alleles[1]
      a2 <- data[[m1]]==alleles[2]
      data[[m1]][a1] <- 1
      data[[m1]][a2] <- 2
      a1 <- data[[m2]]==alleles[1]
      a2 <- data[[m2]]==alleles[2]
      data[[m2]][a1] <- 1
      data[[m2]][a2] <- 2
    }
  }

  if( verbose ) cat( "condGeneP:: Alleles have been recoded.\n" );

  ## nuclify merge the data
  data <- nuclifyMerged( data )

  if( verbose ) cat( "condGeneP:: nuclifyMerged\n" );

  ## Handle trait info
  traitCol <- trait
  if( is.character(trait) )
    traitCol <- which( names(data) == trait )
  if( length(traitCol) != 1 )
    stop( paste( "Trait '", trait, "' was not found.", sep="" ) )

  data_trait <- data[[traitCol]]

  ##print( trait )

  if( traitType=="auto" ) {
    if( trait=="AffectionStatus" || length(unique(data_trait))<=2 ) {
      traitType <- "binary"
    }else{
      traitType <- "continuous"
    }

    if( trait=="AffectionStatus" ) {
      data_trait[ data_trait==0 ] <- NA
      data_trait <- data_trait - 1
      #print( data_trait )
    }
  }

  if( verbose ) cat( "Trait type: ", traitType, "\n", sep="" )

  #######################
  ## Computation piece ##
  #######################



  #############################################
  ## Now handle each of the pairwise fits... ##
  #analyze_allele_index <- 1:length(markerAnalyze) - 1
  #conditional_allele_index <- 1:length(markerCondition) + length(markerAnalyze) - 1
  analyze_allele_index <- 0 ## pairwise, so always moving across this marker...
  conditional_allele_index <- 1:length(markerCondition) + length(analyze_allele_index) - 1

  ##R <- length(analyze_allele_index) + 2*length(conditional_allele_index)
  R <- length(markerAnalyze) + 2*length(markerCondition)*length(markerAnalyze)
  A <- length(markerAnalyze)
  umcMat <- NULL
  umcRMat <- NULL

  CONVERGED <- TRUE
  numInformative <- c()
  for( i in 1:length(markerAnalyze) ) {
    ## Pull out the pairwise pedigree
    #pairwisePed <- data[ , c(  1:6, markerAnalyzeIndex[i*2-c(1,0)], markerConditionIndex  ) ]
    #print( head( pairwisePed ) )

    ## Get the haplotype distribution
    reference <- haplotypeDensity( data, c(markerAnalyzeIndex[i*2-c(1,0)], markerConditionIndex), traitCol, tempPrefix, FBAT=FBAT )

    ##condGeneFBATControl_print( reference ) ## DEBUG
    ##print( condGeneFBATControl_pids( reference ) ) ## DEBUG



    ## Compute the pieces of the statistic
    bc <- NULL
    TOL <- .Machine$double.eps
    umcNuis <- NULL
    ##CENTER_QTL <- FALSE  ## CODETOOLS CLAIM TO FAME
    if( traitType == "binary" ) {
      initial <- rep( 0, 2*length(markerCondition) )
      bcSolve <- NULL
      #try({
      #  require( rootSolve )
      #  bcSolve <- multiroot( f=condGeneFBATControl_uimc_nuis, start=initial, maxiter=100, rtol=TOL, atol=TOL, ctol=TOL, useFortran=FALSE, reference=reference, analyze_allele_index=analyze_allele_index, conditional_allele_index=conditional_allele_index, onlyComputeConditional=TRUE ) ## onlyComputeConditional irrelevant here...
      #  bc <- bcSolve$root
      #}, silent=TRUE )
      ##cat( "multiroot bcSolve" ); print( bc ); bcSolve <- NULL; ## DEBUG ONLY
      if( is.null(bcSolve) ) {
        if( verbose ) cat( "multiroot failed, falling back on slower convergence method.\n" )
        ## Used to be mrroot, changed to mrroot2 for stability (hopefully) of some results
        try( {
          bcSolve <- mrroot2( f=condGeneFBATControl_uimc_nuis, initial=initial, reference=reference, analyze_allele_index=analyze_allele_index, conditional_allele_index=conditional_allele_index, onlyComputeConditional=TRUE, MAXITER=MAXITER, TOL=TOL, verbose=verbose )
          bc <- bcSolve$par
        }, silent=!TRUE );
	##cat( "mrroot bcSolve" ); print( bc ); ## DEBUG ONLY
        if( is.null( bcSolve ) ) {
          #msg <- paste( "Nuisance parameter did not converge for model-based test, analysis marker ", markerAnalyze[i], ".", sep="" )
          cat( msg, "\n" )
          warning( msg )
          CONVERGED <- FALSE ## Results following will need to be erased
          bc <- rep(0,length(conditional_allele_index)+length(analyze_allele_index))
        }
      }
      umcNuis <- condGeneFBATControl_uimc(reference,
                                          rep(0,length(analyze_allele_index)), bc[1:length(conditional_allele_index)], bc[(length(conditional_allele_index)+1):length(bc)],
                                          analyze_allele_index, conditional_allele_index,
                                          FALSE )   ## $conditional0, $conditional1
      umcNuis <- cbind( umcNuis$analyze, umcNuis$conditional1, umcNuis$conditional1 )
      #print( umcNuis )
      #stop()
    }else{ ## traitType == qtl
      ## First center the qtl by the sample mean
      if( is.null( alpha ) ) {
        center <- condGeneFBATControl_centerTrait( reference )
        if( verbose ) cat( "QTL center", center, "\n" )
        #alpha <- 0
        #CENTER_QTL <- TRUE
      }

      ## Then solve for the nuisance parameter bc
      initial <- rep( 0, 2*length(conditional_allele_index) )
      bcSolve <- NULL
      #try( {
      #  require( rootSolve )
      #  bcSolve <- multiroot( f=condGeneFBATControl_contsUimc_nuis, start=initial, maxiter=100, rtol=TOL, atol=TOL, ctol=TOL, useFortran=FALSE, reference=reference, alpha=alpha, sigma=sigma, analyze_allele_index=analyze_allele_index, conditional_allele_index=conditional_allele_index, onlyComputeConditional=TRUE, ignoreBtX=ignoreBtX )
      #  bc <- bcSolve$root
      #}, silent=TRUE )
      #cat( "multiroot bcSolve" ); print( bc ); bcSolve <- NULL;   ## DEBUG ONLY
      if( is.null( bcSolve ) ) {
        #if( verbose ) cat( "multiroot failed, falling back on slower convergence method.\n" )
        ## Used to be mrroot
        try( {
          bcSolve <- mrroot2( f=condGeneFBATControl_contsUimc_nuis, initial=initial, reference=reference, alpha=alpha, sigma=sigma, analyze_allele_index=analyze_allele_index, conditional_allele_index=conditional_allele_index, onlyComputeConditional=TRUE, ignoreBtX=ignoreBtX, MAXITER=MAXITER, TOL=TOL, verbose=verbose )
          bc <- bcSolve$par
        }, silent=TRUE )
        ##cat( "mrroot bcSolve" ); print( bc ); ## DEBUG ONLY
        if( is.null( bcSolve ) ) {
          msg <- paste( "Nuisance parameter did not converge for model-based test, analysis marker ", markerAnalyze[i], ".", sep="" )
          cat( msg, "\n" )
          warning( msg )
          CONVERGED <- FALSE ## Results following will need to be erased
          bc <- rep(0,length(conditional_allele_index)+length(analyze_allele_index))
        }
      }

      ##cat( "mrroot bcSolve" ); print( bc ); ## DEBUG ONLY

      ## only needed for the empirical variance
      umcNuis <- condGeneFBATControl_contsUimc( beta=bc,
                                                reference=reference,
                                                alpha=alpha, sigma=sigma,
                                                analyze_allele_index=analyze_allele_index, conditional_allele_index=conditional_allele_index,
                                                onlyComputeConditional=FALSE )
    }
    if( verbose ) {
      cat( "bc" )
      print( bc )
      cat( "Number of inf fams used: ", condGeneFBATControl_numInfFam(reference), "\n" )
    }
    #cat( "bc " ); print( bc ); ## DEBUG ONLY

    ## potentially needs to be sized...
    numFam <- condGeneFBATControl_numFam( reference )
    if( is.null(umcMat) ) {
      umcMat <- matrix( 0, nrow=numFam, ncol=R )
      umcRMat <- matrix( 0, nrow=numFam, ncol=A )
    }


    #print( dim( umcNuis ) )
    #print( dim( umcMat ) )
    for( col in 1:ncol(umcNuis) ) {
      #cat( "col", col, "\n" )
      #cat( "i+(col-1)*A", i+(col-1)*A, "\n" )
      umcMat[ , i+(col-1)*A ] <- umcNuis[ , col ]
    }

    ## NEW ADDITION FOR THE ROBUST PIECE,
    ##  If dichotomous, see if there is an offset
    if( traitType=="binary" && !is.null(alpha) && alpha!=0 ) {
      #print( "Get your offset on!" ) ## DEBUG ONLY
      #print( alpha )
      ## Ensure that the model-based test doesn't change, just for debug externally

      ## Subtract the offset from the trait
      condGeneFBATControl_centerTrait( reference=reference, center=alpha, mean=FALSE )
    }
    #condGeneFBATControl_print( reference ) ## DEBUG
    #stop()

    ## And the robust statistic piece
    umcRMat[ , i ] <- condGeneFBATControl_robustStat( reference, analyze_allele_index, conditional_allele_index )

    ## NEW! Get the number informative
    numInformative <- c( numInformative, condGeneFBATControl_numInfFam( reference ) )

    ## And finally free the memory
    condGeneFBATControl_free( reference )
  } ## i

  #print( "condGeneP2 umcMat" )
  #print( umcMat )
  #xstop()

  ################################
  ## Compute the test statistic ##
  if( verbose ) cat( "About to compute the test statistic.\n" )

  ## Model-based test using empirical variance
  u <- rep( 0, R )
  uu <- matrix( 0, nrow=R, ncol=R )
  for( k1 in 1:R ) {
    u[k1] <- sum( umcMat[,k1], na.rm=TRUE )
    for( k2 in k1:R ) {
      uu[k1,k2] <- sum( umcMat[,k1] * umcMat[,k2], na.rm=TRUE )
      uu[k2,k1] <- uu[k1,k2]
    }
  }
  u[(A+1):R] <- 0 ## ADDITION -- THEORETICALLY ZERO -- SEMI-QUESTIONABLE???

  numInf <- rep( 0, A )
  for( k1 in 1:A )
    numInf[k1] <- sum( abs(umcMat[,k1]) >= sqrt(.Machine$double.eps), na.rm=TRUE )

  if( verbose ) {
    print( "condGeneP Pairwise u, uu" )
    print( u )
    print( uu )
  }

  #cat( "u\n" ); print( u ); ## DEBUG ONLY
  #cat( "uu\n" ); print( uu ); ## DEBUG ONLY
  #cat( "svd(uu)\n" ); print( solve.svd( svd( uu ) ) ); ## DEBUG ONLY


  pvalue <- 1 ## only computes empirical
  rank <- 1
  try( {
    ## Compute the test statistic
    post <- solve.svd( svd( uu ), u )
    stat <- u %*% post
    #cat( "stat before", stat )  ## DEBUG ONLY
    #require( MASS )
    #stat <- u %*% ginv( uu ) %*% u
    #cat( "stat after", stat, "\n" )  ## DEBUG ONLY

    ## Semi-unstable computing the rank...
    uInv <- solve.svd( svd( uu ) )
    uInvA <- uInv[ 1:length(markerAnalyze), 1:length(markerAnalyze) ]
    uInvAInv <- solve.svd( svd( uInvA ) )
    rank <- attr( uInvAInv, "rank" )

    ## and compute the pvalue
    pvalue <- pchisq( stat, df=rank, lower.tail=FALSE )
  }, silent=!TRUE )

  ## Robust test statistic (empirical variance)
  na <- length(markerAnalyze)
  u <- rep( 0, na )
  uu <- matrix( 0, nrow=na, ncol=na )
  for( k1 in 1:na ) {
    u[k1] <- sum( umcRMat[,k1], na.rm=TRUE )
    for( k2 in k1:na ) {
      uu[k1,k2] <- sum( umcRMat[,k1] * umcRMat[,k2], na.rm=TRUE )
      uu[k2,k1] <- uu[k1,k2]
    }
  }

  numInfR <- rep( 0, A )
  for( k1 in 1:A )
    numInfR[k1] <- sum( abs(umcRMat[,k1]) >= sqrt(.Machine$double.eps), na.rm=TRUE )

  paste( numInfR )

  pvalueR <- 1
  rankR <- 1
  try( {
    ## Compute the test statistic
    post <- solve.svd( svd(uu), u )
    stat <- u %*% post
    ##cat( "stat", stat, "\n" )

    ## Get the rank (should be stable)
    rankR <- attr( post, "rank" )

    ## and compute the pvalue
    pvalueR <- pchisq( stat, df=rankR, lower.tail=FALSE )
  }, silent=TRUE )

  ##stop()

  if( !CONVERGED ) {
    pvalue <- 1
    rank <- 0
  }

  #############
  ## Return! ##
  if( traitType=="continuous" ) {
    offset <- NA
  }else{
    offset <- alpha
  }
  if( is.null(offset) ) offset <- NA
  return( data.frame( trait=trait, traitType=traitType,
                      offset=offset, ignoreBX=ignoreBtX,
#                      numInf=paste(numInformative,collapse=","),
                      pvalue=pvalue, rank=rank, numInf=paste(numInf,collapse=","),
                      pvalueR=pvalueR, rankR=rankR, numInfR=paste(numInfR,collapse=",")  ) )
##  return( data.frame( pvalue=pvalue, rank=rank,
##                      pvalueR=pvalueR, rankR=rankR ) )
}






## Newest -- knockoff of Stijn's paper...
condGeneP3 <- function( ped=NULL, phe=NULL, data=mergePhePed( ped, phe ),
                        trait="AffectionStatus", traitType="auto",
                        markerAnalyze=NULL, markerCondition=NULL,
                        sigma=1, alpha=NULL,
                        ignoreBtX=TRUE,
                        tempPrefix="temp",
                        removeUnphased=FALSE,
                        verbose=FALSE, ## essentially debug mode
                        FBAT="/tmp/th/fbat_2.02c",
                        compVarExpl=FALSE ) {  ## Compute variance explained (takes time) -- keep off for simulations
  ######################
  ## Precompute setup ##
  ######################

  ##solve.svd <- getFromNamespace( "solve.svd", "fbati" ) ## CODETOOLS 03.04.2009

  ## if there wasn't a phenotype, then data is a little messed up... should fix this... -- or is this an R issue?
  if( is.null(data) )
    data <- ped

  ## We want the markerCondition to be the first two markers in marker
  if( is.null(markerAnalyze) || is.null(markerCondition) )
    stop( "markerCondition/markerAnalyze must be non-null." )

  ## Find the marker locations
  allMarkerNames       <- fixMarkerNames( names( ped ) )
  markerAnalyzeIndexTemp   <- match( markerAnalyze,   allMarkerNames ) ## Represents first index only - second is this +1
  markerConditionIndexTemp <- match( markerCondition, allMarkerNames )
  n <- length(markerAnalyzeIndexTemp)*2
  markerAnalyzeIndex <- rep(0,n)
  markerAnalyzeIndex[seq(from=1,to=n,by=2)] <- markerAnalyzeIndexTemp
  markerAnalyzeIndex[seq(from=2,to=n,by=2)] <- markerAnalyzeIndexTemp+1
  nC <- length(markerConditionIndexTemp)*2
  markerConditionIndex <- rep(0,nC)
  markerConditionIndex[seq(from=1,to=nC,by=2)] <- markerConditionIndexTemp
  markerConditionIndex[seq(from=2,to=nC,by=2)] <- markerConditionIndexTemp+1

  if( length(markerAnalyzeIndex)/2 != length(markerAnalyze) )
    stop( "Some markers in markerAnalyze do not exist." )
  if( length(markerConditionIndex)/2 != length(markerCondition) )
    stop( "Some markers in markerCondition do not exist." )

  ## Recode the pedigree alleles (in data) to all be 1/2
  #marker <- c(markerAnalyze, markerCondition)
  marker <- c( markerAnalyzeIndex, markerConditionIndex )
  for( i in seq( from=1, to=length(marker), by=2 ) ) {
    m1 <- marker[i]
    m2 <- marker[i+1]
    alleles <- sort(  setdiff( unique( c(data[[m1]], data[[m2]]) ), 0 )  ) ## sort is really only for debugging purposes -- won't change our simulated data then.
    ##alleles <- setdiff( unique( c(data[[m1]], data[[m2]]) ), 0 )
    if( length(alleles)>2 ) {
      stop( "biallelic loci only" )
    }else if( length(alleles) == 1 ) {
      data[[m1]][ data[[m1]]==alleles ] <- 1
      data[[m2]][ data[[m2]]==alleles ] <- 1
    }else if( length(alleles) == 2 ) {
      ## 1/04/09 addition for consistency when recoding...
      #if( sum(data[[m1]]==alleles[1]) + sum(data[[m2]]==alleles[1]) < sum(data[[m1]]==alleles[2]) + sum(data[[m2]]==alleles[2]) ) { # "DEBUG"
      if( sum(data[[m1]]==alleles[1]) + sum(data[[m2]]==alleles[1]) > sum(data[[m1]]==alleles[2]) + sum(data[[m2]]==alleles[2]) ) { ## Seems to be more powerful
        if( verbose ) cat( "Reverse allele invoked.\n" )
        alleles <- rev(alleles)
      }
      ## The case we expect most of the time
      a1 <- data[[m1]]==alleles[1]
      a2 <- data[[m1]]==alleles[2]
      data[[m1]][a1] <- 1
      data[[m1]][a2] <- 2
      a1 <- data[[m2]]==alleles[1]
      a2 <- data[[m2]]==alleles[2]
      data[[m2]][a1] <- 1
      data[[m2]][a2] <- 2
    }
  }

  if( verbose ) cat( "condGeneP3:: Alleles have been recoded.\n" );

  ## nuclify merge the data
  data <- nuclifyMerged( data )

  if( verbose ) cat( "condGeneP3:: nuclifyMerged\n" );

  ## Handle trait info
  traitCol <- trait
  if( is.character(trait) )
    traitCol <- which( names(data) == trait )
  if( length(traitCol) != 1 )
    stop( paste( "Trait '", trait, "' was not found.", sep="" ) )

  data_trait <- data[[traitCol]]

  ##print( trait )

  if( traitType=="auto" ) {
    if( trait=="AffectionStatus" || length(unique(data_trait))<=2 ) {
      traitType <- "binary"
    }else{
      traitType <- "continuous"
    }

    if( trait=="AffectionStatus" ) {
      data_trait[ data_trait==0 ] <- NA
      data_trait <- data_trait - 1
      #print( data_trait )
    }
  }

  if( traitType=="binary" )
    stop( "Intended for continuous traits, not binary traits." )

  if( verbose ) cat( "Trait type: ", traitType, "\n", sep="" )

  #######################
  ## Computation piece ##
  #######################

  ## The following was a bad idea, before the iteration piece
  ## New, need to compute the offset based on the mean of the 1/1 genotypes across all conditioning alleles
  ##wh <- rep( TRUE, nrow(ped) )
  #traitOffset <- mean( data[[traitCol]], na.rm=TRUE )
  #if( verbose ) cat( "Mean(trait): ", traitOffset, "\n", sep="" )
  #wh <- rep( TRUE, nrow(data) )
  #for( a in 1:(length(markerConditionIndex)/2) )
  #  wh <- wh & ( data[[ markerConditionIndex[a*2-1] ]] == 1 ) & ( data[[ markerConditionIndex[a*2] ]] == 1 )
  #wh <- wh & data$idfath!=0 & data$idmoth!=0  ## should eliminate parents?
  #if( any( wh ) ) {
  #  traitOffsetTemp <- mean( data[[traitCol]][wh], na.rm=TRUE )
  #  if( !is.na(traitOffsetTemp) ) {
  #    traitOffset <- traitOffsetTemp
  #    if( verbose ) cat( "Number used to estimate offset: ", sum(wh), "\n", sep="" )
  #  }else{
  #    cat( "No 1/1 genotypes available, using mean of all traits for offset.\n" )
  #    warning( "No 1/1 genotypes available, using mean of all traits for offset." )
  #  }
  #}else{
  #  cat( "No 1/1 genotypes available, using mean of all traits for offset.\n" )
  #  warning( "No 1/1 genotypes available, using mean of all traits for offset." )
  #}
  #if( verbose ) cat( "Offset: ", traitOffset, "\n", sep="" )

  ## Note that offset must be removed from _each_ of them, as each has it's own copy of the trait (not the best choice)

  #print( "markerCondition" ); print( markerCondition );
  #print( "markerConditionIndex" ); print( markerConditionIndex );
  #print( "markerAnalyze" ); print( markerAnalyze );
  #print( "markerAnalyzeIndex" ); print( markerAnalyzeIndex );
  #print( data[,c(1,2,markerAnalyzeIndex)] )

  ## Compute the density for each markerAnalyze and markerCondition
  analyzeReference <- rep( 0, length(markerAnalyze) )
  for( a in 1:(length(markerAnalyzeIndex)/2) ) {
    ##print( "Before Haplotype density" )
    analyzeReference[a] <- haplotypeDensity( data, markerAnalyzeIndex[a*2-c(0,1)], traitCol, tempPrefix, FBAT=FBAT ) ## ? also links the trait?
    ##print( "After Haplotype density" )
    #condGeneFBATControl_print( reference=analyzeReference[a] )
    #print( tail(ped, n=50) )
    #print( tail(phe, n=50) )
    #stop()
    ## and subtract the offset from the trait if qtl
    #if( traitType=="continuous" ) ## Has to be if in here...
    condGeneFBATControl_centerTrait( reference=analyzeReference[a], center=alpha, mean=TRUE )
    #condGeneFBATControl_centerTrait( reference=analyzeReference[a], center=traitOffset, mean=FALSE )

    #print( "######################")
    #print( "######################")
    #condGeneFBATControl_print( reference=analyzeReference[a] )
  }
  conditionReference <- rep( 0, length(markerCondition) )
  for( a in 1:(length(markerConditionIndex)/2) ) {
    ## Fill in the haplotype density,
    conditionReference[a] <- haplotypeDensity( data, markerConditionIndex[a*2-c(0,1)], traitCol, tempPrefix, FBAT=FBAT )
    ## and subtract the offset from the trait if qtl
    #if( traitType=="continuous" )
    condGeneFBATControl_centerTrait( reference=conditionReference[a], center=alpha, mean=TRUE )
    #condGeneFBATControl_centerTrait( reference=conditionReference[a], center=traitOffset, mean=FALSE )

    #print( "######################")
    #print( "######################")
    #condGeneFBATControl_print( reference=conditionReference[a] )
  }
  ## -- So all markers will be referenced to by allele 0, since we are always just using the first allele

  ## Compute the "explained variance" with both of the parameters, needs to be _before_ nuisance parameter calculations
  NUISANCE_ITER <- 1;
  varExpl <- NULL
  #compVarExpl <- FALSE
  #varExpl <- 1
  ####compVarExpl <- TRUE  ## No, we never want to do it this way now!
  ####verbose <- TRUE
  varMean <- 0
  if( compVarExpl ) {
    ##cat( "COMPVAREXPL BEGIN\n" )

    referenceVarExpl <- c( analyzeReference, conditionReference )

    ## Compute the variance of the mean
    #varMean <- condGeneFBATControl_varContsMean( referenceVarExpl, 0 ) ##betaEst ) -- not really used...
    #cat( "VARMEAN ORIG", varMean, "\n" )

    ## Backup the trait, since we're going to modify it
    condGeneFBATControl_backupTrait( referenceVarExpl );

    ## Get an estimate of beta
    betaEst <- condGeneFBATControl_estEqNuis( referenceVarExpl )
    if( verbose ) { cat( "betaEst for varExpl: " ); print( betaEst ); }
    varModel <- condGeneFBATControl_varContsModel( referenceVarExpl, betaEst )
    varExpl <- 1 - varModel / varMean

    if( verbose ) { cat( "Nuisance iteration ", 0, ", varExpl=", varExpl, ", betaEst=", sep="" ); print(betaEst); }

    for( i in 1:NUISANCE_ITER ) {
      #condGeneFBATControl_estEqNuisUpdate( referenceVarExpl, betaEst )
      condGeneFBATControl_estEqNuisUpdate2( referenceVarExpl, betaEst )
      bc <- condGeneFBATControl_estEqNuis( referenceVarExpl, offset=0.0 )
      #varExpl <- condGeneFBATControl_varExplConts( referenceVarExpl, betaEst )
      varModel <- condGeneFBATControl_varContsModel( referenceVarExpl, betaEst )
      #varModel <- condGeneFBATControl_varContsMean( referenceVarExpl, 0 ) ## DEBUGGING ONLY!!!
      #cat( "varModel after nuis update", varModel, "\n" )
      varExpl <- 1 - varModel / varMean
      if( verbose ) { cat( "Nuisance iteration ", i, ", varExpl=", varExpl, ", betaEst=", sep="" ); print(bc); }
    }

    ## Now restore the trait back for testing
    condGeneFBATControl_restoreTrait( referenceVarExpl );

    ## Compute the variance of the mean
    varMean <- condGeneFBATControl_varContsMean( referenceVarExpl, 0 ) ##betaEst ) -- not really used...
    #cat( "VARMEAN ORIG (SUB)", varMean, "\n" )

		#varExpl <- 1 - varModel / varMean

		varExpl <- varModel

    ##cat( "COMPVAREXPL END\n" )
  }
  #compVarExpl <- TRUE

  ## Compute the nuisance parameters
  #condGeneFBATControl_print( conditionReference[1] )
  bc <- condGeneFBATControl_estEqNuis( conditionReference, offset=0.0 ) ## Offset will need to be 0 for continuous (already subtracted above...)
  if( verbose ) { cat( "condGeneP3 bc " ); print( bc ); }

  for( i in 1:NUISANCE_ITER ) {
#    condGeneFBATControl_estEqNuisUpdate( conditionReference, bc )
    condGeneFBATControl_estEqNuisUpdate2( conditionReference, bc )
    bc <- condGeneFBATControl_estEqNuis( conditionReference, offset=0.0 )
    if( verbose ) { cat( "Nuisance iteration ", i, ", bc=", sep="" ); print(bc); }
  }
  #stop("DEBUG: Updated nuisance parameter.")


  ## DEBUG
  ##varExplCOND <- condGeneFBATControl_varExplConts( conditionReference, bc )
  ##cat( "COND VAR EXPL", varExplCOND, "\n" )
  ## GUBED

  ## Compute the test statistic
  ##  -- offset will need to not be zero for dichotomous case
  res <- condGeneFBATControl_estEq( analyzeReference, conditionReference, bc, offset=0.0, verbose=verbose )
  #cat( "condGene result pvalue=", res$pvalue, ", rank=", res$rank, "\n" )
  ##condGeneFBATControl_print( conditionReference[1] )
  ##stop()

  if( compVarExpl )
    res$varExpl <- varExpl

  ## TWO STEP ADDITION
  ## -- Code now adjusts this...
  #bc <- condGeneFBATControl_estEqNuis( conditionReference, offset=0.0 ) ## Offset will need to be 0 for continuous
  #cat( "condGeneP3 bc (#2) " ); print( bc );
  #res <- condGeneFBATControl_estEq( analyzeReference, conditionReference, bc, offset=0.0 )
  ## NOITIDDA PETS OWT

  ## Can't forget to free the densities then as well...
  for( a in 1:length(analyzeReference) )
    condGeneFBATControl_free( analyzeReference[a] )
  for( a in 1:length(conditionReference) )
    condGeneFBATControl_free( conditionReference[a] )

  ## And return some value
  return( res )
}





## NEWEST -- fixing up the dichotomous case
condGeneP4 <- function( ped=NULL, phe=NULL, data=mergePhePed( ped, phe ),
    trait="AffectionStatus", traitType="auto",
    markerAnalyze=NULL, markerCondition=NULL,
    sigma=1, alpha=NULL,
    ignoreBtX=TRUE,
    tempPrefix="temp",
    removeUnphased=FALSE,
    verbose=FALSE,
    FBAT="/tmp/th/fbat_2.02c" ) { ## essentially debug mode
  ######################
  ## Precompute setup ##
  ######################
  #markerNames <- unique( fixMarkerNames( #names(data)[markerCol] ) )
  if( verbose ) cat( "Entered condGeneP4...\n")

  solve.svd <- getFromNamespace( "solve.svd", "fbati" )

  ## if there wasn't a phenotype, then data is a little messed up... should fix this... -- or is this an R issue?
  if( is.null(data) )
    data <- ped

  ## We want the markerCondition to be the first two markers in marker
  if( is.null(markerAnalyze) || is.null(markerCondition) )
    stop( "markerCondition/markerAnalyze must be non-null." )

  ## Find the marker locations
  allMarkerNames       <- fixMarkerNames( names( ped ) )
  markerAnalyzeIndexTemp   <- match( markerAnalyze,   allMarkerNames ) ## Represents first index only - second is this +1
  markerConditionIndexTemp <- match( markerCondition, allMarkerNames )
  n <- length(markerAnalyzeIndexTemp)*2
  markerAnalyzeIndex <- rep(0,n)
  markerAnalyzeIndex[seq(from=1,to=n,by=2)] <- markerAnalyzeIndexTemp
  markerAnalyzeIndex[seq(from=2,to=n,by=2)] <- markerAnalyzeIndexTemp+1
  nC <- length(markerConditionIndexTemp)*2
  markerConditionIndex <- rep(0,nC)
  markerConditionIndex[seq(from=1,to=nC,by=2)] <- markerConditionIndexTemp
  markerConditionIndex[seq(from=2,to=nC,by=2)] <- markerConditionIndexTemp+1

  if( length(markerAnalyzeIndex)/2 != length(markerAnalyze) )
    stop( "Some markers in markerAnalyze do not exist." )
  if( length(markerConditionIndex)/2 != length(markerCondition) )
    stop( "Some markers in markerCondition do not exist." )

  ## Recode the pedigree alleles (in data) to all be 1/2
  #marker <- c(markerAnalyze, markerCondition)
  marker <- c( markerAnalyzeIndex, markerConditionIndex )
  for( i in seq( from=1, to=length(marker), by=2 ) ) {
    m1 <- marker[i]
    m2 <- marker[i+1]
    alleles <- sort(  setdiff( unique( c(data[[m1]], data[[m2]]) ), 0 )  ) ## sort is really only for debugging purposes -- won't change our simulated data then.
    ##alleles <- setdiff( unique( c(data[[m1]], data[[m2]]) ), 0 )
    if( length(alleles)>2 ) {
      stop( "biallelic loci only" )
    }else if( length(alleles) == 1 ) {
      data[[m1]][ data[[m1]]==alleles ] <- 1
      data[[m2]][ data[[m2]]==alleles ] <- 1
    }else if( length(alleles) == 2 ) {
      ## The case we expect most of the time
      a1 <- data[[m1]]==alleles[1]
      a2 <- data[[m1]]==alleles[2]
      data[[m1]][a1] <- 1
      data[[m1]][a2] <- 2
      a1 <- data[[m2]]==alleles[1]
      a2 <- data[[m2]]==alleles[2]
      data[[m2]][a1] <- 1
      data[[m2]][a2] <- 2
    }
  }

  if( verbose ) cat( "condGeneP:: Alleles have been recoded.\n" );

  ## nuclify merge the data
  data <- nuclifyMerged( data )

  if( verbose ) cat( "condGeneP:: nuclifyMerged\n" );

  ## Handle trait info
  traitCol <- trait
  if( is.character(trait) )
    traitCol <- which( names(data) == trait )
  if( length(traitCol) != 1 )
    stop( paste( "Trait '", trait, "' was not found.", sep="" ) )

  data_trait <- data[[traitCol]]

  ##print( trait )

  if( traitType=="auto" ) {
    if( trait=="AffectionStatus" || length(unique(data_trait))<=2 ) {
      traitType <- "binary"
    }else{
      traitType <- "continuous"
    }

    if( trait=="AffectionStatus" ) {
      data_trait[ data_trait==0 ] <- NA
      data_trait <- data_trait - 1
      #print( data_trait )
    }
  }

  if( verbose ) cat( "Trait type: ", traitType, "\n", sep="" )

  #######################
  ## Computation piece ##
  #######################

  ################################################
  ## First solve for the nuisance parameters... ##

  ## Compute the haplotype density
  conditional_allele_index <- markerConditionIndex
  referenceC <- haplotypeDensity( data, markerConditionIndex, traitCol, tempPrefix, FBAT=FBAT )
  #condGeneFBATControl_print( referenceC ) ## DEBUG
  #print( condGeneFBATControl_pids( referenceC ) ) ## DEBUG
  n <- condGeneFBATControl_numFam( referenceC )

  ## Now solve for the nuisance parameters
  bc <- NULL
  analyze_allele_index <- c()
  conditional_allele_index2 <- 1:length(markerCondition) - 1
  TOL <- .Machine$double.eps
  umcNuis <- NULL
  ##CENTER_QTL <- FALSE ## CODETOOLS 03.04.2009
  if( traitType != "binary" ) stop( "condGeneP4 only for dichotomous traits!" )
  initial <- rep( 0, 2*length(markerCondition) )
  bcSolve <- NULL
  ##bcSolve <- mrroot2( f=condGeneFBATControl_uimc_nuis, initial=initial, reference=reference, analyze_allele_index=analyze_allele_index, conditional_allele_index=conditional_allele_index, onlyComputeConditional=TRUE, MAXITER=MAXITER, TOL=TOL, verbose=verbose )
  try({
        ##require( rootSolve )
        bcSolve <- multiroot( f=condGeneFBATControl_uimc_nuis, start=initial, maxiter=100, rtol=TOL, atol=TOL, ctol=TOL, useFortran=FALSE, reference=referenceC, analyze_allele_index=analyze_allele_index, conditional_allele_index=conditional_allele_index2, onlyComputeConditional=TRUE, verbose=verbose ) ## onlyComputeConditional irrelevant here...
        bc <- bcSolve$root
      }, silent=TRUE )
  if( is.null(bcSolve) ) {
    if( verbose ) cat( "multiroot failed, falling back on slower convergence method.\n" )
    try( {
          bcSolve <- mrroot2( f=condGeneFBATControl_uimc_nuis, initial=initial, reference=referenceC, analyze_allele_index=analyze_allele_index, conditional_allele_index=conditional_allele_index2, onlyComputeConditional=TRUE, verbose=verbose )
          bc <- bcSolve$par
        }, silent=TRUE );
    if( is.null( bcSolve ) ) {
      if( verbose ) cat( "could not solve for nuisance parameter!" )
      return( data.frame( pvalue=1, rank=0,
              pvalueR=1, rankR=0 ) )
    }
  }
  umcNuis <- condGeneFBATControl_uimc( referenceC,
      c(), bc[1:length(conditional_allele_index2)], bc[(length(conditional_allele_index2)+1):length(bc)],
      analyze_allele_index, conditional_allele_index2,
      FALSE )   ## $conditional0, $conditional1
  umcNuis <- cbind( umcNuis$conditional0, umcNuis$conditional1 )

  if( verbose ) {
    cat( "bc" )
    print( bc )
    cat( "Number of inf fams used: ", condGeneFBATControl_numInfFam(referenceC), "\n" )
  }

  #############################################
  ## Now handle each of the pairwise fits... ##
  analyze_allele_index <- 0 ## pairwise, so always moving across this marker...
  conditional_allele_index <- 1:length(markerCondition) + length(analyze_allele_index) - 1

  ##R <- length(markerAnalyze) + 2*length(markerCondition) ## CODETOOLS 03.04.2009
  A <- length(markerAnalyze)
  C <- length(markerCondition)
  C2 <- C*2
  #umcMat <- matrix( 0, nrow=n, ncol=R )
  #umcMat[, (R-ncol(umcNuis)+1):R ] <- umcNuis
  um <- matrix( 0, nrow=n, ncol=A )
  uc <- matrix( 0, nrow=n, ncol=C2 )

  umcRMat <- matrix( 0, nrow=n, ncol=length(markerAnalyze) )

  reference <- rep(0,length(markerAnalyze))
  for( i in 1:length(markerAnalyze) ) {
    ## Get the haplotype distribution
    reference[i] <- haplotypeDensity( data, c(markerAnalyzeIndex[i*2-c(1,0)], markerConditionIndex), traitCol, tempPrefix, FBAT=FBAT )

    ##condGeneFBATControl_print( reference ) ## DEBUG
    ##print( condGeneFBATControl_pids( reference ) ) ## DEBUG

    ## Compute the pieces of the statistic
    bc0 <- bc[1:length(conditional_allele_index)]
    bc1 <- bc[(length(conditional_allele_index)+1):length(bc)]
    umcAddi <- condGeneFBATControl_uimc( reference,
        rep(0,length(analyze_allele_index)), bc0, bc1,
        analyze_allele_index, conditional_allele_index,
        FALSE )
    #print( str( umcAddi ) )
    #print( str( um ) )
    #stop()
    um[, i] <- umcAddi$analyze
    #uc[, i] <- cbind( umcAddi$conditional0, umcAddi$conditional1 )
    ####uc[, i] <- umcAddi$conditional0     ## 11.18.2008
    ####uc[, i+C] <- umcAddi$conditional1   ## 11.18.2008
    ##umcMat[ , i ] <- umcAddi$analyze  ## Will need to adjust later...

    ##print( umcMat )
    ##stop()

    ## And the robust statistic piece
    umcRMat[ , i ] <- condGeneFBATControl_robustStat( reference, analyze_allele_index, conditional_allele_index )
  } ## i

  uc <- umcNuis

  umc_ucc <- condGeneFBATControl_dUdBc( reference, referenceC,
      analyze_allele_index, conditional_allele_index,
      conditional_allele_index2, bc )
  ##print( umc_ucc )
  #stop()
  umc <- umc_ucc$umc
  ucc <- umc_ucc$ucc

  #print( "dim - um, umc, ucc, uc")
  #print(dim(um))
  #print(dim(umc))
  #print(dim(ucc))
  #print(dim(uc))

  ## MAJOR CHANGE MAJOR CHANGE
  ##umcMat <- um  - t(  umc %*% solve.svd( svd(ucc), t(uc) )  )## ??
  umcMat <- um  + t(  umc %*% solve.svd( svd(ucc), t(uc) )  )## ??

  #print( "condGeneP umcMat" )
  #print( umcMat )

  ## Get the new pieces
  ## Then adjust umcMat for the nuisance parameters

  ################################
  ## Compute the test statistic ##
  if( verbose ) cat( "About to compute the test statistic.\n" )

  ## Model-based test using empirical variance
  u <- rep( 0, A )
  uu <- matrix( 0, nrow=A, ncol=A )
  for( k1 in 1:A ) {
    u[k1] <- sum( umcMat[,k1], na.rm=TRUE )
    for( k2 in k1:A ) {
      uu[k1,k2] <- sum( umcMat[,k1] * umcMat[,k2], na.rm=TRUE )
      uu[k2,k1] <- uu[k1,k2]
    }
  }

  ## inf fams added 02.29.2009
  #numInf <- rep( 0, A )
  #for( k1 in 1:A )
  #  numInf[k1] <- sum( abs(umcMat[,k1]) >= sqrt(.Machine$double.eps), na.rm=TRUE )
  numInf <- ""
  for( k1 in 1:A ) {
    temp <- sum(abs(um[,k1])>=sqrt(.Machine$double.eps),na.rm=TRUE)
    if( k1!=A ) {
      numInf <- paste(numInf, temp, sep=",")
    }else{
      numInf <- as.character(temp)
    }
  }
  for( k1 in 1:C ) {
    temp1 <- sum(abs(uc[,k1])>=sqrt(.Machine$double.eps), na.rm=TRUE)
    temp2 <- sum(abs(uc[,k1+C])>=sqrt(.Machine$double.eps), na.rm=TRUE)
    temp <- paste( "(", temp1, ",", temp2, ")", sep="" )
    if( k1==1 ) {
      numInf <- paste( numInf, temp, sep="|" )
    }else{
      numInf <- paste( numInf, temp, sep="," )
    }
  }

  if( verbose ) {
    print( "condGeneP Pairwise u, uu" )
    print( u )
    print( uu )
  }

  pvalue <- 1 ## only computes empirical
  rank <- 1
  try( {
        ## Compute the test statistic
        post <- solve.svd( svd( uu ), u )
        rank <- attr(post,"rank")
        stat <- u %*% post
        ## and compute the pvalue
        pvalue <- pchisq( stat, df=rank, lower.tail=FALSE )
      }, silent=TRUE )
  if( verbose ) cat( 'p-value', pvalue, "\n" )

  ##############################################
  ## Robust test statistic (empirical variance)
  na <- length(markerAnalyze)
  u <- rep( 0, na )
  uu <- matrix( 0, nrow=na, ncol=na )
  for( k1 in 1:na ) {
    u[k1] <- sum( umcRMat[,k1], na.rm=TRUE )
    for( k2 in k1:na ) {
      uu[k1,k2] <- sum( umcRMat[,k1] * umcRMat[,k2], na.rm=TRUE )
      uu[k2,k1] <- uu[k1,k2]
    }
  }

  numInfR <- rep( 0, A )
  for( k1 in 1:A )
    numInfR[k1] <- sum( abs(umcRMat[,k1]) >= sqrt(.Machine$double.eps), na.rm=TRUE )

  pvalueR <- 1
  rankR <- 1
  try( {
        ## Compute the test statistic
        post <- solve.svd( svd(uu), u )
        stat <- u %*% post
        ##cat( "stat", stat, "\n" )

        ## Get the rank (should be stable)
        rankR <- attr( post, "rank" )

        ## and compute the pvalue
        pvalueR <- pchisq( stat, df=rankR, lower.tail=FALSE )
      }, silent=TRUE )

  ##stop()

  ## And finally free the memory
  for( i in 1:length(markerAnalyze) )
    condGeneFBATControl_free( reference[i] )
  condGeneFBATControl_free( referenceC ) ## Free the data here

  if( verbose ) cat( "Going to exit condGeneP4...\n")

  #############
  ## Return! ##
      return( data.frame( pvalue=pvalue, rank=rank,
                      #numInf=paste(numInf,collapse=","),
                      numInf=numInf,
                      pvalueR=pvalueR, rankR=rankR, numInfR=paste(numInfR,collapse=","),
                      stringsAsFactors=FALSE) )
}


## Does pairwise test ONLY (for qtl export in condGene)
condGeneR <- function( ped=NULL, phe=NULL, data=mergePhePed( ped, phe ),
                       trait="AffectionStatus", traitType="auto",
                       markerAnalyze=NULL, markerCondition=NULL,
                       sigma=1, alpha=NULL,
                       ignoreBtX=TRUE,
                       tempPrefix="temp",
                       removeUnphased=FALSE,
                       verbose=FALSE,
                       FBAT="/tmp/th/fbat_2.02c" ) { ## essentially debug mode
  ######################
  ## Precompute setup ##
  ######################

  solve.svd <- getFromNamespace( "solve.svd", "fbati" )

  ## if there wasn't a phenotype, then data is a little messed up... should fix this... -- or is this an R issue?
  if( is.null(data) )
    data <- ped

  ## We want the markerCondition to be the first two markers in marker
  if( is.null(markerAnalyze) || is.null(markerCondition) )
    stop( "markerCondition/markerAnalyze must be non-null." )

  ## Find the marker locations
  allMarkerNames       <- fixMarkerNames( names( ped ) )
  markerAnalyzeIndexTemp   <- match( markerAnalyze,   allMarkerNames ) ## Represents first index only - second is this +1
  markerConditionIndexTemp <- match( markerCondition, allMarkerNames )
  n <- length(markerAnalyzeIndexTemp)*2
  markerAnalyzeIndex <- rep(0,n)
  markerAnalyzeIndex[seq(from=1,to=n,by=2)] <- markerAnalyzeIndexTemp
  markerAnalyzeIndex[seq(from=2,to=n,by=2)] <- markerAnalyzeIndexTemp+1
  nC <- length(markerConditionIndexTemp)*2
  markerConditionIndex <- rep(0,nC)
  markerConditionIndex[seq(from=1,to=nC,by=2)] <- markerConditionIndexTemp
  markerConditionIndex[seq(from=2,to=nC,by=2)] <- markerConditionIndexTemp+1

  if( length(markerAnalyzeIndex)/2 != length(markerAnalyze) )
    stop( "Some markers in markerAnalyze do not exist." )
  if( length(markerConditionIndex)/2 != length(markerCondition) )
    stop( "Some markers in markerCondition do not exist." )

  ## Recode the pedigree alleles (in data) to all be 1/2
  #marker <- c(markerAnalyze, markerCondition)
  marker <- c( markerAnalyzeIndex, markerConditionIndex )
  for( i in seq( from=1, to=length(marker), by=2 ) ) {
    m1 <- marker[i]
    m2 <- marker[i+1]
    alleles <- sort(  setdiff( unique( c(data[[m1]], data[[m2]]) ), 0 )  ) ## sort is really only for debugging purposes -- won't change our simulated data then.
    ##alleles <- setdiff( unique( c(data[[m1]], data[[m2]]) ), 0 )
    if( length(alleles)>2 ) {
      stop( "biallelic loci only" )
    }else if( length(alleles) == 1 ) {
      data[[m1]][ data[[m1]]==alleles ] <- 1
      data[[m2]][ data[[m2]]==alleles ] <- 1
    }else if( length(alleles) == 2 ) {
      ## The case we expect most of the time
      a1 <- data[[m1]]==alleles[1]
      a2 <- data[[m1]]==alleles[2]
      data[[m1]][a1] <- 1
      data[[m1]][a2] <- 2
      a1 <- data[[m2]]==alleles[1]
      a2 <- data[[m2]]==alleles[2]
      data[[m2]][a1] <- 1
      data[[m2]][a2] <- 2
    }
  }

  if( verbose ) cat( "condGeneP:: Alleles have been recoded.\n" );

  ## nuclify merge the data
  data <- nuclifyMerged( data )

  if( verbose ) cat( "condGeneP:: nuclifyMerged\n" );

  ## Handle trait info
  traitCol <- trait
  if( is.character(trait) )
    traitCol <- which( names(data) == trait )
  if( length(traitCol) != 1 )
    stop( paste( "Trait '", trait, "' was not found.", sep="" ) )

  data_trait <- data[[traitCol]]

  ##print( trait )

  if( traitType=="auto" ) {
    if( trait=="AffectionStatus" || length(unique(data_trait))<=2 ) {
      traitType <- "binary"
    }else{
      traitType <- "continuous"
    }

    if( trait=="AffectionStatus" ) {
      data_trait[ data_trait==0 ] <- NA
      data_trait <- data_trait - 1
      #print( data_trait )
    }
  }

  if( verbose ) cat( "Trait type: ", traitType, "\n", sep="" )


  #######################
  ## Computation piece ##
  #######################

  #############################################
  ## Now handle each of the pairwise fits... ##
  #analyze_allele_index <- 1:length(markerAnalyze) - 1
  #conditional_allele_index <- 1:length(markerCondition) + length(markerAnalyze) - 1
  analyze_allele_index <- 0 ## pairwise, so always moving across this marker...
  conditional_allele_index <- 1:length(markerCondition) + length(analyze_allele_index) - 1

  ##R <- length(markerAnalyze) + 2*length(markerCondition)*length(markerAnalyze)  ## CODETOOLS 03.04.2009
  A <- length(markerAnalyze)
  umcRMat <- NULL

  ##CONVERGED <- TRUE ## CODETOOLS 03.04.2009
  numInformative <- c()
  for( i in 1:length(markerAnalyze) ) {
    ## Get the haplotype distribution
    reference <- haplotypeDensity( data, c(markerAnalyzeIndex[i*2-c(1,0)], markerConditionIndex), traitCol, tempPrefix, FBAT=FBAT )

    ## potentially needs to be sized...
    numFam <- condGeneFBATControl_numFam( reference )
    if( is.null(umcRMat) )
      umcRMat <- matrix( 0, nrow=numFam, ncol=A )

    ## NEW ADDITION FOR THE ROBUST PIECE,
    ##  If dichotomous, see if there is an offset
    if( traitType=="binary" && !is.null(alpha) && alpha!=0 ) {
      #print( "Get your offset on!" ) ## DEBUG ONLY
      #print( alpha )
      ## Ensure that the model-based test doesn't change, just for debug externally

      ## Subtract the offset from the trait
      condGeneFBATControl_centerTrait( reference=reference, center=alpha, mean=FALSE )
    }
    #condGeneFBATControl_print( reference ) ## DEBUG
    #stop()

    ## And the robust statistic piece
    umcRMat[ , i ] <- condGeneFBATControl_robustStat( reference, analyze_allele_index, conditional_allele_index )

    ## NEW! Get the number informative
    numInformative <- c( numInformative, condGeneFBATControl_numInfFam( reference ) )

    ## And finally free the memory
    condGeneFBATControl_free( reference )
  } ## i

  ################################
  ## Compute the test statistic ##
  if( verbose ) cat( "About to compute the test statistic.\n" )

  ## Robust test statistic (empirical variance)
  na <- length(markerAnalyze)
  u <- rep( 0, na )
  uu <- matrix( 0, nrow=na, ncol=na )
  for( k1 in 1:na ) {
    u[k1] <- sum( umcRMat[,k1], na.rm=TRUE )
    for( k2 in k1:na ) {
      uu[k1,k2] <- sum( umcRMat[,k1] * umcRMat[,k2], na.rm=TRUE )
      uu[k2,k1] <- uu[k1,k2]
    }
  }

  numInfR <- rep( 0, A )
  for( k1 in 1:A )
    numInfR[k1] <- sum( abs(umcRMat[,k1]) >= sqrt(.Machine$double.eps), na.rm=TRUE )

  ##paste( numInfR )  ## 02.29.2009 -- WTF

  pvalueR <- 1
  rankR <- 0 ## 03.04.2009
  try( {
    ## Compute the test statistic
    post <- solve.svd( svd(uu), u )
    stat <- u %*% post
    ##cat( "stat", stat, "\n" )

    ## Get the rank (should be stable)
    rankR <- attr( post, "rank" )

    ## and compute the pvalue
    pvalueR <- pchisq( stat, df=rankR, lower.tail=FALSE )
  }, silent=TRUE )

  ##if( !CONVERGED ) {
  ##  pvalue <- 1
  ##  rank <- 0
  ##}

  #############
  ## Return! ##
  if( traitType=="continuous" ) {
    offset <- NA
  }else{
    offset <- alpha
  }
  if( is.null(offset) ) offset <- NA
  return( data.frame( trait=trait, traitType=traitType,
                      offset=offset,
                      pvalueR=pvalueR, rankR=rankR, numInfR=as.character(paste(numInfR,collapse=",")),
                      stringsAsFactors=FALSE) )
}


## FOR DEBUG
#ped <- fread.ped( "subped.ped" )
#condGene(ped=ped,marker=c("ser8","ser37"))
