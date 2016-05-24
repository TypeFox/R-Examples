#library( fbati )
#dyn.load( "fbatge.so" )

#solve.svd <- getFromNamespace( "solve.svd", "fbati" )

TRYRR <- 1000 ## 1000

#######################################################
## New gped creation routines from a pedigree object ##
#######################################################

## Always sorts the alleles, so it would be possible to merge two gped objects
as.gped <- function( ped ) {
  gped <- ped[,1:6]
  mnames <- ped.markerNames(ped)
  for( m in 1:length(mnames) ) {
    #if( debug ) cat( "genotypeCode marker", mnames[m], "\n" )

    ## compute locations of the markers in the ped
    c1 <- 6 + m*2 - 1  ## location of first marker
    c2 <- 6 + m*2      ## location of second marker

    ## get the alleles
    alleles <- sort(  setdiff( unique(c(ped[[c1]],ped[[c2]])), 0 )  )
    #if( debug ) { cat( "Alleles " ); print( alleles ); }

    if( length(alleles) > 2 )
      stop( "Only works with biallelic markers." )

    ## Code it in an additive sense
    a <- alleles[1]
    newG <- as.integer(ped[[c1]]==a) + as.integer(ped[[c2]]==a)
    wh <- ped[[c1]]==0 | ped[[c2]]==0
    if( any(wh) )
      newG[wh] <- NA

    #print( newG )

    gped[[ ncol(gped) + 1 ]] <- newG ## ?do we want the as.factor here??
  }

  ## fix the names
  names(gped)[7:ncol(gped)] <- mnames

  ## return the newly coded pedigree object
  return( gped )
}

mergePheGped <- function( gped, phe )
  return( mergePhePed( ped=gped, phe=phe ) )


###########################
## Root solving routines ##
###########################

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

## Improved version always tries a certain number of times
##  to make sure it doesn't converge to a better value... (by default this is turned off)
## Newest improvement first tries to run multiroot (it's much faster)
##  AND the return is just the beta values
mrroot3 <- function( f, ## the function
                     initial, ## the initial value,
                     minRR=rep(-6,length(initial)), maxRR=rep(6,length(initial)),  tryRR=TRYRR, minRRsuccess=1, ## bounds for random restart
                     method="Nelder-Mead",  ##"BFGS",
                     MAXITER=1000, TOL=sqrt(.Machine$double.eps),
                     verbose=FALSE,
                    ... ) {
  ## New 02/18/2009, first try multiroot
  try( {
    ##require( rootSolve )
    res <- multiroot( f=f, start=initial, maxiter=MAXITER, rtol=TOL, atol=TOL, ctol=TOL, useFortran=FALSE, ... )
    return( res$root )
  }, silent=!verbose )

  if( verbose ) print( "multiroot failed." )

  ## for storing the best result
  best <- NULL

  ## No successes
  successes <- 0


  for( t in 1:tryRR ) {
    ##cat( "(",t,") ", sep="" )
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
          return( best$par )
      }else{
        ##cat( "FAILED TO CONVERGE..." )
        #print( res )
        #stop()
      }

      cat( "+" )
    }

    ## And a random start
    ##initial <- runif( n=length(initial), min=minRR, max=maxRR )
    initial <- 3*rnorm( n=length(initial) )
  } ## t
  cat( "\n" )

  if( abs(best$value) > 1 ) {
    return(NULL) ## failure... die...
  }

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

  print( best$par )

  ## and still return the best
  ##return( best )
  return( best$par )
}

##############################
## MMatrix linking routines ##
##############################
cpp_mmatrix_debug <- function() { .C("cpp_mmatrix_debug") }
##cpp_mmatrix_debug(); stop();

############################
## Random number routines ##
############################
cpp_rn_attach <- function() { .C("cpp_rn_attach") }
cpp_rn_detach <- function() { .C("cpp_rn_detach") }
cpp_rn_setNormalSigma <- function( sigma ) {
  cholesky <- chol(sigma)
  print( cholesky )
  p <- nrow(sigma)
  .C("cpp_rn_setNormalSigma", as.double(cholesky), as.integer(p), DUP=TRUE) #DUP=FALSE )
  return( invisible() )
}
cpp_rn_debug <- function() { .C("cpp_rn_debug") }

rn_debug <- function() {
  cpp_rn_attach()

  rho <- 0.5
  sigma <- rbind( c(1,rho,rho), c(rho,1,rho), c(rho,rho,1) )
  cpp_rn_setNormalSigma( sigma )
  cpp_rn_debug()
  cpp_rn_detach()
}
##rn_debug()

#########################################
## GPed & Estimating equation routines ##
#########################################

cpp_gped_clear <- function() { .C("cpp_gped_clear") }
cpp_gped_print <- function( extended=0 ) {
  .C("cpp_gped_print", as.integer(extended) )
  return(invisible())
}
cpp_gped_numCovariates <- function() {
  numCov <- as.integer(0)
  res = .C( "cpp_gped_numCovariates", numCov, DUP=TRUE) #DUP=FALSE )
  #return( numCov )
  return(res[[1]])
}
cpp_gped_setStrategy <- function( strategy="adaptive" ) {
  ## strategy \in geno, pheno, adaptive
  .C( "cpp_gped_setStrategy", as.character(strategy), DUP=TRUE )
  return( invisible() )
}
cpp_gped_estEq <- function( beta ) {
  ## C++ code will crash the whole program if beta is the wrong size
  u <- as.double( rep( 0, length(beta) ) )
  res = .C( "cpp_gped_estEq", as.double(beta), as.integer(length(beta)), u, DUP=TRUE) #DUP=FALSE )
  ####cat( "u" ); print( u );
  #return( u )
  return(res[[3]])
}
cpp_gped_estEq_zeroGE <- function( beta ) {
  uExtra <- cpp_gped_estEq( c( 0, 0, beta ) )
  return( uExtra[3:length(uExtra)] )
}

## Precondition: cpp_gped_setStrategy _MUST_ have been run, with the same strategy!
solveBeta <- function( strategy="adaptive", includeGE=TRUE ) {
  ## Want to not includeGE when solving for the nuisance parameters for the score test
  betaLength <- 4
  if( strategy!="geno" )
    betaLength <- 5 + cpp_gped_numCovariates()
  if( !includeGE )
    betaLength <- betaLength - 2 ## Take off the GE beta parameters

  #cat( "solveNuis debug betaLength:", betaLength, "\n" )

  beta <- rep( 0, betaLength )
  ##print( "solveBeta beta" ); print( beta )
  if( includeGE )
    return( mrroot3( f=cpp_gped_estEq, initial=beta ) )
  return( mrroot3( f=cpp_gped_estEq_zeroGE, initial=beta ) )
}

cpp_gped_statPrecompute <- function( beta ) {
  M <- length(beta)
  if( M == 0 ) return( NULL )
  beta <- c( 0, 0, beta ) ## push in the zeros that were there for the ge stuff
  dbnuis <- as.double( matrix( 0, nrow=M, ncol=M ) )
  res = .C( "cpp_gped_statPrecompute", as.double(beta), as.integer(length(beta)), dbnuis, DUP=TRUE) #DUP=FALSE )
  dbnuis = res[[3]]
  return( matrix( dbnuis, nrow=M, ncol=M ) )
}
cpp_gped_statCompute <- function( dbNuisInv ) {
  stat <- as.double(0.0)
  numInf <- as.integer(0)
  res = .C( "cpp_gped_statCompute", dbNuisInv, stat, numInf, DUP=TRUE )#DUP=FALSE )
  stat = res[[2]]
  numInf = res[[3]]
  pvalue <- pchisq( stat, 2, lower.tail=FALSE )
  attr( pvalue, "numInf" ) <- numInf
  return( pvalue ) ## chi-squared distribution with 2 degrees of freedom
}

testBetaGE <- function( strategy="adaptive" ) {
  ## First set the strategy
  cpp_gped_setStrategy( strategy=strategy )

  ## First solve for the nuisance parameter
  beta <- solveBeta( strategy=strategy, includeGE=FALSE )
  #print( "solved for beta" )
  ####cat( "testBetaGE beta = " ); print( beta );
  if( is.null(beta) ) {
    cat( "Could not solve for nuisance parameters.\n" )
    return( 1 )
  }

  ## Second, compute the test statistic
  ## - first the precompute piece
  dbnuis <- cpp_gped_statPrecompute( beta )
  ####cat( "testBetaGE dbnuis\n" ); print( dbnuis ); ## DEBUG DEBUG
  if( is.null(nrow(dbnuis)) ) {
    cat( "Could not invert nuisance parameters.\n" )
    return( 1 ) ## p-value of 1 safest...
  }
  #cat( "dbnuis\n" ); print( dbnuis )
  ## - then the compute piece, passing back the inverse (not as numerically stable yet...)
  #print( dbnuis )
  #print( solve(dbnuis) )
  #print( dbnuis %*% solve(dbnuis) )
  #pvalue <- cpp_gped_statCompute( solve(dbnuis) )
  ####cat( "testBetaGE dbnuisInv\n" ); print( solve.svd(svd(dbnuis)) ); ## DEBUG DEBUG

  #solve.svd <- getFromNamespace( "solve.svd", "fbati" )
  pvalue <- cpp_gped_statCompute( solve.svd( svd( dbnuis ) ) )

  return( pvalue );
}

## NOW THE ADDITIVE GE ROUTINES ##
cpp_gped_statPrecompute_A<- function( beta ) {
  M <- length(beta)
  if( M == 0 ) return( NULL )
  beta <- c( 0, beta ) ## push in the zeros that were there for the ge stuff
  dbnuis <- as.double( matrix( 0, nrow=M, ncol=M ) )
  res = .C( "cpp_gped_statPrecompute_A", as.double(beta), as.integer(length(beta)), dbnuis, DUP=TRUE ) #DUP=FALSE )
  dbnuis = res[[3]]
  return( matrix( dbnuis, nrow=M, ncol=M ) )
}
cpp_gped_statCompute_A<- function( dbNuisInv ) {
  stat <- as.double(0.0)
  numInf <- as.integer(0)
  res = .C( "cpp_gped_statCompute_A", dbNuisInv, stat, numInf, DUP=TRUE) #DUP=FALSE )
  stat = res[[2]]
  numInf = res[[3]]
  pvalue <- pchisq( stat, 1, lower.tail=FALSE )
  attr( pvalue, "numInf" ) <- numInf
  return( pvalue ) ## chi-squared distribution with 2 degrees of freedom
}

testBetaGE_A <- function( strategy="adaptive" ) {
  ## First set the strategy
  cpp_gped_setStrategy( strategy=strategy )

  ## First solve for the nuisance parameter -- tricky, does not need _A designation!
  beta <- solveBeta( strategy=strategy, includeGE=FALSE )
  ####cat( "testBetaGE_A beta\n" );  print( beta );

  if( is.null(beta) ) {
    cat( "Could not solve for nuisance parameters.\n" )
    return( 1 );
  }

  ## Second, compute the test statistic
  ## - first the precompute piece
  dbnuis <- cpp_gped_statPrecompute_A( beta )
  ####cat( "testBetaGE dbnuis\n" ); print( dbnuis ); ## DEBUG DEBUG
  ####print( dbnuis )
  if( is.null(nrow(dbnuis)) ) {
    cat( "Could not invert nuisance parameters.\n" )
    return( 1 ) ## p-value of 1 safest...
  }
  ##pvalue <- cpp_gped_statCompute_A( solve(dbnuis) )

  #solve.svd <- getFromNamespace( "solve.svd", "fbati" )
  ####cat( "testBetaGE dbnuisInv\n" ); print( solve.svd(svd(dbnuis)) ); ## DEBUG DEBUG
  pvalue <- cpp_gped_statCompute_A( solve.svd( svd( dbnuis ) ) )

  #print( pvalue ) ## DEBUG DEBUG DEBUG DEBUG

  return( pvalue );
}

## They both estimate the same nuisance parameters, so combine them whenever possible
testBetaGE_both <- function( strategy="adaptive" ) {
  ## First set the strategy
  cpp_gped_setStrategy( strategy=strategy ) ## (sets up the permutations)

  ## Solve for the nuisance parameter (same for both)
  beta <- solveBeta( strategy=strategy, includeGE=FALSE )

  ## Compute the test statistic for each
  pvalue <- pvalueA <- 1.0
  ## two df test
  dbnuis <- cpp_gped_statPrecompute( beta )
  if( !is.null(nrow(dbnuis)) )
    pvalue <- cpp_gped_statCompute( solve.svd( svd( dbnuis ) ) )
  ## one df test
  dbnuisA <- cpp_gped_statPrecompute_A( beta )
  if( !is.null(nrow(dbnuisA)) )
    pvalueA <- cpp_gped_statCompute_A( solve.svd( svd( dbnuisA ) ) )

  return( c(pvalue=pvalue, pvalueA=pvalueA) )
}

## And lastly, the routine to link it back with the dataset
cpp_gped_set <- function( pid, id, fathid, mothid,
                         geno, trait, env ) {
  .C( "cpp_gped_set",
     as.integer(pid), as.integer(id), as.integer(fathid), as.integer(mothid),
     as.integer(geno), as.integer(trait), as.double(env), as.integer(length(pid)), DUP=TRUE) #DUP=FALSE )
  return( invisible() )
}
cpp_gped_set_C<- function( pid, id, fathid, mothid,
                          geno, trait, env, cov=NULL ) {
  if( is.null(cov) )
    return( cpp_gped_set( pid, id, fathid, mothid,  geno, trait, env ) )

  ## Because R is an absolute bitch sometimes
  nCov <- ncol(cov)
  if( length(nCov) == 0 ) {
    nCov <- 1
  }

  .C( "cpp_gped_set_C",
     as.integer(pid), as.integer(id), as.integer(fathid), as.integer(mothid),
     as.integer(geno), as.integer(trait), as.double(env), as.double(as.matrix(cov)), as.integer(nCov), as.integer(length(pid)), DUP=TRUE) #DUP=FALSE )
  return( invisible() )
}
fbatge.internal <- function( gped, phe, geno, trait="AffectionStatus", env="env", strategy="adaptive", model="codominant" ) {
  ## Ideally, the ped and phe have already been merged elsewhere
  if( !all(gped$id==phe$id & gped$pid==phe$pid) )
    stop( "fbatge.internal: requires that gped and phe have been merged (and nuclified) previously." )

  ## Find the locations of things
  genoLoc <- which(names(gped)==geno)
  if( length(genoLoc) != 1 ) stop( paste( "fbatge.internal - the specified marker '", geno, "' does not exist in the dataset ped object.", sep="" ) )
  envLoc <- which(names(phe)==env)
  if( length(envLoc) != 1 ) stop( paste( "fbatge.internal - the specified environmental exposure '", env, "' does not exist in the dataset phe object.", sep="" ) )

  ## Trait needs to be handled specially
  traitValue <- NULL
  if( trait=="AffectionStatus" ) {
    traitValue <- gped$AffectionStatus #### phe[[traitLoc]] ## What was I thinking before???
    ## 2=affected, 1=unaffected, 0=missing
    ## FIX, so that 1=affected, 0=unaffected, NA=missing
    traitValue[ traitValue==0 ] <- NA
    traitValue <- traitValue - 1
  }else{
    traitLoc <- which(names(phe)==trait)
    if( length(traitLoc) != 1 ) stop( paste( "fbatge.internal - the specified trait '", trait, "' was not 'AffectionStatus' or one of the value from the dataset phe object.", sep="" ) )

    traitValue <- phe[[traitLoc]]
  }

  ####print( cbind( gped[,c(1:4,genoLoc)], t=traitValue, e=phe[,envLoc] ) )

  ## Eliminate any missing data, do it to the phe, the ped, _and_ the trait
  ## TODO TODO TODO TODO TODO TODO TODO:
  ## --> Ideally, under the 'geno' strategy, we would set missing genotype information
  ##     to be zero, since unaffecteds there aren't used in constructing the test
  ##     statistic, but could still be used to reconstruct S
  ## - ( Missing genotype or missing trait ) and is a child
  kill <- ( is.na(gped[[geno]]) | is.na(traitValue) | is.na(gped[[envLoc]]) ) & gped$idmoth!=0 & gped$idfath!=0
  keep <- !kill
  gped <- gped[keep,]
  phe <- phe[keep,]
  traitValue <- traitValue[keep]

  ## But then for the _parents_ who had missing genotype or phenotype data
  gMiss <- is.na(gped[[geno]])
  if( any(gMiss) ) gped[[geno]][gMiss] <- -1; ## static const int GMISS = -1;
  pMiss <- is.na(traitValue)
  if( any(pMiss) ) traitValue[pMiss] <- -1; ## static const int PMISS = -1;
  eMiss <- is.na(phe[[envLoc]])
  if( any(eMiss) ) phe[[envLoc]][eMiss] <- -999999; ## Won't be used, if it is, this should screw things up sufficiently

  ####cat( "\n\n\nAFTER REMOVING MISSING IN R CODE" )
  ####print( cbind( gped[,c(1:4,genoLoc)], t=traitValue, e=phe[,envLoc] ) )

  ## Now set it
  ## TOCHECK: is it fathid or idfath?
  cpp_gped_set( gped$pid, gped$id, gped$idfath, gped$idmoth,
               gped[[geno]], traitValue, phe[[envLoc]] )
  ####cpp_gped_print() ## DEBUG

  ## Then run the analysis!
  pvalue <- 1
  if( model=="codominant" ) {
    pvalue <- testBetaGE( strategy=strategy )
  }else if( model=="additive" ){
    pvalue <- testBetaGE_A( strategy=strategy )
  }else{
    stop( "model must be 'codominant' or 'additive'." )
  }
  return( pvalue )
}
fbatge.internal_C<- function( gped, phe, geno, cov=NULL, trait="AffectionStatus", env="env", strategy="adaptive", model="codominant" ) {
  ## Ideally, the ped and phe have already been merged elsewhere
  if( !all(gped$id==phe$id & gped$pid==phe$pid) )
    stop( "fbatge.internal_C: requires that gped and phe have been merged (and nuclified) previously." )

  ## Find the locations of things
  genoLoc <- which(names(gped)==geno)
  if( length(genoLoc) != 1 ) stop( paste( "fbatge.internal_C - the specified marker '", geno, "' does not exist in the dataset ped object.", sep="" ) )
  envLoc <- which(names(phe)==env)
  if( length(envLoc) != 1 ) stop( paste( "fbatge.internal_C - the specified environmental exposure '", env, "' does not exist in the dataset phe object.", sep="" ) )
  covLoc <- match( cov, names(phe) )
  if( length(covLoc) != length(cov) ) stop( paste( "fbatge.internal_C - the specified covariates {",paste(cov,collapse=","),"} do not exist in the dataset phe object.", sep="" ) )

  ## Trait needs to be handled specially
  traitValue <- NULL
  if( trait=="AffectionStatus" ) {
    traitValue <- gped$AffectionStatus #### phe[[traitLoc]] ## What was I thinking before???
    ## 2=affected, 1=unaffected, 0=missing
    ## FIX, so that 1=affected, 0=unaffected, NA=missing
    traitValue[ traitValue==0 ] <- NA
    traitValue <- traitValue - 1
  }else{
    traitLoc <- which(names(phe)==trait)
    if( length(traitLoc) != 1 ) stop( paste( "fbatge.internal_C - the specified trait '", trait, "' was not 'AffectionStatus' or one of the value from the dataset phe object.", sep="" ) )

    traitValue <- phe[[traitLoc]]
  }

  ## Eliminate any missing data, do it to the phe, the ped, _and_ the trait
  ## TODO TODO TODO TODO TODO TODO TODO:
  ## --> Ideally, under the 'geno' strategy, we would set missing genotype information
  ##     to be zero, since unaffecteds there aren't used in constructing the test
  ##     statistic, but could still be used to reconstruct S
  ## - ( Missing genotype or missing trait ) and is a child
  is.na.cov <- rep( FALSE, nrow(phe) )
  for( cCol in covLoc )
    is.na.cov <- is.na.cov | is.na(phe[[cCol]])

  kill <- ( is.na(gped[[geno]]) | is.na(traitValue) | is.na(gped[[envLoc]]) | is.na.cov ) & gped$idmoth!=0 & gped$idfath!=0
  keep <- !kill
  gped <- gped[keep,]
  phe <- phe[keep,]
  traitValue <- traitValue[keep]
  is.na.cov <- is.na.cov[keep]

  ## But then for the _parents_ who had missing genotype or phenotype data
  gMiss <- is.na(gped[[geno]])
  if( any(gMiss) ) gped[[geno]][gMiss] <- -1; ## static const int GMISS = -1;
  pMiss <- is.na(traitValue)
  if( any(pMiss) ) traitValue[pMiss] <- -1; ## static const int PMISS = -1;
  eMiss <- is.na(phe[[envLoc]])
  if( any(eMiss) ) phe[[envLoc]][eMiss] <- -999999; ## Won't be used, if it is, this should screw things up sufficiently (only done to parents, where this does nothing)
  if( any(is.na.cov) ) phe[ is.na.cov, covLoc ] <- -999999; ## Again, won't be used

  #cat( "\n\n\nAFTER REMOVING MISSING IN R CODE" )
  #print( cbind( gped[,c(1:4,genoLoc)], t=traitValue, e=phe[,envLoc], phe[,covLoc] ) ); stop(); ## DEBUG ONLY

  ## Now set it
  ## TOCHECK: is it fathid or idfath?
  cpp_gped_set_C( gped$pid, gped$id, gped$idfath, gped$idmoth,
                 gped[[geno]], traitValue, phe[[envLoc]], phe[,covLoc] )
  #cat( "DAMN COVARIATES", cpp_gped_numCovariates() ); stop();
  #cpp_gped_print(); stop(); ## DEBUG

  ## Then run the analysis!
  pvalue <- 1
  if( model=="codominant" ) {
    pvalue <- testBetaGE( strategy=strategy )
  }else if( model=="additive" ){
    pvalue <- testBetaGE_A( strategy=strategy )
  }else{
    stop( "model must be 'codominant' or 'additive'." )
  }
  return( pvalue )
}

####################
## GESim routines ##
####################

cpp_gesim_print <- function() { .C("cpp_gesim_print") }

cpp_gesim_set <- function( numParents=2, numOffspring=1, numFam=500,
                          minAffected=1, maxAffected=1,
                          afreq=0.2, geneticModel="additive",
                          link="log",
                          beta=c(-4, 0, 0.25, 0.5), ## b0 bge bg be
                          env="dichotomous", envCutoff=0.5, rho=0, ## rho is pairwise correlations (symmetric correlation matrix
                          betaCov=0, distCov="normal",
                          phenoCor=0, phenoCutoff=0,
                          markerCor=0, markerAfreq=afreq,
                          phenoOR=1 ) {
  sigma <- matrix( rho, nrow=numOffspring, ncol=numOffspring )
  diag(sigma) <- 1
  cholesky <- chol(sigma)

  if( env=="dichotomous" ) {
    ## Retranslate cutoff based on normal quantiles..
    envCutoff <- qnorm( envCutoff ) ## NOT CORRECT!!!!
  }

  if( phenoCutoff != 0 )
    phenoCutoff <- qnorm( 1-phenoCutoff ) ## Opposite direction as dichotomization of env
  phenoMat <- matrix( markerCor, nrow=numOffspring, ncol=numOffspring )
  diag(phenoMat) <- 1
  pheno.cholesky <- chol(phenoMat)
  #print( pheno.cholesky )
  #print( phenoCutoff )
  #print( markerCor )

  .C("cpp_gesim_set",
     as.integer(numParents), as.integer(numOffspring), as.integer(numFam),
     as.integer(minAffected), as.integer(maxAffected),
     as.double(afreq), as.character(geneticModel),
     as.character(link),
     as.double(beta), as.integer(length(beta)),
     as.character(env), as.double(envCutoff), as.double(cholesky), as.integer(numOffspring),
     as.double(betaCov), as.character(distCov),
     as.double(pheno.cholesky), as.double(phenoCutoff),
     as.double(markerCor), as.double(markerAfreq),
     as.double(phenoOR),
     DUP=TRUE) ## Character vars must be duplicated. OK then.

  return( invisible() )
}

cpp_gesim_draw <- function() { .C("cpp_gesim_draw") }
cpp_gesim_clear <- function() { .C("cpp_gesim_clear") }

cpp_gesim_debug <- function() {
  #cpp_gesim_set()
  #cpp_gesim_print()
  #cpp_gesim_draw()
  #cpp_gped_print()

  #cpp_gesim_clear()
  #cpp_gesim_set( numOffspring=3, minAffected=1, maxAffected=2 )
  #cpp_gesim_print()
  #cpp_gesim_draw()
  #cpp_gped_print()

  ## And how about a mixture
  cpp_gesim_set( numFam=10 );
  cpp_gesim_set( numFam=20, numOffspring=3, numParents=0, minAffected=1, maxAffected=2 )
  cpp_gesim_print()
  cpp_gesim_draw()
  cpp_gped_print()
}
#cpp_gesim_debug()


estEqDebug <- function() {
  #set.seed(13)
  #cpp_gesim_set()
  #cpp_gesim_print()
  #cpp_gesim_draw()
  #cpp_gped_print()
  #cpp_gped_setStrategy("geno")
  #cat( "beta=" )
  #print( solveBeta( strategy="geno" ) )

  b0 <- -4  ## Shouldn't matter
  bge <- 0.5
  bg <- 0.25
  be <- 0.5
  beta <- c(b0,bge,bg,be)

  set.seed(13)
  #cpp_gesim_set(numFam=500,beta=beta) ## Trios
  cpp_gesim_set( beta=beta, numOffspring=2, minAffected=1, maxAffected=1, numParents=0 ) ## DSP's
  NSIM <- 1000
  pvalueLT <- 0
  for( i in 1:NSIM ) {
    cat( i, "" )
    cpp_gesim_draw()
    ####cpp_gped_setStrategy("geno")
    #cat( "beta = " ); print( solveBeta( strategy="geno" ) )
    pvalue <- testBetaGE( strategy="geno" )
    cat( "pvalue =", pvalue, "\n" )  ## Dieing on adaptive, but the parameter was unidentifiable...
    pvalueLT <- pvalueLT + as.numeric(pvalue<0.05)
  }
  cat( "\nEmpirical Alpha/Power = ", pvalueLT/NSIM, "\n" )
}
#estEqDebug()
estEqDebug2 <- function() {
  b0 <- -4  ## Shouldn't matter
  bge <- 0.5  #.5
  bg <- 0.25
  be <- 0.5
  beta <- c(b0,bge,bg,be)

  set.seed(13)
  cpp_gesim_set( beta=beta, numOffspring=2, minAffected=1, maxAffected=1, numParents=0, link="logit", numFam=500 ) ## DSP's
  cpp_gesim_print()
  NSIM <- 1000
  pvalueLT <- 0
  for( i in 1:NSIM ) {
    cat( i, "" )
    cpp_gesim_draw()
    #cpp_gped_setStrategy("pheno")
    #cat( "beta = " ); print( solveBeta( strategy="pheno" ) )
    #stop()
    #pvalue <- testBetaGE( strategy="geno" )
    #pvalue <- testBetaGE( strategy="pheno" )
    pvalue <- testBetaGE( strategy="adaptive" )
    cat( "pvalue =", pvalue, "\n" )  ## Dieing on adaptive, but the parameter was unidentifiable...
    pvalueLT <- pvalueLT + as.numeric(pvalue<0.05)
  }
  cat( "\nEmpirical Alpha/Power = ", pvalueLT/NSIM, "\n" )
}
#estEqDebug2()
