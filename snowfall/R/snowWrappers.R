## Wrappers for Snow function.
##
## The wrappers do the following: decide whether we run in parallel or
## sequential mode.
## In parallel mode the according Snow functions are used.
## In sequential mode, if it makes sense, the sequential counterparts
## of the Snow functions are used.

##****************************************************************************
## Wrapper for: clusterSplit
##****************************************************************************
sfClusterSplit <- function( seq ) {
  sfCheck();

  if( sfParallel() )
    return( clusterSplit( sfGetCluster(), seq ) )
  ## In sequential mode return a list with everything in element 1 (means:
  ## everything is run on one node).
  else
    return( list( seq ) )
}

##****************************************************************************
## Wrapper for: clusterCall
##
## Catches for errors. Return them or stop immidiately.
##****************************************************************************
sfClusterCall <- function( fun, ..., stopOnError=TRUE ) {
  sfCheck();

  if( !checkFunction( fun, stopOnError=FALSE ) ) {
    if( stopOnError )
      stop( "No function or not defined object in sfClusterCall" )
    else {
      warning( "No function or not defined object in sfClusterCall" )
      return( NULL )
    }
  }

  if( sfParallel() ) {
    ## Exec via Snow.
    result <- clusterCall( sfGetCluster(), fun, ... )

    ## Not enough results?
    ## @TODO Check if this test is needed
    if( length( result ) != sfCpus() ) {
      if( stopOnError )
        stop( paste( "Error in sfClusterCall (not all slaves responded).\n",
                     "Call from: ", as.character( sys.call( -1 ) ) ) )
      else {
        message( paste( "Error in sfClusterCall (not all slaves responded).\n",
                        "Call from: ", as.character( sys.call( -1 ) ) ) )
        return( result );
      }
    }

    ## Check if snow throw an exception on any of the slaves.
    if( !all( checkTryErrorAny( result ) ) ) {
      errorsTxt <- sapply( which( inherits( result, "try-error" ) ), function(x) result[[x]] )

      message( "EXCEPTION INFOS:" )
      message( paste( errorsTxt, collapse="\n" ) )
      
      if( stopOnError ) {
        stop( paste( "Error in sfClusterCall (catched TRY-ERROR).\n",
                     "Call from: ", as.character( sys.call( -1 ) ) ) )
      }
      else {
        message( paste( "Error in sfClusterCall (catched TRY-ERROR).\n",
                        "Call from: ", as.character( sys.call( -1 ) ) ) )
        return( result )
      }
    }

    return( result )
  }
  ## Sequential mode.
  else
    return( do.call( fun, list( ... ) ) )
}

##****************************************************************************
## Wrapper for: clusterEvalQ - renamed as indeed "eval" is executed and not
## "evalq".
##****************************************************************************
sfClusterEval <- function( expr, stopOnError=TRUE ) {
  sfCheck();

  if( sfParallel() ) {
    return( sfClusterCall( eval, substitute( expr ), env=globalenv(),
                           stopOnError=stopOnError ) )
  }
  else {
    ## Problems can arise through "enclos", which is default set to parent
    ## and therefore here, too: on this way local variables (higher environments
    ## are visible, which badly are not visible in parallel runs...).
    ## There should be a fix or something.
    return( eval( expr, envir=globalenv(), enclos=parent.frame() ) )
  }
}

## Snows clusterEvalQ uses "eval" and not "evalq", so this wrapper is an alias.
sfClusterEvalQ <- function( expr ) return( sfClusterEval( expr ) )

##****************************************************************************
## Wrapper for: clusterMap.
## Currently not used.
##****************************************************************************
sfClusterMap <- function( fun, ..., MoreArgs=NULL, RECYCLE=TRUE )
  stop( "Currently no wrapper for clusterMap" )

##****************************************************************************
## Wrapper for: clusterApply (snow parallel) - lapply (sequential)
## Adds additional warnings before the execution (esp. in sequential mode,
## where exec works fine but can cause problems runnin in parallel).
##
## PARAMETERS: Parameters like clusterApply
## RETURN:     Result
##****************************************************************************
sfClusterApply <- function( x, fun, ... ) {
  sfCheck();

  checkFunction( fun )

  ## However snow limits list size to cluster nodes in "normal"
  ## execution.
  ## This is a fatal error in parallel mode and a warning in sequential.
  if( length( x ) > sfCpus() ) {
    if( sfParallel() )
      stop( "More list entries as nodes => use sfClusterApplyLB instead. See Snow/Snowfall documentation." )
    else
      warning( "More list entries as nodes => causes error in parallel mode. use sfClusterApplyLB instead." )
  }
  
  if( sfParallel() )
    return( clusterApply( sfGetCluster(), x, fun, ... ) )
  else
    return( lapply( x, fun, ... ) )
}

##****************************************************************************
## Wrapper for: clusterApplyLB (snow parallel) - lapply (sequential)
##
## PARAMETERS: Parameters like clusterApply
## RETURN:     Result
##****************************************************************************
sfClusterApplyLB <- function( x, fun, ... ) {
  sfCheck();

  checkFunction( fun )

  if( sfParallel() )
    return( clusterApplyLB( sfGetCluster(), x, fun, ... ) )
  else
    ## array... korrigieren.
    return( lapply( x, fun, ... ) )
}

##****************************************************************************
## Also snow-Handler handling is hidden to the user.
##
## Wrapper for: parLappy (snow parallel) - lapply (sequential)
##
## As lapply parameters were inkonsitent ("x"/"fun") they were corrected to
## ""x"/"fun".
##
## PARAMETERS: Parameters like lapply
## RETURN:     Result
##****************************************************************************
sfLapply <- function( x, fun, ... ) {
  sfCheck()

  checkFunction( fun )
  
  if( sfParallel() )
    return( parLapply( sfGetCluster(), x, fun, ... ) )
  else
    return( lapply( x, fun, ... ) )
}

##****************************************************************************
## Wrapper for: parSapply (snow parallel) - sapply (sequential)
##
## PARAMETERS: Parameters like sapply
## RETURN:     Result
##****************************************************************************
sfSapply <- function( x, fun, ..., simplify=TRUE, USE.NAMES=TRUE ) {
  sfCheck()

  checkFunction( fun )

  if( sfParallel() )
    return( parSapply( sfGetCluster(), x, fun, ..., simplify=simplify, USE.NAMES=USE.NAMES ) )
  else
    return( sapply( x, fun, ..., simplify=simplify, USE.NAMES=USE.NAMES ) )
}

##****************************************************************************
## Wrapper for: parApply (snow parallel) - apply (sequential)
##
## PARAMETERS: Parameters like apply
## RETURN:     Result
##****************************************************************************
sfApply <- function( x, margin, fun, ... ) {
  sfCheck()

  checkFunction( fun )

  if( sfParallel() )
    return( parApply( sfGetCluster(), x, margin, fun, ... ) )
  else
    return( apply( x, margin, fun, ... ) )
}

sfRapply <- function( x, fun, ... ) {
  stop( "sfRapply does not exists yet. Use Snow's parRapply instead." )
  return( invisible( NULL ) );
}

sfCapply <- function( x, fun, ... ) {
  stop( "sfCapply does not exists yet. Use Snow's parCapply instead." )
  return( invisible( NULL ) );
}

##****************************************************************************
## Wrapper for: parMM (snow parallel) - %*% (sequential)
##
## PARAMETERS: Matrix a, Matrix b
## RETURN:     Result
##****************************************************************************
sfMM <- function( a, b ) {
  sfCheck();

  if( sfParallel() )
    return( parMM( sfGetCluster(), a, b ) )
  else
    return( a %*% b )
}

##****************************************************************************
## Wrappers for the two uniform RNGs used in snow.
## Basically, at the moment these are not used in sequential (means: none
## of the two is included here for sequential execution).
## @TODO Sequential use of the RNGs.
##****************************************************************************
sfClusterSetupSPRNG <- function( seed = round( 2^32 * runif(1) ),
                                 prngkind = "default", para = 0, ... ) {
  sfCheck();

  if( sfParallel() )
    clusterSetupSPRNG( sfGetCluster(), seed, prngkind, para, ... )
  else {
    warning( paste( "Uniform random number streams (currently) not available in serial execution.",
                    "Random numbers may differ in serial & parallel execution." ) )
    set.seed( seed )
  }
}

sfClusterSetupRNGstream <- function( seed=rep( 12345, 6 ), ... ) {
  sfCheck();

  if( sfParallel() )
    clusterSetupRNGstream( sfGetCluster(), seed=seed, ... )
  else {
    warning( paste( "Uniform random number streams (currently) not available in serial execution.",
                    "Random numbers may differ in serial & parallel execution." ) )
    set.seed( seed[1] )
  }
}

sfClusterSetupRNG <- function( type="RNGstream", ... ) {
  sfCheck();

  if( sfParallel() )
    clusterSetupRNG( sfGetCluster(), type=type, ... )
  else {
    warning( paste( "Uniform random number streams (currently) not available in serial execution.",
                    "Random numbers may differ in serial & parallel execution." ) )
  }
}
