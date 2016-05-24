#fixMarkerNames <- getFromNamespace( "fixMarkerNames", "fbati" )
#ADDITIVE <- getFromNamespace( "ADDITIVE", "fbati" )
#DOMINANT <- getFromNamespace( "DOMINANT", "fbati" )
#RECESSIVE <- getFromNamespace( "RECESSIVE", "fbati" )


## Maybe the only thing useful we learned
solve.svd <- function( x, b=NULL, tol=1e-15 ) {
  sub <- x$d >= tol*max(x$d)
  d <- 1/x$d[sub]
  if (is.null(b)) {
    ## Then compute the generalized inverse
    b <- x$v[,sub] %*% (d * t(x$u[,sub]))
  } else {
    ## Give the least squares solution
    b <- as.matrix(b)
    b <- x$v[,sub] %*% (d * (t(x$u[,sub]) %*% b))
  }
  attr( b, 'rank' ) <- length(d)
  attr( b, 'rcond' ) <- min(d) / max(d)
  return( b )
}

## The joint routine
fbatj <- function( ped=NULL, phe=NULL,
                   data=mergePhePed(ped,phe),
                   marker = NULL,
                   trait = "AffectionStatus",
                   env = "env",
                   model = "additive",
                   mode = "univariate",
                   fixNames = TRUE,
                   verbose = FALSE ) {
  if( is.null(data) )
    return(fbatjGUI())

  if( mode != "univariate" )
    stop( "Only 'univariate' is supported for mode." );

  if( is.na(trait) || is.null(trait) )
    stop( "Trait must be specified (e.g. AffectionStatus)." )

  ##################################
  ## TAKEN DIRECTLY FROM 'ibat.R' ##

  ## New preprocessing 01/14/2008 safety precaution
  data <- data[ order(data[,1]), ] ## order it so pids are grouped together
  ## NEW preprocessing - always nuclify the data...
  data <- nuclifyMerged( data )

  ## get the environment column
  envCols <- match( env, names(data) )
  if( length(envCols) < 1 )
    stop( "At least one environmental exposure must be selected (each will be tested seperately)." )
  ##print( "envCols" )
  #print( envCols )
  #print( names(data) )

  ## get the markers
  if( is.null(marker) ) {
    if( is.null(phe) )
      stop( "If using 'data' as the option, then you must specify the marker locations in _pairs_." )

    marker <- 7:ncol(ped)
  }


  markerNames <- NULL
  if( is.character(marker) ) {
    ## Then its strings of which markers
    potentialMarkerNames <- fixMarkerNames( names(data) ) ## some of these aren't markers...
    m <- match( marker, potentialMarkerNames )
    if( length(m)==0 )
      stop( paste( "Marker names do not exist, choices are:", marker, collapse=" " ) )
    n <- length(m)*2
    marker <- rep(0,n)
    marker[seq(1,n,by=2)] <- m
    marker[seq(2,n,by=2)] <- m+1

    markerNames <- potentialMarkerNames[m]
  }else{
    ## fix up the marker names
    markerNames <- names(data)[ marker[seq(from=1,to=length(marker),by=2)] ]
    if( fixNames )
      markerNames <- fixMarkerNames( markerNames )
  }

  if( length(marker) %% 2 == 1 )
    stop("Markers must be specified in _pairs_.")

  if( mode!="multivariate" && length(marker)>2 ) {
    ## Go into univariate!
    #print( "Univariate instigated." )
    ares <- NULL
    for( markerN in markerNames ) {
      res <- fbatj( data=data, marker=markerN, trait=trait, env=env, model=model, fixNames=fixNames, verbose=verbose )
      if( is.null(ares) ) {
        ares <- res
      }else{
        ##ares <- rbind( res, ares )
        ares <- rbind( ares, res )
      }
    }
    return( ares )
  }

  ## Convert model into constant format
  model <- c(ADDITIVE,DOMINANT,RECESSIVE)[match( model, c("additive","dominant","recessive") )]
  if( length(model)==0 )
    stop( "model must be 'additive', 'dominant', or 'recessive'" )

  ## TAKEN DIRECTLY FROM 'ibat.R' ##
  ##################################

  ## Update the trait columns...
  traitCols <- match( trait, names(data) )
  if( length(traitCols) < 1 )
    stop( "The trait must be specified." )

  ## handle the recoding, which is necessary for my part of the code... hmm... this test only makes sense under an additive or a genotype coding...
  for( m in seq( from=1, to=length(marker), by=2 ) ){
    m0pos <- marker[m]
    m1pos <- marker[m+1]
    alleles <- setdiff( unique( c(data[,m0pos],data[,m1pos]) ), 0 ) ## elts excluding zero
    if( length(alleles) > 2 )
      stop( paste( "Currently only biallelic markers are supported, marker", markerNames[m], "values are not supported." ) )

    ## recode the dataset
    wh1 <- data[,m0pos]==alleles[1]
    wh2 <- data[,m0pos]==alleles[2]
    data[wh1,m0pos] <- 1
    if( !all(is.na(wh2)) )
      data[wh2,m0pos] <- 2

    wh1 <- data[,m1pos]==alleles[1]
    wh2 <- data[,m1pos]==alleles[2]
    data[wh1,m1pos] <- 1
    if( !all(is.na(wh2)) )
      data[wh2,m1pos] <- 2
  }

  if( trait=="AffectionStatus" ) {
    ##print( "AffectionStatus recoded." )
    #data[ data[,traitCols]==0, traitCols ] <- NA
    ##data[ data[,traitCols]==1, traitCols ] <- 0
    ##data[ data[,traitCols]==2, traitCols ] <- 1
    temp <- rep(NA,nrow(data))
    temp[ data[[traitCols]]==1 ] <- 0
    temp[ data[[traitCols]]==2 ] <- 1
    data[[traitCols]] <- temp
  }

  ## Can't be handled yet
  if( length(traitCols) > 1 ) stop( "Can only handle a single trait." )
  if( length(envCols) > 1 ) stop( "Can only handle a single environmental exposure." )

  ## Create the output pieces
  NN <- length(marker)
  a <- rep( as.double(0), NN )
  b <- matrix( as.double(0), nrow=NN, ncol=NN )
  numInf <- rep( as.double(0), 1 )   ## This doesn't work with multi-marker mode yet...
  as.matrix( data )
  results <- REXP_joint( as.matrix(data), marker-1, traitCols-1, envCols-1, model, a, b, numInf )

  a <- results$a
  b <- matrix( results$b, nrow=NN, ncol=NN )
  numInf <- results$numInf
  #print( results )

  ## Main effect would almost come for free (debugging -- it is correct at least!)
  #print( "length(marker)" )
  #print( length(marker) )
  #print( "ME Result = " )
  #print( a[1]^2 / b[1,1] )

  #print( "IDEALLY WE WANT TO DO SOMETHING LIKE df=rank(b)." )
  ##rank <- qr(b)$rank
  #print( rank )

  if( verbose ) {
    print( "A" )
    print( a )
    print( "B" )
    print( b )
  }

  ## Give a try to removing the bloody 0 rows...
  kill <- c()
  for( rr in 1:nrow(b) ) {
    if( all( b[rr,]==0 ) )
      kill <- c(kill, rr)
  }
  #print( kill )


  ##stat <- a %*% solve(b) %*% a
  ##stat <- a %*% solve( b, a )

  #print( solve( b, a ) )
  #print( solve.svd( svd(b), a ) )
  #return(NULL)

  #b_svd <- svd(b)
  #genInv <- b_svd$v %*% diag(1/b_svd$d) %*% b_svd$u
  #range <- 1:min( 5, ncol(genInv) )
  #print( genInv[range,range] )

  #genInv2 <- solve.svd(b_svd)
  #print( genInv2[range,range] )

  #print( solve( b ) )
  #stop()
  ##stat <- a %*%  b_svd$v %*% diag(1/b_svd$d) %*% b_svd$u  %*% a
  ##print( stat )
  #stat <- a %*% genInv2 %*% a
  #print( "stat" )
  #print( stat )

  post <- solve.svd( svd(b), a )
  rank <- attr(post,"rank")
  stat <- a %*% post
  #print( stat )
  #print( rank )

  ## DEBUG BEGIN (usual emp var)
  ##stat2 <- a[1]^2 / b[1,1]
  ##pval2 <- pchisq( stat2, df=1, lower.tail=FALSE )
  ##cat( "Univariate pvalue", pval2, "\n" )
  ## DEBUG END

  pvalue <- pchisq( stat, df=rank, lower.tail=FALSE )
  if( length(markerNames) > 1 )
    return( data.frame( stat=stat, pvalue=pvalue, rank=rank ) )
  return( data.frame( marker=markerNames, numInf=numInf, stat=stat, pvalue=pvalue, rank=rank ) )
}