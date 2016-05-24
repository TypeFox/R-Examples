## Function auto-generated...
REXP_fbatmeev <- function(  data, marker, trait, model, RET_stat, RET_numInf ) {
  rexport_result <- .C( 'eREXP_fbatmeev',  as.double(data), as.integer(dim(data)), as.double(marker), as.integer(length(marker)), as.double(trait), as.double(model), as.double(RET_stat), as.integer(length(RET_stat)), as.double(RET_numInf), as.integer(length(RET_numInf)), dup=FALSE, NAOK=TRUE )
  list(   stat = rexport_result[[7]],  numInf = rexport_result[[9]] )
}

## The fbat main effects empirical variance routine
fbatmeev <- function( ped=NULL, phe=NULL,
                      data=mergePhePed(ped,phe),
                      marker=NULL,
                      trait="AffectionStatus",
                      model="additive",
                      fixNames=TRUE,
                      verbose = FALSE ) {





  ##################################
  ## TAKEN DIRECTLY FROM 'ibat.R' ##

  ## New preprocessing 01/14/2008 safety precaution
  data <- data[ order(data[,1]), ] ## order it so pids are grouped together
  ## NEW preprocessing - always nuclify the data...
  data <- nuclifyMerged( data )

  ## NO ENVIRONMENTAL COLUMN --> DELETED

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
  afreq <- rep(NA,length(marker)/2) ## new
  for( m in seq( from=1, to=length(marker), by=2 ) ){
    m0pos <- marker[m]
    m1pos <- marker[m+1]
    alleles <- setdiff( unique( c(data[,m0pos],data[,m1pos]) ), 0 ) ## elts excluding zero
    if( length(alleles) > 2 )
      stop( paste( "Currently only biallelic markers are supported, marker", markerNames[m], "values are not supported." ) )

    ## recode the dataset

    sumwh1 <- 0
    sumwh2 <- 0

    wh1 <- data[,m0pos]==alleles[1]
    wh2 <- data[,m0pos]==alleles[2]
    data[wh1,m0pos] <- 1
    if( !all(is.na(wh2)) )
      data[wh2,m0pos] <- 2

    sumwh1 <- sumwh1 + sum(wh1,na.rm=TRUE)
    sumwh2 <- sumwh2 + sum(wh2,na.rm=TRUE)

    wh1 <- data[,m1pos]==alleles[1]
    wh2 <- data[,m1pos]==alleles[2]
    data[wh1,m1pos] <- 1
    if( !all(is.na(wh2)) )
      data[wh2,m1pos] <- 2

    sumwh1 <- sumwh1 + sum(wh1,na.rm=TRUE)
    sumwh2 <- sumwh2 + sum(wh2,na.rm=TRUE)

    afreqNew <- sumwh1 / (sumwh1 + sumwh2)
    if( afreqNew > 0.5 )
      afreqNew <- 1 - afreqNew
    afreq[(m+1)/2] <- afreqNew
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

  ## Create the output pieces
  NN <- length(marker)
  stats <- rep( as.double(0), NN/2 )
  numInf <- rep( as.double(0), NN/2 )
  res <- REXP_fbatmeev( as.matrix(data), marker-1, traitCols-1, model, stats, numInf )
  pvalues <- pchisq( res$stat, df=1, lower.tail=FALSE )  ## Damn you and your stupid lower.tail...
  try( names(pvalues) <- ped.markerNames(ped) )

  return( data.frame(marker=names(pvalues), afreq=afreq, numInf=res$numInf, pvalue=pvalues) )
}













## Function auto-generated...
REXP_fbatme <- function(  data, marker, trait, model, RET_stat, RET_numInf ) {
  rexport_result <- .C( 'eREXP_fbatme',  as.double(data), as.integer(dim(data)), as.double(marker), as.integer(length(marker)), as.double(trait), as.double(model), as.double(RET_stat), as.integer(length(RET_stat)), as.double(RET_numInf), as.integer(length(RET_numInf)), dup=FALSE, NAOK=TRUE )
  list(   stat = rexport_result[[7]],  numInf = rexport_result[[9]] )
}

## The fbat main effects empirical variance routine
fbatme <- function( ped=NULL, phe=NULL,
                    data=mergePhePed(ped,phe),
                    marker=NULL,
                    trait="AffectionStatus",
                    model="additive",
                    fixNames=TRUE,
                    verbose = FALSE ) {





  ##################################
  ## TAKEN DIRECTLY FROM 'ibat.R' ##

  ## New preprocessing 01/14/2008 safety precaution
  data <- data[ order(data[,1]), ] ## order it so pids are grouped together
  ## NEW preprocessing - always nuclify the data...
  data <- nuclifyMerged( data )

  ## NO ENVIRONMENTAL COLUMN --> DELETED

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
  afreq <- rep(NA,length(marker)/2) ## new
  for( m in seq( from=1, to=length(marker), by=2 ) ){
    m0pos <- marker[m]
    m1pos <- marker[m+1]
    alleles <- setdiff( unique( c(data[,m0pos],data[,m1pos]) ), 0 ) ## elts excluding zero
    if( length(alleles) > 2 )
      stop( paste( "Currently only biallelic markers are supported, marker", markerNames[m], "values are not supported." ) )

    ## recode the dataset

    sumwh1 <- 0
    sumwh2 <- 0

    wh1 <- data[,m0pos]==alleles[1]
    wh2 <- data[,m0pos]==alleles[2]
    data[wh1,m0pos] <- 1
    if( !all(is.na(wh2)) )
      data[wh2,m0pos] <- 2

    sumwh1 <- sumwh1 + sum(wh1,na.rm=TRUE)
    sumwh2 <- sumwh2 + sum(wh2,na.rm=TRUE)

    wh1 <- data[,m1pos]==alleles[1]
    wh2 <- data[,m1pos]==alleles[2]
    data[wh1,m1pos] <- 1
    if( !all(is.na(wh2)) )
      data[wh2,m1pos] <- 2

    sumwh1 <- sumwh1 + sum(wh1,na.rm=TRUE)
    sumwh2 <- sumwh2 + sum(wh2,na.rm=TRUE)

    afreqNew <- sumwh1 / (sumwh1 + sumwh2)
    if( afreqNew > 0.5 )
      afreqNew <- 1 - afreqNew
    afreq[(m+1)/2] <- afreqNew
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

  ## Create the output pieces
  NN <- length(marker)
  stats <- rep( as.double(0), NN/2 )
  numInf <- rep( as.double(0), NN/2 )
  res <- REXP_fbatme( as.matrix(data), marker-1, traitCols-1, model, stats, numInf )
  pvalues <- pchisq( res$stat, df=1, lower.tail=FALSE )  ## Damn you and your stupid lower.tail...
  try( names(pvalues) <- ped.markerNames(ped) )

  return( data.frame(marker=names(pvalues), afreq=afreq, numInf=res$numInf, pvalue=pvalues) )
}
