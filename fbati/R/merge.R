## merging datasets
mergePhePed <- function( ped, phe ) {
  if( is.null(phe) || is.null(ped) ) {
    ## Should we warn?
    return( NULL )
  }
  data <- merge( ped, phe, all.x=TRUE )
  data <- data[ order( data$pid, data$id ), ]
  return( data )
}

## data is a merged pedigree and phenotype object
## assumes the data is sorted
## This probably isn't very fast, but it'll get the job done.
#nuclifyMerged <- function( data, debug=FALSE ) {
#  res <- NULL;
#  cpid <- 1; ## will be renumbered...
#
#  ## go accross the families
#  cat( "nuclifyMerged family completion (by pid): " )
#  for( pid in unique(data$pid) ) {
#    family <- data[ data$pid==pid, ]
#    for( idfath in unique(family$idfath) ) {
#      ffamily <- family[family$idfath==idfath,]
#      for( idmoth in unique(family$idmoth) ) {
#        mfamily <- ffamily[ffamily$idmoth==idmoth,]
#
#        ## this one should be a nuclearized family...
#        if( nrow(mfamily)> 0 ) {
#
#          ## here's where we can add in some empty clauses
#          ##  in case the mother / father is completely missing
#          resAddi <- rbind( mfamily,
#                            family[ family$id==idfath, ],
#                            family[ family$id==idmoth, ] )
#
#          if( debug )
#            resAddi$oldPid <- resAddi$pid ## so I can make sure it chopped right
#          resAddi$pid <- cpid
#          cpid <- cpid+1
#
#          if( is.null(res) ) {
#            res <- resAddi
#          }else{
#            res <- rbind( res, resAddi )
#          }
#        }
#      }
#    }
#
#    cat( pid, "" )
#  }
#  cat( "\n" )
#
#  return(res)
#}

## New c++ code does this a lot faster
nuclifyMerged <- function( data, OUT_MULT=2 ) {
  OUT_MAXSIZE <- nrow(data)*OUT_MULT ## Actually should be fine for _all_ datasets?

  dataOut <- matrix(as.double(0),nrow=OUT_MAXSIZE,ncol=ncol(data))
  ##dataOutDim <- as.integer(dim(dataOut))
  dataOutDim <- as.integer( c( dim(dataOut)[1], dim(dataOut)[2] ) ) ## Bug in R?, see strataReduce
  failure <- as.integer(0)
  res = .C("nuclify",
     as.double(as.matrix(data)), as.integer(dim(data)),
     dataOut, dataOutDim, failure,
     DUP=TRUE, NAOK=TRUE) #DUP=FALSE, NAOK=TRUE )
  dataOut = res[[3]]
  dataOutDim = res[[4]]
  failure = res[[5]]

  if(failure == 1) {
    ## output isn't big enough! not the best implementation...
    if(OUT_MULT > 100)
      cat("Probable error in nuclification unless there are really strange nuclear families.\n")

    return(nuclifyMerged(data=data, OUT_MULT=OUT_MULT*2))
  }

  #print( "nuclifyMerged" )
  #print( dataOutDim )

  if( dataOutDim[1] == 0 )
    return( NULL )

  dataOut <- data.frame( dataOut[ 1:dataOutDim[1], 1:dataOutDim[2] ] )
  names(dataOut) <- names(data)
  return(dataOut)
}
## returns a list of a pedigree and a phenotype object
##  that have been nuclified
nuclify <- function( ped, phe ) {
  ped <- sort( ped )
  phe <- sort( phe )

  ## nuclify the merged data
  data <- mergePhePed( ped, phe )
  ndata <- nuclifyMerged( data )

  ## and unmerge back into a pedigree and phenotype object
  nped <- ndata[, 1:ncol(ped)]
  class(nped) <- c("ped","data.frame")
  nphe <- ndata[, c(1,2,(ncol(ped)+1):ncol(ndata))]
  class(nphe) <- c("phe","data.frame")
  return( list( ped=nped,
                phe=nphe ) )
}


## Returns a more strata friendly family
strataReduce <- function( data, envCol, m0, m1=m0+1, maxSib=3 ) {
  OUT_MAXSIZE <- nrow(data) ## here, can't be any larger

  #print( "MAXSIB" )
  #print( maxSib )

  ## NASTY little new introduction of R!!! ARGH!!!
  dataOut <- matrix( as.double(0), nrow=OUT_MAXSIZE, ncol=ncol(data) )
  ##dataOutDim <- as.integer(dim(dataOut))
  ##dataOutDim <- as.integer(c(5,9))
  dataOutDim <- as.integer( c( dim(dataOut)[1], dim(dataOut)[2] ) ) ## Bug in R?

  res = .C( "strataReduce",
      as.double(as.matrix(data)), as.integer(dim(data)),
      dataOut, dataOutDim,
      as.integer(envCol-1),
      as.integer(m0-1), as.integer(m1-1),
      as.integer(maxSib),
      DUP=TRUE, NAOK=TRUE) #DUP=FALSE, NAOK=TRUE )
  dataOut = res[[3]]
  dataOutDim = res[[4]]

  if( dataOutDim[1] == 0 )
    return(NULL)

  dataOut <- data.frame( dataOut[ 1:dataOutDim[1], 1:dataOutDim[2] ] )

  names(dataOut) <- names(data)
  return( dataOut )
}

###################################
## DEBUG ONLY (for strataReduce) ##
debug.merge <- function() {
  ## Leaving off with needing to debug by using the srFam(...) routine
  ##  and then we need to implement it in the general analysis routines
  #dyn.load( "../src/nuclify.so" )

  ## Creates a family with the specified markers and statuses
  ## NOTE: affection is false/true, whereas it is coded 1/2 in the ped file
  createFam <- function( pa=c(0,0), pb=c(0,0), ca, cb, caffected=rep(TRUE,length(ca)), env=1:length(ca) ) {
    ## pid, id, idfath, idmoth, sex, affection, m0a, m0b
    numC <- length(ca)
    return( data.frame( pid=rep(1,2+numC), id=1:(2+numC), idfath=c(0,0,rep(1,numC)), idmoth=c(0,0,rep(2,numC)), sex=c(2,1,rep(0,numC)), affection=c(0,0,as.integer(caffected)+1), m0.a=c(pa,ca), m0.b=c(pb,cb), env=c(NA,NA,env) ) )
  }
  ## Tests/Exemplifies the strataReduce(...) routine
  srFam <- function( ... ) {
    data <- createFam( ... )
    data2 <- strataReduce( data=data, envCol=9, m0=7 )
    cat( "Original data:\n" )
    print( data )
    cat( "Reduced stratification data:\n" )
    print( data2 )
  }

  ## Basic sib test
  srFam( ca=c(1,1,2), cb=c(1,2,2) )

  ## Basic trio test
  srFam( ca=c(1,1,2), cb=c(1,2,2), pa=c(1,1), pb=c(2,2) )

  ## a fairly comprehensive test here
  ## The affected should always be one of the first three,
  ##  the unaffected could be one the first eight
  for( i in 1:10 )
    srFam( ca=c(1:8,0,0), cb=c(1:8,0,0), pa=c(1,1), caffected=c(rep(TRUE,6),rep(FALSE,4)), env=c(1:3,rep(NA,7)) )

  ## Now just to make sure, a full pedigree, rather than just one family
  data <- createFam( ca=1:2, cb=1:2 )
  for( i in 2:10 )
    data <- rbind( data, createFam( ca=1:2, cb=1:2 ) )
  cat( "Original data (full pedigree):\n" )
  print( data )
  cat( "Reduced stratification data (full pedigree)\n" )
  print( strataReduce( data=data, envCol=9, m0=7 ) )
}
