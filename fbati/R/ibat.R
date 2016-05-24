fixMarkerNames <- function( names )
{
  pop <- function(x) return( x[1:(length(x)-1)] )
  ## pops off the .a/.b
  for( i in 1:length(names) )
    names[i] <- paste( pop(unlist(strsplit(names[i],"\\."))), collapse="." )
  return(names)
}

computeAfreq <- function( idfath, idmoth, m0, m1, allele=2 ) {
  ## they are a parent if they have no parent...
  ##wh <- idfath==0 & idmoth==0 & m0!=0 & m1!=0
  wh <- rep(TRUE, length(m0))
  #print( wh )
  #print( m0 )
  #print( m1 )

  ## 05/28/2008 update
  ##afreq <- ( sum(m0[wh]==allele) + sum(m1[wh]==allele) ) / ( sum(m0[wh]!=0) + sum(m1[wh]!=0) )
  afreq <- ( sum(m0[wh]==allele,na.rm=TRUE) + sum(m1[wh]==allele,na.rm=TRUE) ) / ( sum(m0[wh]!=0,na.rm=TRUE) + sum(m1[wh]!=0,na.rm=TRUE) )
  #cat( "afreq num", ( sum(m0[wh]==allele,na.rm=TRUE) + sum(m1[wh]==allele,na.rm=TRUE) ),
  #     "afreq den", ( sum(m0[wh]!=0,na.rm=TRUE) + sum(m1[wh]!=0,na.rm=TRUE) ),
  #     "\n" )
  return( afreq )
}

## The main exported function
## - needs to be modified to handle 'sym' arguments...
## Right now it only takes in bi-allelic markers...
## - We will need code that translates it to this...
##
## Only uses the first affected (traces to dataComputeGroupG C function),
##  so this is true for the LL method, DR could be different
fbati <- function( ped=NULL, phe=NULL,
                   data=mergePhePed(ped,phe),
                   marker=NULL, ## pairs...
                   env=NULL,
                   method="fbati",
                   model="additive",
                   iter=10000,
                   seed=7,
                   maxSib=3,
                   fixNames=TRUE,
                   debug=FALSE ) {
  ## NEW! if no args, then start the GUI! nifty!
  if( is.null(data) )
    return( fbatiGUI() )

  ## New preprocessing 01/14/2008 safety precaution
  data <- data[ order(data[,1]), ] ## order it so pids are grouped together

  ## NEW preprocessing - always nuclify the data...
  data <- nuclifyMerged( data )

  ## get the environment column
  #envCols <- which(names(data)==env)
  #if( length(envCols) != 1 ) {
  #  stop( "No environment specified." )
  #}
  envCols <- match( env, names(data) )
  if( length(envCols) < 1 )
    stop( "At least one environmental exposure must be selected (each will be tested seperately)." )

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

  return( ibat2( data, marker, markerNames, envCols, env, method, model, iter, seed, maxSib, debug ) )
}

ibat2 <- function( data, marker, markerNames, envCols, envColsNames, method, model, iter, seed, maxSib=2, debug=FALSE ) {
  NUMTESTS <- (length(marker)/2)
  results <- NULL

  for( envCol in envCols ) {
    for( i in 1:NUMTESTS ) {
      m0pos <- marker[2*i-1]
      m1pos <- marker[2*i]

      ## Would be fastest to _subset_ the data at this point
      ##  -- otherwise with large datasets, it keeps copying the whole dataset in!
      s.data <- data[, c(1:6, m0pos,m1pos, envCol)]
      s.m0pos <- 7
      s.m1pos <- 8
      s.envCol <- 9

      ## Reduce the strata
      dataStrataFix <- NULL
      if( maxSib!=0 ) {
        if( !is.null(seed) && !is.na(seed) )
          set.seed(seed)
        dataStrataFix <- strataReduce( data=s.data, envCol=s.envCol, m0=s.m0pos, m1=s.m1pos, maxSib=maxSib )
      }else{
        dataStrataFix <- s.data
      }
      if( debug ) {
        cat( "*** dataStrataFix (recoded + strataReduce(d).) ***\n" )
        print( dataStrataFix )
      }

      ## go accross each of the alleles
      alleles <- setdiff( unique(  c( dataStrataFix[,s.m0pos], dataStrataFix[,s.m1pos] )  ), 0 ) ## elts excluding zero
      #print( alleles )
      if( length(alleles) > 2 )
        stop( paste("Currently only biallelic loci are supported - this is not the case at marker ", markerNames[i] ) )
      for( a in 1:length(alleles) ) {
        #print( alleles )
        ## Recode the dataset
        dataRecode <- dataStrataFix
        for( aa in 1:length(alleles) ) {
          wh <- dataStrataFix[,s.m0pos]==alleles[aa]
          dataRecode[wh,s.m0pos] <- aa

          wh <- dataStrataFix[,s.m1pos]==alleles[aa]
          dataRecode[wh,s.m1pos] <- aa
        }

        ## test allele 2 on the recoded dataset
        pvalue <- NA
        numInf <- NA
        strataSum <- NULL
        res <- dataComputeGroupG( dataRecode, s.m0pos, s.m1pos )
        #print( head( res$groupsG ) )
        #print( res$affectedIndex )
        if( method=="fbati" ) {
          if( debug )
            cat( "***", markerNames[i], "***\n" )

          if( !is.null(dataRecode) ) {
            res_fbati <- fbati.calc( dataRecode, s.m0pos, s.m1pos, res$groupsG, res$affectedIndex, s.envCol, model, iter,  debug )
            pvalue <- res_fbati$pvalue
            numInf <- res_fbati$numInf
            strataSum <- res_fbati$strataSum
          }
        }

        ## how about the allele frequency?
        aafreq <-computeAfreq( dataRecode$idfath, dataRecode$idmoth, dataRecode[,s.m0pos], dataRecode[,s.m1pos] ) ## s., 5/28/2008

        ## also nice would be the number of informative families...

        ## and put the results together
        resultsAddi <- data.frame(environment=names(s.data)[s.envCol],
                                  marker=markerNames[i],
                                  allele=alleles[2],
                                  afreq=aafreq,
                                  numInf=numInf,
                                  model=c("additive","dominant","recessive")[model+1],
                                  pvalue=pvalue )
        if( length(strataSum)>=1 ) {  ## 05.28.2007
          ## If strataSum is empty, then we can't tack this on!
          resultsAddi <- cbind( resultsAddi, strataSum )
        }

        if( is.null(results) ) {
          results <- resultsAddi
        }else{
          results <- mbind( results, resultsAddi ) ##
        }

        ## lastly cycle the alleles for the next iteration
        alleles <- alleles[ c(2:length(alleles),1) ]
      }

    }
  }

  ## Lastly eliminate empty groups
  colKeep <- rep(TRUE,ncol(results))
  for( cc in 8:ncol(results) ) {
    if( sum( results[,cc], na.rm=TRUE ) == 0 )
      colKeep[cc] <- FALSE
  }
  results <- results[,colKeep]


  ## And more lastly, rename the stratas!
  if( ncol(results) > 7 ) {
    for( cc in 8:ncol(results) )
      names(results)[cc] <- strataNameDehash( names(results)[cc] )

    ## And, yes, really the last one now, reorder the strata informative
    resNames <- names(results)
    resNames <- resNames[8:length(resNames)]
    results <- results[, c(1:7, 7+order(nchar(resNames))) ]
  }

  ## And most lastly, sort the results
  results <- results[order(results$marker),]

  return(results)
}

## rbind-ing two data.frame objects, but inserting
##  NA's if the columns have different names
mbind <- function( group1, group2 ) {
  if( ncol(group1)==ncol(group2)
      && sum(sort(names(group1))==sort(names(group2))) == length(group1) )
    return( rbind(group1,group2) );

  names2 <- names(group2);
  for( name in names(group1) )
    if( sum(name==names2)==0 )
      group2[name] = NA;

  mbind( group2, group1 );
}


## DEBUG
ibat.R.debug <- function()
{
  ##library( pbatR )
  dyn.load( "src/ibat.so" ) ## unix way
  source( "R/fbati.R" )
  source( "R/datamatrix.R" )
  source( "R/merge.R" )
  source( "R/defines.R" )
  ped <- read.ped( "test", sym=FALSE ) ## coded 1/2
  ped[ped$id==3,]$AffectionStatus <- 2
  set.seed(0)
  phe <- as.phe( data.frame( pid=ped$pid, id=ped$id, env=rnorm(nrow(ped)) ) )

  print( fbati( ped=ped, phe=phe, env="env" ) )
}
#ibat.R.debug()