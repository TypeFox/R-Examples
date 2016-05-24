## Exported routine for fbatc
fbatc <- function( ped=NULL, phe=NULL, data=mergePhePed( ped, phe),
                   trait="AffectionStatus", traitType="auto",
                   markerAnalyze=NULL, markerCondition=NULL,
                   offset=NULL,
                   tempPrefix="temp",
                   MAXITER=1000, TOL=sqrt(.Machine$double.eps),
                   verbose=FALSE )
{
  ## Do we need to run the GUI version
  if( is.null( ped ) && is.null(data) )
    return(  fbatcGUI()  )

  ## ignoreBtX is a leftover -- passed from several functions, and  so now a necessary placeholder, but it isn't actually used...
  ignoreBtX <- TRUE

  if( is.null(markerCondition) || length(markerCondition)==0 ) {
    ## Possibly need to call this then on every markerAnalyze
    if( is.null(markerAnalyze) || length(markerAnalyze)<2 )
      stop( "Either markerCondition must be specified, or length(markerAnalyze) must be at least 2, in which case each marker will be conditioned on in turn." )

    ## Recurse, calling on every markerAnalyze
    res <- NULL
    for( i in 1:length(markerAnalyze) ) {
      markerAnalyze2 <- markerAnalyze[-i]
      markerCondition2 <- markerAnalyze[i]
      #print( "recurse" )
      #print( markerAnalyze2 )
      #print( markerCondition2 )
      resPlus <- fbatc( ped=ped, phe=phe, data=data,
                        trait=trait, traitType=traitType,
                        markerAnalyze=markerAnalyze2, markerCondition=markerCondition2,
                        offset=offset,
                        tempPrefix=tempPrefix,
                        MAXITER=MAXITER, TOL=TOL,
                        verbose=verbose )
      if( is.null(res) ) {
        res <- resPlus
      }else{
        res <- rbind( res, resPlus )
      }
    }

    ## And return the bound results
    return( res )
  }

  ## Prerequisites
  fbat.install()
  sh.install()

  ## Overwrite this information
  ##assign( "FBAT", fbat.exename(), inherits=TRUE )
  ##cat( "FBAT", FBAT, "\n" )
  ##stop()
  FBAT <- fbat.exename()

  if( traitType=="auto" )
    traitType <- resolveTraitType( phe=phe, trait=trait )

  ## 04.28.2009 NEW, safety precautions for linking the vile trait
  ## - first merge together the pedigree and phenotype object (taken largely from 'nuclify' code)
  data <- mergePhePed( ped, phe )
  ncolped <- ncol(ped)
  ped <- data[, 1:ncolped]
  phe <- data[, c(1,2,(ncolped+1):ncol(data))]
  ## - Find the analyzing and conditioning markers
  markers <- sort(match( c(paste(c(markerAnalyze,markerCondition),".a",sep=""),
                           paste(c(markerAnalyze,markerCondition),".b",sep="")),
                         names(ped)))
  ####print( markers )
  ## -- pull out only the useful markers into the pedigree object
  ped <- ped[ , c(1:6, markers) ]
  keep <- rep(TRUE,nrow(ped))
  for( i in 7:ncol(ped) )
    keep <- keep & ped[,i]!=0
  ped <- ped[ keep, ]
  phe <- phe[ keep, ]
  ## Now sort by 1) pid, 2) is.na(trait) / trait==2 or 1
  if( trait=="AffectionStatus" ) {
    ord <- order( ped$pid, ped$AffectionStatus, na.last=TRUE, decreasing=c(FALSE,TRUE) )
  }else{
    ord <- order( ped$pid, phe[[trait]], na.last=TRUE, decreasing=c(FALSE,TRUE) )
  }
  ped <- ped[ ord, ]
  phe <- phe[ ord, ]
  data <- mergePhePed( ped=ped, phe=phe )
  ## WEN 04.28.2009
  
  ## Now call the routines
  if( traitType=="binary" ) {
    res <- condGeneP4( ped=ped, phe=phe, data=data,
                      trait=trait, traitType=traitType,
                      markerAnalyze=markerAnalyze, markerCondition=markerCondition,
                      alpha=offset,
                      ignoreBtX=ignoreBtX,
                      tempPrefix=tempPrefix,
                      #MAXITER=MAXITER, TOL=TOL,
                      verbose=verbose,
                      FBAT=FBAT )  ## The new addition that fixes this..

    res <- cbind( markerAnalyze=paste(markerAnalyze,collapse=","),
                  markerCondition=paste(markerCondition,collapse=","),
                  res, stringsAsFactors=FALSE )
    return( res ) ## 01/03/09
  }

  ## Otherwise traitType=="continuous"

  res <- condGeneR( ped=ped, phe=phe, data=data,
                    trait=trait, traitType=traitType,
                    markerAnalyze=markerAnalyze, markerCondition=markerCondition,
                    alpha=offset,
                    tempPrefix=tempPrefix,
                    verbose=verbose,
                    FBAT=FBAT )

  ## New -- add on the new condGene Routine
  compVarExpl <- FALSE ## NEVER want to do this, right?
  #res2 <- condGeneP2( ped=ped, phe=phe, data=data,
  #                    trait=trait, traitType=traitType,
  #                    markerAnalyze=markerAnalyze, markerCondition=markerCondition,
  #                    alpha=offset,
  #                    tempPrefix=tempPrefix,
  #                    verbose=verbose,
  #                    FBAT=FBAT )
  res2 <- condGeneP3( ped=ped, phe=phe, data=data,
                      trait=trait, traitType=traitType,
                      markerAnalyze=markerAnalyze, markerCondition=markerCondition,
                      alpha=offset,
                      tempPrefix=tempPrefix,
                      verbose=verbose,
                      FBAT=FBAT,
                      compVarExpl=compVarExpl )
  res <- cbind( res, pvalue=res2$pvalue, rank=res2$rank, numInf=res2$numInf, stringsAsFactors=FALSE )
  if( compVarExpl )
    res <- cbind( res, varExpl=res2$varExpl, stringsAsFactors=FALSE )

  return( res )
}

fbatcFunc <- function(  ped, phe, trait, traitType="auto", markerAnalyze=markerAnalyze, markerCondition=markerCondition, offset=NULL, tempPrefix="temp", MAXITER=1000, TOL=sqrt(.Machine$double.eps), write_results ) {
  #print( "offset" )
  #print( offset )
  #print( "offset" )

  if( is.null(ped) || is.na(ped) || (is.character(ped) & nchar(ped)==0) )
    return( "A 'pedigree' file must be specified." )
  #if( is.null(phe) || is.na(phe) || (is.character(phe) & nchar(phe)==0) )
  #  return( "A 'phenotype' file must be specified." )

  if( is.na(trait) || is.null(trait) )
    return( "A 'trait' must be specified." )

  if( is.null(markerAnalyze) || is.na(markerAnalyze[1]) ) {
    if( is.null(markerCondition) || is.na(markerCondition) ) {
      markerAnalyze <- as.vector( guiGetSafe("fbati_internal_possibleMarkers") )
      if( is.null(markerAnalyze) || is.na(markerAnalyze) )
        return( "The pedigree file must have at least two markers!" )
      markerCondition <- NULL
    }else{
      return( "A marker to be analyzed must be specified." )
    }
  }
  if( is.null(markerCondition) || is.na(markerCondition[1]) ){
    if( length(markerAnalyze) < 2 )
      return( "A marker to be conditioned on must be specified, or at least two markers to analyze, in which case each marker will be conditioned on in turn." )
    markerCondition <- NULL
  }
  if( length(intersect(markerAnalyze,markerCondition)) != 0 )
    return( "There can be no overlap in the set of markers to be analyzed and the set of markers to be conditioned on." )

  if( !is.null(offset) && !is.na(offset) && is.character(offset) && nchar(offset)==0 )
    offset <- NULL
  if( !is.null(offset) && is.na(offset) )
    offset <- NULL

  fped <- fread.ped(ped,lowercase=FALSE)

  fphe <- NULL
  if( is.character(phe) & nchar(phe)!=0 )
    try( fphe <- fread.phe( phe, lowercase=FALSE ), silent=TRUE )

  #print( "markerAnalyze, markerCondition" )
  #print( markerAnalyze )
  #print( str( markerAnalyze ) )
  #print( markerCondition )

  fbatcRes <- fbatc( ped=fped, phe=fphe, trait=trait, traitType=traitType, markerAnalyze=markerAnalyze, markerCondition=markerCondition, offset=offset, MAXITER=MAXITER, TOL=TOL, tempPrefix=tempPrefix )

  guiSet( "fbatcRes", fbatcRes )
  #print( fbatcRes )

  return( "Processed." )
}

updateFbatcGUI <- function( arg ) {
  if( arg=="ped" ) {
    pedfile <- guiGetValue(1)
    ped <- read.ped(pedfile)

    possibleMarkers <- names(ped)[7:length(names(ped))]
    setListElements( "markerAnalyze", sort(possibleMarkers) )
    setListElements( "markerCondition", sort(possibleMarkers) )
    guiSet( "fbati_internal_possibleMarkers", sort(as.vector(possibleMarkers)) )

    file.strip.extension <- getFromNamespace( "file.strip.extension", "pbatR" )
    phename <- paste( file.strip.extension(pedfile), ".phe", sep="" )
    guiSetValue( 2, phename )
    arg <- "phe"
  }
  if( arg=="phe" ) {
    phe <- NULL
    try( phe <- read.phe(guiGetValue(2)), silent=TRUE )  ## potentially try to auto-load
    if( !is.null(phe) ) {
      possibleEnv <- names(phe)[-c(1,2)]
      setListElements( "trait", c("AffectionStatus",possibleEnv) )
    }else{
      setListElements( "trait", c("AffectionStatus") )
    }
  }
}

writeFbatcGUI <- function() {
  res <- guiGetSafe( "fbatcRes" )

  if( !is.null(res) && !is.na(res[1]) ) {
    defaultFile <- "results"
    outStr <- tclvalue(tkgetSaveFile(title="Write FBAT-C Results", filetypes="{{CSV (spreadsheet)} {.csv}}", initialfile=defaultFile))
    if( nchar(outStr) > 0 ) {
      outStr <- getFromNamespace( "str.file.extension", "pbatR" )( outStr, "csv" )
      write.csv( res, outStr, row.names=FALSE )
      cat( "Results written to disk.\n" )
    }else{
      print( outStr )
      tkmessageBox( message="Could not write file to disk.", title="Write Failure" )
    }
  }else{
    tkmessageBox( message="There are no results to write to disk.", title="No Results" )
  }
}

fbatcGUI <- function() {
  gui( fbatcFunc,
       argFilename=list(ped=NULL,phe=NULL),
       argFilter=list(ped="{{Pedigree file} {.ped}}", phe="{{Phenotype file} {.phe}}"),
       argOption=list(traitType=c("auto","binary","continuous")),
       argList=list(trait=c("AffectionStatus"), markerAnalyze=NULL, markerCondition=NULL),
       argCommand=list(write_results=writeFbatcGUI),
       callback=updateFbatcGUI,
       helpsFunc="fbatc",
       title="FBAT-C GUI",
       argText=list(ped="Pedigree File...", phe="Phenotype File...", markerAnalyze="Marker Analyze", markerCondition="Marker Condition", ignoreBtX="Ignore BX", tempPrefix="Temp Prefix") )
  return( guiGetSafe( "fbatcRes" ) )
}


## Might be better to instead write a separate linkTrait function that
##  1) verify's the trait is linked based on the trait _and_ the genotype at that trait
##  2) just verify that the trait is the same across all genotypes --> I.e. all references end up linking to the same trait...
## Safety precaution in case condGeneFBATControl_linkTrait is misbehaving
## The idea in that routine is that those individuals who had missing genotypes should not be linked
#safety_linkTrait <- function( ped, phe ) {
#  ## Merge phe and ped files
#  ## Sort s.t. those with traits are listed _first_ in the pedigree (for each pedigree)
#  
#  data <- mergePhePed( ped, phe )
#  
#  
#}