## Might want to put a message in the fbati function suggesting that the user use this methodology instead...

#source( "fbatge.R" )
#library( fgui )

## Exported function
fbatge <- function( ped=NULL, phe=NULL, env=NULL, cov=NULL, trait="AffectionStatus", geno=NULL, strategy="hybrid", model="additive" ) {
  if( is.null(ped) && is.null(phe) && is.null(env) )
    fbatgeGUI()

  ## Me and my stupid naming conventions
  strategy <- tolower(strategy)
  model <- tolower(model)
  if( strategy=="hybrid" ) strategy <- "adaptive"
  if( strategy=="rr" ) strategy <- "geno"
  if( strategy=="clr" ) strategy <- "pheno"

  ## nuclify and merge it
  data <- nuclify( ped=ped, phe=phe )
  ped <- data$ped
  phe <- data$phe
  ##print( head( ped ) )  ## DEBUG ONLY
  ##print( head( phe ) )  ## DEBUG ONLY

  ## convert the ped to a gped
  gped <- as.gped( ped )
  ##print( head( gped ) )

  ##print( head(gped) )
  ##stop()

  ## if geno is NULL, then set it to be all of the markers
  if( is.null(geno) )
    geno <- names(gped)[7:ncol(gped)]

  ## Safety checks are handled in fbatge.internal_C

  ## Finally go ahead and run the analysis
  ret <- NULL
  for( g in geno ) {
    ## Try wrapper added 06.19.2009
    try( {
      res <- fbatge.internal_C( gped=gped, phe=phe, geno=g, cov=cov, trait=trait, env=env, strategy=strategy, model=model )
      retaddi <- c( strategy=strategy, model=model, marker=g, trait=trait, env=env, cov=paste(cov,collapse=","), numInf=attr(res,"numInf"), pvalue=res )
      if( is.null(ret) ) {
        ret <- as.data.frame( t(retaddi), stringsAsFactors=FALSE )
      }else{
        ret <- rbind( ret, retaddi )
      }
    } )
  }

  ## Umm... do we really want to do this?
  ret$strategy[ret$strategy=="adaptive"] <- "Hybrid"
  ret$strategy[ret$strategy=="geno"] <- "RR"
  ret$strategy[ret$strategy=="pheno"] <- "CLR"

  return( ret )
}

## structure modeled after fbatc stuff (cgFbatRun.R)
fbatgeFunc <- function( ped, phe, env, cov=NULL, trait="AffectionStatus", geno=NULL, strategy="Hybrid", model="additive", print_results, write_results ) {
  ## Make sure a pedigree file was specified
  if( is.null(ped) || is.na(ped) || (is.character(ped) && nchar(ped)==0) )
    return( "A 'pedigree' file must be specified." )

  ## This routine cannot be without a phenotype file
  if( is.null(phe) || is.na(phe) || (is.character(phe) && nchar(phe)==0) )
    return( "A 'phenotype' file must be specified." )

  ## Make sure an environment was specified
  if( is.na(env) || is.null(env) )
    return( "An environmental exposure must be specified." )

  ## Make sure a trait was specified
  if( is.na(trait) || is.null(trait) )
    return( "A 'trait' must be specified." )

  ## Fix up the covariate if empty
  if( is.na(cov) || is.null(cov) )
    cov <- NULL

  ## Fix up geno if empty
  if( is.na(geno) || is.null(geno) )
    geno <- NULL

  ## DEBUG
  #cat( "ped",ped,"\n", "phe",phe,"\n", "env",env,"\n", "cov",cov,"\n", "trait",trait,"\n", "geno",geno,"\n", "strategy",strategy,"\n", "model",model,"\n" )

  fped <- fread.ped( ped )
  fphe <- fread.phe( phe )

  fbatgeRes <- fbatge( ped=fped, phe=fphe, env=env, cov=cov, trait=trait, geno=geno, strategy=strategy, model=model )

  guiSet( "fbatgeRes", fbatgeRes )

  return( "Processed." )
}

updateFbatgeGUI <- function( arg ) {
  if( arg=="ped" ) {
    pedfile <- guiGetValue(1)
    ped <- read.ped( pedfile )

    possibleMarkers <- ped.markerNames( ped )
    setListElements( "geno", sort(possibleMarkers) )

    file.strip.extension <- getFromNamespace( "file.strip.extension", "pbatR" )
    phename <- paste( file.strip.extension(pedfile), ".phe", sep="" )
    guiSetValue( 2, phename )
    arg <- "phe"
  }
  if( arg=="phe" ) {
    phe <- NULL
    try( phe <- read.phe( guiGetValue(2) ), silent=TRUE ) ## potentially auto-loaded
    if( !is.null(phe) ) {
      nm <- names(phe)[-c(1,2)]
      setListElements( "trait", c("AffectionStatus", nm) )
      setListElements( "env", nm )
      setListElements( "cov", nm )
    }
  }
}

writeFbatgeGUI <- function() {
  res <- guiGetSafe( "fbatgeRes" )

  if( !is.null(res) && !is.na(res[1]) ) {
    defaultFile <- "results_fbatge"
    outStr <- tclvalue(tkgetSaveFile(title="Write FBAT Gene Environment Interaction Results", filetypes="{{CSV (spreadsheet)} {.csv}}", initialfile=defaultFile))
    if( nchar(outStr) > 0 ) {
      outStr <- getFromNamespace( "str.file.extension", "pbatR" )( outStr, "csv" )
      write.csv( res, outStr, row.names=FALSE )
      cat( "Results written to disk.\n" ) ## semi-debug
    }else{
      print( outStr )
      tkmessageBox( message="Could not write file to disk.", title="Write Failure" )
    }
  }else{
    tkmessageBox( message="There are no results to write to disk.", title="No Results" )
  }
}

printFbatgeGUI <- function() {
  res <- guiGetSafe( "fbatgeRes" )

  if( !is.null(res) && !is.na(res[1]) ) {
    print( res )
  }else{
    cat( "There are no results...\n" )
  }
}

fbatgeGUI <- function() {
  gui( fbatgeFunc,
      argFilename=list(ped=NULL,phe=NULL),
      argFilter=list(ped="{{Pedigree file} {.ped}}", phe="{{Phenotype file} {.phe}}"),
      argOption=list(strategy=c("Hybrid", "RR", "CLR"),
                     model=c("Additive","Codominant")),
      argList=list(env=NULL,
                   cov=NULL,
                   trait=c("AffectionStatus"),
                   geno=NULL),
      argCommand=list(print_results=printFbatgeGUI, write_results=writeFbatgeGUI),
      callback=updateFbatgeGUI,
      helpsFunc="fbatge",
      title="FBAT Gene Environment Interaction Tests",
      argText=list( ped="Pedigree File...",
                    phe="Phenotype File...",
                    env="Environmental Exposure",
                    cov="Covariates",
                    trait="Trait",
                    geno="SNPs",
                    strategy="Strategy",
                    model="Genetic Model",
                    print_results="Print Results",
                    write_results="Write Results" ),
      exec="Process...",
      argGridOrder=c(1,2,3,3,3,3,4,5,6,6)
      )
  return( guiGetSafe( "fbatgeRes" ) )
}

  