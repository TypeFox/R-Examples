## DEBUG, must be left out of package... sigh...
#library( fgui )
#library( ibat )

fbatiFunc <- function( ped, phe, env, marker, model="additive", iter=10000, seed=7, maxSib=3, write_results ) {
  ## 05/28/2008
  if( is.null( ped ) || is.na( ped ) || (is.character(ped) & nchar(ped)==0) )
    return( "A 'pedigree' must be specified." )

  if( is.null( phe ) || is.na( phe ) || (is.character(phe) & nchar(phe)==0) )
    return( "The 'phenotype' must be specified." )

  fped <- fread.ped(ped,lowercase=FALSE)  ## lowercase added 05/26/2008
  fphe <- fread.phe(phe,lowercase=FALSE)

  if( is.null( env ) || is.na( env ) )
    return( "Environment must be specified." )

  if( !is.null(marker) && is.na(marker[1]) ) marker <- NULL

  fbatiRes <- fbati( ped=fped, phe=fphe, env=env, marker=marker, model=model, seed=seed, iter=iter, maxSib=2 )

  #print( fbatiRes ) ## debugging
  guiSet( "fbatiRes", fbatiRes ) ## Set the results for the 'write' function

  return( "Processed." )
}

updateFbatiGUI <- function(arg) {
  #print( paste( "update", arg ) )

  ## Do I need to tryCatch(...) any errors here???

  if( arg=="ped" ) {
    pedfile <- guiGetValue(1)
    ped <- read.ped(pedfile) ## symbolic for now
    ## Set the list of possible markers
    possibleMarkers <- names(ped)[7:length(names(ped))]
    setListElements( "marker", sort(possibleMarkers) )

    ## Set the phenotype file, and then arg to be that so it falls through...
    file.strip.extension <- getFromNamespace( "file.strip.extension", "pbatR" )
    phename <- paste( file.strip.extension(pedfile), ".phe", sep="" )
    guiSetValue( 2, phename )
    arg <- "phe"
  }
  if( arg=="phe" ) { ## CANNOT BE "ELSE IF"!
    phe <- NULL
    try( phe <- read.phe(guiGetValue(2)) ) ## symbolic for now
    if( !is.null(phe) ) {
      ## Set the list of possible environmental variables
      possibleEnv <- names(phe)[-c(1,2)]
      setListElements( "env", sort(possibleEnv) )
    }
  }
}

writeFbatiGUI <- function() {
  res <- guiGetSafe( "fbatiRes" )

  if( !is.null(res) ) {
    #pedfile <- guiGetValue(1)
    #file.strip.extension <- getFromNamespace( "file.strip.extension", "pbatR" )
    #defaultFile <- paste( file.strip.extension(pedfile), ".csv", sep="" )
    defaultFile <- "results"

    ## Prompt user for a filename
    outStr <- tclvalue(tkgetSaveFile(title="Write FBAT-I Results",filetypes="{{CSV (spreadsheet)} {.csv}}", initialfile=defaultFile))
    if( nchar(outStr)>0 ) {
      outStr <- getFromNamespace( "str.file.extension", "pbatR" )( outStr, "csv" ) ## 05/28/2008
      write.csv( res, outStr, row.names=FALSE )
      cat( "Results written to disk.\n" )
    }else{
      print( outStr )
      tkmessageBox( message="Could not write file to disk.", title="Write failure" )
    }
  }
}

## *** EXPORT ***
fbatiGUI <- function() {
  gui( fbatiFunc,
       argFilename=list(ped=NULL,phe=NULL),
       argFilter=list(ped="{{Pedigree file} {.ped}}", phe="{{Phenotype file} {.phe}}"),
       argOption=list(model=c("additive","dominant","recessive"),
                      strataFix=c("TRUE","FALSE")),
       argList=list(env=NULL,marker=NULL),
       argCommand=list(write_results=writeFbatiGUI),
       callback=updateFbatiGUI,
       helpsFunc="fbati",
       title="FBAT-I GUI",
       argText=list(ped="Pedigree File ...", phe="Phenotype file ...", env="Environmental Exposure", marker="Choose genetic marker (or all will be done)", write_results="Write data results in csv format (press after running)", iter="Monte-Carlo iterations", seed="Random Seed", maxSib="Maximum sibs"),
       verbose=FALSE )
  return( guiGetSafe( "fbatiRes" ) )
}

#fbatiGUI() ## DEBUG