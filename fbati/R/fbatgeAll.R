## Provides a wrapper that does all of the GxE tests, the joint tests, and the ME tests

# uni <- function( ped, phe, alpha=1, trait="trait", debug=TRUE, tempPrefix="temp_", FBATEXE="~/bin/fbat"  ) {
#   univariateFilter <- getFromNamespace( "univariateFilter", "fbati" )
#   #fbatShell <- getFromNamespace( "fbatShell", "fbati" )
#
#   numMarkers <- (ncol(ped) - 6) / 2
#   markers <- ped.markerNames(ped)
#   res <- univariateFilter( ped=ped, phe=phe, markers=markers, trait=trait, alpha=alpha/numMarkers, tempPrefix=tempPrefix, FBATEXE=FBATEXE )
#   #print( res )
#   #stop()
#   return( res$filteredMarkers )
# }

fbatgeAll <- function( ped=NULL, phe=NULL, env=NULL, trait="AffectionStatus" ) {
  #require( "fbati" ) ## TAKE OUT IF PACKAGED

  if( is.null(ped) )
    return( fbatgeAllGUI() )

  model <- "additive"

  ## The main effect
  # #   me <- NULL
  # #   try( {
  # #     ## - Get helper functions from fbati (control FBAT)
  # #     univariateFilter <- getFromNamespace( "univariateFilter", "fbati" )
  # #     fbat.install <- getFromNamespace( "fbat.install", "fbati" )
  # #     sh.install <- getFromNamespace( "sh.install", "fbati" )
  # #     fbat.exename <- getFromNamespace( "fbat.exename", "fbati" )
  # #     fbat.install()
  # #     sh.install()
  # #     FBAT <- fbat.exename()
  # #     ## - and then get the results
  # #     me <- univariateFilter( ped=ped, phe=phe, markers=NULL, trait=trait, alpha=1, tempPrefix="temp_", FBATEXE=FBAT )
  # #     me <- me$univariateResults
  # #   } )
  # #   if( is.null(me) )
  # #     cat( "Failed to properly install FBAT2.02c; main effects results will not be shown.\n" )

  me <- fbatme( ped=ped, phe=phe, trait=trait, model=model )
  ##print( head( me ) )

  ## The joint test
  joint <- fbatj( ped=ped, phe=phe, marker=NULL, trait=trait, env=env, model=model )

  ## The GxE tests
  #model <- "codominant" ## PERSONAL HACK, MUST BE TAKEN OUT!!!
  rr <- fbatge( ped=ped, phe=phe, env=env, trait=trait, geno=NULL, strategy="rr", model=model )
  clr <- fbatge( ped=ped, phe=phe, env=env, trait=trait, geno=NULL, strategy="clr", model=model )
  hybrid <- fbatge( ped=ped, phe=phe, env=env, trait=trait, geno=NULL, strategy="hybrid", model=model )

  ## Some cleanup...
  cat( "\n" )

  ## Now paste all the results together
  #print( head(me) )
  #print( head(joint) )
  #print( head(rr) )
  #print( head(clr) )
  #print( head(hybrid) )

  # # names( me ) <- c("Marker","Afreq","G#","G")
  names( me ) <- c("Marker","Afreq","G#","G")
  joint <- joint[,c("marker","numInf","pvalue")];  names( joint ) <- c("Marker","G-GE#","G-GE")
  rr <- rr[,c("marker","numInf","pvalue")];  names( rr ) <- c("Marker","GE.RR#","GE.RR")
  clr <- clr[,c("marker","numInf","pvalue")];  names( clr ) <- c("Marker","GE.CLR#","GE.CLR")
  hybrid <- hybrid[,c("marker","numInf","pvalue")];  names( hybrid ) <- c("Marker","GE.Hybrid#","GE.Hybrid")

  res <- me
  res <- merge( res, joint, all=TRUE )
  res <- merge( res, rr, all=TRUE )
  res <- merge( res, clr, all=TRUE )
  res <- merge( res, hybrid, all=TRUE )
  res <- res[ order(res$"G-GE"), ]
  return( res )
}




#############################
## For Graphical Interface ##
#############################

## structure modeled after fbatc stuff (cgFbatRun.R)
fbatgeAllFunc <- function( ped, phe, env, trait="AffectionStatus", print_results, write_results ) {
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

  ## DEBUG
  #cat( "ped",ped,"\n", "phe",phe,"\n", "env",env,"\n", "cov",cov,"\n", "trait",trait,"\n", "geno",geno,"\n", "strategy",strategy,"\n", "model",model,"\n" )

  fped <- fread.ped( ped )
  fphe <- fread.phe( phe )

  fbatgeAllRes <- fbatgeAll( ped=fped, phe=fphe, env=env, trait=trait )

  guiSet( "fbatgeAllRes", fbatgeAllRes )

  return( "Processed." )
}

updateFbatgeAllGUI <- function( arg ) {
  if( arg=="ped" ) {
    pedfile <- guiGetValue(1)
    ## ped <- read.ped( pedfile )  ## 02/20/2014 codetools

    #possibleMarkers <- ped.markerNames( ped )
    #setListElements( "geno", sort(possibleMarkers) )

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
    }
  }
}

writeFbatgeAllGUI <- function() {
  res <- guiGetSafe( "fbatgeAllRes" )

  if( !is.null(res) && !is.na(res[1]) ) {
    defaultFile <- "results_fbatgeAll"
    outStr <- tclvalue(tkgetSaveFile(title="Write Extended FBAT Gene Environment Interaction Results", filetypes="{{CSV (spreadsheet)} {.csv}}", initialfile=defaultFile))
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

printFbatgeAllGUI <- function() {
  res <- guiGetSafe( "fbatgeAllRes" )

  if( !is.null(res) && !is.na(res[1]) ) {
    print( res )
  }else{
    cat( "There are no results...\n" )
  }
}

fbatgeAllGUI <- function() {
  gui( fbatgeAllFunc,
      argFilename=list(ped=NULL,phe=NULL),
      argFilter=list(ped="{{Pedigree file} {.ped}}", phe="{{Phenotype file} {.phe}}"),
      argOption=list(strategy=c("Hybrid", "RR", "CLR"),
                     model=c("Additive","Codominant")),
      argList=list(env=NULL,
                   trait=c("AffectionStatus")),
      argCommand=list(print_results=printFbatgeAllGUI, write_results=writeFbatgeAllGUI),
      callback=updateFbatgeAllGUI,
      helpsFunc="fbatge",
      title="Additive GxE (fbatge more itx options)",
      argText=list( ped="Pedigree File...",
                    phe="Phenotype File...",
                    env="Environmental Exposure",
                    trait="Trait",
                    print_results="Print Results (After 'Process...')",
                    write_results="Write Results (After 'Process...')" ),
      exec="Process...",
      argGridOrder=c(1,2,3,3,4,4)
      )
  return( guiGetSafe( "fbatgeAllRes" ) )
}
