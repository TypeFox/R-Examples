## Wrapper to programControl for fbat
fbatControl <- function( commands, FBATEXE="~/bin/fbat_2.02c", ... )
  programControl( FBATEXE, c(commands,"quit"), ... )

############################################
## Main fbat routine                      ##
## model = 'a', 'd', 'r'                  ##
## Updated to handle multivariate traits! ##
############################################
fbatShell <- function( ped, phe=NULL, trait=NULL, model="a", markers=c(), fbatCmd="fbat", tempPrefix="temp_", FBATEXE="~/bin/fbat_2.02c", debug=FALSE, justLines=FALSE ) {
  ## When taking this out of the context of these sims,
  ## just add the following line to the beginning of this function
  ##  seed <- ""

  ## Are we using affection status?
  usingAffectionStatusOnly <- FALSE
  if( is.null(trait) ) trait <- "AffectionStatus"
  if( all(trait=="AffectionStatus") )
    usingAffectionStatusOnly <- TRUE  ## Really _only_ using affection status...

  ## write the datasets to disk
  ## write out data to disk
  pedFilename <- paste( tempPrefix, "_fbat.ped", sep="" )
  pheFilename <- paste( tempPrefix, "_fbat.phe", sep="" )
  write.ped( pedFilename, ped )
  if( !usingAffectionStatusOnly )
    write.phe( pheFilename, phe )

  wh <- which( trait=="AffectionStatus" )
  if( sum(wh) > 0 )
    trait[wh] <- "affection"

  pheAddi <- c()
  if( !usingAffectionStatusOnly )
    pheAddi <- paste("load ",pheFilename,sep="")

  logfile <- paste( tempPrefix, "_log.log", sep="" )
  ## figure out the fbat commands
  cmds <- c( "minsize 0",
             paste( "log", logfile ),  ## 03.09.2009
             paste("load ",pedFilename,sep=""),
             pheAddi,
             paste( "model", model ), ## i.e. a=additive, d=dominant, ...
             paste( "trait", paste(trait,collapse=" ") ),
             paste( c( fbatCmd, markers ), collapse=" " )  )

  #print( cmds )
  #stop()

  #print( "just before fbat control" )
  controlFilename <- paste( tempPrefix, "_control.sh", sep="" )
  ##lines <- fbatControl( cmds, FBATEXE=FBATEXE, filename=controlFilename ) ## 03.09.2009
  programControl( program=FBATEXE, commands=c(cmds,"quit"),  filename=controlFilename, intern=FALSE, output="/dev/null" )
  lines <- readLines( logfile )
  #print( lines )

  ## New, just return the lines when something particularly special
  if( justLines )
    return( lines )

  if( debug )
    print( lines )

  results <- fbatDataframe( lines, length(trait), debug=debug )

  ## CLEAN-UP
  # No longer done at this level - done at the end of the sim..
  #file.remove( "/tmp/th" )  ## Somewhat superfluous
  ## UNIX ONLY
  #system( paste( "rm -fr", thdir ) )

  return( results )
}


## LEFT OFF HERE, BEFORE THE DATAFRAME PIECE!!!

## takes results of "lunetta", e.g.
fbatDataframe <- function( fbatLines, numTraits, debug=FALSE ) {
  ## get the lines that have commands on them
  cmdLines <- which( substr( fbatLines, 1, 2 ) == ">>" )
  nCmdLines <- length(cmdLines)

  ## figure out which lines should be the dataframe...
  dfLines <- (cmdLines[nCmdLines-1]+2 + 2):(cmdLines[nCmdLines]-2)
  if( numTraits>1 ) dfLines <- dfLines + 1
  headerLine <- dfLines[1]-1 - 1;

  if( fbatLines[headerLine] == "" ) {
    ## fucking FBAT
    headerLine <- dfLines[1]-1 - 1 + 1;
    dfLines <- dfLines[-1]
  }

  if( debug ) {
    print( "dfLines" )
    print( dfLines )
    print( "header line" )
    print( fbatLines[headerLine] )
    print( "data lines" )
    print( fbatLines[dfLines] )
  }

  dfNames <- unlist( strsplit( fbatLines[headerLine], split="  *" ) )
  #print( "dfNames" )
  #print( dfNames )
  ## KLUDGE FOR THE NAMES
  if( dfNames[1] == "MarkerAllele" )
    dfNames <- c( "Marker", "Allele", dfNames[2:length(dfNames)] )
  df <- strsplit( fbatLines[dfLines], split="  *" )
  df <- data.frame( t( data.frame( df, stringsAsFactors=FALSE ) ), stringsAsFactors=FALSE )
  row.names(df) <- NULL
  if( ncol(df) == length(dfNames)-1 )
    dfNames <- dfNames[-6] ## Stupid, stupid, wretched, evil, FBAT! AHHHHHH!!!
  names( df ) <- dfNames

  ## Old code -- we get factors here...
  #df <- data.frame( matrix( unlist( strsplit( fbatLines[dfLines], split="  *" ) ), nrow=length(dfLines), byrow=TRUE ) )
  #names( df ) <- dfNames

  return( df )
}

fbatSingleLine <- function( fbatLines, strPrefix ) {
  locs <- grep( paste("^",strPrefix,sep=""), fbatLines )
  return( fbatLines[locs] )
}

fbatShellMM <- function( ... ) {
  fbatLines <- fbatShell( ..., fbatCmd="fbat -m", justLines=TRUE )
  line <- fbatSingleLine( fbatLines, "FBAT_MM" )
  #print( line )
  #stop()
  if( length(line) != 1 ) {
    print( "fbatMM -- error" )
    return( 1.0 )
  }
  toks <- unlist( strsplit(line, " ") )
  return( as.numeric( toks[length(toks)] ) )
}

fbatShellCorrelation <- function( fbatCmd="hapfreq -r", ... ) {
  fbatLines <- fbatShell( ..., fbatCmd=fbatCmd, justLines=TRUE )
  #print( fbatLines )

  upperLoc <- grep( "^Pairwise LD", fbatLines ) + 1
  lowerLoc <- grep( "^>>", fbatLines ) - 1
  lowerLoc <- lowerLoc[length(lowerLoc)]
  header <- fbatLines[upperLoc]
  #print( header )
  matlines <- fbatLines[ (upperLoc+1):(lowerLoc-1) ]
  #print( matlines )

  dfNames <- unlist( strsplit( header, split="  *" ) )
  #print( "dfNames" ); print( dfNames );
  df <- strsplit( matlines, split="  *" )
  #print( "df 1" ); print( df );

  M <- length( df[[length(df)]] ) #- 1
  corr <- matrix( 1, nrow=M, ncol=M )
  for( r in 2:M ) {
    for( c in 1:(r-1) ) {
      ## Currently LD(r^2), but we want the R^2 values!
      toks <- unlist( strsplit( df[[r-1]][c+1], "\\(|\\)" ) )
      corr[r,c] <- as.numeric(toks[2])
      #corr[r,c]  <- df[[r-1]][c+1]
    }
  }
  for( r in 1:(M-1) )
    for( c in (r+1):M )
      corr[r,c] <- corr[c,r]

  corrNames <- c( dfNames[2:length(dfNames)], df[[M-1]][1] )
  #print( corrNames )
  rownames(corr) <- corrNames
  colnames(corr) <- corrNames
  #print( corr )
  #stop()

  return( corr )
}

## DEBUG
fbatShell.debug <- function() {
  source( "../cgFbat/cgFbatCode.R" ) ## for programControl
  ## library( pbatR )
  ped <- fread.ped( "../results/analyze/fbatDatasets/data/CAMP" )
  phe <- fread.phe( "../results/analyze/fbatDatasets/data/CAMPZ" )
  print( "Pedigree names" )
  print( names( ped ) )
  #print( ped$AffectionStatus )
  print( "Phenotype names" )
  print( names( phe ) )

  DEBUG <- TRUE
  #print( fbatShell( ped=ped, phe=phe, debug=DEBUG ) )
  #print( fbatShell( ped=ped, phe=phe, trait="zposfevp", debug=DEBUG ) )

  #print( fbatShellMM( ped=ped, phe=phe, debug=TRUE ) )

  #print( fbatShellCorrelation( ped=ped ) )

  print( fbatShell( ped=ped, phe=phe, markers=c("m47","m654"), debug=DEBUG ) )
}
#fbatShell.debug()