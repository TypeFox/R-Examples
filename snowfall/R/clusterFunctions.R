##*****************************************************************************
## Functions which extend Snow usage or implement some higher level usage.
## Wrappers for Snow functions are in snowWrappers.R
##
## Functions:
##    sfLoadLib   - Load library in cluster (path conversion, flags...)
##    sfSource    - Load source (path conversion...)
##
##    sfExport    - Export local and global objects to cluster
##    sfExportAll - Export all global objects (with given exception list)
##    sfRemove    - Remove objects from nodes
##    sfRemoveAll - Remove all objects from nodes (with given excpetion list)
##
##    sfCat       - Cat something on cluster
##    sfPrint     - Print something on cluster
##*****************************************************************************

##*****************************************************************************
## Load a library depending on sequential or parallel execution.
##
## Should behave most likely the build-in library() function.
##
## Running in sequential mode: normal "library()" command.
## Running in parallel mode: the library is loaded on ALL nodes.
##
## PARAMETERS: 
## RETURN:     Logical TRUE/FALSE on success. On failure, if noStopOnError is not
##             set, stop immidiately.
##*****************************************************************************
sfLibrary <- function( package,
                       pos = 2,
                       lib.loc = NULL,
                       character.only = FALSE,
                       warn.conflicts = TRUE,
# keep.source is removed in R3, but kept here for back-compatible API
#                       keep.source = getOption("keep.source.pkgs"),
                       keep.source = NULL,
                       verbose = getOption("verbose"),
                       version,
                       stopOnError = TRUE ) {
  sfCheck();
  
  ## Generate (global) names list with all parameters.
#  setVar( ".sfPars", list() )
  sfPars <- list()

  ## Help does not make sense.
##  if( !missing( help ) )
##    stop( "Help is not allowed in sfLibrary. Use 'library' instead." )
  
  if( !missing( package ) ) {
    if( character.only ) {
      if( is.character( package ) )
        sfPars$package <- package
      else
        stop( paste( "Package", package, "is no character string." ) )
    }
    else
      sfPars$package <- deparse( substitute( package ) )
  }

  ## package is now a string in any case.
  sfPars$character.only <- TRUE

  sfPars$pos            <- pos
  sfPars$lib.loc        <- lib.loc
  sfPars$warn.conflicts <- warn.conflicts

  # Raw quickfix for R 3: simply remove argument if we are on version 3.
  # But keep it to not break R 2.x (although default argument changed
  # there).
  if( as.integer(R.version$major) < 3 ) {
    # As new default for keep source is "NULL" (from 1.84-2), we need
    # to rebuild the old default behavior of the functions arguments.
    if( is.null( keep.source ) )
      keep.source = getOption("keep.source.pkgs")

    sfPars$keep.source    <- keep.source
  }

  sfPars$verbose        <- verbose

  ## All libraries are loaded internally with logical.return.
  sfPars$logical.return <- TRUE
  
  if( !missing( version ) )
    sfPars$version <- version

  if( sfParallel() ) {
    ## On Nodes load location with absolute path.
    if( !is.null( sfPars$lib.loc ) )
      sfPars$lib.loc <- absFilePath( sfPars$lib.loc )

    ## Export to namespace.
    setVar( ".sfPars", sfPars )

    ## Weird enough ".sfPars" need to be exorted (else it would not be found
    ## on slave, although it is a parameter)
    sfExport( ".sfPars", local=FALSE, namespace="snowfall" )

    ## Load libs using require as Exception on nodes doesn't help us here.
    ## @todo Check on correct execution via logical.return
    ## @todo Exporting of .sfPars needed?
    ## CHANGE FOR R-3: attribute 'keep-source' is removed.
    result <- try( sfClusterEval( do.call( "library",
                                          .sfPars ) ) )
#    result <- try( sfClusterEval( do.call( "library", .sfPars ) ) )

    if( inherits( result, "try-error" ) ||
        ( length( result ) != sfCpus() ) ||
        !all( checkTryErrorAny( result ) ) ||
        !all( unlist( result ) ) ) {
      if( stopOnError )
        stop( paste( "Stop: error loading library on slave(s):",
                     .sfPars$package ) )
      else {
        warning( paste( "Error loading library on slave(s):", package ) )
        return( invisible( FALSE ) )
      }
    }
    else {
      ## Load message in slave logs.
      sfCat( paste( "Library", .sfPars$package, "loaded.\n" ) )

      ## Message in masterlog.
      message( paste( "Library", .sfPars$package, "loaded in cluster.\n" ) )
    }
  }

  result <- try( do.call( "library", sfPars ) )

  ## Remove global var from cluster (and local).
  ## Do before exception checks, as there might by a stop.
  sfRemove( ".sfPars" )

  if( inherits( result, "try-error" ) || !result ) {
    if( stopOnError ) {
      warning( paste( "Unable to load library:", package ) )
      return( invisible( FALSE ) )
    }
    else
      stop( paste( "Unable to load library:", package ) )
  }
  else {
    if( verbose )
      message( paste( "Library", package, "loaded.\n" ) )

    ## If logical return is requested here it comes.
    ## In clustermode the programm immidiately stops it a library couldn't be
    ## load on a slave. In sequentially mode it behaves like library().
    return( invisible( TRUE ) )
  }
}

##*****************************************************************************
## Include a source file.
## @todo Include complete source() parameter list.
##
## Should behave most likely the build-in source() function.
##
## Running in sequential mode: normal "source()" command (relative path)
## Running in parallel mode: the sourcefile is loaded on ALL nodes
## (abs. path and no echo).
##
## PARAMETERS: filename
## RETURN:     Logical true on success, false else. is stopOnError=TRUE => stop
##*****************************************************************************
sfSource <- function( file,
                      encoding = getOption("encoding"),
                      stopOnError = TRUE
                    ) {
  sfCheck();

  absFile <- absFilePath( file )

  if( file.exists( absFile ) ) {
    if( sfParallel() ) {
      ## Load source on all nodes (with globalized arguments)
      success <- sfClusterCall( source, file=absFile, encoding=encoding,
                                echo=FALSE, local=FALSE, chdir=FALSE )

      if( inherits( success, "try-error" ) ) {
        if( stopOnError )
          stop( paste( "Try error on cluster source call: 'source ",
                       absFile, "'", sep="" ) )
        else {
          message( paste( "Try error on cluster source call: 'source ",
                          absFile, "'", sep="" ) )
          return( FALSE )
        }
      }
      else {
        sfCat( paste( "Source", file, "loaded.\n" ) )
        message( paste( "Source", file, "loaded in cluster.\n" ) )
      }
    }

    ## Same as in sfLibrary(): include file on master as well (with
    ## original filename and echo setting).
    res <- try( source( file=absFile, encoding=encoding,
                        echo=TRUE, local=FALSE, chdir=FALSE ) )

    if( inherits( res, "try-error" ) ) {
      if( stopOnError )
        stop( paste( "Try error loading on master: '", file, "'", sep="" ) )
      else {
        message(  paste( "Try error loading on master: '", file, "'", sep="" ) )
        return( invisible( FALSE ) )
      }
    }

    return( invisible( TRUE ) )
  }
  ## File not found?
  else {
    if( stopOnError )
      stop( paste( "File does not exist:", file ) )
    else {
      message( paste( "File does not exist:", file ) )
      return( invisible( FALSE ) )
    }
  }
}

##****************************************************************************
## Export a single variable.
## On slaves, ALL exported variables are global!
## From the master, either global or local variables can be exported.
##
## PARAMETERS: ...   - Names or Variables
##             local - if TRUE, local vars will be exported
##             namespace - also exports from namespace
##             debug - with local=TRUE, prints where vars are found
##             list  - List of variable names (like snows clusterExport)
## RETURN:     -
##****************************************************************************
sfExport <- function( ..., list=NULL, local=TRUE, namespace=NULL,
                      debug=FALSE, stopOnError=TRUE ) {
  sfCheck();

  ## Export is only done if running parallel, on master all vars are visible.
  ## @TODO: Although all vars are visible, they are not surely in global env!
  ##        => export to global env.
  ##        Test if global object of this name exists. If yes: warning,
  ##        if not: assign in global space.
  if( !sfParallel() ) {
    warning( "sfExport() writes to global environment in sequential mode.\n" )
##    return( invisible( TRUE ) )
  }

  ## List of given names in dot arguments.
  names <- fetchNames( ... )

  ## If extra list is given (calling style of clusterExport()), just add this
  ## list.
  if( !is.null( list ) ) {
    ## Test from rm, see fetchNames for details.
    if( !length( list ) ||
        !all( sapply( list, function(x) is.symbol(x) || is.character(x) ) ) ) {
      if( stopOnError )
        stop( "'list' must contain names or character strings" )
      else {
        warning( "Error in sfExport: 'list' must contain names or character strings" )
        return( invisible( FALSE ) )
      }
    }

    names <- c( names, list )
  }

  for( name in names ) {
    ## Also examine namespace (from snowfall package?). Only needed for internal
    ## functions.
    if( !is.null( namespace ) && is.character( namespace ) ) {
      ## On some strange behavior, this only works with given error
      ## function. Else errors on not found objects are not caught.
      val <- tryCatch( getFromNamespace( name, namespace ),   ##, pos=-1
                       error = function(x) { NULL } )

      if( !is.null( val ) && !inherits( val, "try-error" ) ) {
        res <- sfClusterCall( assign, name, val, env = globalenv(),
                              stopOnError = FALSE )

        ## Error on export?
        if( is.null( res ) || !all( checkTryErrorAny( res ) ) ) {
          if( stopOnError )
            stop( paste( "Error exporting '", name, "': ",
                         geterrmessage(), sep="" ) )
          else {
            warning( paste( "Error exporting '", name, "': ",
                            geterrmessage(), sep="" ) )
            return( invisible(FALSE) )
          }
        }

        ## Skip local tests.
        next
      }
    }
    
    ## Check if exists before exporting.
    if( local ) {
      found <- FALSE

      # Traverse back through scopes (get() with inherit only finds
      # global presence of variables).
      # sys.nframe() at least 1 at this point, globalenv() to check last.
      for( pframe in seq( 1, sys.nframe() ) ) {
        ## Var exists in this frame?
        if( exists( name, inherits=FALSE, envir=sys.frame( -pframe ) ) ) {
          found <- TRUE

          ## If debug Messages are wanted, print these (especially for local
          ## mode with nested functions important, to locate probably
          ## overwriting variables.
          if( debug ) {
            definedIn <- gsub( "\n", "|", as.character( sys.call( -pframe ) ) )

            cat( "Export '", name, "' defined in '", definedIn, "'",
                 "\n", sep="" )
            print( get( name, envir=sys.frame( -pframe ) ) )
          }

          ## Export it.
          ## Direct call to assign is far slower as a call to a function
          ## doing it (however...)
#          res <- sfClusterCall( simpleAssign, name, 
#                                get( name,
#                                     envir=sys.frame( -pframe ) ),
#                                stopOnError = FALSE )
          ## <= 1.70
          res <- sfClusterCall( assign, name,
                                get( name, envir=sys.frame( -pframe ) ),
                                env = globalenv(),
                                stopOnError = FALSE )

          ## Error on export?
          ## 1.84: object can be null if source variable was null, too.
          if( ( is.null( res ) && !is.null( get( name, envir=sys.frame( -pframe ) ) ) )
              || !all( checkTryErrorAny( res ) ) ) {
            if( stopOnError )
              stop( paste( "Error exporting '", name, "': ",
                           geterrmessage(), sep="" ) )
            else {
              message( paste( "Error exporting '", name, "': ",
                              geterrmessage(), sep="" ) )
              return( invisible(FALSE) )
            }
          }

          break
        }
      }

      ## If variable to export is not found.
      if( !found ) {
        if( stopOnError )
          stop( paste( "Unknown/unfound variable ", name,
                       " in export. (local=", local, ")", sep="" ) )
        else {
          message( paste( "Unknown/unfound variable ", name,
                          " in export. (local=", local, ")", sep="" ) )
          return( invisible( FALSE ) )
        }
      }
    }
    ## Global export only.
    else {
      ## 1.84-3 typo
      if( exists( name, inherits=FALSE, envir=globalenv() ) ) {
#        res <- sfClusterCall( simpleAssign, name, 
#                                get( name, inherit=FALSE,
#                                     envir=globalenv() ),
#                                stopOnError = FALSE )
        ## <= 1.70
        ## 1.84-3 typo
        res <- sfClusterCall( assign, name,
                              get( name, inherits=FALSE, envir=globalenv() ),
                              env = globalenv(), stopOnError = FALSE  )

        if( is.null( res ) || !all( checkTryErrorAny( res ) ) ) {
          if( stopOnError )
            stop( paste( "Error exporting global '", name, "': ",
                         geterrmessage(), sep="" ) )
          else {
            warning( paste( "Error exporting global '", name, "': ",
                            geterrmessage(), sep="" ) )
            return( invisible( TRUE ) )
          }
        }
      }
      else {
        if( stopOnError )
          stop( paste( "Unknown variable ", name, " in export." ) )
        else {
          warning( paste( "Unknown variable ", name, " in export." ) )
          return( invisible( TRUE ) )
        }
      }
    }
  }

  invisible( TRUE )
}

##****************************************************************************
## Export all GLOBAL variables and functions to the whole cluster.
## Aware of memory usage.
##
## PARAMETERS: [Vector/List Names or variables NOT to export]
## RETURN:     Logical Success
##****************************************************************************
sfExportAll <- function( except=NULL, debug=FALSE ) {
  sfCheck();

  if( sfParallel() ) {
    ## Vector with all global variables.
    expList <- as.list( objects( pos = globalenv() ) )

    ## Now remove all variables which are listed in list except from
    ## parameters.
    if( !is.null( except ) ) {
      if( is.list( except ) )
        except <- unlist( except )

      if( !is.vector( except ) ) {
        warning( "sfExportAll: except is not a vector.\n" )
        return( invisible( FALSE ) )
      }

      ## Remove those elements which are included in except.
      ## Reverse matches for correct indices after removing of
      ## single elements (start removing on right border).
##      for( i in rev( match( except, expList ) ) )
##        expList <- expList[-i]

      ## Nicer version, proposal Greggory Jefferis (1.7.2)
      ## na.omit is not a must here though.
      expList <- expList[-na.omit(match(except, expList))]
    }

    ## Exporting mode with explicit global mode.
    sfExport( list=expList, local=FALSE )

    if( debug ) {
      message( "sfExportAll: Following variables are exported:" )
      message( paste( expList, collapse=", " ) )
    }
  }
  else {
    message( "sfExportAll() ignored in sequential mode.\n" )
    return( invisible( TRUE ) )
  }

  invisible( TRUE )
}

##****************************************************************************
## Remove objects from global environment (on whole cluster cluster)
## or at least from master (sequentially mode).
##
## PARAMETERS: List with Names(!) of the variables.
## RETURN:     -
##****************************************************************************
sfRemove <- function( ..., list=NULL, master=FALSE, debug=FALSE ) {
  sfCheck();

  ## List of given names in dot arguments.
  .sfNames <- fetchNames( ... )

  ## If extra list is given (calling style of clusterExport()), just add this
  ## list.
  if( !is.null( list ) ) {
    ## Test from rm, see fetchNames for details.
    if( !length( list ) ||
        !all( sapply( list, function(x) is.symbol(x) || is.character(x) ) ) )
      stop( "list must contain names or character strings" )

    .sfNames <- c( .sfNames, list )
  }

  ## If running parallel, remove objects from slaves.
  if( sfParallel() ) {
    if( debug )
      for( name in .sfNames )
        cat( "REMOVE:", name, "\n" )

    sfExport( ".sfNames", local=TRUE )
    sfClusterEval( rm( list=.sfNames, pos=globalenv() ) )
    sfClusterEval( rm( .sfNames, pos=globalenv() ) )
  }

  ## Remove on master as well?
  if( master )
    rm( list=.sfNames, pos=globalenv() )

  invisible( NULL )
}

##****************************************************************************
## Remove all variables from nodes (important: only global vars from nodes
## - NOT the master R process - are deleted).
## To delete on master as well, use sfRemove().
##
## PARAMETERS: [Vector/List Names of objects NOT to remove].
## RETURN:     Boolean Success (invisible)
##****************************************************************************
sfRemoveAll <- function( except=NULL, debug=FALSE, hidden=TRUE ) {
  sfCheck();

  if( sfParallel() ) {
    ## @TODO Also hidden vars?
    if( hidden )
      sfTmpAll <- sfClusterEval( ls( pos=globalenv(), all.names=TRUE ) )
    else
      sfTmpAll <- sfClusterEval( ls( pos=globalenv(), all.names=FALSE ) )

    if( length( sfTmpAll ) == 0 ) {
      message( "sfRemoveAll: problems fetching variables from nodes (or none existant)...\n" )
      return( invisible( FALSE ) )
    }

    ## Only take result from one node.
    ## We assume all nodes have exactly the same variables in global space.
    ## It may be the case, that there are different variables on each node
    ## (like a node-routine writes different vars on different cases).
    ## Take that node with the most variables in object space.
    ## @todo: Merge result lists from all nodes.
    sfTmp <- sfTmpAll[[which.max(sapply(sfTmpAll,length))]]

    ## If there are any variables on nodes.
    if( length( sfTmp ) > 0 ) {
      ## Now remove all variables which are listed in list except from
      ## parameters.
      if( !is.null( except ) ) {
        if( is.list( except ) )
          except <- unlist( except )

        if( !is.vector( except ) ) {
          warning( "sfRemoveAll: except is not a vector.\n" )
          return( invisible( FALSE ) )
        }

        ## Remove those elements which are included in except.
        ## Not very elegant... However.
        ## Bugfix see sfExportAll
        for( i in match( except, sfTmp ) )
          sfTmp[i] <- NA     ## sfTmp <- sfTmp[-i] would fail for multiple removals

        ## Remove NAs
        sfTmp <- sort( sfTmp, na.last = NA )
      }

      ## Create a new namespace vector (temporary).
#      setVar( ".sfTmpList", sfTmp )

      if( debug ) {
        message( "sfRemoveAll: Remove variables from nodes:" )
        message( paste( sfTmp, collapse=", " ) )
      }

      ## Export the list to cluster.
      sfExport( "sfTmp", local=TRUE )
     
      ## Delete all variables in the list.
      sfClusterEval( rm( list=sfTmp, pos=globalenv() ) )
      sfClusterEval( rm( "sfTmp", pos=globalenv() ) )
    }
    else {
      message( "sfRemoveAll: no variables on nodes.\n" )
      return( invisible( FALSE ) )
    }

    return( invisible( TRUE ) )
  }
  ## In sequential mode nothing is done and it counts as success ;)
  else {
    message( "sfRemoveAll() ignored in sequential mode.\n" )
    return( invisible( TRUE ) )
  }
}

##****************************************************************************
## Fast messages on the cluster (all nodes!).
## For testing proposes mainly.
##
## PARAMETER: Vector  x      Objects to print,
##            String  sep    Separator
##            Boolean master Print on master as well
##****************************************************************************
sfCat <- function( ..., sep=" ", master=TRUE ) {
  sfCheck();

  .sfTmpX   <- c( ... )
  .sfTmpSEP <- sep

  if( length( .sfTmpX ) == 0 )
    return( invisible( NULL ) )

  ## ...it's unbelievable...
  if( sfParallel() ) {
    sfExport( ".sfTmpSEP", ".sfTmpX", local=TRUE )

    sfClusterCall( cat, .sfTmpX, sep=.sfTmpSEP )

    sfRemove( ".sfTmpX", ".sfTmpSEP", master=FALSE )
  }

  ## Master&sequential mode.
  if( master )
    cat( .sfTmpX, sep=.sfTmpSEP )

  invisible( TRUE )
}

##****************************************************************************
## Mainly an parallised lapply with intermediate result savings and restore
## on later callback.
## Resultfiles are saved on each step, where a step is defined by the amount
## of CPUs given (e.g. 4 cpus, 100 steps => 25 savings).
##
## Characteristic of the called functions: they are not allowed to return NULL
## values, as these indicate uncalculated potions in the result file. Please
## use NA or any other marker for undefined values.
##
## Files are saved under directory .sfOption$RESTDIR with the form:
##    SAVE_file_name
##
## where  file : name of current R-file or "DEFAULT" in interactive mode
##        name : usergiven name for this current calculation (default: "default")
##               If a program uses more than one call to sfClusterApplySR(),
##               then name MUST be set!
##
## As this function itself calls sfLappy(), there is no explicit sequential form
## here.
##
## To disable printing of progress, set perupdate to 100.
##
## PARAMETERS: List x,         \ Like lapply
##             Function fun,   |
##             ...             /
##             [String name Name for this call of sfClusterApplySR],
##             [perupdate int Percent Update frequency for process report],
##             [Logical restore: restore previous results or don't restore]
## RETURN:     List
##****************************************************************************
sfClusterApplySR <- function( x, fun,
                              ...,
                              name="default",
                              perUpdate=NULL,
                              restore=sfRestore()
                             ) {
  sfCheck();

  checkFunction( fun )

  ## If none or no regular update frequency is given.
  if( is.null( perUpdate ) || !is.numeric( perUpdate ) ||
      ( perUpdate < 0 ) || ( perUpdate > 100 ) )
    perUpdate <- .sfOption$RESTOREUPDATE

  ## Ensure destination directory's existing, if not: create.
  if( !file.exists( .sfOption$RESTDIR ) )
    dirCreateStop( .sfOption$RESTDIR )

  ## No R-file given?
  if( is.null( .sfOption$CURRENT ) )
    setOption( "CURRENT", "DEFAULT" )

  ## Abs. file path
  file <- file.path( .sfOption$RESTDIR,
                     paste( "SAVE_", .sfOption$CURRENT, "_", name, sep="" ) )

  ## Mark this file for deletion on (regular) cluster stop - even if the file
  ## itself does not exist atm.
  addRestoreFile( file )
  
  ## Resultfile is present: try to load it, check if variable result is included
  ## and check how many results are present in the file.
  ## If it seems that the results are ok, take them and continue at their end.
  if( file.exists( file ) && restore ) {
    ## Temp global var for saving possible loading errors (in namespace).
    setVar( ".sfLoadError", "" )

    ## Load in current environment.
    tryCatch( load( file ), error=function( x ) { setVar( ".sfLoadError", x ) } )

    if( .sfLoadError != "" )
      stop( paste( "Loading error:", .sfLoadError ) )

    cat( "Restoring previous made results from file:", file, "\n" )
    errMsg <- "\nPlease remove file manually.\n"

    ## First check the contents of the file.
    ## If these don't match, do NOT remove file or overwrite it automatically,
    ## as (due to the weak filenames) mistakes could be done.
    ## Commit removal to user (with message).
    
    ## Variable "result" is loaded?
    if( length( ls( pattern="^result$" ) ) == 0 )
      stop( paste( "Result variable not found in datafile:", file, errMsg ) )

    ## Check if variable result is present.
    if( !is.list( result ) )
      stop( paste( "Variable result is no list in datafile:", file, errMsg ) )

    ## Check if variable result has correct length.
    if( length( result ) != length( x ) )
      stop( paste( "Variable result from resultfile has different length to data:",
                   length( result ), "<->", length( x ), errMsg ) )
   
    ## Set marker to NA.
    startIndex <- NA

    ## Fetch the last non-NULL value in result (which declares the last result
    ## value which does not have to be recalculated).
    for( index in seq( length( result ), 1 ) ) {
      if( !is.null( result[[index]] ) ) {
        ## Flip to first NULL value => means the first unprocessed value.
        startIndex <- index + 1
        break
      }
    }

    ## Complete unprocessed resultset found? Then start at first element.
    if( is.na( startIndex ) ) {
      startIndex <- 1
      perCent    <- 0
    }
    ## At least some parts in the resultset are given.
    else {
      ## Complete processed result? Can happen in programs with more than one
      ## parallised call.
      if( startIndex >= length( result ) ) {
        return( result )
      }

      ## Message for user where restore begins.
      perCent <- ( ( startIndex - 1 ) * 100 ) / length( result )
      cat( "Starting calculation at ", round( perCent, 1 ),
           "% (", startIndex, "/", length( result ), ")\n" )
    }
  }
  ## No resultfile given/present: generate clear result with all NULL fields.
  else {
    if( !restore )
      message( "Restore is not active! No results are loaded." )

    message( paste( "Saving results to: ", file, "\n" ) )
  
    ## Resultlist: init with NULL in any element.
    result <- lapply( 1:length( x ), function( x ) NULL )

    startIndex <- 1        # Start at the beginning
    perCent    <- 0        # Nothing done yet
  }

  lastPrintPercent <- 0
  
  ## Calculating list parts in cluster.
  for( sIndex in seq( startIndex, length( x ), by=sfCpus() ) ) {
    ## Endindex.
    eIndex <- sIndex + sfCpus() - 1

    ## End out of bounds?
    if( eIndex > length( x ) )
      eIndex <- length( x )

##    cat( "Calculating Indizes: ", sIndex, eIndex, "\n" )

    newResult <- sfLapply( x[sIndex:eIndex], fun, ... )

    ## Fill cells with new results.
    result[sIndex:eIndex] <- newResult[1:length( newResult )]

    ## Intermediate save of current results.
    save( result, file=file )
    
    ## Calculated percentage.
    perCent <- eIndex * 100 / length( result )

    ## If message about process is wanted, print it (also is a connector for
    ## sfCluster to show the calculation process).
    ## Also catch the case where mod rarely matches (if amount of CPUs is
    ## bigger than perUpdate).
    if( ( ( round( perCent, 0 ) -
            round( sIndex * 100 / length( result ), 0 ) ) >= perUpdate )
        ||
        ( ( ( round( perCent, 0 ) %% perUpdate ) == 0 ) &&
          ( ( round( perCent, 0 ) - lastPrintPercent ) >= perUpdate ) ) ) {
      cat( "SR '", name, "' processed: ",
           round( perCent, 1 ), "%\n", sep="" )
      lastPrintPercent <- round( perCent, 0 )
    }

##    cat( "Finished Indizes: ", sIndex, eIndex, "\n" )
  }
   
  return( result )
}

##****************************************************************************
## Complete "unit test" or most of the buildin functions.
## Mainly integrated for development, but can be used for testing the
## R functionality on all nodes, too.
##
## PARAMETER: -
## RETURN:    Int amount of errors (0: everything is ok).
##****************************************************************************
sfTest <- function() {
  sfCheck();

  if( !sfParallel() ) {
    message( "Tests only work in parallel mode." )
    return( invisible( FALSE ) )
  }

  ##***************************************************************************
  ## Basic checks for Calls/Evals.
  ##***************************************************************************
  checkResultBasic <- function( result ) {
    if( is.null( result ) )
      return( c( FALSE, "Result was NULL" ) )

    if( !is.list( result ) )
      return( c( FALSE, "No proper return type (no list)." ) )

    if( length( result ) != sfCpus() )
      return( c( FALSE, "No proper return type (wrong length)." ) )

    if( inherits( result, "try-error" ) )
      return( c( FALSE, "TRY-ERROR raised on result." ) )

    if( !all( sapply( result,
                      function( x )
                        if( inherits( x, "try-error" ) )
                          return( FALSE )
                        else
                          return( TRUE )
                    ) ) )
      return( c( FALSE, "Result elements raised TRY-ERROR(s)." ) )

    return( c( TRUE, "" ) )
  }

  ##***************************************************************************
  ## Checks if each element of a given list is equal to a certain value !=
  ## NA/NULL).
  ##***************************************************************************
  checkAllEqual <- function( result, equal ) {
    if( !all( sapply( unlist( result ), function( x ) return( x == equal ) ) ) )
      return( FALSE )

    return( TRUE )
  }

  ##***************************************************************************
  ## Test a list of lists against a vector (each sublist must be equal to the
  ## given list).
  ##***************************************************************************
  checkAllEqualList <- function( result, equal ) {
    if( is.list( equal ) )
      equal <- sort( unlist( equal ) )
    else
      equal <- sort( equal )

    for( res in result ) {
##      res <- sort( res )

      if( ( length( res ) != length( equal ) ) ||
          ( length( which( sort( res ) == equal ) ) < length( equal ) ) )
        return( FALSE )
      
##      i <- 1
##     while( i <= length( res ) ) {
##        if( res[i] != equal[i] )
##          return( FALSE )
##        i <- i + 1
##      }
    }

    return( TRUE )
  }

  ##***************************************************************************
  ## Compare vectors.
  ##***************************************************************************
  checkVecCmp <- function( x, y ) {
    if( length( x ) != length( y ) )
      return( FALSE )

    for( i in seq( 1, length( x ) ) ) {
      ## If NULL or NA, nothing is to compare. But only if both vals are NA/NULL
      ## If not, compare (throws exception/stop)
      if( ( is.na( x[i] ) && is.na( y[i] ) ) ||
          ( is.null( x[i] ) && is.null( y[i] ) ) )
        next

      if( x[i] != y[i] )
        return( FALSE )
    }

    return( TRUE )
  }

  ##***************************************************************************
  ## Testing sfLibrary.
  ##***************************************************************************
  testLib <- function() {
    ## Package always be installed.
    if( !sfLibrary( "boot", character.only=TRUE, stopOnError=FALSE ) )
      return( c( FALSE, "Unable to load library 'tools'" ) )

    ## calcium is a dataframe.
    result <- sfClusterEval( as.matrix( get('calcium') )[,2] )

    ## Compare if all nodes delivered the same data for variable "calcium"
    ## get needed to avoid R CMD check warnings.
    for( res in result )
      if( !checkVecCmp( res, as.matrix( get( "calcium" ) )[,2] ) )
        return( c( FALSE, "Wrong data delivered..." ) )

    ## Load surely uninstalled package to test if lib call fail safely.
    if( try( sfLibrary( "xxxyyyzzz", character.only=TRUE, stopOnError=FALSE ),
             silent=TRUE ) )
      return( c( FALSE, "Irregular return on loading inexisting library." ) )

    return( c( TRUE, "ok" ) )
  }

  ##***************************************************************************
  ## testing sfSource.
  ##***************************************************************************
  testSource <- function() {
    sfRemoveAll()

    ## Find path of the installed snowfall Package.
    res <- NULL
    res <- try( find.package( "snowfall" ) )   ## CHG 131712 from .find.package

    if( inherits( res, "try-error" ) )
      return( c( FALSE, paste( "Exception: cannot locate package snowfall.",
                               geterrmessage() ) ) )

    if( is.null( res ) )
      return( c( FALSE, "Cannot locate package snowfall." ) )

    res <- file.path( res, "data", "test.R" )

    cat( "PACKAGE...: ", res, "\n" )

    con <- file( res, "r", blocking=FALSE )
    a <- readLines( con, n=-1 )
    debug( "test.R content:" )
    debug( a )

    result <- sfSource( res, stopOnError=FALSE )

    if( inherits( result, "try-error" ) )
      return( c( FALSE, paste( "Exception: cannot source on slaves.",
                               geterrmessage() ) ) )

    ## get to satisfy R CMD check
    result <- sfClusterEval( get("f1")(), stopOnError=FALSE )

    resBasic <- checkResultBasic( result )

    if( resBasic[1] == FALSE )
      return( resBasic )

    if( !checkAllEqual( result, 999 ) )
      return( c( FALSE, "Wrong results on sourced function f1." ) )

    ## get to satisfy R CMD check
    result <- sfClusterEval( get("f2")( 99, 1 ), stopOnError=FALSE )

    resBasic <- checkResultBasic( result )

    if( resBasic[1] == FALSE )
      return( resBasic )

    if( !checkAllEqual( result, 100 ) )
      return( c( FALSE, "Wrong results on sourced function f2." ) )
    
    return( c( TRUE, "ok" ) )
  }
  
  ##***************************************************************************
  ## Testing sfClusterCall. Allways first test.
  ##***************************************************************************
  testCall <- function() {
    # Test 1 on Call
    result <- sfClusterCall( paste, "a", "b", "c", sep="", stopOnError=FALSE )

    resBasic <- checkResultBasic( result )

    if( resBasic[1] == FALSE )
      return( resBasic )

    if( !checkAllEqual( result, "abc" ) )
      return( c( FALSE, "Wrong results on paste." ) )

    # Test 2 on Call
    sums <- c( 99, 7, 3.4 )
    result <- sfClusterCall( sum, sums, stopOnError=FALSE )

    if( !checkAllEqual( result, sum( sums ) ) )
      return( c( FALSE, "Wrong result on sum." ) )

    return( c( TRUE, "ok" ) )
  }

  ##***************************************************************************
  ## Testing sfClusterEval
  ##***************************************************************************
  testEval <- function() {
    # Test 1 on Eval
    result <- sfClusterEval( sum( sapply( 1:10, exp ) ) )
   
    resBasic <- checkResultBasic( result )

    if( resBasic[1] == FALSE )
      return( resBasic )

    if( !checkAllEqual( result, sum( sapply( 1:10, exp ) ) ) )
      return( c( FALSE, "Wrong results on sum." ) )

    return( c( TRUE, "ok" ) )
  }

  ##***************************************************************************
  ## Testing Export Funktion
  ##***************************************************************************
##  testExport <- function() {
##    ## Needed to have a clean comparison global env.
##    sfRemoveAll( hidden=TRUE )
##
##    vars <- sfClusterEval( ls( all.names=TRUE, envir=globalenv() ) )
##
##    print( vars )
##    
##     if( length( vars ) != 0 )
##       if( !all( sapply( vars,
##                         function( x ) return( length( x ) == 0 ) ) ) )
##         return( c( FALSE, "sfRemoveAll() didn't kill everything" ) )

##     ## Setting global variable via assign, as <<- invokes warnings on
##     ## package check.
##     assign( "var1", 99, pos=globalenv() )
##     assign( "var2", 101, pos=globalenv() )
## #    var1 <<- 99    # Global
## #    var2 <<- 101
##     var3 <-  103   # Local
##     var4 <-  7

##     ## Setting var in namespace ("snowfall").
##     setVar( ".sfTestVar5", 77 )

##     if( getVar( ".sfTestVar5" ) != 77 )
##       return( c( FALSE, "Access to namespace failed." ) )
    
##     iTest <- function() {
##       var3 <- 88

##       res <- FALSE

##       res <- sfExport( "var1", "var2", ".sfTestVar5",
##                        list=list( "var3", "var4" ),
##                        local=TRUE, namespace="snowfall", stopOnError=FALSE )

##       if( inherits( res, "try-error" ) )
##         return( c( FALSE, "Exception on export." ) )

##       if( !res )
##         return( c( FALSE, "Unexpected Exception on export." ) )

##       print( "GLOBALENV..." )
##       print( sfClusterCall( ls, envir=globalenv() ) )
      
##       if( !checkAllEqualList( sfClusterCall( ls, all.names=TRUE,
##                                              envir=globalenv() ),
##                               c( "var1", "var3", "var2", "var4",
##                                  ".sfTestVar5" ) ) )
##         return( c( FALSE, "Not all vars exported." ) )

##       ## get to satisfy R CMD check
##       if( !checkAllEqual( sfClusterEval( get("var1") ), 99 ) ||
##           !checkAllEqual( sfClusterEval( get("var2") ), 101 ) )
##         return( c( FALSE, "Error exporting global var." ) )

##       ## get to satisfy R CMD check
##       if( !checkAllEqual( sfClusterEval( get("var3") ), 88 ) ||
##           !checkAllEqual( sfClusterEval( get("var4") ), 7 ) )
##         return( c( FALSE, "Error exporting local var." ) )

##       if( !checkAllEqual( sfClusterEval( get(".sfTestVar5") ), 77 ) )
##         return( c( FALSE, "Error exporting namespace var." ) )
      
##       ## Test removeAll with Exception-List
##       sfRemoveAll( except=list( "var2", "var3" ) )

##       if( !checkAllEqualList( sfClusterCall( ls, envir=globalenv() ),
##                               list( "var2", "var3" ) ) )
##         return( c( FALSE, "Error on removeAll except-list." ) )

##       sfRemoveAll()

##       return( c( TRUE, "ok" ) )
##     }   
    
##     return( iTest() )
##   }    


  ##***************************************************************************
  ## Testing Calculation Function Part 1
  ##***************************************************************************
  testCalc1 <- function() {
    size <- 50
    mat  <- matrix( 0, size, size )

    for( var in 1:nrow( mat ) ) mat[var,] = runif( nrow( mat ) )

    rSum <- function( row, mat ) {
      s <- 0
      for( col in 1:ncol( mat ) )
        s <- s + mat[row,col]
      return( s )
    }

    cmp <- unlist( lapply( seq( 1, ncol( mat ) ), rSum, mat ) )

    sfExport( "rSum", local=TRUE, debug=TRUE )

    # Test 1 on Eval
    result <- sfLapply( seq( 1, ncol( mat ) ), rSum, mat )

##    cat( "FINISHED...\n" )
##    print( result )
    
    if( !checkVecCmp( unlist( result ), cmp ) )
      return( c( FALSE, "Wrong results on sfLapply." ) )

    ## Testing sfClusterApplyLB
    result <- sfClusterApplyLB( seq( 1, ncol( mat ) ), rSum, mat )

    if( !checkVecCmp( unlist( result ), cmp ) )
      return( c( FALSE, "Wrong results on sfClusterApplyLB." ) )

    ## Testing sfClusterApplySR
    result <- sfClusterApplySR( seq( 1, ncol( mat ) ), rSum, mat,
                                name="TEST",
                                restore=FALSE, perUpdate=100 )

    if( !checkVecCmp( unlist( result ), cmp ) )
      return( c( FALSE, "Wrong results on sfClusterApplySR." ) )

    ## As clusterApply only works for #nodes samples, reduce data size depending
    ## on it.
    result <- sfClusterApply( seq( 1, min( sfCpus(), ncol( mat ) ) ), rSum, mat )

    if( !checkVecCmp( unlist( result ),
                      unlist( lapply( seq( 1, min( sfCpus(), ncol( mat ) ) ), rSum, mat ) ) ) )
      return( c( FALSE, "Wrong results on sfClusterLapply." ) )
    
    return( c( TRUE, "ok" ) )
  }

  ##***************************************************************************
  ## Testing Calculation Function Part 2
  ## Further snow Wrappers.
  ##***************************************************************************
  testCalc2 <- function() {
    size <- 50
    mat1 <- matrix( 0, size, size )
    mat2 <- matrix( 0, size, size )

    for( var in 1:nrow( mat1 ) ) mat1[var,] = runif( nrow( mat1 ) )
    for( var in 1:nrow( mat2 ) ) mat2[var,] = runif( nrow( mat2 ) )

    matRes1 <- sfMM( mat1, mat2 )
    matRes2 <- mat1 %*% mat2

    for( row in seq( 1, size ) )
      if( !checkVecCmp( matRes1[row,], matRes2[row,] ) )
        return( c( FALSE, "Wrong results on sfParMM." ) )
   
    return( c( TRUE, "ok" ) )
  }

  
  ##***************************************************************************
  ## Run single test (with given functionname)
  ##***************************************************************************
  runTest <- function( fun ) {
    cat( "Run test: ", fun, "\n" )

    res <- c( NA, "" )
    res <- try( do.call( fun, list() ) )

    if( inherits( res, "try-error" ) )
      return( c( FALSE, paste( "TRY-ERROR on testCall()", geterrmessage() ) ) )

    if( is.null( res ) || !is.vector( res ) || is.na( res[1] ) )
      return( c( FALSE, paste( "Hidden exception on test.", geterrmessage() ) ) )
   
    return( res )
  }

  complete <- list( errors=0, warnings=0 )

  ## @todo - Bibliotheken / Source
  ## @todo - anderen Applies / parMM
  ## @todo - exportAll
  ## testExport removed because of a R 3.0.0 warning (not error!)
  tests <- c( "testCall",
              "testEval",
##              "testExport",
              "testCalc1",
              "testCalc2",
              "testLib",
              "testSource" )

  ## Run tests.
  for( test in tests )
    complete[[test]] <- runTest( test )

  ## Print results.
  cat( "\n\nRESULTS ON TEST:\n\n" )
  errors <- 0

  for( test in tests ) {
    ## Align names to same length.
    if( as.logical( complete[[test]][1] ) )
      cat( test, sapply( seq( 1, 13 - nchar( test ) ), function( x ) return( " " ) ),
           ": ok", "\n", sep="" )
    else {
      cat( test, sapply( seq( 1, 13 - nchar( test ) ), function( x ) return( " " ) ),
           ": FAILED! (", complete[[test]][2], ")\n", sep="" )

      errors <- errors + 1
    }
  }

  cat( "\n----------------------------\n", errors, "tests failed.\n\n" )

  return( invisible( errors ) )
}
