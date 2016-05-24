##*****************************************************************************
## Function for initialisizing the Cluster.
##
## Also the predefinition of (internal) global variables is done here,
## mainly because of the code check of R.
##
## Compability issue: snowfall needs to know if sfCluster is working (also
## old versions of sfCluster).
## --session could be set by other solutions.
## So, --lockfile is decided to be the sfCluster indicator, as this most
## likely will have no use if not used with sfCluster.
## Therefore setting of --lockfile (LOCKFILE) can cause troubles.
##*****************************************************************************

## These variables are used by internal function. Need to declared for
## the compiler warnings. As these cannot be altered directly, setOption() does
## this (in the namespace). So no global objects are set, internal state is
## kept inside the namespace.
DEBUG          <- FALSE     ## Static switch for debugging messages
.sfOption      <- list()    ## Configuration of this master
.sfPresetCPUs  <- 0         ## Presetted CPU amount (max. allocatable).

## Some vars needed at specific points, which can handled correctly otherwise,
## but would raise R CMD check warnings if not defined before.
.sfPars        <- ''        ## Tmp. var for sfLibrary.
.sfLoadError   <- ''        ## Tmp. var for loading.
.sfTestVar5    <- 0         ## Exporting test in sfTest().

##*****************************************************************************
## Function for initialisizing the Cluster.
##
## Attention: this package does nasty things with explicit sideeffects (not
## only using "require" and "options" etc.)...
##
## "Nodes" and "CPUs" are used identical (unbeautiful).
##
## PARAMETER:  [Boolean parallel - overwrite CMDline settings],
##             [Int nodes        - overwrites commandline settings / DEPRECATED]
##             [Int cpus         - overwrites commandline settings]
##             [Boolean nostart  - If set, no cluster start will be run.
##                                 Needed for nested usage of snowfall]
##             [Boolean restore  - Globally set restore]
##             [String type      - {'MPI','SOCK', 'PVM', 'NWS'}
##             [Vector socketHosts - List of all hosts used in socketmode]
##             [slaveOutfile     - filename for output on slaves]
##             [useRscript       - Startup via R Script or shellscript. Only snow>0.3]
## RETURN:     Boolean TRUE
##*****************************************************************************
sfInit <- function( parallel=NULL,
                    cpus=NULL,
                    type=NULL,
                    socketHosts=NULL,
                    restore=NULL,
                    slaveOutfile=NULL,
                    nostart=FALSE,
                    useRscript=FALSE      ## snow: Default is TRUE.
                   ) {
  ## Flag for detection of reconnect (means: non-first calls to sfInit())
  reconnect <- FALSE

  ## Get rid of that stupid data load to global env.
  initEnv <- new.env()
  
  ## Saves users from many own if-clauses probably.
  if( nostart ) return( TRUE )

  ## Are options setted?
  if( length( .sfOption ) == 0 ) {
    debug( "Setup sfOption..." )

    ## Add 1.62: list from sysdata cleared and created again
    setOption( "parallel", FALSE )
    setOption( "session", NULL )
    setOption( "priority", 1 )
    setOption( "nodes", 1 )
    setOption( "stopped", FALSE )
    setOption( "init", FALSE )

    ## Load configuration file: delivered with package and changeable by user.
    data( "config", package="snowfall", envir=initEnv )
    configM <- as.matrix( t( config ) )
    config  <- as.list( configM )
    names( config ) <- dimnames( configM )[[2]]

    ## Node count are limited in snowfall as well (as it is useable without
    ## sfCluster) and you probably don't want an arbitrary amount of CPUs
    ## requested by a DAU.
    ## If changed preset exists, take this number.
    if( .sfPresetCPUs > 0 )
      setOption( "MAXNODES", .sfPresetCPUs )
    else
      setOption( "MAXNODES", as.numeric( config[["MAXNODES"]] ) )

    ## Startup lockfile (only coming from sfCluster and if available
    ## signalling that snowfall is started through sfCluster).
    ## LOCKFILE can only be set through commandline --lockfile
    setOption( "LOCKFILE", "" )

    ## Temporary directory (for logfiles, esp. on the slaves)
    ## Only if set, if not, take default.
    if( as.character( config[["TMPDIR"]] ) != "-" )
      setOption( "TMPDIR", path.expand( as.character( config[["TMPDIR"]] ) ) )
    else {
      ## Default tempdir on Unix systems is R session tempdir
      if( .Platform$OS.type == "unix" )
        setOption( "TMPDIR", file.path( Sys.getenv( "R_SESSION_TMPDIR" ), "sfCluster" ) )
      ## On any non *nix system: take local dir (R_SESSION_TMPDIR unset on Win)
      else
        setOption( "TMPDIR", "" )
    }

    ## Addition variables for save/restore (only used in sfClusterApplySR).
    setOption( "RESTOREFILES", NULL )   ## List with restore files (for cleanup)
    setOption( "RESTOREUPDATE", 5 )     ## Updates percent output any 5%
    setOption( "RESTORE", FALSE )       ## Restore previous results?
    setOption( "CURRENT", NULL )        ## Currently executed R-File

    ## Default cluster type (unchangeable by config to ensure runnability
    ## of a specific code in any setting).
    setOption( "type", "SOCK" )
    setOption( "sockHosts", NULL )
    
    ## Restore file directory (for saved intermediate results) - not neccessary
    ## under/in TMPDIR.
    ## (As log files prob woul be set in global dir, restore files should be
    ## stored under users home - as they don't contain a session-ID or something
    ## generic unique thing to differ them.
    if( as.character( config[["RESTDIR"]] ) != "-" )
      setOption( "RESTDIR", path.expand( as.character( config[["RESTDIR"]] ) ) )
    else
      setOption( "RESTDIR", file.path( Sys.getenv( "HOME" ), ".sfCluster", "restore" ) )

    ## Remove config (as data() writes it as global variable).
    rm( config, envir=initEnv )
    #pos=globalenv() )
  }
  ## If .sfOption exists, sfInit() was called before: restart.
  ## (sfCluster should be able to handle this - although slaves are iterated and
  ## slave killing is only done through snow).
  else {
    reconnect <- TRUE

    if( .sfOption$stopped && !.sfOption$init )
      debug( "Irregluar init state (error on previous init)..." )

    ## If not stopped, but initialised.
    if( !.sfOption$stopped && .sfOption$init ) {
      message( "Explicit sfStop() is missing: stop now." )
      sfStop()
    }
  }

  ##**************************************************************************
  ## Values for parallel/session can be in the commandline or the environment.
  ## Function parameters overwrite commandline.
  ##**************************************************************************
  searchCommandline( parallel, cpus=cpus, type=type, socketHosts=socketHosts,
                     restore=restore )

  if( getOption( 'verbose' ) && !reconnect )
    print( .sfOption )

  ## If given restore-directory does not exist, create it.
  ## if( !file.exists( .sfOption$RESTDIR ) ) {
  ##   ## 1.62: removed
  ##   ##    .sfOption$RESTDIR <<- path.expand( "~/.sfCluster/restore" )
  ##   dirCreateStop( .sfOption$RESTDIR )
  ## }

  ## Running in parallel mode? That means: Cluster setup.
  ## Will be blocked if argument "nostart" is set (for usage of snowfall
  ## inside of packages).
  if( .sfOption$parallel && !nostart ) {
    ## Internal stopper. Running in parallel mode a session-ID is needed.
    ## For testing purposes can be anything (mainly used for pathnames
    ## of logfiles).
    if( startedWithSfCluster() && is.null( .sfOption$session ) )
      stop( "No session-ID but parallel run with sfCluster (something went wrong here?)..." )
    ## @TODO regenerate session id if missing.

    ## If amount of nodes not set via commandline, then it will be 2
    if( is.null( .sfOption$nodes ) || is.na( as.numeric( .sfOption$nodes ) ) )
      setOption( "nodes", 2 )
    else
      setOption( "nodes", as.numeric( .sfOption$nodes ) )

    ## Preload required libraries if needed (as an extended error check).
    libList <- list( "PVM"="rpvm", "MPI"="Rmpi", "NWS"="nws", "SOCK"="" )

    if( libList[[.sfOption$type]] != "" ) {
      if( !require( libList[[.sfOption$type]], character.only=TRUE ) ) {
        message( paste( "Failed to load required library:", libList[[.sfOption$type]],
                        "for parallel mode", .sfOption$type, "\nFallback to sequential execution" ) )

        ## Fallback to sequential mode.
        return( sfInit( parallel=FALSE ) )
      }
      else
        message( paste( "Library", libList[[.sfOption$type]], "loaded." ) )
    }
    
    ## In any parallel mode, load snow if needed.
    ## CHG 131217: not needed because of package depends.
   ## if( !require( snow ) ) {
   ##   message( paste( "Failed to load library 'snow' required for parallel mode.\n",
   ##                   "Switching to sequential mode (1 cpu only)!." ) );

   ##   ## Fallback to sequential mode.
   ##   return( sfInit( parallel=FALSE ) )
   ## }

    ## Chg. 1.62
    ## Temporary file for output.
    ## If sfCluster is running (LOCKFILE given): session is taken.
    ## If sfCluster not running but user setted slaveOutfile option: take arg.
    ## Else (default): no slave outfiles (writing to /dev/null|nul).
    if( startedWithSfCluster() ) {
      tmp <- file.path( .sfOption$TMPDIR,
                        paste( "rout_", .sfOption$session, sep="" ) )

      ## Only create temporary directory once and if needed.
      ## Only needed if running with sfCluster. If user sets it's own
      ## slaveOutfile, he has to ensure himself about existing pathes.
      ## If needed create temporary path. Problem: this is executed only on
      ## master, not on slaves. The clusterstarter needs to manage this.
      if( !reconnect )
        dirCreateStop( .sfOption$TMPDIR )
    }
    else
      tmp <- ifelse( is.null( slaveOutfile ), '/dev/null', slaveOutfile )

    ## @TODO Exception handler.
    ## @TODO Timeout on init.
    ## Ebenso: Timeout - das ist extrem hässlich, wenn das Cluster nicht
    ## korrekt startet und hängen bleibt (z.B. wenn zuviele CPUs für das
    ## Cluster angefordert werden - was PVM schluckt, macht MPI anscheinend
    ## Kopfzerbrechen).
    setDefaultClusterOptions( type = .sfOption$type )
    setDefaultClusterOptions( homogenous = FALSE )

    ## On socket connections the list of hosts needs to be given.
    ## If no is set, use localhost with default R.
    if( .sfOption$type == "SOCK" ) {
      ## No host information given: use localhost with wished CPUs.
      ## Else: host settings overwrite wished CPUs (important for error checks!).
      if( is.null( .sfOption$sockHosts ) || ( length( .sfOption$sockHosts ) == 0 ) )
        setOption( "sockHosts", c( rep( "localhost", .sfOption$nodes ) ) )
      else
        setOption( "nodes", length( .sfOption$sockHosts ) )

      setOption( "cluster", try( makeCluster( .sfOption$sockHosts,
                                              type = "SOCK",
                                              outfile = tmp,
                                              homogenous = TRUE
                                            ) ) )
    }
    # PVM cluster
    else if( .sfOption$type == "PVM" ) {
      setOption( "cluster", try( makeCluster( .sfOption$nodes,
                                              outfile = tmp ) ) )
    }
    # Network Spaces
    else if( .sfOption$type == "NWS" ) {
      if( is.null( .sfOption$sockHosts ) || ( length( .sfOption$sockHosts ) == 0 ) )
        setOption( "sockHosts", c( rep( "localhost", .sfOption$nodes ) ) )
      else
        setOption( "nodes", length( .sfOption$sockHosts ) )

        ## Patch Markus Schmidberger (Mail 11/25/2008).
        setOption( "cluster", try( makeNWScluster(
                                   .sfOption$sockHosts[1:.sfOption$nodes],
                                   type = "NWS",
                                   outfile = tmp
                                 ) ) )
    }
    # MPI cluster (also default for irregular type).
    else {
      ## 1.81: useRScript must be FALSE. Else sfCluster wont work
      ##       with snow > 0.3 (on older snow Versions this option
      ##       is ignored. Also homogenous is always on.
      ## 1.83: But for non-sfCluster usage at least it has to be modifyable.
      setOption( "cluster", try( makeMPIcluster( .sfOption$nodes,
                                                 outfile = tmp,
                                                 homogenous = TRUE,
                                                 useRscript = useRscript
                                               ) ) )
    }

    ## Startup successfull? If not: stop.
    if( is.null( .sfOption$cluster ) ||
        inherits( .sfOption$cluster, "try-error" ) )
      stop( paste( "Starting of snow cluster failed!",
                   geterrmessage(), .sfOption$cluster ) )

    ## Cluster setup finished. Set flag (used in error handlers and stop).
    ## Also: no function can be called if init is not set.
    setOption( "init", TRUE )
    setOption( "stopped", FALSE )

    if( !reconnect ) {
      ## As Snow Init spawn all the requires R-processes, the proprietary
      ## lockfile can be deleted now (if it exists).
      ## Problem: now all R procs are spawned, but the observer most
      ## likely didn't catch them until the next time of his observing
      ## loop.
      if( !is.null( .sfOption$LOCKFILE ) && file.exists( .sfOption$LOCKFILE ) ) {
        if( unlink( .sfOption$LOCKFILE ) != 0 )
          warning( "Unable to remove startup lockfile: ", .sfOption$LOCKFILE )
        else
          message( "Startup Lockfile removed: ", .sfOption$LOCKFILE )
      }
    
      if( getOption( 'verbose' ) ) {
        if( tmp == '/dev/null' )
          message( "Slave output suppressed. Use 'slaveOutfile' to activate." )
        else
          message( paste( "Temporary log for STDOUT/STDERR (on each node): ", tmp, "\n",
                          "Cluster started with", .sfOption$nodes, "CPUs.", "\n" ) )
      }
      else
        debug( paste( "Temporary log for STDOUT/STDERR (on each node): ", tmp, "\n",
                      "Cluster started with", .sfOption$nodes, "CPUs.", "\n" ) )
   
      ## Write R-Version and Time in (slave-)logfiles.
      .startInfo <- strsplit( Sys.info(), "\n" );
      .startMsg <- paste( sep="",
                          "JOB STARTED AT ", date(),      # Global Var!
                          " ON ", .startInfo$nodename, " (OS", .startInfo$sysname,
                          ") ", .startInfo$release, "\n" )

      sfExport( ".sfOption", ".startMsg", local=TRUE, namespace="snowfall", debug=DEBUG )
      sfCat( .startMsg, "\n", master=FALSE )    ## No master
      sfCat( paste( "R Version: ", R.version$version.string, "\n\n" ) )

      ## Remove starting message.
      sfRemove( ".startMsg" )
    }
    ## @TODO Checken, ob dieser Export wirklich noch benötigt ist (jan 10)
    else
      sfExport( ".sfOption", local=FALSE, namespace="snowfall" )
  }
  ## Sequential mode or option "nostart":
  ## init will be set. If someone calls sfInit with nostart and aims
  ## it to be started, it's his or her problem.
  else {
    ## Cluster setup finished. Set flag (used in error handlers and stop).
    ## Also: no function can be called if init is not set.
    setOption( "init", TRUE )
    setOption( "stopped", FALSE )

    setOption( "cluster", NULL )
  }
    
  ## Print init Message (esp. print State of parallel and snowfall
  ## version.
  if( sfParallel() ) {
    message( paste( "snowfall ", packageDescription( "snowfall" )$Version,
                    " initialized (using snow ", packageDescription( "snow" )$Version,
                    "): parallel execution on ", sfCpus(), " CPUs.\n", sep="" ) );
  }
  else {
    message( paste( "snowfall", packageDescription( "snowfall" )$Version,
                    "initialized: sequential execution, one CPU.\n" ) );
  }

  return( invisible( TRUE ) )
}

##*****************************************************************************
## Check if sfInit() was called.
## This function is called before any function which need initialised cluster.
##
## Previous it stops with error, now it calls sfInit() without parameters,
## so sfInit() does not have to be called explicitely (requested from Harald).
##
## (Not exported to namespace).
##*****************************************************************************
sfCheck <- function() {
  if( !sfIsRunning() ) {
    message( paste( "Calling a snowfall function without calling 'sfInit'",
                    "first or after sfStop().\n'sfInit()' is called now." ) )
    return( invisible( sfInit() ) )
  }

  return( invisible( TRUE ) )
}

##*****************************************************************************
## Exported as userfunction.
## Give the user information if sfInit() was called and cluster is not stopped.
## (Maybe helpful inside of packages etc.).
##*****************************************************************************
sfIsRunning <- function() {
  ## Add 1.62: stopped as argument
  if( ( length( .sfOption ) == 0 ) || !.sfOption$init || .sfOption$stopped )
    return( FALSE )
  else
    return( TRUE )
}

##*****************************************************************************
## Stop the (snow)-Cluster. Just calls Snows stopCluster.
##
## PARAMETER: [Boolean nostop: don't stop]
##*****************************************************************************
sfStop <- function( nostop=FALSE ) {
  ## Saves users from many own if-clauses probably.
  if( nostop ) return( TRUE );

  if( exists( ".sfOption" )  && ( length( .sfOption ) > 0 ) ) {
    ## Only stop if initialisized and running parallel.
    if( !.sfOption$stopped && .sfOption$init && .sfOption$parallel ) {
      message( "\nStopping cluster\n" )

      ## Stopping snow cluster.
      ## NO call to sfGetCluster() here, as sfGetCluster sfCheck()s again.
      stopCluster( .sfOption$cluster )
    }

    ## Reset default values.
    ##.sfOption$init     <<- FALSE
    setOption( "stopped", TRUE )
    setOption( "parallel", FALSE )

    ## Delete probably stored resultfiles (can also be used in sequential mode!)
    deleteRestoreFiles()
  }

  invisible( NULL )
}

##*****************************************************************************
## Is programm running parallel? Wrapper for internal Optionblock (therefore
## exported of course).
## Also: get cluster Handler (prob. not exported in the final).
##
## RETURN: Boolean Running in parallel mode
##*****************************************************************************
sfParallel <- function() {
  sfCheck()

  return( .sfOption$parallel )
}

##*****************************************************************************
## Shall sfClusterApplySR restore results?
##*****************************************************************************
sfRestore <- function() {
  sfCheck()

  return( .sfOption$RESTORE )
}

##*****************************************************************************
## Receive snow cluster handler (for direct calls to snow functions).
##*****************************************************************************
sfGetCluster <- function() {
  sfCheck()

  return( .sfOption$cluster )
}

##*****************************************************************************
## Receive amount of currently used CPUs (sequential: 1).
##*****************************************************************************
sfCpus <- function() {
  sfCheck()

  return( .sfOption$nodes )
}

## getter for amount of nodes. Wrapper for sfCPUs.
sfNodes <- function() return( sfCpus() )

##*****************************************************************************
## Receive type of current cluster.
##*****************************************************************************
sfType <- function() {
  sfCheck()

  if( sfParallel() )
    return( .sfOption$type )
  else
    return( "- sequential -" )
}

##*****************************************************************************
## Receive list with all socket hosts.
##*****************************************************************************
sfSocketHosts <- function() {
  if( sfType() == "SOCK" ) {
    sfCheck()

    return( .sfOption$sockHosts )
  }
  else {
    warning( paste( "No socket cluster used:", sfType() ) )

    return( invisible( NULL ) )
  }
}

##*****************************************************************************
## getter for session-ID.
##*****************************************************************************
sfSession <- function() {
  sfCheck();

  return( .sfOption$session )
}

##*****************************************************************************
## Increase max. numbers of CPUs used per process.
## No check for sensefull values (if user wants 1000, you get 1000 :)).
##*****************************************************************************
sfSetMaxCPUs <- function( number=32 ) {
  setVar( ".sfPresetCPUs", number )
}

##*****************************************************************************
## Internal function:
##
## Search commandline arguments for Parallel and Session values.
## If there are arguments on function call, these overwrites the values on the
## commandline.
##
## Basically the arguments on the commandline come from sfCluster, but of
## course set manually or via another load- or sessionmanager.
##
## Commandline arguments: --parallel(=[01])*
##                        --session=\d{8}
##                        --nodes=\d{1,2}
##                        --tmpdir=\/[a-z_].*
##                        --hosts=((\s+:\d+))+
##                        --restoreDir=\/[a-z_].*
##                        --restoreSR
##                        --lockfile
## Results will be saved in options .parallel (bool) and .session (8 chars)
##*****************************************************************************
searchCommandline <- function( parallel=NULL, cpus=NULL,
                               socketHosts=NULL, type=NULL,
                               restore=NULL ) {
#  if( !exists( ".sfOption", envir=globalenv() ) )
#    stop( "Global options missing. Internal error." )

  ## If set, copy to sfCluster data structure.
  if( !is.null( cpus ) ) {
    setOption( "nodes", max( 1, cpus ) )

    ## For socket/NWS clusters: force rebuild of hostlist (as probably changed).
    ## (If not overwritten later by users own arguments).
    setOption( "sockHosts", NULL )
    
    ## If more than one CPU is wanted, parallel mode is forced.
    ## Probably this is not an intended behavior.
#    if( .sfOption$nodes > 1 ) {
#      ## Potential misuse of argument: inform user.
#      if( !is.null( parallel ) && ( parallel == FALSE ) )
#        warning( "Explicit parallel=FALSE, but required >1 CPUs ==> parallel mode forced." )
#
#      parallel = TRUE
#    }
  }

  ## Defaults come from calling arguments on sfInitCluster.
  if( !is.null( parallel ) ) {
    setOption( "parallel", parallel )

    if( parallel ) {
      ## There is a slightly problem: as many users can use sfCluster without
      ## session-ID, the session number "XXXXXXXX" is not good enough.
      ## Problem: we need the filename on clusterinit so we cannot use cluster
      ## here.
      ## Win: USERNAME, *nix: LOGNAME
      ## LOGNAME/USER ist not set under Windows (tried Win Server 2003)
      uname <- ifelse( Sys.getenv( "LOGNAME" ) != "", Sys.getenv( "LOGNAME" ),
                                                      Sys.getenv( "USERNAME" ) )

      if( uname == "" )
        uname <- "___"

      ## Add R for RunSnowMode heterogenous mode.
      ## XXX Check R version and fill in correct version.
      setOption( "session", paste( sep="_",
                                   "XXXXXXXXR",
                                   uname,
                                   format( Sys.time(), "%H%M%S_%m%d%y" ) ) )

##      message( "Forced parallel. Using session: ", .sfOption$session, " \n" )
    }
    ## Sequential mode: reduce to one CPU.
    else {
      setOption( "nodes", 1 )

##      message( "Forced to sequential mode.\n" )
    }
  }

  ## If socket hosts are set, take them.
  if( !is.null( socketHosts ) || is.vector( socketHosts ) )
    setOption( "sockHosts", socketHosts )

  ## Type of the cluster ({SOCK|PVM|MPI|NWS} are allowed).
  if( !is.null( type ) ) {
    if( length( grep( "PVM|MPI|SOCK|NWS", type ) ) > 0 )
      setOption( "type", type )
    else {
      warning( paste( "Unknown cluster type:", type, "Allowed are: {PVM,MPI,SOCK,NWS}. Fallback to SOCKet." ) )
      setOption( "type", "SOCK" )
    }
  }
  ## Default value: socket cluster.
  else
    setOption( "type", "SOCK" )

  ## Global restore setting (for sfClusterApplySR).
  if( !is.null( restore ) )
    setOption( "RESTORE", restore )
  
  arguments <- commandArgs()

  ## Search for currently executed R-file (if there is any). Detected by
  ## argument followed to option "-f" ("R CMD BATCH" adds -f implicitely).
  ## Save filename for options (for save/restore)
  ## @todo Find a better way to detect R-file (is there any?)
  ## Last argument to be ignored (as no follow-up exists).
  if( length( arguments ) >= 2 ) {
    for( entry in seq( 1, length( arguments ) - 1 ) ) {
      if( !is.null( arguments[entry] ) && ( arguments[entry] == '-f' ) ) {
        ## Switch to next entry and check if this is valid.
        entry <- entry + 1;

        ## If yes, take it as filename.
        if( !is.null( arguments[entry] ) && ( arguments[entry] != "" ) ) {
          setOption( "CURRENT", arguments[entry] )
          break
        }
      }
    }
  }

  ## No R-file given: set to DEFAULT filename (always occurs in interactive
  ## mode).
  if( is.null( .sfOption$CURRENT ) )
    setOption( "CURRENT", "DEFAULT" )

  ## Go through all arguments from commandline.
  for( arg in arguments ) {
    ## Non sfCluster-like argument? Skip.
    ## (Only empty argument are '--parallel' and '--restoreSR')
    if( ( length( grep( "=", arg ) ) == 0 ) &&
        !( ( arg == "--parallel" ) || ( arg == "--restoreSR" ) || ( arg == "--restore" ) ) )
      next;

    ## Arguments in form "--name=value"
    args <- strsplit( arg, "=" )

    ## Marker for parallel execution.
    ## If parallel was set via function arguments, commandline is ignored.
    if( args[[1]][1] == "--parallel" ) {
      if( !is.null( args[[1]][2] ) && !is.na( as.numeric( args[[1]][2] ) ) )
        cmdParallel <- ifelse( ( as.numeric( args[[1]][2] ) > 0 ), TRUE, FALSE )
      ## --parallel is allowed to use without value (means: true).
      else
        cmdParallel <- TRUE

      ## Ask here, instead there will be a warning if used with commandline arg
      ## --parallel and sfInit( parallel=TRUE ).
      ## Rise warning if command arguments are overwritten by sfInit() arguments.
      if( is.null( parallel ) )
        setOption( "parallel", cmdParallel )
      else if( parallel != cmdParallel )
        warning( paste( "Commandline argument --parallel",
                        "overwritten with sfInit argument parallel=", parallel ) )
    }
    ## Marker for general restore (only used in sfClusterApplySR).
    ## Both --restoreSR/--restore are allowed.
    else if( ( args[[1]][1] == "--restoreSR" ) || ( args[[1]][1] == "--restore" ) ) {
      if( is.null( restore ) )
        setOption( "RESTORE", TRUE )
      else if( !restore )
        warning( "Commandline argument --parallel",
                 "overwritten with sfInit argument restore=TRUE" )
    }
    ## Marker for Session-ID.
    else if( args[[1]][1] == "--session" ) {
      ## Session-ID is allways 8 Chars long.
      ## Not anymore since sfCluster >=0.23
      if( !is.null( args[[1]][2] ) ) { ##&& ( nchar( args[[1]][2] ) == 8 ) ) {
        setOption( "session", args[[1]][2] )
      }
      else
        warning( paste( "Empty or irregular Session-ID: '", args[[1]][2], "'\n" ) )
    }
    ## Amount of CPUs (formerly called "nodes", kept for backward
    ## compatibility).
    ## If set via function arguments, commandline is ignored.
    else if( ( args[[1]][1] == "--nodes" ) || ( args[[1]][1] == "--cpus" ) ) {
      nodes <- try( as.numeric( args[[1]][2] ) )

      if( !is.null( nodes ) && !is.na( nodes ) ) {
        if( nodes > .sfOption$MAXNODES ) {
          stop( paste( "Too much CPUs allocated:", nodes, "Max.:",
                       .sfOption$MAXNODES,
                       "\n - Call sfSetMaxCPUs() before sfInit() if you need more." ) )
        }
        else
          nodes <- max( 1, nodes )

        ## Really set amount of CPUs? Rise overwrite warning if needed.
        if( is.null( cpus ) )
          setOption( "nodes", nodes )
        else if( cpus != nodes )
          warning( paste( "Commandline --cpus=", nodes,
                          " overwritten by sfInit() argument cpus=", cpus, sep="" ) )
      }
      else
        warning( paste( "Empty or irregular nodes amount: '", nodes, "'\n" ) )
    }
    ## Type of the network.
    else if( args[[1]][1] == "--type" ) {
      if( !is.null( args[[1]][2] ) && ( nchar( args[[1]][2] ) > 0 ) ) {
        if( length( grep( "PVM|MPI|SOCK|NWS", args[[1]][2] ) ) > 0 ) {
          if( is.null( type ) )
            setOption( "type", args[[1]][2] )
          else if( type != args[[1]][2] )
            warning( paste( "Commandline --type=", args[[1]][2],
                            " overwritten by sfInit() argument type=", type, sep="" ) )
        }
        else {
          warning( paste( "Unknown cluster type on commandline:", args[[1]][2],
                          "Allowed are: {PVM,MPI,SOCK,NWS}" ) )
        }
      }
      else
        warning( "No cluster-type is given as value for argument --type" )
    }
    ## Hosts for socket mode.
    ## Arguments come in format:
    ##   nodename:cpus  ->  On node X are Y cpus used.
    ##   nodename       ->  On node X one cpu is used.
    ## Any entries are comma seperated (no whitespace allowed!):
    ##  node1:3,node2,node3:2
    else if( args[[1]][1] == "--hosts" ) {
      if( !is.null( args[[1]][2] ) && ( nchar( args[[1]][2] ) > 0 ) ) {
        cmdHosts <- c()

        hosts = unlist( strsplit( args[[1]][2], "," ) )

        ## Examine single host
        for( host in hosts ) {
          info <- unlist( strsplit( host, ":" ) )

          ## No CPU amount given: assume 1.
          if( is.null( info[2] ) || is.na( info[2] ) )
            info[2] <- 1

          offset <- as.integer( info[2] )

          if( offset <= 0 )
            offset <- 1
          
          if( !is.numeric( offset ) )
            stop( paste( "NOT NUMERIC: '", offset, "'", sep="" ) )

          len <- length( cmdHosts ) + 1

          ## Insert Host n-times where n is amount of CPUs
          ## (required for snows argument format).
          cmdHosts[seq(len,len+offset-1)] <- rep( as.character( info[1] ), offset )
        }

        if( is.null( socketHosts ) )
          setOption( "sockHosts", cmdHosts )
        else if( paste( cmdHosts, collapse="" ) != paste( socketHosts, collapse="" ) ) {
          warning( paste( "Commandline --hosts=", args[[1]][2],
                          " overwritten by sfInit() argument hosts=", paste( socketHosts, collapse="," ),
                          sep="" ) )
        }
      }
      else
        warning( "No hosts are given as value for --hosts" )
    }
    ## Temporary directory: slave logs.
    else if( args[[1]][1] == "--tmpdir" ) {
      if( !is.null( args[[1]][2] ) && ( nchar( args[[1]][2] ) > 0 ) )
        setOption( "TMPDIR", args[[1]][2] )
      else
        warning( "No temporary directory given as value for --tmpdir" )
    }
    ## Restore directory: intermediate results are lawn here.
    else if( args[[1]][1] == "--restdir" ) {
      if( !is.null( args[[1]][2] ) && ( nchar( args[[1]][2] ) > 0 ) )
        setOption( "RESTDIR", args[[1]][2] )
      else
        warning( "No restore/result directory given as value for --restdir" )
    }
    ## Startup lock.
    ## Add 1.62:
    ## should only used from sfCluster => is the marker snowfall is started
    ## though sfCluster!
    else if( args[[1]][1] == "--lockfile" ) {
      if( !is.null( args[[1]][2] ) && ( nchar( args[[1]][2] ) > 0 ) )
        setOption( "LOCKFILE", args[[1]][2] )
      else
        warning( "No lockfile given as value for --lockfile" )
    }
    ## Unknown option
    ## Add 1.62
    ## Add 1.72: patch from Michael Siegel for Mac OS X, which sets --gui.
    else if( args[[1]][1] != "--gui" ) {
      warning( paste( "Unknown option on commandline:", args[[1]][1] ) )
    }
  }

  invisible( NULL )
}
