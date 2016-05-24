## Thomas Hoffmann
## Created:  10/15/2007
## Modified: 10/21/2007
##
## Pieces of this file were taken from my R package pbatR,
##  as I never want to have to create a GUI from scratch
##  again. Once you do one though, everyone wants one.
## The idea is to have a general creation routine,
##  that allows rapid creation of a gui interface, and
##  that is flexible enough to allow the user to use this
##  to create additional tcltk objects on afterward if desired
##  for an even more complicated interface if they
##  so desire.
## Most internal routines will be made public then, in that idea.

#library( tcltk )
#source( "parseHelp.R" )

## TODO:
## - any way to stop the windows focus stealing? NO
## - Group horizontal widgets onto one frame?
##   - ? or just put the output in a separate frame than the
##     other widgets?
## - Make R package, further debugging
##   - Create some examples -- could do the standard normal
##     sample size e.g.
##   - An example of loading in a file (to the multi-edit)
##     -- Could even create a function to display text/an editor

## User: Your code should throw errors if you are using the text output
##        via 'stop', 'warning' will only appear in the console
## -- they will be caught by the gui widget with an error message dialog!

## Things I can think of that you might want...
## - Text entry for values (the default - at it's base, it just does this)
## - Slider (the pbat power package) - for ranges
## - Input for filenames
## - Options box (e.g. TRUE/FALSE)
## - Selection box (selecting a group of individuals)
## - Command button (calls another function -- this allows nesting!?)

#SLIDER_LENGTH <- 500
#ENTRYWIDTH <- 40
#LISTHEIGHT <- 15
#LISTWIDTH <- 15
#TEXTWIDTH <- 65
#TEXTHEIGHT <- 5

#########################
## Namespace functions ##
#########################
guiEnv <- new.env()
## *** EXPORT *** (allows setting constants, getting at internal things)
guiSet <- function( x, value )
  assign( x, value, envir=guiEnv );
guiGet <- function( x, mode="any", ifnotfound=NA )
  get( x, envir=guiEnv, mode=mode, inherits=FALSE )
guiGetSafe <- function( x, ifnotfound=NA )
  mget( x, envir=guiEnv, ifnotfound=ifnotfound )[[1]]

###############
## Constants ##
###############
guiSet( "SLIDER_LENGTH", 500 )
guiSet( "ENTRY_WIDTH", 40 ) ## ?10
guiSet( "LIST_HEIGHT", 15 )
guiSet( "LIST_WIDTH", 15 )
guiSet( "EDIT_WIDTH", 65 )
guiSet( "EDIT_HEIGHT", 5 )

#getPowerV <- function( ... )
#  return( tclvalue( getPower( ... ) ) )
#getPowerN <- function( ... )
#  return( as.numeric(tclvalue( getPower( ... ) ) ) )
#getPowerT <- function( ... )
#  return( tclvalue( getPower( ... ) ) == TRUE )

clearEnvironment <- function() {
  ## not sure if this is possible... (don't know how to traverse, use 'rm')
}

loadTclTkOrDie <- function()
{
  if (!require("tcltk"))
    stop("You need to have loaded tcl/tk to use this function {e.g. library(tcltk)}.")
}

numericFix <- function( v ) {
  ## try to make it numeric
  vn <- NA
  tryCatch( {vn <- as.numeric(v)}, error=function(e){}, warning=function(e){} )
  if( !is.na(vn[1]) ) {
    return(vn)
  }

  ## it is a string
  if( length(v)==1 && !is.na(v) ) {
    if( v=="TRUE" ) return( 1 )  ## TRUE
    if( v=="FALSE" ) return( 0 ) ## FALSE
  }
  return(v)
}

elt <- function(list,str,die=TRUE) {
  wh <- which( names(list)==str )
  if( length(wh)!=1 ) {
    if( die ) {
      stop( "Element fetching does not work." )
    }else{
      return( NULL )
    }
  }
  return( list[[wh]] )
}
elti <- function( list, str ) {
  res <- which( names(list)==str )
  if( length(res) == 0 )
    return( length(list)+1 )
  return( res )
}

fixName <- function( str ) {
  for( ch in c('\\.', '\\_') ) {
    s <- unlist( strsplit( str, ch ) )
    str <- paste( s, collapse=" " )
  }
  return( str )
}

## Allow nesting!
## Returns an object to be passed to unnest when exiting
##  modal stage (and being destroyed)
nest <- function() {
  ## Backup:
  ##  argList, argListSelected, func,
  ##  object, guiObject, fargTypes, fargNames, farg, output
  return( list( argList=guiGetSafe("argList"),
                argListSelected=guiGetSafe("argListSelected"),
                func=guiGetSafe("func"),
                object=guiGetSafe("object"),
                guiObject=guiGetSafe("guiObject"),
                fargTypes=guiGetSafe("fargTypes"),
                fargNames=guiGetSafe("fargNames"),
                farg=guiGetSafe("farg"),
                output=guiGetSafe("output"),
                MAIN=guiGetSafe("MAIN")  ) )
}
unnest <- function( n ) {
  nNames <- names(n)
  for( i in 1:length(n) )
    guiSet( nNames[i], n[[i]] )

  #print( "NESTED FUNCTION" )
  #print( guiGet("func") )
}

## *** EXPORT ***
## Aside: mention to the user the power of getFromNamespace
##
## func - the function that should be called upon execution
##
## argType - list (unspecified is auto-detected)
##  (unspecified): auto-detect
##       - if in options, assumes is an option box
##       - if in selections, assumes is a selection box
##       - if an integer,
##  't': text entry
##  's': slider
##  'f': input for filenames
##  'o': options box (options are put in argOption)
##  'l': list box (lists are put in argList, which is 'set', and can be modified by user)
##  'c': command button
##  'm': multi-line text entry
##  'i': _ignore_ -- really not necessary anymore...
## argOption    - list of options vectors (names should be the same as args)
## argList      - list of strings for lists (names should be the same as args)
## argSlider    - list of slider ranges (names should be the same as args)
## argCommand   - list of functions to execute on command
## argEdit      - list of (width,height) both optional, NULL/NA/missing for default
## argFilter    - list of file filters (empty for all files)
##
## exec - name of string to use when user should press button to have them
##         call your function
##      - empty indicates it should not be drawn
##
## callback - name of function to handle callbacks, takes one parameter,
##            which is a string for the arg that was updated
##
## output - one of the above, 't', 's', 'm', or NULL; will try
##           to auto-choose this as well. If not 'm', then
##           an initial value will be set by running
##           the default parameters
##
## helps - 'auto' indicates it will try to load in the help from the package help,
##          if possible
##       - otherwise this can be a list of strings for help
## helpsFunc - optional name of the function (string) of where to take the help info from
##
## grid - whether to grid the objects or not (otherwise, just let the user do it)
## modal - lock input away from R
## nameFix - boolean, tries to fix names (replaces '_' & '.' with ' ').
##
## clearEnv - TRUE for the most part
##          - FALSE added, for the idea of nesting GUI objects!
gui <- function( func,
                 argOption=NULL, argFilename=NULL, argList=NULL, argSlider=NULL, argCommand=NULL, argEdit=NULL, argFilter=NULL,
                 argText=NULL, argType=NULL,
                 argGridOrder=1:length(formals(func)),
                  argGridSticky=rep("a",length(formals(func))),
                  argGridFrame=rep("f",length(formals(func))),
                 title=NULL,
                 exec="OK",
                 closeOnExec=is.null(output), cancelButton=TRUE,
                 callback=NULL,
                 output='m',
                 helps='auto', helpsFunc=NULL,
                 grid=TRUE, modal=NULL, nameFix=TRUE, getFix=TRUE,
                 verbose=FALSE ) {
  require( tcltk )

  if( verbose ) {
    print( "argGrid..." )
    print( argGridOrder )
    print( argGridSticky )
    print( argGridFrame )
  }
  #stop()

  ## Store the getFix (05/09/2008)
  guiSet( "GUIINTERNALS_getFix", getFix )

  ## 01/23/2009
  modalDefault <- FALSE
  if( is.null(modal) ) {
    modalDefault <- TRUE ## special message to developer
    modal <- TRUE ## and set the default to be true
  }

  ## Allow nesting
  if( modal )
    n <- nest() ## allow nesting

  ## parses the passed in function to widgets!
  farg <- formals( func )
  farg <- farg[ names(farg) != "..." ] ## (05/09/2008)
  fargNames <- names(farg)
  fargTypes <- rep('t',length(farg))

  #print( "farg" )
  #print( farg )
  #print( farg[[1]] )

  ## Fix up farg a little more... infernal boxplot.default... (05/09/2008)
  for( i in 1:length(farg) ) {
    #if( !is.element( class(farg[[i]]), c("numeric","logical","character") ) || is.na(farg[[i]] ) )
    #  farg[[i]] <- ""
    if( !is.element( class(farg[[i]]), c("numeric","logical","character","NULL") ) )
      farg[[i]] <- ""

    ## 01.23.2009
    #cat( str(farg[[i]]), "before\n" )
    ####try(  if( is.nan( farg[[i]] ) ) farg[[i]] <- "NaN",  silent=TRUE  ) ## A little confused between NA and NaN??
    ####try(  if( is.na( farg[[i]] ) ) farg[[i]] <- "NA",  silent=TRUE  )
    tryCatch(  if( is.nan( farg[[i]] ) ) farg[[i]] <- "NaN",  warning=function(e){},  error=function(e){}  ) ## A little confused between NA and NaN??
    tryCatch(  if( is.na( farg[[i]] ) ) farg[[i]] <- "NA",  warning=function(e){}, error=function(e){}  )
    #try(  if( is.null( farg[[i]] ) ) farg[[i]] <- "NULL" )  ## Not Necessary?
    #cat( str(farg[[i]]), "after\n" )
  }

  ## store in the text (in case user wants a different name)
  fargText <- fargNames
  textNames <- names(argText)
  for( i in 1:length(fargText) ) {
    if( is.element(fargText[i],textNames) ) {
      fargText[i] <- elt(argText,fargText[i])
    }else if( nameFix ) {
      fargText[i] <- fixName( fargText[i] )
    }
  }
  if( verbose ) {
    cat( "Arg text:\n" )
    print( fargText )
  }

  ## Get the function name
  call <- match.call(expand.dots=FALSE)
  funcName <- call[[match("func",names(call))]]
  funcName <- as.character( as.expression( funcName ) )
  if( is.null(title) ) title <- funcName

  ## Try to get help on that function
  if( is.character(helps) && helps=='auto' ) {
    helps <- NULL
    if( !is.null(helpsFunc) ) {
      try( helps <- parseHelp( helpsFunc ), silent=TRUE )
    }else{
      try( helps <- parseHelp( funcName ), silent=TRUE )
    }
    if( verbose ) {
      cat( "Help parsing:\n" )
      print( helps )
    }
  }
  if( verbose ) cat( "Finished help parsing...\n" )

  ## prereqs
  ## - need to run a few things - e.g. do the lengths of argType and function match?

  ## Clear the environment (doesn't really do anything right now)
  clearEnvironment()

  ## Set the lists (can be modified by user later...)
  guiSet( "argList", argList )
  ##KILL##, guiGet( "argList" )
  argListSelected <- list()
  guiSet( "argListSelected", argListSelected )

  guiSet( "func", func )

  if( !is.function(callback) )
    callback <- function(str){}; ## do nothing...
  ##guiSet( "callback", callback ) ## not needed?

  ############################
  ## auto-detect the types! ##
  ############################

  ## - get lists of names to search
  typeList <- names(argType)
  optionList <- names(argOption)
  filenameList <- names(argFilename)
  listList <- names(argList)
  sliderList <- names(argSlider)
  commandList <- names(argCommand)
  editList <- names(argEdit)
  filterList <- names(argFilter)

  ## - do the search
  gv <- function( i, default="" ) {
    if( i > length(farg) ) return( default )

    #print( farg[[i]] )
    if( class(farg[[i]]) != "name" ) {
      if( !is.null(farg[[i]]) && !is.na(farg[[i]]) )
        return( farg[[i]] )
    }
    return( default )
  }
  for( i in 1:length(farg) ) {
    if( verbose ) cat( "Determining widget type ", i, "\n" )
    curArg <- fargNames[i]

    if( is.element( curArg, typeList ) ) {
      ## The user specified a type
      fargTypes[i] <- elt(argType,curArg)
    }else if( is.element( curArg, optionList ) ) {
      fargTypes[i] <- 'o'
    }else if( is.element( curArg, sliderList ) ) {
      fargTypes[i] <- 's'
    }else if( is.element( curArg, listList ) ) {
      fargTypes[i] <- 'l'
    }else if( is.element( curArg, filenameList ) || is.element( curArg, filterList ) ) {
      fargTypes[i] <- 'f'
    }else if( is.element( curArg, commandList ) ) {
      fargTypes[i] <- 'c'
    }else if( is.element( curArg, editList ) ) {
      fargTypes[i] <- 'm'
    }else{
      fargTypes[i] <- 't'
    }
  }
  if( verbose ) {
    cat( "Widget types determined...\n" )
    print( fargTypes )
  }

  ## Add a command button on the end if exec isn't empty
  okButtonGridOrder <- NULL
  if( !is.null(exec) && nchar(exec)>0 ) {
    if( verbose ) cat( "Adding command button...\n" )

    pos <- length(fargTypes)+1
    fargTypes[pos] <- 'c'
    fargNames[pos] <- exec
    fargText[pos] <- exec ## "Execute ..."  (05/09/2008)
    if( is.null(argCommand) ) argCommand <- list()
    pos2 <- length(argCommand)+1
    argCommand[[pos2]] <- guiExec
    if( closeOnExec ) { ## 01.27.2009
      argCommand[[pos2]] <- function() {
        res <- guiExec()
        tkdestroy(guiGet("MAIN"))
        guiSet( "FGUI_INTERNAL_last_result", res ) ## 02.02.2009
        return(res)
      }
    }
    names(argCommand)[pos2] <- exec

    commandList <- names(argCommand)
    argGridOrder[pos] <- max(argGridOrder)+1
    argGridSticky[pos] <- "ns"
    argGridFrame[pos] <- "f"
    okButtonGridOrder <- argGridOrder[pos]
  }
  if( verbose ) cat( "Determined output format 1...\n" )

  ## New, add a cancel button
  guiSet( "GUIINTERNALS_cancelled", FALSE )
  if( cancelButton ) {
    if( verbose ) cat( "Adding cancel button...\n" )

    pos <- length(fargTypes)+1
    fargTypes[pos] <- 'c'
    fargNames[pos] <- "cancel"
    fargText[pos] <- "Cancel"
    if( is.null(argCommand) ) argCommand <- list()
    pos2 <- length(argCommand)+1
    argCommand[[pos2]] <-
      function() {
        cat( "Cancelled.\n" )
        guiSet( "GUIINTERNALS_cancelled", TRUE )
        tkdestroy(main)
      }
    names(argCommand)[pos2] <- "cancel"
    commandList <- names(argCommand)
    if( !is.null(okButtonGridOrder) ) {
      argGridOrder[pos] <- okButtonGridOrder  ## Right next to the OK button
      argGridSticky[pos] <- "ns"
      argGridFrame[pos] <- "f"
    }else{
      argGridOrder[pos] <- max(argGridOrder)+1
      argGridSticky[pos] <- "ns"
      argGridFrame[pos] <- "f"
    }
  }

  ## add in a status box
  #if( statusBox && (output!='m'&&output!='t') ) {
  #  pos <- length(fargTypes)+1
  #  fargTypes[pos] <- 't'
  #  fargNames[pos] <- 'status'
  #  if( is.null(argEdit) ) argEdit <- list()
  #  pos2 <- length(argEdit)+1
  #  argEdit[[pos2]] <- NA
  #  names(argEdit)[pos2] <- 'status'
  #
  #  editList <- names(argEdit)
  #  argGridOrder[pos] <- max(argGridOrder)+1
  #}

  ## add in the output
  ## - first auto-choose if possible
  curArg <- "output"
  if( is.element( curArg, sliderList ) ) {
    output <- 's'
  }else if( is.element( curArg, editList ) ) {
    output <- 'm'
  }
  if( verbose ) cat( "Determined output format 2...\n" )
  #else{
  #  fargTypes[i] <- 't'
  #}
  ## - then add it into the lists
  if( !is.null(output) && (output=='m' || output=='t' || output=='s') ) {
    pos <- length(fargTypes)+1
    fargNames[pos] <- 'output'
    fargText[pos] <- "Output:"
    if( is.element('output',textNames) )
      fargText[pos] <- elt(argText,'output')
    if( output=='m' ) {
      fargTypes[pos] <- 'm'
      if( is.null(argEdit) ) argEdit <- list()
      pos2 <- length(argEdit)+1
      argEdit[[pos2]] <- NA
      names(argEdit)[pos2] <- 'output'

      editList <- names(argEdit)
    }else if( output=='s' ) {
      fargTypes[pos] <- 's'
      if( is.null(argSlider) ) argSlider <- list()
      pos2 <- length(argSlider)+1
      argSlider[[pos2]] <- NA
      names(argSlider[[pos2]]) <- 'output'

      sliderList <- names(argSlider)
    }else if( output=='t' ) {
      fargTypes[pos] <- 't'
    }
    argGridOrder[pos] <- max(argGridOrder)+1
    argGridSticky[pos] <- "news"
    argGridFrame[pos] <- "g"
  }
  if( verbose ) cat( "Output added to lists...\n" )

  ##print( argCommand )

  ########################
  ## create the widgets ##
  ########################

  ## Create the main window
  require( tcltk )
  main <- tktoplevel()
  tkwm.title( main, title )
  if( verbose ) cat( "Created main window, proceeding to create all widgets...\n" )

  guiSet( "MAIN", main ) ## for cursor

  ## Create all the widgets
  ##for( i in 1:length(fargNames) ) {
  object <- list()
  guiObject <- list()
  i <- 1
  for( ago in unique(argGridOrder) ) {
    #print( "ago" )
    #print( ago )
    #print( "argGridOrder" )
    #print( argGridOrder )

    agos <- which( argGridOrder == ago )

    sframe <- main
    #print( "argGridFrame" )
    #print( argGridFrame )
    #print( "agos" )
    #print( agos )
    if( argGridFrame[agos][1] == "f" ) {
      sframe <- guiFrame( main, borderwidth=0 )
      tkgrid( sframe )
      ##sticky <- "news" ## codeTools -- 02/11/2009
      if( argGridSticky[agos][1] != 'a' ) {
        tkgrid.configure( sframe, sticky=argGridSticky[agos][1] )
      }else if( fargTypes[agos][1] == 't' || fargTypes[agos][1] == 'f' ) {
        tkgrid.configure( sframe, sticky="nes" )
      }else{
        tkgrid.configure( sframe, sticky="nws" )
      }
      #if( fargTypes[agos][1] == 't' ) {
      #  if( argGridSticky[agos][1] == 'a' ) {
      #    tkgrid.configure( sframe, sticky="nes" ) ## default placement hack
      #  }else{
      #    tkgrid.configure( sframe, sticky=argGridSticky[agos][1] )
      #  }
      #}else{
      #  ## default (essentially nws, really)
      #  tkgrid.configure( sframe, sticky="news" )
      #}
    }

    ## Create an object for each
    for( a in agos ) {
      if( verbose ) cat( "Creating widget ", i, "\n" )
      helpsi <- elt(helps,fargNames[i],die=FALSE)
      res <- NULL

      if( fargTypes[i] == 't' ) {
        ## Text entry
        #var <-  tclVar(gv(i))
        #res <- guiTextEntry( var, fargText[i], main, helps=helpsi )
        res <- guiTextEntry( sframe=sframe, text=fargText[i], default=gv(i), helps=helpsi )
        if( argGridSticky[i] == 'a' ) {
          argGridSticky[i] <- "nes"
          ##print( "WTF?" )
        }
      }else if( fargTypes[i] == 's' ) {
        ## Slider
        range <- elt( argSlider, fargNames[i] )
        if( length(range) < 2 )
          stop( "Slider ranges must be a vector of two numeric values, with an optional 3rd argument stepsize." )
        if( length(range) == 2 )
          range[3] <- (range[2]-range[1])/100
        value <- gv( i, range[1] )
        res <- guiSlider( text=fargText[i], default=value, min=range[1], max=range[2], step=range[3], sframe=sframe, update=guiExecUpdateFunc(fargNames[i],callback), helps=helpsi )
        if( argGridSticky[i] == 'a' )
          argGridSticky[i] <- "nws"
      }else if( fargTypes[i] == 'f' ) {
        ## Input filenames
        value <- gv( i )
        #callback <- NULL
        filter <- elt( argFilter, fargNames[i], die=FALSE )
        if( is.null(filter) ) filter <- "{{All files} {.*}}"
        if( argGridSticky[i] == 'a' )
          argGridSticky[i] <- "nws"
        res <- guiFilename( sframe, text=fargText[i], default=value, filter=filter, callback=guiExecUpdateFunc(fargNames[i],callback), helps=helpsi )
      }else if( fargTypes[i] == 'o' ) {
        ## options box
        res <- guiOption( sframe, fargText[i], elt(argOption,fargNames[i]), update=guiExecUpdateFunc(fargNames[i],callback), helps=helpsi )
        if( argGridSticky[i] == 'a' )
          argGridSticky[i] <- "nws"
      }else if( fargTypes[i] == 'l' ) {
        ## list box
        res <- guiList( sframe, fargText[i], fargNames[i], update=guiExecUpdateFunc(fargNames[i],callback), helps=helpsi )
        if( argGridSticky[i] == 'a' )
          argGridSticky[i] <- "nws"
      }else if( fargTypes[i] == 'c' ) {
        ## command button
        res <- list()
        res$object <- "no object"
        #res$guiObject <- tkbutton( main, text=fargNames[i], command=elt(argCommand,fargNames[i]) )
        res$guiObject <- tkbutton( sframe, text=fargText[i], command=guiExecUpdateFunc(fargNames[i],callback,elt(argCommand,fargNames[i],die=FALSE)) )
        if( argGridSticky[i] == 'a' )
          argGridSticky[i] <- "nws"
      }else if( fargTypes[i] == 'm' ) {
        ## multi-line edit
        values <- elt(argEdit,fargNames[i])
        #print( values )
        width <- guiGet("EDIT_WIDTH")
        height <- guiGet("EDIT_HEIGHT")
        readonly <- TRUE
        if( length(values) > 0 && !is.na(values[1]) ) width <- values[1]
        if( length(values) > 1 && !is.na(values[2]) ) height <- values[2]
        if( length(values) > 2 && !is.na(values[3]) ) readonly <- (values[3]==TRUE)
        res <- guiEdit( sframe, fargText[i], gv(i), width, height, readonly, helps=helpsi )
        if( argGridSticky[i] == 'a' )
          argGridSticky[i] <- "news"
      }else if( fargTypes[i] == 'i' ) {
        res <- list(object=NULL,guiObject=NULL)
      }

      #i <- length(object) + 1
      object[[i]] <- res$object
      guiObject[[i]] <- res$guiObject
      i <- i + 1
    }## a

    ## Now grid the objects
    ##argGridSticky[agos] <- "nws"
    if( verbose ) cat( "Gridding widgets ", paste(agos,collapse="-"), "\n", sep="" )
    ##guiGrid( guiObject, sticky=argGridSticky[agos] )
    ###guiGrid( subset( guiObject, argGridOrder==ago ), sticky="nws" )
    guiGrid( subset( guiObject, argGridOrder==ago ), sticky=argGridSticky[agos] )

  }## unique(argGridOrder)
  ##}

  ## Draw all the objects
  #if( grid ) {
  #  for( ago in unique(argGridOrder) ) {
  #    if( verbose ) cat( "Gridding widget",ago,"\n" )
  #    guiGrid( subset( guiObject, argGridOrder==ago ), sticky=sticky )
  #  }
  #}

  ## Set some things in the globs, so user can do it, etc.
  guiSet( "object", object )
  guiSet( "guiObject", guiObject )
  guiSet( "fargTypes", fargTypes )
  guiSet( "fargNames", fargNames )
  guiSet( "farg", farg )
  guiSet( "output", output )

  ## Update output if necessary
  if( !is.null(output) && output!='m' ) ##&& !closeOnExec ) ## Not sure if I want this last piece
    guiExec()

  ## lastly, should we go modal?
  if( modalDefault ) {
    ## Notify developer it can be set
    cat( "Close '", title, "' window to allow entering commands in the R console; that window has gone modal. Note that the '", title, "' window may be hidden from view (esp. in windows), and you may have to find it in the taskbar and close it. To the developer, the modality of a window can be set.\n" )
    tkwait.window(main)
  }else if( modal ) {
    cat( "Close '", title, "' window to allow entering commands in the R console; that window has gone modal. Note that the '", title, "' window may be hidden from view (esp. in windows), and you may have to find it in the taskbar and close it.\n" )
    #cat( "Close '", title, "' window to return to R (sometimes R will steal the focus, esp. in windows).\n", sep="" )
    tkwait.window(main)
  }else{
    ## Not modal
    cat( "The '", title, "' window has been launched. Note that the '", title, "' window may be hidden from view (esp. in windows) and you may need to click on it in the taskbar.\n" )
  }

  ## If cancelled, don't return a value
  if( guiGet( "GUIINTERNALS_cancelled" ) == TRUE )
    return( NULL )

  ## And return all the values...
  allValues <- guiGetAllValues()
  if( modal ) unnest(n) ## RECOVER FROM NESTING
  return( allValues )
}

## *** EXPORT ***
## getting values (of objects) from the gui
guiGetValue <- function( i ) {
  #cat( "guiGetValue", i, "\n" )
  object <- guiGet( "object" )
  fargTypes <- guiGet( "fargTypes" )
  fargNames <- guiGet( "fargNames" )

  #cat( "fargType =", fargTypes[i] )

  ## handle if it is a string to find
  if( is.character(i) ) {
    wh <- which( i==fargNames )
    if( length(wh)==0 )
      warning( paste( "guiGetValue(), '",i,"' not found.", sep="" ) )
    return( guiGetValue(wh[1]) )
  }

  ## Do the real value finding
  value <- NA
  if( fargTypes[i] == 't' ) {
    ## Text entry
    value <- tclvalue(object[[i]])
  }else if( fargTypes[i] == 's' ) {
    ## Slider
    value <- tclvalue(object[[i]])
  }else if( fargTypes[i] == 'f' ) {
    ## Input filenames
    value <- tclvalue(object[[i]])
  }else if( fargTypes[i] == 'o' ) {
    ## options box
    value <- tclvalue(object[[i]])
  }else if( fargTypes[i] == 'l' ) {
    ## list box
    value <- getSelectedListElements(fargNames[i])
  }else if( fargTypes[i] == 'c' ) {
    ## command button
    value <- guiGetSafe( fargNames[i] ) ## a _glorious_ hack
  }else if( fargTypes[i] == 'm' ) {
    ## multi-line edit
    value <- "Cannot get edit values on closing."
    try( value <- tclvalue(tkget(object[[i]],"0.0","end")), silent=TRUE )
  }else if( fargTypes[i] == 'i' ) {
    ## Use the default value
    farg <- guiGet("farg")
    value <- farg[[i]]
  }

  ## really want to convert to numeric only _if_ numeric
  if( length(value) == 0 ) value <- NA
  value <- numericFix( value )

  return( value )
}
guiSetValue <- function( i, value ) {
  object <- guiGet( "object" )
  fargTypes <- guiGet( "fargTypes" )
  fargNames <- guiGet( "fargNames" )

  ## handle if it is a string to find
  if( is.character(i) ) {
    wh <- which( i==fargNames )
    if( length(wh)==0 )
      warning( paste( "guiGetValue(), '",i,"' not found.", sep="" ) )
    i <- wh[1]
  }

  ## Do the real value finding
  if( fargTypes[i] == 't' ) {
    ## Text entry
    tclvalue(object[[i]]) <- value
  }else if( fargTypes[i] == 's' ) {
    ## Slider
    tclvalue(object[[i]]) <- value
  }else if( fargTypes[i] == 'f' ) {
    ## Input filenames
    tclvalue(object[[i]]) <- value
  }else if( fargTypes[i] == 'o' ) {
    ## options box
    tclvalue(object[[i]]) <- value
  }else if( fargTypes[i] == 'l' ) {
    ## list box
    stop( "Cannot set 'a' single value in a list, try instead setListElements(...) instead." )
  }else if( fargTypes[i] == 'c' ) {
    ## command button
    stop( "Cannot set a value for a command." )
  }else if( fargTypes[i] == 'm' ) {
    ## multi-line edit
    ##tclvalue(tkget(object[[i]],"0.0","end")) <- value
    stop( "Cannot set a multi-line edit currently." )
  }else if( fargTypes[i] == 'i' ) {
    ## Use the default value
    stop( "Cannot set a value for an invisible component!" )
  }
}

## *** EXPORT ***
guiGetAllValues <- function() {
  farg <- guiGet( "farg" )
  #print( "farg guiGetAllValues" )
  #print( farg )

  value <- list()
  for( i in 1:length(farg) ) {
    value[[i]] <- guiGetValue(i)
    #if( is.na(value[[i]]) ) value[[i]] <- NULL
  }
  names(value) <- names( farg )

  return(value)
}

## Handle callback to the function
guiExec <- function( lastTouched=NULL ) {
  cancelled <- guiGetSafe("GUIINTERNALS_cancelled")
  if( !is.na(cancelled) && !is.null(cancelled) && cancelled==TRUE )
    return(NULL)

  #print( "callback to gui function." );

  object <- guiGet( "object" )
  ## CODETOOLS ## fargTypes <- guiGet( "fargTypes" )
  ## CODETOOLS ## fargNames <- guiGet( "fargNames" )
  ## CODETOOLS ## farg <- guiGet( "farg" )
  func <- guiGet( "func" )
  output <- guiGet( "output" )
  main <- guiGet( "MAIN" )

  ## first need to parse out all of the options...
  ##print( "about to get all the options" )
  value <- guiGetAllValues()
  #print( "got the options" )
  ## debug for now..
  #print( value )

  ## Set the cursor to business
  try( tkconfigure( main, cursor="watch" ) )

  ## RUN THE FUNCTION!
  #if( is.list(value) ) {
  #  for( i in 1:length(value) )
  #    if( !is.na(value) )
  #
  #}
  if( !is.list(value) ) {
    value2 <- list;
    for( i in 1:length(value) )
      value2[[i]] <- value[i]
    value <- value2
  }else{
    #print( "is list" )
  }
  ## Below really should work just fine
  ## But some funniness with R deciding lists
  ## Are just vectors, to piss us off
  #formals(func) <- value
  ## (05/09/2008)

  namesValue <- names(value)
  namesFormalsFunc <- names(formals(func))
  for( v in 1:length(value) ) {                      ## v indices value
    f <- which( namesFormalsFunc == namesValue[v] )  ## f indices formals(func)

    ##cat( "names", namesValue[[v]], namesFormalsFunc[[f]], "\n" )

    ## precaution
    if( names(value)[1] == "..." ) {
      stop("guiExec: '...' should have been caught earlier!")
    }

    if( is.character(value[[v]]) && nchar(value[[v]])==0 ) {
      ## Then it's an empty character, don't modify anything
    }else{
      formals(func)[[f]] <- value[[v]]
      if( guiGet("GUIINTERNALS_getFix") == TRUE )
        try( formals(func)[[f]] <- eval(parse(text=value[[v]])), silent=TRUE )
    }
  }

  #cat("JUST BEFORE EVALUATING...\n")

  res <- NULL;
  tryCatch( {res <- func()},
            error=function(e){gui_errorMessage(e$message)} )

  #cat("JUST AFTER EVALUATING...\n")

  ## Unset the business
  try( tkconfigure( main, cursor="arrow") )

  if( !is.null(res) ) {
    if( is.null(output) || ( output!='m' && output!='t' && output!='s' ) ) { ## 01/27/2009
      return(res) ## If no output window, then return the value
    }else if( output=='m' ) {
      gui_tkInsertText( object[[length(object)]], as.character(res) )    ## 05/09/2008 -- as.character
      gui_tkInsertText( object[[length(object)]], "\n" )
    }else if( output=='t' ) {
      ##gui_tkSetText( object[[length(object)]], res )
      tclvalue(object[[length(object)]]) <- as.character(res)  ## 05/09/2008 -- as.character
    }else if( output=='s' ) {
      tclvalue(object[[length(object)]]) <- as.character(res)  ## 05/09/2008
    }

    ## 05/16/2008 addition!
    winExists <- guiGetSafe( "INTERNALMENU_USED" )
    if( is.na( winExists ) || !winExists ) {
      ## oops --- nothing really! had it backwards
    }else{
      fguiWindowPrint( as.character(res) )
    }
  }
}
## Returns a function to handle executing the
##  update (for immediate, NULL if shouldn't
##  be immediate), and sets what was last touched
##  beforehand, in case this is needed (would have
##  been needed in my more complicated power interface)
guiExecUpdateFunc <- function(argLastTouched, callback, callback2=NULL) {
  argLastTouched
  callback2
  return(
    function() {
      callback(argLastTouched)
      if( is.function( callback2 ) )
        callback2()
    } );
}

## *** EXPORT *** (optional)
## Handle gridding
## - guiObject is a list of objects to be gridded together
guiGrid <- function( guiObject, sticky="nws" ) {
  n <- length(guiObject)
  if( n==1 ) {
    if( !is.null(guiObject[[1]]) )
      tkgrid( guiObject[[1]] )  ## columnspan=2
  }else if( n==2 ) {
    tkgrid( guiObject[[1]], guiObject[[2]]  )
  }else if( n==3 ) {
    tkgrid( guiObject[[1]], guiObject[[2]], guiObject[[3]] )
  }else if( n==4 ) {
    tkgrid( guiObject[[1]], guiObject[[2]], guiObject[[3]], guiObject[[4]] )
  }else if( n==5 ) {
    tkgrid( guiObject[[1]], guiObject[[2]], guiObject[[3]], guiObject[[4]], guiObject[[5]] )
  }else if( n==6 ) {
    tkgrid( guiObject[[1]], guiObject[[2]], guiObject[[3]], guiObject[[4]], guiObject[[5]], guiObject[[6]] )
  }

  ## Allow sticky to be a vector now for different stickinesses
  #print( "guiGrid sticky" )
  #print( sticky )
  if( length(sticky) != n )
    sticky <- rep( sticky[1], n )

  for( i in 1:n )
    if( !is.null(guiObject[[i]]) )
      tkgrid.configure( guiObject[[i]], sticky=sticky[i] )
}

gui_tkClearText <- function( obj ) {
  loadTclTkOrDie()  ## has to be here to pass the check

  tkdelete( obj, 0, 9999999 );
  #tkdelete( obj, '@0,0', 'end' )
}
gui_tkClearTextM <- function( obj ) {
  loadTclTkOrDie()  ## has to be here to pass the check

  tkdelete( obj, '@0,0', 'end' )
}
gui_tkSetText <- function( obj, text, READONLY=TRUE ) {
  loadTclTkOrDie()  ## has to be here to pass the check

  if( READONLY ) try( tkconfigure( obj, state="normal" ) );  ## try
  gui_tkClearText( obj );
  tkinsert( obj, "end", text );
  if( READONLY ) try( tkconfigure( obj, state="readonly" ) );  ## try
}
gui_tkInsertText <- function( obj, text ) {
  tkinsert( obj, "end", text )
}
gui_populateList <- function( lst, items ) {
  loadTclTkOrDie()  ## has to be here to pass the check

  if( length(items) > 0 ) { ## semi-bug fix 01/20/2006
    for( i in 1:length(items) )
      tkinsert( lst, "end", items[i] );
  }
}
gui_errorMessage <- function( message ) {
  loadTclTkOrDie()  ## has to be here to pass the check

  tkmessageBox( title="ERROR",
                message=message,
                icon="error", type="ok" );
}

## All the GUI routines return a list of
## - the object
## - the guiObject to be gridded via tkgrid (for layouts)

## the frame (basic requirement of most)

## *** EXPORT ***
guiFrame <- function( sframe, grid=FALSE, relief="groove", borderwidth=2, sticky="nws" ) {
  frame <- tkframe( sframe, relief=relief, borderwidth=borderwidth );
  if( !grid ) return(frame);
  tkgrid( frame );
  tkgrid.configure( frame, sticky=sticky );
  return(frame)
}

## NEW: A help button, if the user wants to know more about the argument.
helpButton <- function( sframe, helps, title ) {
  helps  ## suspect will be necessary for scoping
  if( !is.null(helps) && !is.na(helps) && nchar(helps)>0 ) {
    but <- tkbutton( sframe, text="?",
                     command=function(){ tkmessageBox( title=title, message=helps, icon="info", type="ok" ) }
                     )
    return( but )
  }
  return( NULL )
}

## *** EXPORT ***
## Text entry
guiTextEntry <- function( sframe, text, default, width=NULL, helps=NULL ) {
  tVar <- tclVar(default)

  if( is.null(width) ) width <- guiGet("ENTRY_WIDTH")
  #frame <- guiFrame( borderwidth=3, sframe=sframe );
  frame <- guiFrame( borderwidth=0, sframe=sframe, relief="flat" );
  tkgrid.configure( frame, sticky="news" )
  lab <- tklabel( frame, text=text );
  entry <- tkentry( frame, width=width, textvariable=tVar );
  butHelp <- helpButton( frame, helps, text )
  tkgrid( lab, entry, butHelp )
  #tkgrid( entry, butHelp, lab )
  tkgrid.configure( lab, sticky="nes" )
  tkgrid.configure( entry, sticky="nes" )
  #tkgrid.configure( butHelp, sticky="nes" )
  #tkgrid.configure( lab, sticky="nes" )

  return( list( object=tVar, guiObject=frame ) )
}


## *** EXPORT ***
## Creates a slider, etc.
guiSlider <- function( sframe, text, default, min, max, step=(max-min)/100, update=NULL, state="enabled", helps=NULL )
{
  frame <- guiFrame( sframe )
  sliderVal <- tclVar( default )
  sliderValLabel <- tkentry( frame, width=5, textvariable=sliderVal )
  slider <- tkscale( frame, from=min, to=max, resolution=step,
                     showvalue=FALSE, variable=sliderVal, orient="horizontal",
                     length=guiGet("SLIDER_LENGTH") )
  butHelp <- helpButton( sframe, helps, text )
  #tkgrid( slider,
  #        tklabel(frame, text=text),
  #        sliderValLabel,
  #        butHelp )
  tkgrid( slider, sliderValLabel, tklabel(frame,text=text), butHelp ) ## 02/01/2009
  tkconfigure( sliderValLabel,textvariable=sliderVal )
  if( state=="disabled" ) {
    try( {
      tkconfigure( slider, state="disabled" )
      tkconfigure( sliderValLabel, state="disbled" )
    } )
  }

  ## bind the update firing
  if( !is.null(update) ) {
    tkbind(sliderValLabel, "<Return>", update)
    tkbind(slider, "<ButtonRelease-1>", update)
  }

  #tkgrid.configure( slider, sticky="news" )

  #return( sliderVal )  ## return the tclvariable corresponding to the new object

  return( list( object=sliderVal, guiObject=frame ) )
}

## *** EXPORT ***
## e.g. "{{Un/compressed pedigree files} {.ped .pped}}"
##                                      ^-- note this space here is required!
guiFilename <- function( sframe, text="Filename ...", default="", title="", filter="{{All files} {.*}}", callback=NULL, helps=NULL ) {
  frame <- tkframe( sframe )

  if( is.null(default) || is.na(default) ) default <- ""
  var <- tclVar( default )
  te <- tkentry( frame, width=guiGet("ENTRY_WIDTH"), textvariable=var )

  callback ## needs to be here to scope it correctly... odd...
  filter
  cmd <- function(){
    #loadTclTkOrDie()  ## has to be here to pass the check

    # get the globals
    #globs <- getPbatGUI( "globs" );

    # See if we can get the filename, return if cancelled
    tempstr <- tclvalue(tkgetOpenFile(title=title,filetypes=filter))
    if( !nchar(tempstr) ) return()

    gui_tkSetText( te, tempstr )

    ## And a callback
    if( !is.null(callback) ) callback()
  }

  but <- tkbutton( frame, text=text, command=cmd )
  try( tkconfigure( te, state="readonly" ) ) ## try, told this fails on Red Hat old school
  butHelp <- helpButton( frame, helps, title )
  tkgrid( but, te, butHelp )
  tkgrid.configure( but, sticky="ew" )

  return( list( object=var, guiObject=frame ) )
}

## *** EXPORT ***
## Options
guiOption <- function( sframe, text, choices, defaultChoice=1, update=NULL, helps=NULL ){
  ## Create a frame for everything
  frame <- guiFrame( sframe, borderwidth=1, relief="ridge" ); ## groove, sunken, flat, ridge
  ######tkgrid( frame );
  ##tkgrid.configure( frame, sticky="news" );
  ######tkgrid.configure( frame, sticky="nws" );

  ## The label
  lbl <- tklabel( frame, text=paste(text,": ",sep="") );

  ## The choices
  rb <- list();
  rbLbl <- list();
  rbSubframe <- list();

  ## Create the subframe for each choice
  for( i in 1:length(choices) )
    rbSubframe[[i]] <- tkframe( frame, relief='groove', borderwidth=0 ); #1 );

  ## help button
  butHelp <- helpButton( frame, helps, text )

  ## grid the subframes
  if( length(choices)==2 ){
    tkgrid( lbl, rbSubframe[[1]], rbSubframe[[2]], butHelp );
  }else if( length(choices)==3 ){
    tkgrid( lbl, rbSubframe[[1]], rbSubframe[[2]], rbSubframe[[3]], butHelp );
  }else if( length(choices)==4 ){
    tkgrid( lbl, rbSubframe[[1]], rbSubframe[[2]], rbSubframe[[3]], rbSubframe[[4]], butHelp );
  }else if( length(choices)==5 ){
    tkgrid( lbl, rbSubframe[[1]], rbSubframe[[2]], rbSubframe[[3]], rbSubframe[[4]], rbSubframe[[5]], butHelp );
  }else if( length(choices)==6 ){
    tkgrid( lbl, rbSubframe[[1]], rbSubframe[[2]], rbSubframe[[3]], rbSubframe[[4]], rbSubframe[[5]], rbSubframe[[6]], butHelp );
  }else{
    stop( "Number of options cannot be greater than 6, try using a list instead." )
  }

  ## Create the variable
  rbVal <- tclVar( choices[defaultChoice] );

  ## Create each choice
  for( i in 1:length(choices) ) {
    rb[[i]] <- tkradiobutton( rbSubframe[[i]] );
    tkconfigure( rb[[i]], variable=rbVal, value=choices[i]  );
    rbLbl[[i]] <- tklabel( rbSubframe[[i]], text=choices[i] );
    tkgrid( rb[[i]], rbLbl[[i]] );

    ## how about binding the fire event?
    if( !is.null( update ) )
      tkbind( rb[[i]], "<ButtonRelease-1>", update )
  }

  ## Return the variable
  #return( rbVal );
  return( list( object=rbVal, guiObject=frame ) )
}

## *** EXPORT ***
## Creates a button that launches the appropriate listForm function
guiList <- function( sframe, text, name=text, update=NULL, helps=NULL ) {
  text ## Really, really, unbelievably strange, and rather messed up, but this _NEEDS_ to be here! fine then!
  name
  update
  helps
  cmd <- function() {
    listForm( name, helps )
    if( is.function(update) )
      update()
  }

  button <- tkbutton( sframe, text=text, command=cmd )
  return( list( object="no object", guiObject=button ) )
}

listForm <- function( name, helps=NULL ){
  loadTclTkOrDie()  ## has to be here to pass the check

  ## get the globals
  argList <- guiGet( "argList" )
  list <- elt( argList, name )
  argListSelected <- guiGet( "argListSelected" )
  listSelected <- elt( argListSelected, name, die=FALSE )

  ## create a modal dialog
  form <- tktoplevel()
  tkwm.deiconify(form)
  tkgrab.set(form) # make it modal
  tkfocus(form)
  tkwm.title( form, name )

  ## Draw a list box for the phenos stuff ( no order here... )
  scr <- tkscrollbar( form, repeatinterval=5,
                      command=function(...)tkyview(lst,...) );
  lst <- tklistbox( form, height=guiGet("LIST_HEIGHT"), width=guiGet("LIST_WIDTH"),
                    selectmode="multiple", background="white",
                    yscrollcommand=function(...)tkset(scr,...) );
  tkgrid( lst, scr );
  tkgrid.configure( scr, sticky="ns" );

  ## populate the list
  gui_populateList( lst, list );

  ## And select some of them which might have been selected
  if( length(listSelected) > 0 ) {
    ## and select them!
    for( i in 1:length(listSelected) )
      tkselection.set( lst, listSelected[i] );
  }

  on.exit <- function() {
    ## get the globals  -- is this necessary?
    ## CODETOOLS ## argList <- guiGet( "argList" )
    ## CODETOOLS ## list <- elt( argList, name )
    argListSelected <- guiGet( "argListSelected" )
    listSelected <- elt( argListSelected, name, die=FALSE )

    ## Need to translate any selected phenotypes to globs$phenos
    listSelected <- as.numeric(tkcurselection(lst))
    if( length(listSelected) > 0 ) {
      pos <- elti(argListSelected,name)
      argListSelected[[ pos ]] <- listSelected
      names(argListSelected)[pos] <- name
    }else{
      ## need to clear
      pos <- elti(argListSelected,name)
      argListSelected[[ pos ]] <- NULL
    }

    ## set the global variables
    guiSet( "argListSelected", argListSelected );
  }

  ## Bind mouse presses to on.exit(), since we can't capture closing the form!
  tkbind( lst, "<ButtonRelease>", on.exit );

  ## lastly, an OK button
  cmdOK <- function() {
    on.exit();
    tkdestroy(form);
  }
  but.ok <- tkbutton( form, text="   OK   ", command=cmdOK );
  butHelp <- helpButton( form, helps, name )
  tkgrid( but.ok, butHelp )

  ## Now make the form go modal
  tkfocus( form );
  tkwait.window( form ); # this makes it go modal
  ##on.exit(); # run after going modal --> doesn't work anymore!
}
## *** EXPORT ***
getSelectedListElements <- function( name ) {
  argList <- guiGet( "argList" )
  list <- elt( argList, name )

  argListSelected <- guiGet( "argListSelected" )
  listSelected <- elt( argListSelected, name, die=FALSE )

  if( length(listSelected)>0 )
    return( list[listSelected+1] )
  return( NULL )
}
## *** EXPORT ***
setListElements <- function( name, elements ) {
  argList <- guiGet( "argList" )
  pos <- which( names(argList) == name )
  if( length(pos)!=1 ) stop( paste("Cannot set list elements, list'",name,"' does not exist.",sep="") )
  argList[[pos]] <- elements
  guiSet( "argList", argList )
  return(invisible())
}

## *** EXPORT ***
guiEdit <- function( sframe, text="", default="", width=NULL, height=NULL, readonly=FALSE, helps=NULL ) {
  if( is.null(width) ) width <- guiGet("EDIT_WIDTH" )
  if( is.null(height) ) height <- guiGet("EDIT_HEIGHT" )

  frame <- tkframe( sframe, relief="groove", borderwidth=2 )
  #tkgrid( frame )

  text
  lbl <- NULL
  if( text != "" ) lbl <- tklabel( frame, text=text )
  butHelp <- helpButton( frame, helps, text )
  if( !is.null(lbl) ) {
    tkgrid( lbl, butHelp )
    if( !is.null(butHelp) ) tkgrid.configure( butHelp, sticky="sw" )
    tkgrid.configure( lbl, sticky="sw" )
  }else if( !is.null(butHelp) ) {
    tkgrid( butHelp )
    tkgrid.configure( butHelp, sticky="sw" )
  }

  yscr <- tkscrollbar( frame, repeatinterval=5,
                       command=function(...)tkyview(statusText,...))
  statusText <- tktext( frame, height=height, width=width,
                        yscrollcommand=function(...)tkset(yscr,...) )
  tkgrid(statusText, yscr)
  tkgrid.configure(yscr, sticky="nse")
  tkgrid.configure(statusText, sticky="nse" )
  #if( readonly==TRUE ) try( tkconfigure( statusText, state="disabled" ) );  ## try
  #tkinsert( statusText, "end", "Status of run (last line, enter 'pbat.status(n=0)' from the command prompt for more):\n" )

  return( list( object=statusText, guiObject=frame ) );
}

## *** EXPORT ***
## nestArg - the name of the argument in the other function
guiNestedF <- function( func, nestArg, title=nestArg, exec=NULL, output=NULL, ... ) {
  func
  nestArg
  return( function(){
    nestArgVal <- guiGetSafe(nestArg)
    if( length(nestArgVal)>1 || !is.na(nestArgVal) ) formals(func) <- nestArgVal

    res <- gui( func, title=title, exec=exec, output=output, ... )
    guiSet( nestArg, res )
    return( res )
    #print( guiGet( nestArg ) )
    #return( guiGet( nestArg ) )
  } )
}

## Because R is being insane with transforming lists into vectors!!!! why?????
guiFormals <- function( func, object ) {
  if( is.list(object) ) {
    for( i in 1:length(object) )
      if( !is.null(object[[i]]) && !is.na(object[[i]]) )
        formals(func)[[i]] <- object[[i]]
  }else{
    for( i in 1:length(object) )
      if( !is.null(object[i]) && !is.na(object[i]) )
        formals(func)[[i]] <- object[i]
  }
  return( func )
}

## 01.23.2009 -- Additions for the reviewers

## Helper function for guiv, to return the evaluation of the last function
guiv_last <- function() {
  return( guiGetSafe( "FGUI_INTERNAL_last_result" ) ) ## so handles 'rivers'

  #func <- guiGetSafe( "FGUI_INTERNAL_guiv_function", ifnotfound=0 )
  #if( is.numeric(func) )
  #  stop("No output available from the last gui function ran." )
  #
  #formals( func ) <- guiGetAllValues()
  #return( func() )
}
## Instead of returning the list of parameters to the function, this evaluates the function at those parameters and returns those values
guiv <- function( func=NULL, output=NULL, modal=TRUE, title=NULL, ... ) {
  if( class(func) != "function" )
    return( guiv_last() )

  ## hack -- we need the function name!
  call <- match.call(expand.dots=FALSE)
  funcName <- call[[match("func",names(call))]]
  funcName <- as.character( as.expression( funcName ) )
  print( funcName )
  if( is.null(title) ) title <- funcName

  ##ret <- gui( func=func, modal=modal, output=output, title=title, ... ) ## codetools 02/11/09
  gui( func=func, modal=modal, output=output, title=title, ... ) ## codetools 02/11/09
  if( modal ) {
    cancelled <- guiGetSafe("GUIINTERNALS_cancelled")
    if( !is.na(cancelled) && !is.null(cancelled) && cancelled==TRUE )
      return(NULL)

    ##formals( func ) <- ret
    ##return( func() )
    return( guiv_last() )
  }

  ## Otherwise it is not modal
  guiSet( "FGUI_INTERNAL_guiv_function", func )
  return( invisible() )
}

## 01.23.2009 -- End of additions for the reviewers


################
## DEBUG ONLY ##
################
#testFunc <- function( a=0.0, b=1.0, c=12.0, f=13.5, g ) {
#  print("hello")
#  if( a==666 ) stop( "daemon..." )
#  return( a + b + c + f )
#}
#callback <- function( str ) {
#  print( paste( "Callback (", str, ").", sep="" ) )
#}
#
#helps <- list(a="alpha",b="beta",f="fool",c="corr",g="gummi")
#gui( testFunc,
#     argGridOrder=c(1,1,2,2,2),
#     output='s',
#     argSlider=list(output=c(0,100)),
#     helps=helps,
#     verbose=TRUE,
#     argText=list(a="alpha.a",c="cali")
#     )

#gui( testFunc,
#     argSlider=list(b=c(0,69)),
#     argOption=list(c=c("hi","bye","what")),
#     argFilename=list(g=""),
#     argCommand=list(a=function(){print("hello")}),
#     argList=list(f=c("hi","low","med")),
#     callback=callback,
#     helps=helps,
#     verbose=TRUE,
#     argFilter=list(g="{{R files} {.R .R~}}"),
#     argText=list(a="alpha.a",c="cali",b="butch",f="fairies",g="goodness_gracious") )

#gui( testFunc, argEdit=list(a=c(NA,10,TRUE), b=c(5,5,FALSE)), helps=helps )

## How to do a power interface, e.g.
#gui( testFunc, argSlider=list(a=c(1,100,0.001),b=c(1,100), c=c(1,100), f=c(1,100), output=c(1,400)), callback=guiExec, exec="", verbose=TRUE, helps=helps )


#guiNested <- function( func, nestArg, ... ) {
#  nestArgVal <- guiGetSafe(nestArg)
#  if( !is.na(nestArg) ) formals(func) <- nestArg
#
#  res <- gui( func, ... )
#  guiSet( nestArg, res )
#}
