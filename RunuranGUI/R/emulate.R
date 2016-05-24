#############################################################################
##
##  Emulate parameter settings and mouse clicks in GUI
##
#############################################################################
##
##  Synopsis:
##
##  emul.stage1 ( main, delay, button, distr, method )
##  emul.stage2 ( main, delay, button, distr, method, action )
##
##  Arguments:
##    main ... main window
##    delay ... sleeping time until button is pressed
##    button ... kind of button
##    distr  ... list of parameters for distribution
##    method ... list of parameters for generation method
##    action ... list of parameters for performing tasks
##
##  See below for more details
##
#############################################################################
## 
##  These functions are currently not exported!
## 
#############################################################################


## --------------------------------------------------------------------------

emul.stage1 <- function(main, delay = 0,
                        button = c("none","ok","cancel","new"),
                        distr = list(), method = list())
  ## ........................................................................
  ## Emulate mouse clicks for stage 1
  ##
  ##   main   ... "Runuran" window in stage 1
  ##   delay  ... delay before button is clicked
  ##   button ... clicked button; one of
  ##                "none", "ok", "cancel", or "new"
  ##   distr  ... new settings for widgets in distribution frame;
  ##              accepted variables in list:
  ##                 'type', 'howdef', 'cont', 'discr'
  ##              values are either:
  ##               - index of position in widget (type='numeric')
  ##               - or a (sub)string of a valid entry (type='character')
  ##   method ... new settings for widgets in method frame;
  ##              accepted variables in list:
  ##                 'type', 'cont', 'discr'
  ##              values are either:
  ##               - index of position in widget (type='numeric')
  ##               - or a (sub)string of a valid entry (type='character')
  ## ........................................................................
{
  ## check arguments
  if( id(main) != "Runuran" || is.null(tag(main,"stage1")) )
    stop ("Invalid argument 'main'")
  if(! is.list(distr))
    stop ("Invalid argument 'distr'")
  if(! is.list(method))
    stop ("Invalid argument 'method'")

  ## wait
  Sys.sleep(delay)

  ## set parameters for distribtions
  emul.set.select(main, "distr.type.rbx",   distr$type,   DISTRIBUTIONS.TYPE)
  emul.set.select(main, "distr.howdef.rbx", distr$howdef, DISTRIBUTIONS.HOWDEF)
  emul.set.select(main, "distr.cont.cbb",   distr$cont,   DISTRIBUTIONS[["continuous"]]["name",])
  emul.set.select(main, "distr.discr.cbb",  distr$discr,  DISTRIBUTIONS[["discrete"]]  ["name",])

  ## set parameters for methods
  emul.set.select(main, "method.type.rbx",  method$type,  METHODS.TYPE)
  emul.set.select(main, "method.cont.cbb",  method$cont,  METHODS[["continuous"]]["name",])
  emul.set.select(main, "method.discr.cbb", method$discr, METHODS[["discrete"]]  ["name",])
  
  ## emulate click on button
  button <- match.arg(button)
  switch (button,
          "none"   = { FALSE },
          "ok"     = { type.evaluate(main) },  
          "cancel" = { dispose(main) },
          "new"    = { type.clearup(main); stage1(main) }
          )

  return(invisible(NULL))
}


## --------------------------------------------------------------------------

emul.stage2 <- function(main, delay=0,
                        button=c("none","ok","cancel","new"),
                        distr=list(), method=list(), action=list() )
  ## ........................................................................
  ## Emulate mouse clicks for stage 2
  ##   main   ... "Runuran" window in stage 1
  ##   delay  ... delay before button is clicked
  ##   button ... clicked button; one of
  ##                "none", "ok", "cancel", or "new"
  ##   distr  ... settings for widgets in distribution frame;
  ##              accepted variables in list depend on the
  ##              choosen distribution.
  ##              keys 'VAR1', 'VAR2', ... refer to the
  ##              first, second, ... element in the list of parameter widgets.
  ##   method ... settings for widgets in method frame;
  ##              accepted variables in list depend on the
  ##              choosen method.
  ##   action ... settings for widgets in perform frame;
  ##              accepted variables in list:
  ##                 'distr = c(cbx,var) = c(logical, character)',
  ##                 'gen = c(cbx,var) = c(logical, character)',
  ##                 'sample = c(cbx,var,size) = c(logical, character, integer)',
  ##                 'notebook = cbx = logical'
  ## ........................................................................
{
  ## check arguments
  if( id(main) != "Runuran" || is.null(tag(main,"stage2")) )
    stop ("Invalid argument 'main'")

  ## wait
  Sys.sleep(delay)

  ## get data
  distr.params  <- tag(main,"distr.params.gwl")
  method.params <- tag(main,"method.params.gwl")
  action.params <- tag(main,"action.params.gwl")

  ## set parameters for distribtions
  if (length(distr)) {
    for (key in names(distr))
      emul.set.edit(main,distr.params,key,distr[[key]])
  }

  ## set parameters for methods
  if (length(method)) {
    for (key in names(method))
      emul.set.edit(main,method.params,key,method[[key]])
  }

  ## set parameters for performing tasks
  action.items <- c("cbx","var","size")
  if (length(action)) {
    for (key in names(action)) {
      for (pos in 1:length(action[[key]])) {
        var <- paste(key,".",action.items[pos], sep="")
        emul.set.edit(main,action.params,var,action[[key]][pos])
      }
    }
  }

  ## emulate click on button
  button <- match.arg(button)
  switch (button,
          "none"   = { FALSE },
          "ok"     = { param.evaluate(main) },
          "cancel" = { dispose(main) },
          "new"    = { param.clearup(main); stage1(main) }
          )

  return(invisible(NULL))
}


## --------------------------------------------------------------------------
##
## Auxiliary functions: Set parameters
##
## --------------------------------------------------------------------------

emul.set.select <- function(main,name,val,vtable)
  ## ........................................................................
  ## Set parameters in radio group and droplist widgets.
  ##
  ##   main   ... "Runuran" window
  ##   name   ... name of widget stored in 'main'
  ##   val    ... new value for widget:
  ##              val is either
  ##               - an index (type='numeric'),
  ##               - or a (sub)string of a valid entry (type='character')
  ##   vtable ... table of arguments (possible strings) for 'val'
  ##              (only required if 'val' is of type 'character')
  ## ........................................................................
{
  ## check argument
  if (is.null(val))
    return()  ## nothing to do

  ## get widget
  obj <- tag(main,name)
  if (is.null(obj)) {
    warning(paste("unkown tag '",name,"'",sep=""))
    return ()
  }
  
  ## check type of value
  if ( typeof(val) == "character" ) {
    if (missing(vtable))
      stop ("table of arguments missing")
    ## get position in table
    idx <- pmatch(val,vtable)
    if (is.na(idx))
      stop (paste("invalid argument 'val': \"",
                  val,"\" does not match a valid entry", sep=""))
    val <- idx
  }
  
  ## set new value
  svalue(obj,index=TRUE) <- val
}  


## --------------------------------------------------------------------------

emul.set.edit <- function(main,objs,name,val)
  ## ........................................................................
  ## Set parameters in gedit or glabel widgets.
  ##
  ##   main ... "Runuran" window
  ##   objs ... list of widget
  ##   name ... name of widget stored in 'objs';
  ##            names 'VAR1', 'VAR2', ... refer to the
  ##            first, second, ... element in the list of parameter widgets.
  ##   val  ... new value for widget
  ## ........................................................................
{
  ## check argument
  if (is.null(val))
    return()  ## nothing to do

  ## get widget
  if (regexpr("^VAR(\\d+)$", name, perl=TRUE)==1) {
    ## name is of type 'VAR?'
    key <- as.integer(sub("^VAR(\\d+)$", "\\1", name, perl=TRUE))
    if (key > length(objs) || key < 1) {
      warning(paste("variable '",name,"' out of bounds",sep=""))
      return()
    }
    obj <- objs[[key]]
  }

  else {
    obj <- objs[[name]]
    if (is.null(obj)) {
      warning(paste("unknown widget '",name,"'",sep=""))
      return()
    }
  }

  ## set new value
  svalue(obj) <- val
}  

## --------------------------------------------------------------------------

