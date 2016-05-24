
#############################################################################
##
##  Stage 2: Insert parameters for distribution and generation method
##
#############################################################################
##
##  Synopsis:
##
##    stage2 ( main )
##
##  Arguments:
##    main ... main window. It must contain all information created in stage 1.
##
##  Frames:
##    'Distributions'
##    'Generation Methods'
##    'Action'
##
##  All data are stored in 'main'.
##
##  Actions (buttons)
##
##  'ok'
##     Compute UNU.RAN objects and proceed to stage 3 (if requested)
## 
##  'cancel'
##     Quit and destroy window.
##
##  'new'
##     Restart with stage 1.
##
##
#############################################################################

## --------------------------------------------------------------------------

stage2 <- function(main) {

  ## pack everything in a group
  group <- ggroup(horizontal=FALSE, spacing=10, container=main)
  tag(main,"stage2") <- group

  ## Remark: we an create a scrollbar using
  ##   group <- ggroup(horizontal=FALSE, spacing=10, use.scrollwindow=TRUE, container=main)
  
  ## glabel("Stage 2:  Insert parameters for distribution and generation method",
  ##        container=ggroup(horizontal=TRUE,container=group))

  param.distr(main,group)
  param.method(main,group)
  param.action(main,group)
  param.buttons(main,group)
}


## --------------------------------------------------------------------------
## Insert parameter for distribution

param.distr <- function(main, group) {

  ## get data about distribution
  distr  <- tag(main,"distribution")
  type   <- distr$type
  howdef <- distr$howdef
  name   <- distr$name
  constr <- distr$constr
  
  ## the frame
  text <- paste(sub("^(\\w)"," \\U\\1",type,perl=TRUE),
                "Distribution ")
  frame <- gframe(text, horizontal=FALSE, spacing=5, container=group)

  ## print name of distribution and add handler to call help page
  help.grp <- ggroup(horizontal=TRUE, container=frame)
  glabel(name, container=help.grp)
  addSpring(help.grp)
  glabel(paste("built-in function:  ",constr),
         handler=show.help, action=list(constr), container=help.grp)
  gimage("info",dirname="stock", container=help.grp,
         handler=show.help, action=list(constr) )

  ## check for constructor
  if (! exists(constr, mode="function",
               envir=as.environment("package:Runuran"), inherits=FALSE) )
    internal.error(paste("Function '",constr,"' does not exist!",sep=""))
  
  ## two cases
  switch(howdef,
         "built-in" = param.distr.builtin(main,frame),
         "user-defined" = param.distr.userdef(main,frame),
         internal.error())
  
  ## store all widgets in main window
  tag(main,"distr.frame") <- frame
}

## ..........................................................................
## Insert parameter for distribution: built-in

param.distr.builtin <- function(main, group) {

  ## get data about distribution
  distr  <- tag(main,"distribution")
  constr <- distr$constr
  symbol <- distr$symbol

  ## create frame for parameters
  frame <- gframe(" Parameters ",
                  horizontal=FALSE, container=group)

  ## Do we have a special function?
  param.distr.special <- paste("param.distr.",toupper(symbol),sep="")
  if ( exists(param.distr.special, mode="function") ) {
    ## use special function
    params.gwl <- do.call(param.distr.special, c(main,frame))   
  }

  else {
    ## use generic function
    params.gwl <- param.distr.generic(main,frame)    
  }

  ## domain (for truncated distributions)
  domain.gwl <-
    param.distr.domain(main, group," Domain of (truncated) distribution ", params.gwl)

  ## store widgets for parameters in main window
  tag(main,"distr.params.gwl") <- c(params.gwl, domain.gwl)
  tag(main,"distr.pm.gwl") <- params.gwl
}

## ..........................................................................
## Insert parameter for distribution: user-defined

param.distr.userdef <- function(main, group) {

  ## get data about distribution
  distr  <- tag(main,"distribution")
  method <- tag(main,"method")
  type   <- distr$type
  name   <- distr$name
  symbol <- method$symbol
  constr <- distr$constr
  
  ## create frame for parameters
  frame <- gframe(" Parameters ",
                  horizontal=FALSE, container=group)

  ## get name of function for creating widgets 
  param.method.special <- paste("param.distr.",(symbol),sep="")
  if ( !exists(param.method.special, mode="function") )
    internal.error(paste("Function '",param.method.special,"' does not exist!",sep=""))

  ## call function
  params.gwl <- do.call(param.method.special, list(main=main,group=frame))

  ## domain 
  domain.gwl <-
    param.distr.domain(main, group," Domain of distribution ", params.gwl)

  ## store widgets for parameters in main window
  tag(main,"distr.params.gwl") <- c(params.gwl, domain.gwl)
}

## ..........................................................................
## Insert domain for distribution

param.distr.domain <- function(main, group, text, params.gwl=NULL) {

  ## get data about distribution
  distr  <- tag(main,"distribution")
  type   <- distr$type
  constr <- distr$constr

  ## parse arguments of constructor
  args <- function.args(constr)

  ## there must be arguments 'lb' and 'ub' 
  if ( is.null(args[["lb"]]) || is.null(args[["ub"]]) ) {
    internal.error("Could not find arguments 'lb' and 'ub' for domain!")
  }

  ## create frame
  ##   frame <- gframe(text, horizontal=TRUE, container=group)
  ## well it seems to be better to use a group instead:
  glabel(text,anchor=c(-1,0),container=group)
  frame <- ggroup(horizontal=TRUE, container=group)

  ## store widgets for domain boundary in a list
  domain.gwl <- list()

  ## left boundary of domain
  addSpace(frame,10)

  if (args[["lb"]] != "NA")
    lb <- args[["lb"]]
  else
    lb <- if (type == "continuous") "-Inf" else "0"

  glabel("lower bound ", container=frame)
  domain.gwl[["lb"]] <- gedit(text=lb, width=20,
                              coerce.with=as.numeric, container=frame)
  tag(domain.gwl[["lb"]],"label") <- "lower bound"

  if (! (is.null(params.gwl) || is.null(params.gwl[[lb]])) ) {
    val <- svalue(params.gwl[[lb]])
    if (!is.na(val)) svalue(domain.gwl[["lb"]]) <- val
    addHandlerKeystroke(params.gwl[[lb]], handler=function(h,...) {
      svalue(domain.gwl[["lb"]]) <- svalue(h$obj) })
  }

  ## right boundary of domain
  addSpace(frame, 20)

  ub <- if (args[["ub"]] != "NA") args[["ub"]] else "Inf"

  glabel("upper bound ", container=frame)
  domain.gwl[["ub"]] <- gedit(text=ub, width=20,
                              coerce.with=as.numeric, container=frame)
  tag(domain.gwl[["ub"]],"label") <- "upper bound"

  if (! (is.null(params.gwl) || is.null(params.gwl[[ub]])) ) {
    val <- svalue(params.gwl[[ub]])
    if (!is.na(val))  svalue(domain.gwl[["ub"]]) <- val
    addHandlerKeystroke(params.gwl[[ub]], handler=function(h,...) {
      svalue(domain.gwl[["ub"]]) <- svalue(h$obj) })
  }
      
  ## return list of widgets for domain
  return (domain.gwl)
}


## --------------------------------------------------------------------------
## Insert parameters for generation method

param.method <- function(main,group) {

  ## get data about method
  method <- tag(main,"method")
  type   <- method$type
  name   <- method$name
  symbol <- method$symbol
  constr <- method$constr

  ## the frame
  frame <- gframe(text=" Generation Method ",
                  horizontal=FALSE, spacing=5, container=group)

  ## print name of distribution and add handler to call help page
  help.grp <- ggroup(horizontal=TRUE, container=frame)
  glabel(name,container=help.grp)
  addSpring(help.grp)
  glabel(paste("built-in function:  ",constr),
         handler=show.help, action=list(constr), container=help.grp)
  gimage("info",dirname="stock", container=help.grp,
         handler=show.help, action=list(constr) )

  ## get name of function for creating widgets 
  param.method.special <- paste("param.method.",(symbol),sep="")
  if ( !exists(param.method.special, mode="function") )
    internal.error(paste("Function '",param.method.special,"' does not exist!",sep=""))

  ## call function
  params.gwl <- do.call(param.method.special, list(main=main,group=frame))

  ## store all widgets in main window
  tag(main,"method.frame")      <- frame
  tag(main,"method.params.gwl") <- params.gwl
}


## --------------------------------------------------------------------------
## Tasks to be performed

param.action <- function(main,group) {
  
  ## the frame
  frame <- gframe(text=" Perform ",
                  horizontal=FALSE, spacing=5, container=group)

  ## store widgets for parameters in a list
  params.gwl <- list()
  
  ## create widgets (use a tabular layout)
  tbl <- glayout(container=frame)
  
  ## assign to variables
  tbl[1,1:2,anchor=c(-1,0)] <- " Assign ..."

  ## we create variables names using the symbol of the distribution
  symbol <- tag(main,"distribution")[["symbol"]]
  if (is.na(symbol)) symbol <- "userdef"
  
  ## distribution object
  varname <- paste("distr.", symbol, sep="")
  cet <- gcheckedit(label="distribution object to", checked=FALSE,
                    text=varname, container=tbl)
  tbl[2,2,anchor=c(-1,0)] <- params.gwl[["distr.cbx"]] <- cet$cbx
  tbl[2,3,anchor=c(-1,0)] <- params.gwl[["distr.var"]] <- cet$edt


  ## generator object
  varname <- paste("gen.", symbol, sep="")
  cet <- gcheckedit(label="generator object to", checked=TRUE,
                    text=varname, container=tbl)
  tbl[3,2,anchor=c(-1,0)] <- params.gwl[["gen.cbx"]] <- cet$cbx
  tbl[3,3,anchor=c(-1,0)] <- params.gwl[["gen.var"]] <- cet$edt
  
  ## sample
  cet <- gcheckedit(label="random sample to", checked=FALSE,
                    text="x", container=tbl)
  params.gwl[["sample.size"]]<- gedit(text="100",width=10, coerce.with=as.integer,container=tbl)
  tbl[4,2,anchor=c(-1,0)] <- params.gwl[["sample.cbx"]] <- cet$cbx
  tbl[4,3,anchor=c(-1,0)] <- params.gwl[["sample.var"]] <- cet$edt
  tbl[4,4,anchor=c(-1,0)] <- "   with samplesize"
  tbl[4,5,anchor=c(-1,0)] <- params.gwl[["sample.size"]]
  
  ## show code ?
  tbl[5,1:6] <- gseparator(container=tbl)
  params.gwl[["notebook.cbx"]] <-
    gcheckbox("show R code and validate generator", checked=FALSE, container=tbl)
  tbl[6,2:3,anchor=c(-1,0)] <- params.gwl[["notebook.cbx"]]

  ## store all widgets in main window
  tag(main,"action.frame")      <- frame
  tag(main,"action.params.gwl") <- params.gwl
}


## --------------------------------------------------------------------------
## Buttons: Insert

param.buttons <- function(main,group) {

  ## the group
  buttons.grp <- ggroup(horizontal=TRUE, spacing=15, container=group)
  
  gbutton(action=gaction(label="Restart", icon="new",
            handler=function(h,...){param.clearup(main); stage1(main)}),
          container=buttons.grp)

  gbutton("cancel", container=buttons.grp,
          handler = function(h,...){dispose(main)})
  
  gbutton("help", container=buttons.grp,
          handler=show.help,
          action=list("stage2", tag(main,"distribution")$constr,
            tag(main,"method")$constr))

  gbutton("ok", container=buttons.grp,
          handler=function(h,...){param.evaluate(main)})
  
}


## --------------------------------------------------------------------------
## Evaluate input, create UNU.RAN objects and proceed to stage 3

param.evaluate <- function(main) {

  ## list of parameters in 'action' frame
  action.gwl  <- tag(main,"action.params.gwl")
  
  ## names for R variables to which values are assigned to
  distr.var   <- svalue(action.gwl[["distr.var"]])
  gen.var     <- svalue(action.gwl[["gen.var"]])
  sample.var  <- svalue(action.gwl[["sample.var"]])
  sample.size <- svalue(action.gwl[["sample.size"]])
  
  ## check boxes
  distr.assign  <- svalue(action.gwl[["distr.cbx"]])
  gen.assign    <- svalue(action.gwl[["gen.cbx"]])
  sample.assign <- svalue(action.gwl[["sample.cbx"]])
  start.notebook<- svalue(action.gwl[["notebook.cbx"]])

  ## validate given parameters
  if (!isTRUE(param.check(main))) return()
  
  ## compile R code
  Rcode <- param.get.Rcode(main)
  Rcode.show <- param.get.Rcode(main, internal=FALSE)

  ## create UNU.RAN distribution object
  ## and assign to variable in requested environment
  distr.obj <- param.unuran.run(Rcode$distr, "UNU.RAN Distribution Object")

  if (is.null(distr.obj)) return()
  
  if (isTRUE(distr.assign)) {
    cat(">",Rcode.show$distr,"\n")
    assign(distr.var, val=distr.obj, envir=tag(main,"envir"))
  }

  ## create UNU.RAN generator object 
  ## and assign to variable in requested environment
  gen.obj <- param.unuran.run(Rcode$gen, "UNU.RAN Generator Object", distr.obj=distr.obj)

  if (is.null(gen.obj)) return()

  if (isTRUE(gen.assign)) {
    cat(">",Rcode.show$gen,"\n")
    assign(gen.var, val=gen.obj, envir=tag(main,"envir"))
  }

  ## create random sample
  ## and assign to variable in requested environment
  if (isTRUE(sample.assign)) {
    cat(">",Rcode.show$sample,"\n")
    x <- param.unuran.run(Rcode$sample, "UNU.RAN Generator Object",
                          gen.obj=gen.obj, n=sample.size)
    assign(sample.var, val=x, envir=tag(main,"envir"))
  }

  ## start stage 3
  if (isTRUE(start.notebook))
    stage3(main, gen.obj) 

  ## or pop up a message
  else
    mygmessage("UNU.RAN generator object successfully generated.",
             title="UNU.RAN - message", icon="info")
}

## --------------------------------------------------------------------------
## check parameters

param.isvalid.val <- function(params.gwl,type) {
  ## check parameters given in a list of widgets.
  ## parameters must not be NAs!
  
  ## parameters values
  params.val <- lapply(params.gwl,svalue)

  ## are there any NAs?
  if (any(is.na(params.val))) {

    ## name/label of parameter
    idx <- which(is.na(params.val))[1]
    param <- names(params.gwl)[idx]
    label <- tag(params.gwl[[param]],"label")
    
    if (is.null(label))
      name <- paste("'", param, "'", sep="")
    else
      name <- paste("'", label, "' ('", param, "')", sep="")

    ## print error message
    error.message(paste("Parameter for",type,"missing or invalid:\n", name))
    return (FALSE)
  }

  return (TRUE)
}

## ..........................................................................

param.isvalid.varname <- function(param, varname.edt) {

  ## get character string
  varname <- svalue(varname.edt)

  ## check for empty string
  if ((gsub("\\s+", "", varname, perl=TRUE) == "") ) {
    error.message(paste("Variable name for",param,"missing!"))
    return (FALSE)
  }

  ## check for non-word characters
  if ((gsub("^[a-zA-Z]+[a-zA-Z0-9\\.]*", "", varname, perl=TRUE) != "")) {
    error.message(paste("Variable name for ",param," contains invalid characters:\n'",
                        varname, "'", sep=""))
    return (FALSE)
  }
  
  return(TRUE)
}

## ..........................................................................

param.check <- function(main) {

  ## parameters for distributions
  if (!isTRUE( param.isvalid.val(tag(main,"distr.params.gwl"), "distribution") ))
    return (FALSE)
  
  ## parameters for methods
  if (!isTRUE(param.isvalid.val(tag(main,"method.params.gwl"), "generation method") ))
    return (FALSE)

  ## parameters for actions
  action.gwl  <- tag(main,"action.params.gwl")
  if (! (param.isvalid.varname("distribution object", action.gwl[["distr.var"]]) &&
         param.isvalid.varname("generator object", action.gwl[["gen.var"]])      &&
         param.isvalid.varname("random sample", action.gwl[["sample.var"]]) ))
    return (FALSE)

  if (is.na( svalue(action.gwl[["sample.size"]]) )) {
    error.message("Samplesize missing or invalid!")
    return (FALSE)
  }

  ## o.k.
  return (TRUE)
}


## --------------------------------------------------------------------------
## run UNU.RAN code and catch errors

param.unuran.run <- function(Rcode, title, distr.obj=NULL, gen.obj=NULL, n=NULL ) {

  ## we have to capture messages from UNU.RAN.
  ## since we use Rprintf() statements, we need text connections.
  ## [remark: warning() in the UNU.RAN code instead of Rprintf() does not work.]
  unuran.message <- ""  ## keep R CMD check happy
  unuran.con <- textConnection("unuran.message", "w", local = TRUE)
  sink(unuran.con)

  ## evaluate given Runuran code
  robj <-
    tryCatch(eval(parse(text=Rcode)),
             error=function(e){
               msg.text <- paste(paste(unuran.message, collapse="\n"),
                                 e$message, sep="\n")
               error.message(msg.text, title=title)
               NULL })

  ## close text connection again
  sink()
  close(unuran.con)

  ## return result
  return(robj)
}


## --------------------------------------------------------------------------
## compose R code

param.get.Rcode <- function(main, internal=TRUE)
  ## internal ... if TRUE internal variable names are used
  ##              if FALSE the names given by user are used
{
  
  ## list of parameters in 'action' frame
  action.gwl  <- tag(main,"action.params.gwl")
  
  ## names for R variables to which values are assigned to
  distr.var   <- svalue(action.gwl[["distr.var"]])
  gen.var     <- svalue(action.gwl[["gen.var"]])
  sample.var  <- svalue(action.gwl[["sample.var"]])
  sample.size <- svalue(action.gwl[["sample.size"]])
  
  ## compile code for creating distribution object
  constr     <- tag(main,"distribution")$constr
  params.gwl <- tag(main,"distr.params.gwl")
  if (length(params.gwl))
    args <- paste(names(params.gwl),lapply(params.gwl,svalue), sep="=",collapse=", ")
  else
    args <- ""

  if (isTRUE(internal))
    distr.code <- paste(constr,"(",args,")", sep="")
  else
    distr.code <- paste(distr.var," <- ",constr,"(",args,")", sep="")

  ## compile code for creating generator object
  constr     <- tag(main,"method")$constr
  params.gwl <- tag(main,"method.params.gwl")
  if (length(params.gwl)) {
    args <- paste(names(params.gwl),lapply(params.gwl,svalue), sep="=",collapse=", ")
    args <- paste(",",args)
  }
  else
    args <- ""

  if (isTRUE(internal))
    gen.code <- paste(constr,"(distr=distr.obj",args,")", sep="")
  else
    gen.code <- paste(gen.var," <- ",constr,"(distr=",distr.var,args,")", sep="")

  ## compile code for drawing sample
  if (isTRUE(internal))
    sample.code <- paste("ur(gen.obj, n=",sample.size,")",sep="")
  else
    sample.code <- paste(sample.var," <- ur(",gen.var,", n=",sample.size,")",sep="")

  return (list(distr=distr.code, gen=gen.code, sample=sample.code))
}


## --------------------------------------------------------------------------
## common calls for all frames in stage 2

param.clearup <- function(main) {
  
  ## delete all widgets in window
  delete(main, tag(main,"stage2"))
  tag(main,"stage2")            <- NULL
  
  ## remove data on distribution
  tag(main,"distr.frame")       <- NULL
  tag(main,"distr.params.gwl")  <- NULL
  
  ## remove data on method
  tag(main,"method.frame")      <- NULL
  tag(main,"method.params.gwl") <- NULL
  
  ## remove actions
  tag(main,"action.frame")      <- NULL
  tag(main,"action.params.gwl") <- NULL
}

## --------------------------------------------------------------------------
