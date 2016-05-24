##' @include guiComponents.R

## A widget to provide a simple command line


##################################################

##' Class for a commandline widget
setClass("gCommandline",
         contains="guiComponent",
         prototype=prototype(new("guiComponent"))
         )

##' constructor of widget for use as a command line
##'
##' Basically a giant hack
##' @exports
##' @param command Initial command to evaluation
##' @param assignto Character or NULL. assign value to this name if non-NULL
##' @param useGUI use the GUI
##' @param useConsole use the console
##' @param prompt what prompt to use
##' @param width width in pixels
##' @param height height in pixels
##' @param container parent container
##' @param ... passed to containers \code{add} method
##' @param toolkit toolkit
##' @return an object of class \code{gCommandline} with methods:
##'
##' \enumerate{
##'
##' \item{svalue<-} Character. An expression to evaluate. If it has a names attribute this is used for assignment
##'
##' \item{[} Lists history
##'
##' }
gcommandline =function(
  command = "", assignto = NULL, useGUI = TRUE, useConsole = FALSE,
  prompt = getOption("prompt"), width = 500, height = 0.6 * width,
  container = NULL, ... ,
  toolkit=guiToolkit()){
  widget =  .gcommandline (toolkit,
    command=command, assignto=assignto,
    useGUI = useGUI, useConsole=useConsole,
    prompt=prompt, width=width, height=height, container=container, ...
    )
  obj = new( 'guiComponent',widget=widget,toolkit=toolkit) 
  return(obj)
}


##' generic for toolkit dispatch
##' @alias gcommandline
setGeneric( '.gcommandline' ,
           function(toolkit,
                    command = "", assignto = NULL,
                    useGUI = TRUE, useConsole = FALSE,
                    prompt = getOption("prompt"), width = 500,
                    height = 0.6 * width, container = NULL, ... )
           standardGeneric( '.gcommandline' ))


##' gcommandline implementation for any toolkit
setClass("gCommandlineANY",
         representation=representation("gComponentANY",
           commandBox="guiWidgetOrNULL",
           outputBox="guiWidgetOrNULL",
           histList="guiWidgetOrNULL",
           useGUI = "logical",
           useConsole = "logical"
           ),
         contains="gComponentANY",
         prototype=prototype(new("gComponentANY"))
         )


##' gcommandline constructor for ANY toolkit
##' @alias gcommandline
setMethod(".gcommandline",
          signature(toolkit="ANY"),
          function(toolkit,
                   command = "", assignto=NULL,
                   useGUI = TRUE, 
                   useConsole = FALSE,
                   prompt = getOption("prompt"),
                   width = 500, height = 0.6 * width,
                   container = NULL,
                   ...) { 

            force(toolkit)

            .history = character()

            
            ## adjust command if need be
            if(nchar(command) > 0 && !is.null(assignto))
              command = addAssignto(command, assignto)


            if(is(container,"logical") && container)
              container = gwindow("gCommandLIne")
          ##
            if(useGUI == TRUE) {
              
              pg <- gpanedgroup(container=container, horizontal=FALSE, expand=TRUE)
              g <- ggroup(horizontal=FALSE, expand=FALSE, container=pg)
              bg <- ggroup(container=g, horizontal=TRUE)
              addSpring(bg)
              history <- gcheckbox("history", use.togglebutton=TRUE, container=bg)
              clear <- gbutton("clear", container=bg)
              evalCmd <- gbutton("evaluate", container=bg)
              
              commandBox <- gtext("## Type commands here:\n", container=g, expand=TRUE)
              size(commandBox) <- c(width, floor(height/3))
              
              outputBox <- gtext("", container=pg, expand=TRUE, fill="both")
              svalue(pg) <- 0.33
              
              ## History
              histWindow <- gwindow("History", parent=container, visible=FALSE)
              addHandlerUnrealize(histWindow, handler=function(h,...) {
                showHistory(FALSE)
                FALSE                                 # don't close just hide
              })
              histg <- ggroup(horizontal=FALSE, container=histWindow)
              glabel("History of recent commands", container=histg)
              histList <- gtable(data.frame(commands=character(0), stringsAsFactors=FALSE), container=histg, expand=TRUE)
              histbg <- ggroup(container=histg)
              addSpring(histbg)
              gbutton("cancel", container=histbg, handler=function(h,...) {
                svalue(history) <- FALSE
                visible(histWindow) <- FALSE # for tcltk
                focus(commandBox) <- TRUE
              })
              
              showHistory <- function(bool) {
                visible(histWindow) <- bool
              }
              
              addHandlerDoubleclick(histList, handler=function(h, ...) {
                svalue(commandBox) <- svalue(h$obj)
                focus(commandBox) <- TRUE
              })
              
              
              ## handlers
              addHandlerChanged(history, handler=function(h, ...) {
                showHistory(svalue(h$obj))
              })
              
              
              ## clear
              addHandlerChanged(clear, handler=function(h,...) {
                svalue(commandBox) <- ""
                focus(commandBox) <- TRUE
              })
              
              ## evaluate
              evaluateCommand <- function(...) {
                cmd <- svalue(commandBox)
                ## update history list
                tmp <- histList[,, drop=TRUE]
                tmp <- c(cmd, tmp)
                histList[] <- data.frame(commands=tmp, stringsAsFactors=FALSE)
                ## evaluate command, update outputbox
                out <- evalChunkReturnOutput(cmd)
                if(!out$error) {
                  insert(outputBox, out$output)
                } else {
                  ## an errow, now what
                  insert(outputBox, out$output, font.attr=c(color="red"))
                }
                ## clearm move focus to commandBox
                dispose(commandBox) 
                focus(commandBox) <- TRUE
              }
              
              addHandlerChanged(evalCmd, handler=evaluateCommand)
              
              
              
              
              obj = new("gCommandlineANY",
                block=pg,
                widget = pg,
                toolkit=toolkit,
                ID=getNewID(),
                commandBox=commandBox,
                outputBox=outputBox,
                histList=histList,
                useGUI = useGUI,
                useConsole = useConsole
                )
              
          ## initialize history
          tag(obj,"history")  <- c()
        } else {
          
          obj = new("gCommandlineANY",
            block=container,
                widget = container,
                toolkit=toolkit,
                ID=getNewID(),
                commandBox = container,
                outputBox = container,
                useGUI = useGUI,
                useConsole = useConsole
                )
            }

            if(nchar(command) > 0 )
              svalue(obj) <- command
            
            return(obj)
              
          })
          


##' return all previous, or just the index most recent
setMethod(".svalue",
          signature(toolkit="ANY",obj="gCommandlineANY"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            theArgs = list(...);

            histList <- obj@histList
            commandHistory <- histList[,1,drop=TRUE]
            
            if(length(commandHistory) == 0)
              return(c())
            if(is.null(index)) {
              return(commandHistory)
            } else {
              n = length(commandHistory)
              m = max(1, n - index + 1)
              return(rev(commandHistory[m:n]))
            }
          })

## evaluate command, store in history, swqp out widgets
setReplaceMethod(".svalue",
                 signature(toolkit="ANY",obj="gCommandlineANY"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   ## get commane
                   command = value;
                   assignto = names(value)
                   if(!is.null(assignto)) {
                     command = addAssignto(command, assignto)
                   }
                   if(obj@useGUI)
                     svalue(obj@commandBox,font.attr = "monospace") <-  command 

                   ## add to history
                   tag(obj, "history", replace=FALSE) <- command

                   retVal = evalChunkReturnOutput(command)

                   if(obj@useGUI)
                     add(obj@outputBox, retVal$output)

                   if(obj@useConsole) {
                     cat(retVal$output)
                   }
                   
                   return(obj)
                 })

## history function
setMethod("[",
          signature(x="gCommandlineANY"),
          function(x, i, j, ..., drop=TRUE) {
            .leftBracket(x, x@toolkit, i, j, ..., drop=drop)
          })

setMethod(".leftBracket",
          signature(toolkit="ANY",x="gCommandlineANY"),
          function(x, toolkit, i, j, ..., drop=TRUE) {
            
            histList <- x@histList
            commandHistory <- histList[,1,drop=TRUE]

            if(missing(i))
              return(commandHistory)
            else
              commandHistory(i)
          })



##################################################
## Helpers
## taken from Sweave
## takes a chunk, interweaves command and output
evalChunkReturnOutput = function(chunk, prompt = getOption("prompt")) {

  output = ""
  addToOutput = function(...) 
    output <<- paste(output,..., sep=" ", collapse="\n")
  
  chunkexps <- try(parse(text=chunk), silent=TRUE)
  if(inherits(chunkexps,"try-error")) {
   addToOutput(chunkexps)
   cat("Houston, we have a problem with:\n",chunkexps,"\n")
   return(list(output=output, error=TRUE))
 }

  if(length(chunkexps) == 0)
    return(list(output=output, error=FALSE))

  for(nce in 1:length(chunkexps)) {
    ce <- chunkexps[[nce]]
    dce <- deparse(ce, width.cutoff=0.75*getOption("width"))
    command = paste(prompt,
      paste(dce,collapse=paste("\n", getOption("continue"), sep="")),
      sep="", collapse=""
      )

    addToOutput(command,"\n")

    ## is there output?
    tmpcon <- file()
    sink(file=tmpcon)
    err <- RweaveEvalWithOpt(ce, list(eval=TRUE,print=FALSE,term=TRUE,visible=FALSE))
    cat("\n") # make sure final line is complete
    sink()
    theOutput <- readLines(tmpcon)
    close(tmpcon)
    ## delete empty output
    if(length(theOutput)==1 & theOutput[1]=="") theOutput <- NULL
    
    if(inherits(err, "try-error")) {
      addToOutput(err,"\n")
    } else {
      if(!is.null(theOutput)) {
        addToOutput(paste(theOutput,sep="",collapse="\n"))
      }
    }
  }

  return(list(output = output, error=FALSE))
}


### working functions


## parse command(s) and make assingment on last one.
addAssignto = function(command,assignto) {
  assignto = make.names(assignto)
  tmp = unlist(strsplit(command, ";"))
  if(length(tmp)>1) {
    command = paste(tmp[-length(tmp)], Paste(assignto,"<-",tmp[length(tmp)]), collapse=";", sep=";")
  } else {
    command =  Paste(assignto,"<-", command)
  }
  return(command)
}


