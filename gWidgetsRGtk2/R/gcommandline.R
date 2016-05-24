## command line widget
## toggles between gtext() instances containing text to edit, and output to display.

setClass("gCommandlineRGtk",
         representation=representation("gComponentRGtk",
           textGroup="guiWidget",
           editText="guiWidget",
           showText="guiWidget",
           textGroupState="character",
           editButton="guiWidget",
           clearButton="guiWidget",
           runButton="guiWidget",
           historyButton="guiWidget",
           width="numeric",
           height="numeric",
           prompt="character",
           useGUI = "logical",
           useConsole="logical"),
         contains="gComponentRGtk",
         prototype=prototype(new("gComponentRGtk"))
         )


## constructor
setMethod(".gcommandline",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   command = "", assignto=NULL,
                   useGUI = TRUE, 
                   useConsole = FALSE,
                   prompt = getOption("prompt"),
                   width = 500, height = .6*width,
                   container = NULL,
                   ...) { 

            force(toolkit)
            
            ## adjust command if need be
            if(nchar(command) > 0 && !is.null(assignto))
              command = addAssignto(command, assignto)


          ##
          if(useGUI == FALSE) {
            container = NULL
          }
            
            ## the main widgets
            group = ggroup(horizontal = FALSE, container = container, ...)
            toolbarGroup = ggroup(container=group, spacing = 0)
            textGroup = ggroup()                  # holds editText or showText
            add(group, textGroup, expand=TRUE)
            editText = gtext()
            showText = gtext()
            
            
            ## set up widget,
            
            ## toolbar
            ## the handlers
            openFile = function(h,...) {
              icl = h$action
              gfile("Select a file to read into command line",
                    type="open",
                    action = icl,
                    handler = function(h,...) {
                      file = h$file
                      tmp = tag(icl,"editText")
                      svalue(tmp) <- readLines(file)
                      if(tag(icl,"textGroupState") != "edit") {
                        delete(tag(icl,"textGroup"), tag(icl,"showText"))
                        add(tag(icl,"textGroup"), tag(icl,"editText"), expand=TRUE)
                        ## adjust buttons
                        enabled(runButton) <- TRUE
                        enabled(historyButton) <- TRUE
                        enabled(clearButton) <- TRUE
                        enabled(editButton) <- FALSE
                        ## set focus on editText??
                      }
                    })
            }
            saveFile = function(h,...) {
              icl = h$action
              win = gwindow("Save buffer contents", toolkit=toolkit)
              group = ggroup(horizontal=FALSE, container=win)
              saveFileName = gfilebrowse("",type="save", container = group)

              gp = ggroup(container=group)
              glabel("Save which values?", container=gp)
              saveType = gradio(c("commands","output"), index=FALSE, container=gp)
              gseparator(container=group)

              buttonGroup = ggroup(container=group)
              addSpring(buttonGroup)
              gbutton("save",handler=function(h,...) {
                filename = svalue(saveFileName)
                if(is.empty(filename)) {
                  cat(gettext("Need file to save to\n"))
                  return()
                }
                if(svalue(saveType) == "commands")
                  values = svalue(tag(icl,"editText"))
                else
                  values = svalue(tag(icl,"showText"))
                ## strop quotes off filename
                filename = gsub("^'","", filename)
                filename = gsub("'$","", filename)
                filename = gsub('^"',"", filename)
                filename = gsub('"$',"", filename)
                err = gtktry(writeLines(values, filename),silent=TRUE)
                if(inherits(err, "try-error")) {
                  cat(sprintf("Saving gave an error: %s\n",err))
                }
                dispose(win)
              }, container=buttonGroup)
            }
            editCode = function(h,...) {
              icl = h$action
              ## switch widgets
              delete(tag(icl,"textGroup"), tag(icl,"showText"))
              add(tag(icl,"textGroup"), tag(icl,"editText"), expand=TRUE)
              tag(icl,"textGroupState") <- "edit"

              enabled(runButton) <- TRUE
              enabled(historyButton) <- TRUE
              enabled(clearButton) <- TRUE
              enabled(editButton) <- FALSE

              ## set focus on editText?
            }
            runCode = function(h,...) {
              icl = h$action
              chunk = svalue(tag(icl,"editText"))
              svalue(icl) <- chunk
            }
            
            selectHistory = function(h,...) {
              previous = svalue(h$action, index=25)
              if(length(previous) == 0) {
                cat(gettext("No previous commandline history\n"))
                return()
              }
              win = gwindow("Select a previous value", visible=TRUE)
              group = ggroup(horizontal = FALSE, container = win)
              add(group, glabel("double click selection"))
              theHistory = gtable(previous, action = h$action,
                handler = function(h,...) {
                  newcommand = svalue(h$obj)
                  icl = h$action
                  svalue(tag(icl,"editText"), font.attr = c(style="monospace")) <- newcommand
                  ## set focus on editText?
                  dispose(win)
                })
              add(group, theHistory, expand=TRUE)
              buttonGroup = ggroup(container=group)
              addSpring(buttonGroup)
              add(buttonGroup, gbutton("cancel",handler = function(h,...) dispose(win)))
            }
            
  ## pack into widget
            add(textGroup, editText, expand=TRUE)
            ## toolbars
            sourceButton = gbutton("open", container=toolbarGroup)
            saveButton = gbutton("save", container = toolbarGroup)
            editButton = gbutton("edit",  container = toolbarGroup)
            clearButton = gbutton("clear",container = toolbarGroup)
            runButton = gbutton("evaluate",  container = toolbarGroup)
            historyButton = gbutton("history",  container= toolbarGroup)

            obj = new("gCommandlineRGtk",
              block=group,
              widget = group,
              toolkit=toolkit,
              textGroup = textGroup, editText = editText, showText = showText,
              textGroupState = "edit",
              editButton = editButton, clearButton=clearButton,
              runButton = runButton, historyButton = historyButton, 
              width=width, height=height,
              prompt =prompt,
              useGUI = useGUI,
              useConsole = useConsole)
            
            tag(obj,"showText")<-showText
            tag(obj,"editText")<-editText # delete doesn't work if it makes copis using @ slot
            tag(obj,"textGroup") <- textGroup
            tag(obj,"textGroupState") <- "edit"
            ## add handlers to buttons
            addhandlerclicked(sourceButton,  handler = openFile, action=obj)
            addhandlerclicked(saveButton,  handler = saveFile, action=obj)
            addhandlerclicked(editButton,  handler = editCode, action=obj)
            addhandlerclicked(clearButton, action=obj, function(h,...)
                              dispose(h$action@editText))
            addhandlerclicked(runButton,  handler = runCode, action=obj)
            addhandlerclicked(historyButton,  handler = selectHistory, action=obj)
            ## initialize history
            tag(obj,"history")  <- c()
            ## initialize state: used to check if swap is needed
            tag(obj,"textGroupState") <- "edit"
            
            ## which text widget?
            if(command == "") {
              enabled(editButton) <- TRUE
            } else {
              #svalue(editText) <- command
              svalue(obj) <- command
            }
            
            return(obj)
            
          })
          

### Methods
## return all previous, or just the index most recent
setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCommandlineRGtk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            theArgs = list(...);
            
            commandHistory = tag(obj,"history")
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
                 signature(toolkit="guiWidgetsToolkitRGtk2",obj="gCommandlineRGtk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   ## get commane
                   command = value;
                   assignto = names(value)
                   if(!is.null(assignto)) {
                     command = addAssignto(command, assignto)
                   }
                   if(obj@useGUI)
                     svalue(obj@editText,font.attr = c(style="monospace")) <-  command 

                   ## add to history
                   tag(obj, "history", replace=FALSE) <- command

                   evalChunk(command, obj@showText, obj@prompt, obj@useConsole, obj@useGUI)

                   ## switch widgets -- if not correct
                   textGroupState = tag(obj,"textGroupState")
                   if(!is.null(textGroupState) && textGroupState == "edit") {
                     delete(tag(obj,"textGroup"), tag(obj,"editText"))
                     add(tag(obj,"textGroup"), tag(obj,"showText"), expand=TRUE)
                   }

                   tag(obj,"textGroupState") <- "text"
                   enabled(tag(obj,"showText")) <- FALSE        # no editing of this display
                   enabled(obj@runButton) <- FALSE
                   enabled(obj@historyButton) <- FALSE
                   enabled(obj@editButton) <- TRUE
                   
                   enabled(obj@clearButton) <- FALSE
                   return(obj)
                 })

## history function
setMethod("[",
          signature(x="gCommandlineRGtk"),
          function(x, i, j, ..., drop=TRUE) {
            .leftBracket(x, x@toolkit, i, j, ..., drop=drop)
          })

setMethod(".leftBracket",
          signature(toolkit="guiWidgetsToolkitRGtk2",x="gCommandlineRGtk"),
          function(x, toolkit, i, j, ..., drop=TRUE) {
            history = tag(x, "history")

            if(missing(i))
              return(history)
            else
              history(i)
          })

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




## taken from Sweave
## takes a chunk, iterweaves command and output
evalChunk = function(chunk, widget, prompt = getOption("prompt"),
  useConsole=FALSE, useGUI = TRUE) {
  svalue(widget) <- ""                 # clear out

  chunkexps <- gtktry(parse(text=chunk), silent=TRUE)
  if(inherits(chunkexps,"try-error")) {
    if(useGUI)
      add(widget, chunkexps, font.attr = c(style="monospace"))
#    addTextWidget(widget, chunkexps)
    cat(sprintf("Houston, we have a problem with: %s\n",chunk))
    return(c())
  }
  if(length(chunkexps) == 0)
    return(c())
#  output = c()

  for(nce in 1:length(chunkexps)) {
    ce <- chunkexps[[nce]]
    dce <- deparse(ce, width.cutoff=0.75*getOption("width"))
    command = Paste(prompt,
      paste(dce,collapse=paste("\n", getOption("continue"), sep=""))
      )
    if(useGUI)
      add(widget, command, font.attr = c(style="monospace",color="red",weight="italic"))

    if(useConsole)
      cat(command,"\n")
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
      if(useGUI)
        add(widget, err, font.attr=c(style="monospace",color="red",weight="bold"))
      if(useConsole)
        cat(err,"\n")
    } else {
      if(!is.null(theOutput)) {
        if(useGUI)
          add(widget, theOutput, font.attr = c(style="monospace"))
        if(useConsole) 
          cat(paste(theOutput,sep="",collapse="\n"),"\n")
      }
    }
  }
}


