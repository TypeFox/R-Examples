## file chooser dialog: creates gfile and gfilebrowser
setMethod(".gfile",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   text="",
                   type=c("open","save","selectdir"),
                   initialfilename = NULL,
                   filter =  list(
                     "All files"=list(
                       patterns=c("*")
                       ),
                     "R files"=list(
                       patterns=c("*.R","*.Rdata")
                       ),
                     "text files"=list(
                       mime.types=c("text/plain")
                       )
                     ),
                   handler = NULL,
                   action = NULL,                     # 
                   ...
                   ) {
            
            force(toolkit)
            
            args = list(...)

            ## this will be in the API, for now we pass in through ...
            multiple <- getWithDefault(args$multiple, FALSE)
            ## pass in initial dir information, or get from filename, or get from option, or setwd
            initialdir <- args$initialdir
            if(is.null(initialdir) && !is.null(initialfilename))
              initialdir <- dirname(initialfilename)
            if(is.null(initialdir))
              initialdir <- getOption("gWidgetstcltk::gfile_initialdir")
            ## may still be NULL that is okay
            options("gWidgetstcltk::gfile_initialdir"=initialdir) # store
            
            type = match.arg(type)

            ## different things depending on type
            if(type == "open") {

              
              theFilter = ""
              if(!is.null(filter)) {

                ## turn named character vector into list of patterns
                if(is.character(filter)) {
                  filter <- sapply(names(filter), function(nm) {
                    list(patterns=paste(".", filter[nm], sep=""))
                  }, simplify=FALSE)
                  filter[['All files']]$patterns = "*"
                } 
                
                for(i in names(filter)) {
                  pats = filter[[i]]$patterns
                  if(!is.null(pats)) {
                    theFilter = paste(theFilter,"{{",
                      i,"} ",
                      if(length(pats) > 1)
                       paste("{",paste(filter[[i]]$patterns,collapse=" "),
                             "}} ", sep="",collapse="")
                      else
                      paste(pats,"} ",sep="",collapse=""),
                      sep="",collapse="")
                  }
                }
              } else {
                theFilter = "{{All files} *}"
              }

              l <- list(title=text, filetypes=theFilter, multiple=multiple)
              if(!is.null(initialfilename))
                l$initialfile=initialfilename
              if(!is.null(initialdir))
                l$initialdir=initialdir
              
              val <- do.call("tkgetOpenFile", l)
              if(multiple) {
                val <- as.character(val) # empty = character(0)
                if(length(val) == 0)
                  val <- NA
              } else {
                val <- tclvalue(val)    # empty=""
                if(val == "")
                  val <- NA
              }
              ## save initialdir information
              if(!is.na(val[1]) && nchar(val[1]) > 0)
              options("gWidgetstcltk::gfile_initialdir"=dirname(val[1]))

            } else if(type == "save") {

              l <- list(title=text)
              l$initialfile <- initialfilename
              val <- do.call(tkgetSaveFile, l)
              val <- tclvalue(val)
              
            } else if(type == "selectdir") {

              val <- tkchooseDirectory()
              val <- tclvalue(val)
              
            }


            
            if (length(val) > 1 || nchar(val) > 0) {
              h = list(obj = NULL, action=action, file=val)
              if(!is.null(handler)) 
                handler(h)
              
              ## how to return filename?
              return(val)
            } else {
              ## cancel
              return(NA)
            }
            
            
          })


##################################################
## gfilebrowser is not modal, like gfile
setClass("gFilebrowsetcltk",
         contains="gEdittcltk",
         prototype=prototype(new("gEdittcltk"))
         )


## create a browse button -- put value into text box
setMethod(".gfilebrowse",
          signature(toolkit="guiWidgetsToolkittcltk"),
          function(toolkit,
                   text="Select a file...", type="open",  quote=TRUE,
                   container=NULL, ...) {

            if(is(container,"logical") && container)
              container = gwindow()
            if(!is(container,"guiWidget")) {
              warning("Container is not correct. No NULL containers possible\n" )
              return()
            }


            
            group = ggroup(horizontal=TRUE, container=container)
            entry = gedit(text=text, container=group, ...)
            browseButton = gbutton("browse",container=group)

            args <- list(...)
            filter <- args$filter
            initialfilename <- args$initialfilename
            
            file.cb = function(h,...) {
              ## called when button is clicked
              
              ## in this h is gFile object, not gBrowse object
              gfile(text=text,
                    type = type,
                    handler = function(h,...) svalue(entry) <- h$file,
                    quote = TRUE,
                    filter=filter, initialfilename=initialfilename
                    )
            }
            addhandlerclicked(browseButton,handler=file.cb)


            ## put entry as widget to pick up gEdit methods
            obj = new("gFilebrowsetcltk",
#              block=group, widget=entry@widget@widget, toolkit=toolkit,ID=getNewID())
              block=group, widget=entry@widget,
              toolkit=toolkit,ID=getNewID(),e = new.env())

            tag(obj,"entry") <- entry
            
            invisible(obj)
          })


setMethod(".svalue",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gFilebrowsetcltk"),
          function(obj, toolkit, index=NULL, drop=NULL, ...) {
            entry <- tag(obj,"entry")
            svalue(entry,index,drop,...)
          })

## svalue<-
setReplaceMethod(".svalue",
                 signature(toolkit="guiWidgetsToolkittcltk",
                           obj="gFilebrowsetcltk"),
                 function(obj, toolkit, index=NULL, ..., value) {
                   entry <- tag(obj,"entry")
                   svalue(entry, index, ...) <- value
                   return(obj)
          })

## Pass down to entry -- id must good for entry though XXX could be fixed
setMethod(".addhandlerchanged",
          signature(toolkit="guiWidgetsToolkittcltk",obj="gFilebrowsetcltk"),
          function(obj, toolkit, handler, action=NULL, ...) {
            entry <- tag(obj, "entry")
            addHandlerChanged(entry, handler, action, ...)
          })
