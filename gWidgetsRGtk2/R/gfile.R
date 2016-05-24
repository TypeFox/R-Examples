## file chooser dialog: creates gfile and gfilebrowser
setMethod(".gfile",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
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
                   multi=FALSE, ## XXX uncomment at some point
                   handler = NULL,
                   action = NULL,                     # 
                   ...
                   ) {
            
            force(toolkit)

            args = list(...)
            
            type = match.arg(type)

            availTypes = c(
              "open"="open",
              "save"="save",
              "selectdir"="select-folder",
              "createdir"="create-folder"
              )
            
            actiontype = GtkFileChooserAction[availTypes[type]]

            
            buttonWithId = list(
              "ok"= c("gtk-ok",GtkResponseType["ok"]),
              "cancel" = c("gtk-cancel",GtkResponseType["cancel"])
              )
            
            whichButtons = switch(type,
              "save"=c("ok","cancel"),
              "open"=c("ok","cancel"),
              "selectdir"=c("ok","cancel")
              )
            
            
            okhandler.default = function(h,...) {
              if(is.gComponent(h$action)) {
                if(!is.null(args$quote))
                  svalue(h$action) <- Paste("'",h$file,"'")
                else
                  svalue(h$action) <- h$file
              } else {
                do.call(h$action,list(h$file))
              }
            }
            ## give a default of printing.
            if(is.null(handler)) {
              handler = okhandler.default
              if(is.null(action))
                action="print"
            }
            
            cancelhandler = function(h,...) {
              dispose(h$obj)
              return(NA)
            }


            
            filechooser = gtkFileChooserDialogNew(title=text, action=actiontype)
            filechooser$setSelectMultiple(multi)
            
            for(i in whichButtons) 
              filechooser$AddButton(buttonWithId[[i]][1],buttonWithId[[i]][2])
            
            ## add a filter
            if(!is.null(filter) && type %in% c("open","save")) {
              if(is.character(filter)) {
                ## make alist
                filter <- sapply(names(filter), function(nm) {
                  list(patterns=paste("*.", filter[nm], sep=""))
                }, simplify=FALSE)
                filter[[gettext("All files")]]$patterns <- "*.*"
              }

              for(i in names(filter)) {
                filefilter = gtkFileFilterNew()
                filefilter$SetName(i)
                if(!is.null(filter[[i]]$patterns)) {
                  for(pattern in filter[[i]]$patterns)
                    filefilter$AddPattern(pattern)
                }
                if(!is.null(filter[[i]]$mime.types)) {
                  for(mime.type in filter[[i]]$mime.types)
                    filefilter$AddMimeType(mime.type)
                }
                filechooser$AddFilter(filefilter)
              }
            }
            
            
            ## initialize
            if(!is.null(initialfilename)) {
              if(type == "open") {
                filechooser$SetFilename(Paste(getwd(),.Platform$file.sep,initialfilename))
              } else if(type == "save") {
                filechooser$setCurrentFolder(getwd())
                filechooser$setCurrentName(initialfilename)
              }
            }
            
            ## this makes it modal
            response = filechooser$Run()
            
#            file=filechooser$GetFilename()
            
            ## return a vector of chars for multi select - TT
            file=unlist(filechooser$GetFilenames())
            Encoding(file) <- "UTF-8"
            
            h = list(obj=filechooser,action=action,file=file)
            if(response == GtkResponseType["cancel"]) {
              ## just close
              filechooser$Destroy()
              return(NA)
            } else if(response == GtkResponseType["ok"]) {
              filechooser$Destroy()
              if(!is.null(handler)) handler(h)
              if(!is.null(args$quote) && as.logical(args$quote))
                return(paste("'",file,"'",sep=""))
              else
                return(file)
            } else {
              filechooser$Destroy()
              return(NA)
            }
          })


##################################################
## gfilebrowse is not modal, like gfile
setClass("gFilebrowseRGtk",
         contains="gEditRGtk",
         prototype=prototype(new("gEditRGtk"))
         )


## create a browse button -- put value into text box
setMethod(".gfilebrowse",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   text="Select or drag name here...", type="open",
                   quote=TRUE,
                   container=NULL, ...) {


            theArgs = list(...)
            if(!is.null(theArgs$expand) && as.logical(theArgs$expand))
              group = ggroup(horizontal=FALSE, container=container, expand=TRUE)
            else
              group = ggroup(horizontal=FALSE, container=container)

            g <- ggroup(container =group, horizontal=TRUE, expand=FALSE)   # stop growth of button
            gedit.local <- function(...) gedit(..., expand=TRUE)
            entry = gedit.local(text=text, container=g, ...)

            file.cb = function(h,...) {
              ## called when button is clicked
              ## pop up a calendar, when date selected, copy to entry
              
              ## in this h is gFile object, not gBrowse object
              val = gfile(text=text,
                type = type,
                quote = quote,          
                filter = theArgs$filter,
                action = function(...) invisible(...) # don't print
                )
              if(!is.na(val))           # check return, set if not NA
                svalue(entry) <- val
            }


            gbutton("browse", handler=file.cb, container=g, expand=FALSE)

            ## put entry as widget to pick up gEdit methods
            obj = new("gFilebrowseRGtk",
              block=group, widget=getWidget(entry), toolkit=toolkit)

            invisible(obj)
          })

