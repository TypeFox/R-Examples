# Taken from John Verzani's gfile
## file chooser dialog: creates gfile and gfilebrowser
my.gfile = function( text="",
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
                        
            filechooser = gtkFileChooserDialogNew(title=text, action=actiontype)
            filechooser$setSelectMultiple(multi)
            
            for(i in whichButtons) 
              filechooser$AddButton(buttonWithId[[i]][1],buttonWithId[[i]][2])
            
            ## add a filter
            if(!is.null(filter) && type %in% c("open","save")) {
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
                filechooser$SetFilename(paste(getwd(),.Platform$file.sep,initialfilename))
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
              if(!is.null(args$quote) && as.logical(args$quote))
                return(paste("'",file,"'",sep=""))
              else
                return(file)
            } else {
              filechooser$Destroy()
              return(NA)
            }
          }


