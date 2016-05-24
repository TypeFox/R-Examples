## creates a notebook interface tohandle plots

setClass("gGraphicsNotebookGtk",
         representation=representation(
           width="numeric",height="numeric"
           ),
         contains="gNotebookRGtk",
         prototype=prototype(new("gNotebookRGtk"))
         )
setMethod(".ggraphicsnotebook",
          signature(toolkit="guiWidgetsToolkitRGtk2"),
          function(toolkit,
                   width=dpi*6, height=dpi*6,dpi=75,
                   container = NULL,
                   ...) {
            ## ... passed onto gnotebook
            force(toolkit)
            
            group = ggroup(horizontal = FALSE, container=container, ...)
            
            ## make toolbar
            toolbargroup = ggroup(horizontal=TRUE, container=group)
            notebook = gnotebook(closebuttons = TRUE,dontclosethese=1)
            add(group, notebook, expand=TRUE)
            
            ## store both gnotebook
            obj = new("gGraphicsNotebookGtk",
              block=group, widget=getWidget(notebook), # gWidgetNotebookRGtk
              toolkit=toolkit, width=width, height=height)

  
            ##   ## shove in width, height into notebook
            ##   notebook$width=width; notebook$height=height
            
            toolbar = list()
            if(is.gWindow(container)) {
              toolbar$Quit$handler = function(h,...) dispose(container)
              toolbar$Quit$icon = "quit"
              toolbar$tmp1$separator = TRUE
            }
            
            toolbar$New$handler = function(h,...) {
              addNewPage(obj)
            }
            toolbar$New$icon = "new"
            
            toolbar$Close$handler = function(h,...) {
              ## we need to close the current device -- this was fixed??
              ##currentDevice = names(obj)[svalue(obj)]
              ##currentDevice = as.numeric(gsub("^dev:","",currentDevice))
              dispose(obj)
              ##dev.off(currentDevice) 
            }
            toolbar$Close$icon = "close"
            toolbar$tmp2$separator = TRUE
            
            toolbar$Print$handler = function(h,...) {
              printCurrentPage(obj)
            }
            toolbar$Print$icon="print"
            
            toolbar$Save$handler = function(h,...) {
              saveCurrentPage(obj)
            }
            toolbar$Save$icon="save"

            toolbar$tmp3$separator = TRUE

            ## toolbar$Record$handler = function(h,...) {
            ##   recordCurrentPage(obj)
            ## }
            ## toolbar$Record$icon = "media-record"
            
            toolbar$Replay$handler = function(h,...) {
              replayAPlot(obj)
            }
            toolbar$Replay$icon = "media-play"
            
            gtoolbar(toolbar, container=toolbargroup)
            

            ## start with a plot, then add handler!
            addNewPage(obj)
            

            ## add handler to raise device when page is changed
            addhandler(obj, signal="switch-page", function(h,...) {
              ## we need the new page not the old page
              ## oldPage = svalue(obj)

              theArgs = list(...)
              newPage = theArgs[[3]] + 1 ## gtk thingy

              
              plotWidget = notebook[newPage] ## notebook from scope
              devNo = tag(plotWidget, "device")



              
              availDevs = dev.list()
              availDevs = availDevs[names(availDevs) == "Cairo"]
              
              if(!is.null(devNo) && is.numeric(devNo) && devNo %in% availDevs)
                dev.set(devNo)

              return(TRUE)
            })
            

            
            return(obj)
            
          })
          
### Two key handlers for when page is raised or lowered
## remove plot device when it is unrealized
## unrealizePage = function(h,...) {
##  remove from list of devices
##  in the new cairoDevice ths gives errors
##     theDevice = h$action$device
##     if(!is.null(theDevice) && theDevice > 1) {
##     if(.Platform$OS != "windows")
##       try(dev.off(theDevice), silent=TRUE)
##     ## doesn't work! we get problems with devices here big time
##     else
##       if(theDevice %in% dev.list()) dev.off(theDevice)
##     }
##   return(TRUE)
## }

  
## when a page is entered, we set the plot device
## enterPage = function(h,...) {
##    devNo = h$action$device

##    if(!is.null(theDevice)) {
##      dev.set(theDevice)
##    }
##   return(TRUE)
## }

## value is ignored
addNewPage = function(obj,  ...) {

  theArgs = list(...)                      # ignored for now

  width = obj@width;height=obj@height


  plotwindow = ggraphics(width,height)

  
  ## These two are now made obsolete
##  addhandlerexpose(plotwindow,
##             handler=enterPage, action=list(device=dev.cur()))

##  addhandlerunrealize(plotwindow,
##             handler=unrealizePage, action=list(device=dev.cur()))

  noPlots = tag(obj,"noPlots")
  if(is.null(noPlots)) noPlots = 0
  
  label = paste("plot:",noPlots + 1,sep="",collapse="")
  tag(obj,"noPlots") = noPlots+1
  tag(plotwindow@widget,"device") <- dev.cur()

  
  ## add to notebook
  add(obj, plotwindow, label=label)

  
}


printCurrentPage = function(obj,...) {
  ## do a confirmation dialog first
  gconfirm("Really print this graph?",
           ok.handler= function(...) dev.print()
           )
}

## recordCurrentPage = function(obj,...) {
##   ## get variable name, record to this
##   ginput(message="Enter a variable to assign recorded plot to:",
##          icon="question",
##          handler = function(h,...) {
##            assign(h$input,recordPlot(),envir=.GlobalEnv)
##          })
## }

replayAPlot = function(...) {

  recordedPlots =c()
  for(i in ls(envir=.GlobalEnv)) {
    if(inherits(getObjectFromString(i),"recordedplot")) {
      recordedPlots = c(recordedPlots, i)
    }
  }
  
  if(length(recordedPlots) == 0) {
    gmessage("No recorded plots were found")
  } else {
    
    win = gwindow("Replay plot", visible=TRUE)
    group = ggroup(horizontal=FALSE, container=win)
    glabel("Select variable name with recorded plot", container=group)
    dl = gdroplist(recordedPlots, container = group)
    buttonGroup = ggroup(container=group)
    addSpring(buttonGroup)
    okButton = gbutton("ok",handler = function(h,...) {
      thePlot = svalue(dl)
      if(nchar(thePlot) > 0) {
        replayPlot(getObjectFromString(thePlot))
        dispose(win)
      } else {
        warning("Select a plot")
      }
    },
      container=buttonGroup)
    
    cancelButton = gbutton("cancel",
      handler = function(h,...) {
        dispose(win)
      },
      container=buttonGroup)
  }
}

## dialog to save current page
saveCurrentPage = function(obj) {
  win = gwindow("Save current plot", visible=TRUE)
  group = ggroup(horizontal=FALSE, container=win)

  warningMessage = glabel("Saving a plot is currently kind of flaky", container=group)
  gseparator(container=group)

  tbl = glayout()

  knownFileTypes = c("ps","eps","pdf","jpg","jpeg","png")
  knownFileTypes = c("png")
  
  filetype = gdroplist(knownFileTypes)
  filename = gfilebrowse(".png")
  
  saveButton = gbutton("save",handler = function(h,...) {
    theFileName = svalue(filename)
    if(nchar(theFileName) == 0) {
      gmessage("No filename selected")
    } else {
      ## get filetype from filename -- not from filetype
      tmp = unlist(strsplit(theFileName,"\\."))
      
      if(nchar(tmp[1]) == 0) {
        gmessage("No filename given")
      } else if(length(tmp) == 1) {
        gmessage("No filetype selected (based on extension)")
      } else {
        theFileType = tmp[length(tmp)]


        if(theFileType %in% knownFileTypes) {
          ## use undocumented part of gnotebook
#          drawarea = getNotebookPageWidget(obj$notebook)
          drawarea = getNotebookPageWidget(obj)
          newobj = as.gGd(drawarea)
          svalue(newobj) <- list(file=theFileName, extension=theFileType)
        } else {
          cat(sprintf("***\n Don't know this extension: %s\n",theFileType))
        }
        dispose(win)
      }
    }
  })

    addhandlerchanged(filetype, handler = function(h,...) {
    curFileName = svalue(filename)
    newFileType = svalue(filetype)
    tmp = unlist(strsplit(curFileName,"\\."))
    if(length(tmp) > 1)
      tmp[length(tmp)] = newFileType
    else
      tmp = c(tmp, newFileType)
    svalue(filename)<- paste(tmp,collapse=".")
#    focus(filename)<-TRUE
  })
    
  
  tbl[1,1] = glabel("filetype:")
  tbl[1,2] = filetype
  tbl[2,1] = glabel("filename:")
  tbl[2,2:4] = filename
  
  visible(tbl, set=TRUE)
  add(group, tbl, expand=TRUE)
  buttonGroup = ggroup(container=group)
  addSpring(buttonGroup)
  add(buttonGroup,saveButton)
  add(buttonGroup, gbutton("cancel",handler=function(h,...) dispose(win)))
}



