## TODO
## * manage the open windows somehow
## * left drop targetes need fixing
## * idlehandler for tree

### globals
## main windows

pmg.helpBrowser.window = NULL
pmg.plotnotebook.window = NULL
pmg.dialogs.window = NULL
pmg.cli.window = NULL

## for interactions
pmg.menuBar = NULL
pmg.toolBar = NULL
pmg.dialog.notebook = NULL
pmg.cli = NULL
pmg.statusBar = NULL

##
pmg.prompt = getOption("prompt")

pmg.window = NULL

###################################################################
##
## Functions to add to pmg 
pmg.help = function(h,...) {
  ## what to call for help page
  ## h$action contains help topic

  ## open helpBrowser if not yet
  ## else we deal with GUI
  if(is.null(pmg.helpBrowser.window) ||
     !is.gWindow(pmg.helpBrowser.window) ||
     is.invalid(pmg.helpBrowser.window)
     ) {

##      pmg.helpBrowser.window <<- ghelpbrowser()
    assignInNamespace("pmg.helpBrowser.window", ghelpbrowser(),"pmg")
  } else {
    focus(pmg.helpBrowser.window) <- TRUE
  }
  ## open page
  add(pmg.helpBrowser.window,label=h$action)
}

## function for generic widget usage
pmg.gw = function(lst, label=NULL) {
  if(!is.list(lst) || is.null(lst$variableTypeExtra)) {
    widget = ggenericwidget(lst, container=NULL, cli=pmg.cli,help.cb = pmg.help)
  } else {
    argList = list(lst=lst,cli = pmg.cli,helphandler=pmg.help, container=NULL)
    tmp = lst$variableTypeExtra ## a list
    argList[[tmp$name]] <- tmp$value
    widget = do.call("ggenericwidget",argList)
  }
  if(is.null(label)) {
    if(is.list(lst))
      label = lst$title
    else
      label = Paste(lst,"()")                         # a character string,
  }

  g = ggroup(use.scrollwindow=TRUE)
  add(g,widget, expand=TRUE)
  add(pmg.dialog.notebook, g, label=label, pageno = 3) # add near beginnign
}

### Add to the dialog notebook
pmg.add = function(widget, label) {
  add(pmg.dialog.notebook, widget, label=label, pageno=3) # add near beginning
}

### add to the menu bar
pmg.addMenubar = function(menulist) {
  add(pmg.menuBar, menulist)
}

pmg.eval = function(command, assignto=NULL) {
  if(!is.null(assignto)) names(command) <- assignto
  svalue(pmg.cli) <- command
}

### -- Not working right now
pmg.closeAll = function() {
  for(i in c(
             "pmg.cli.window","pmg.helpBrowser.window",
             "pmg.plotnotebook.window","pmg.dialogs.window")
      ) {
    window = getFromNamespace(i,"pmg")
    try(dispose(window), silent=TRUE)
    assignInNamespace(i,NULL, "pmg")
  }
}


##################################################
## call with "console" to use console, defaults to GUI
pmg = function(cliType="console", width=1000, height=.75*width,
  guiToolkit="RGtk2") {                 # getOption("guiToolkit")
  if(!interactive()) {
    cat("PMG requires an interactive environment\n")
    return()                            # no sense to have GUI if not
  }


  ## sizes
  rightWidth = width*.6                 # notebook, command area
  mainHeight = height*.8                # height without menu, tool bars

  ### which toolkit to load. If there is a gWidgets, then do that, else try pmggWidgetsRGtk
  if(guiToolkit != "RGtk2") {
    cat("pmg uses gWidgets and gWidgetsRGtk2, overriding choice of toolkit\n")
  }
  
  options("guiToolkit"="RGtk2")         # must have RGtk2 here
  
  ## what type of cli
  if(cliType != "console")
    cliType = "GUI"

  ## Make a window for pmg.gw to load into
  if(is.null(pmg.dialogs.window) ||
     !is.gWindow(pmg.dialogs.window) ||
     is.invalid(pmg.dialogs.window)
     ) {
    assignInNamespace("pmg.dialogs.window", pmgWC$new("P M G Dialogs", visible=FALSE), "pmg")
    size(pmg.dialogs.window) <- c(width, height)
  } else {
    ## raise window, exit
    return()
  }

  ## Define the main widgets
  assignInNamespace("pmg.menuBar", gmenu(pmg.menu, container=NULL), "pmg")
  assignInNamespace("pmg.dialog.notebook", gnotebook(closebuttons = TRUE,
                                                     dontCloseThese = 1, # was 1:2 before commands area moved
                                                     tearable = FALSE),
                    "pmg"
                    )
  assignInNamespace("pmg.statusBar", gstatusbar("", container=NULL),"pmg")


  ## Main layout
  mainGroup = ggroup(horizontal = FALSE, spacing=0, container=pmg.dialogs.window, expand=TRUE)

  add(mainGroup, pmg.menuBar)
  ## optional menu for user The user menu is a named list, the
  ## top-level names yield the name across the menubar
  ## if(exists("pmg.user.menu")) {
  ##     pmg.user.menu <- get("pmg.user.menu")
  ##     for(i in names(pmg.user.menu)) {
  ##         userMenu = gmenu(pmg.user.menu[[i]], name=i)
  ##         add(pmg.menubar, userMenu)
  ##     }
  ## }

  helpMenu = gmenu(help.menu, name="Help")
  add(pmg.menuBar, helpMenu)
  
  buttonBar = ggroup(spacing=0)
  add(mainGroup, buttonBar)             # toolbar

  bottomGroup = ggroup(horizontal=TRUE)
  add(mainGroup, bottomGroup, expand=TRUE)
  pmg.droparea = ggroup(horizontal=FALSE, container=bottomGroup)
  pmg.varbrowser = gvarbrowser(
    handler = function(h,...) {         # double click handler calls pmgSummary
      value = svalue(pmg.varbrowser)
      add(pmg.dialog.notebook, pmgSummary(value),
          label=Paste("Summary of ",svalue(h$obj)))
    }
    )


  ### How to layout the notebook?
  ### Try with command area below
  ## pg = gpanedgroup(pmg.varbrowser, pmg.dialog.notebook)
  ## put commands on bottom able to be expanded

  commandGroup = gexpandgroup("Command area")
  visible(commandGroup) <- TRUE
  rightPanedGroup = gpanedgroup(pmg.dialog.notebook,commandGroup,horizontal=FALSE)
  pg = gpanedgroup(pmg.varbrowser, rightPanedGroup)
  size(pmg.dialog.notebook) <- c(rightWidth,mainHeight*.67)
  
  add(bottomGroup, pg, expand=TRUE)

  add(mainGroup, pmg.statusBar)
  
  ## add buttons to buttonbar
  ## define list structure
  toolbar = list()
  ## quit
  toolbar$quit$handler = function(h,...) {
    dispose(pmg.dialogs.window)
    assignInNamespace("pmg.dialogs.window", NULL,"pmg")
#    pmgWC$closeAll()
    pmg.closeAll()
  }
  toolbar$quit$icon = "quit"
  ##
  toolbar$tmp1$separator = TRUE         # add line

  ## save workspace
  toolbar$save$handler = function(h,...) {
    gfile("Save workspace",type="save", action="save.image")
  }
  toolbar$save$icon = "save"

  ## plot notebook
  ### XXX This is an issue: cairoDevice needs to be 
    toolbar$plotnotebook$handler = function(h,...) {
      if(is.null(pmg.plotnotebook.window) ||
         !is.gWindow(pmg.plotnotebook.window) ||
         is.invalid(pmg.plotnotebook.window)
         ) {
        assignInNamespace("pmg.plotnotebook.window", pmgWC$new("P M G plot notebook", visible=TRUE ),"pmg")
        add(pmg.plotnotebook.window, ggraphicsnotebook())
      } else {
        focus(pmg.plotnotebook.window) <- TRUE
      }
    }
    toolbar$plotnotebook$icon = "plot"
##  }

  toolbar$tmp2$separator = TRUE

  
  ## fill these in
##   toolbar$print$handler = function(h,...) print("print")
##   toolbar$print$icon = "print"

  ## help
##   toolbar$help$handler = function(h,...) {
##     if(class(pmg.helpBrowser.window) != "pmgHelpBrowser")  {
## ##      pmg.helpBrowser.window <<- ghelpbrowser()
##       assignInNamespace("pmg.helpBrowser.window", ghelpbrowser(),"pmg")
##     } else {
##       ## raise window pmg.helpBrowser.window
##       focus(pmg.helpBrowser.window) <- TRUE
##     }
##   }
##   toolbar$help$icon = "help"
  
  ## make the toolbar
  tmp = gtoolbar(toolbar)
  assignInNamespace("pmg.toolBar",tmp,"pmg")
  add(buttonBar, pmg.toolBar, expand=TRUE)

  ##################################################
  ## add drop targets to left side
  ## for quick actions from varbrowser
  editDrop = gimage("edit",dirname="stock",container=pmg.droparea)
  addSpace(pmg.droparea,10);add(pmg.droparea,gseparator());addSpace(pmg.droparea,10)
  plotDrop = gimage("plot",dirname="stock",container=pmg.droparea)
  addSpace(pmg.droparea,10);add(pmg.droparea,gseparator());addSpace(pmg.droparea,10)
  summaryDrop = gimage("info",
    dirname="stock",container=pmg.droparea)
  addSpace(pmg.droparea,10);add(pmg.droparea,gseparator());addSpace(pmg.droparea,10)
  removeDrop = gimage("delete",dirname="stock",container=pmg.droparea)
  addSpring(pmg.droparea)
  
  ## add handlers
  adddroptarget(summaryDrop,handler=function(h,...) {
    svalue(pmg.cli) <-  Paste("summary(",list(h$dropdata),")")
  })
  adddroptarget(plotDrop,handler=function(h,...) {
    svalue(pmg.cli) <- Paste("plot(",list(h$dropdata),")")
  })
  adddroptarget(editDrop,handler=function(h,...) {
    ## don't do this in CLI
    rpel(Paste("fix(",list(h$dropdata),")"))
  })
  adddroptarget(removeDrop,handler=function(h,...) {
    svalue(pmg.cli) <- Paste("rm(",list(h$dropdata),")")
  })


  ## add big page to notebook to give instructions, and fix size of notebook
  useConsole = ifelse(cliType == "console",TRUE, FALSE)
  assignInNamespace("pmg.cli",
                    gcommandline("",width=rightWidth,height=mainHeight*.33,
                                              useConsole=useConsole),"pmg")
  ## put CLI and editing sheet here
  ## This adds to notebook
##   add(pmg.dialog.notebook,pmg.cli,label = "Commands",
##       pageno = 1, override.closebutton = TRUE,
##       tearable = FALSE
##       )
  add(commandGroup, pmg.cli, expand=TRUE)
  
  ## add notebook page for editing data and cli
  ## the hack keeps charaacter not factor
  x = as.numeric(NA);df=data.frame(X1=x)

  pmg.data.frame.viewer.nb = gdfnotebook(tab.pos=1, dontCloseThese=1)
  ## for some reason this gives error message
  ##  add(pmg.data.frame.viewer.nb, gdf(df,do.subset=TRUE), label= "*scratch:1*")

  add(pmg.dialog.notebook, pmg.data.frame.viewer.nb, label = "Data",
      pageno = 2, override.closebutton = TRUE
      )
  
  ## some blurb
  add(pmg.dialog.notebook,pmg.about(),label="About PMG")

  ## Finally, draw window
  visible(pmg.dialogs.window)<-TRUE
}

