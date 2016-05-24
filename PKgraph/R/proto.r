#############################################################################################
## Project: PKgraph
## File: proto.R
## Author: Xiaoyong Sun
## Date: 08/19/2009
## Goal: PKgraph
## Notes:
##      - class style from proto
#############################################################################################




###############################################################################
## Main sub-GUI
###############################################################################
BasicGUI = proto(
  new = function(., message = "Basic GUI",...)
  {
     .$proto(message=message,...)
  },

  ## method to check if window has been drawn or destroyed
  isVisible = function(.,win = .$window)
  {
    if(!is.null(win) && is(win,"guiWidget") && isExtant(win)) return(TRUE)
    return(FALSE)
  },

  show = function(.,...)
  {
    ## ... passed to gwindow
    ## check if window is already there
    if(.$isVisible()) return()

    .$window <- gwindow(title=.$message, parent=parent.win, height=getSubHeight(), width=getSubWidth(),...)
    g = ggroup(horizontal=FALSE, cont=.$window, use.scrollwindow = TRUE)

    ## group for toolbar and menubar
    if(!is.null(.$menubarList) || !is.null(.$toolbarList))
    {
      g2 = ggroup(horizontal=FALSE, container=g, expand=FALSE)
      if(!is.null(.$menubarList)) .$menubar <- gmenu(.$menubarList, cont=g2)
      if(!is.null(.$toolbarList)) .$toolbar <- gtoolbar(.$toolbarList, style="icons",cont=g2)
      gseparator(cont=g2)
    }

    ## container for body --e xpand = TRUE
    g1 = ggroup(horizontal=FALSE, container=g, expand=TRUE) # expand
    .$makeBody(container = g1)
    .$makeSave(container = g1)
    .$makeButtons(container = g1)
    ## command area
    #commandGroup = gexpandgroup("Command area", expand=TRUE)
    commandGroup = ggroup(expand=TRUE)
    gline = gcommandline("", cont=commandGroup)
    global.data <- getCurrentData()
    svalue(gline) <- "#data name: global.data"
    add(g1, commandGroup)

    ########################################

    leftpane =gpanedgroup(g1)
    gim = ggroup(horizontal=TRUE, expand=TRUE) # expand
    #addSpace(gim, 50)
    pk.dialog.image1 <- gimage("save", dirname="stock", cont=gim, size="large_toolbar")
    gseparator(horizontal=FALSE, cont=gim)
    pk.dialog.image2 <- gimage("ts", dirname="stock", cont=gim, size="large_toolbar")
    gseparator(horizontal=FALSE, cont=gim)
    pk.dialog.image3 <- gimage("dataframe", dirname="stock", cont=gim, size="large_toolbar")
    gseparator(horizontal=FALSE, cont=gim)
    pk.dialog.image31 <- gimage("save-as", dirname="stock", cont=gim, size="large_toolbar")    
    gseparator(horizontal=FALSE, cont=gim)    
    pk.dialog.image4 <- gimage("close", dirname="stock", cont=gim, size="large_toolbar")
    addhandlerclicked(pk.dialog.image1, handler=.$saveImageHandler)
    addhandlerclicked(pk.dialog.image2, handler=.$ggobiImageHandler)
    addhandlerclicked(pk.dialog.image3, handler=.$dataframeImageHandler)
    addhandlerclicked(pk.dialog.image31, handler=.$saveasImageHandler)    
    addhandlerclicked(pk.dialog.image4, handler=.$closeImageHandler)
    
    pk.dialog.notebook <<- gnotebook(closebuttons = TRUE, tearable = FALSE)
    size(pk.dialog.notebook) <<- c(getSubHeight()*0.3, getSubHeight()*0.5)
    
    rightpane = gpanedgroup(gim, pk.dialog.notebook, horizontal = FALSE)

    part.pg = gpanedgroup(leftpane, rightpane)

    add(.$window, part.pg)
  },
  
  makeBody = function(., container)
  {
    glabel(.$message, cont=container)
    if(length(.$widgetList) > 0)
    {
      tbl <- glayout(cont=container)
      ctr = 1;
      for(i in names(.$widgetList))
      {
        tmp = .$widgetList[[i]]
        FUN = tmp[[1]]
        tmp[[1]] <- NULL
        tbl[ctr,1, anchor = c(-1,-1)] = i
        tbl[ctr,2] <-
              (.$widgets[[i]] <- do.call(FUN, c(tmp, container = tbl)))
     
        ctr = ctr + 1
      }
      visible(tbl) <- TRUE
    }
  },

  makeSave = function(., container)
  {
    gseparator(cont=container)
    glabel("General Configure", cont=container)
    if(length(.$saveList) > 0)
    {
      tbl <- glayout(cont=container)
   
      ctr = 1;
      for(i in names(.$saveList))
      {
        tmp = .$saveList[[i]]
        FUN = tmp[[1]]
        tmp[[1]] <- NULL
        tbl[ctr,1, anchor = c(-1,0)] = i
        tbl[ctr,2, anchor = c(-1,0)] <-
          (.$savewd[[i]] <- do.call(FUN, c(tmp, container = tbl)))
        ctr = ctr + 1
      }
      visible(tbl) <- TRUE
    }
  },

  makeButtons = function(., container)
  {
    ## add buttons help, cancel, ok (if xxxButtonHandler is not NULL)
    gseparator(cont=container)
    bg = ggroup(cont=container)
    if(!is.null(.$helpButtonHandler))
      helpButton = gbutton("help", cont=bg,
        action = list(self=., super=.super),
        handler = .$helpButtonHandler)
     addSpace(bg, .pk$getSubHeight()/10, horizontal=TRUE)

    if(!is.null(.$cancelButtonHandler))
      cancelButton = gbutton("cancel", cont=bg,
        action = list(self=., super=.super),
        handler = .$cancelButtonHandler)
    if(!is.null(.$cleanFigureButtonHandler))
      cleanFigureButton = gbutton("Clean Figures", cont=bg,
        action = list(self=., super=.super),
        handler = .$cleanFigureButtonHandler)        
    if(!is.null(.$clearButtonHandler))
      clearButton = gbutton("clear", cont=bg,
        action = list(self=., super=.super),
        handler = .$clearButtonHandler)
    if(!is.null(.$okButtonHandler))
      okButton = gbutton("ok", cont=bg,
        action = list(self=., super=.super),
        handler = .$okButtonHandler)
  },

  ## Notice, the signature includes the initial "."
  helpButtonHandler = NULL,             # make a handler if interested
  cancelButtonHandler = NULL,           # make non-NULL handler
  clearButtonHandler = NULL,           # make non-NULL handler
  drawButtonHandler = NULL,
  saveImageHandler = NULL,
  ggobiImageHandler = NULL,
  closeImageHandler = NULL,
  okButtonHandler = NULL,
  dataframeImageHandler = NULL,

  cleanFigureButtonHandler = function(.,h,...) #1031
  {
  },
    
  saveImageHandler = function(.,h,...)
  {
  },

  closeImageHandler = function(.,h,...)
  {
      if (!is.null(ggobi_get())) close(ggobi_get())
      else
      {
          ErrorMessage("There is no ggobi instance to close!")
      }
  },

  dataframeImageHandler = function(.,h,...)
  {
      if (!is.null(ggobi_get()))
      {
          ggobi.data <- ggobi_get()[1]
          ggobi.select <- glyph_colour(ggobi.data)
          ## not equal default value: 1
          mydata <- ggobi.data[ggobi.select != 1,]
          
          child.win <- gwindow(parent=.$window)
          child.gt <-gtable(mydata, cont=child.win, expand=TRUE)
      }
      else
      {
          ErrorMessage("You have to start ggobi using second ICON first!")
      }
  },

  saveasImageHandler = function(.,h,...)   # 1110
  {
    if (!is.null(ggobi_get()))
    {
            ggobi.data <- ggobi_get()[1]
            ggobi.select <- glyph_colour(ggobi.data)
            ## not equal default value: 1
            mydata <- ggobi.data[ggobi.select != 1,]
            
            gfile("Save file",type="save", handler = function(h,...)
              {
                filename <- paste(h$file, "csv", sep=".")
                write.csv(mydata, file=filename)
                focus(.$window) <- TRUE
              })          
  
     }
     else
     {
            ErrorMessage("You have to start ggobi using second ICON first!")
     }
  
  
  },
  
  cancelButtonHandler = function(.,h,...)
  {
      dispose(.$window)
      ## others?
  },

  ## menubar
  menubarList = NULL,                   # non-null to have menubar
  menubar = NULL,
  getMenubar = function(.) return(.$menubar),
  setMenubar = function(.,lst) svalue(.$menubar) <- lst,
  ## toolbar
  toolbarList = NULL,                   # non-null to have toolbar
  toolbar = NULL,
  getToolbar = function(.) return(.$toolbar),
  setToolbar = function(.,lst) svalue(.$toolbar) <- lst,
  ## statusbar
  statusbarText = NULL,                 # non-null for statusbar
  statusbar = NULL,
  getStatusbar = function(.) return(.$statusbar),
  setStatusbar = function(.,value) svalue(.$statusbar) <- value,
  ## gwindow stuff
  window = NULL,                      # top-level gwindow
  parent.win = NULL,
  pk.dialog.notebook = NULL,

  ## properties
  message = "Basic widget",
  props = list(),                     # for storing properties of widgets
  ## for generic use
  widgetList =  list(),
  widgets = list(),
  
  saveList = list(),
  savewd = list(),
  
  # check save group figure or single figure
  figureGroup = 0,
  datagroup = 0
  
)

################################################################################
simple.okButtonHandler = function(.,h,...)
{
    
    myName <- colnames(getCurrentData())
    termDf <- data.frame(Varname=rep(NA, length(myName)), TermName=rep(NA, length(myName)))
    sapply(1:length(myName), function(i) termDf[i,] <<- c(myName[[i]], svalue(.$widgets[[myName[i]]])) )
    colnames(termDf) <- c("VarName", "TermName")
    setTerm(termDf)

    svalue(pmg.statusBar) <- "Data is loaded AND configured successfully."
    
    dispose(.$window)
}

SimpleGUI = proto(
  new = function(., message = "Simple GUI", note = NULL,...)
  {
     .$proto(message=message, note=note,...)
  },

  ## method to check if window has been drawn or destroyed
  isVisible = function(.,win = .$window)
  {
    if(!is.null(win) && is(win,"guiWidget") && isExtant(win)) return(TRUE)
    return(FALSE)
  },

  show = function(.,...)
  {
    if(.$isVisible()) return()

    .$window <- gwindow(title=.$message, parent=parent.win, height=getSubHeight()*0.7, width=getSubWidth()*0.45,...)
    g = ggroup(horizontal=FALSE, cont=.$window, use.scrollwindow = TRUE)

    ## container for body --e xpand = TRUE
    g1 = ggroup(horizontal=FALSE, container=g, expand=TRUE) # expand

    .$makeBody(container = g1)

    .$makeButtons(container = g1)

  },

  makeBody = function(., container)
  {
    gseparator(cont=container)
    glabel("General Configure", cont=container)
    if(length(.$widgetList) > 0)
    {
      tbl <- glayout(cont=container)
      ctr = 1;
      j = 1;
      tline = 1
      for(i in names(.$widgetList))
      {

        tmp = .$widgetList[[i]]
        FUN = tmp[[1]]
        tmp[[1]] <- NULL
        
        if (j > 4) j <- 1
        tbl[ctr,j, anchor=c(-1,-1)] = i
        tbl[ctr,j+1] <-
          (.$widgets[[i]] <- do.call(FUN, c(tmp, container = tbl)))
        j = j + 2

        if ((tline%%2)==0) ctr = ctr + 1
        tline <- tline + 1
      }
      visible(tbl) <- TRUE
    }
  },
  makeButtons = function(., container)
  {
    ## add buttons help, cancel, ok (if xxxButtonHandler is not NULL)
    gseparator(cont=container)
    bg0 = ggroup(cont=container, horizontal=FALSE) 
    bg = ggroup(cont=bg0) 

     addSpace(bg, getSubWidth()/10, horizontal=TRUE)

    if(!is.null(.$cancelButtonHandler))
      cancelButton = gbutton("cancel", cont=bg,
        action = list(self=., super=.super),
        handler = .$cancelButtonHandler)

    addSpace(bg, getSubWidth()/10, horizontal=TRUE)
    if(!is.null(.$okButtonHandler))
      okButton = gbutton("ok", cont=bg,
        action = list(self=., super=.super),
        handler = .$okButtonHandler)
        
    if (!is.null(.$note))  
    {
        gf.note = ggroup(cont=bg0)
        glabel(.$note, cont=gf.note)
    }        
  },

  ## Notice, the signature includes the initial "."

  cancelButtonHandler = NULL,           # make non-NULL handler
  clearButtonHandler = NULL,           # make non-NULL handler
  okButtonHandler = NULL,

  cancelButtonHandler = function(.,h,...)
  {
      dispose(.$window)
      ## others?
  },

  ## menubar
  menubarList = NULL,                   # non-null to have menubar
  menubar = NULL,
  getMenubar = function(.) return(.$menubar),
  setMenubar = function(.,lst) svalue(.$menubar) <- lst,
  ## toolbar
  toolbarList = NULL,                   # non-null to have toolbar
  toolbar = NULL,
  getToolbar = function(.) return(.$toolbar),
  setToolbar = function(.,lst) svalue(.$toolbar) <- lst,
  ## statusbar
  statusbarText = NULL,                 # non-null for statusbar
  statusbar = NULL,
  getStatusbar = function(.) return(.$statusbar),
  setStatusbar = function(.,value) svalue(.$statusbar) <- value,
  ## gwindow stuff
  window = NULL,                      # top-level gwindow
  parent.win = NULL,
  pk.dialog.notebook = NULL,

  ## properties
  message = "Basic widget",
  note = NULL,                         
  props = list(),                     # for storing properties of widgets
  ## for generic use
  widgetList =  list(),
  widgets = list(),

  saveList = list(),
  savewd = list()
  
)

################################################################################
GGobiGUI = proto(
  new = function(., message = "Basic GUI",...)
  {
     .$proto(message=message,...)
  },

  ## method to check if window has been drawn or destroyed
  isVisible = function(.,win = .$window)
  {
    if(!is.null(win) && is(win,"guiWidget") && isExtant(win)) return(TRUE)
    return(FALSE)
  },

  show = function(.,...)
  {
    ## ... passed to gwindow
    ## check if window is already there
    if(.$isVisible()) return()

    .$window <- gwindow(title=.$message, parent=parent.win, height=getSubHeight(), width=getSubWidth(),...)

    ## container for body --e xpand = TRUE
    g1 = ggroup(horizontal=FALSE, expand=FALSE) # expand
    glabel("Step 2: Map data", cont=g1)
    gf = gframe(horizontal=FALSE, expand=FALSE, cont=g1) # expand
    .$makeLeftBody(container = gf)
    .$makeLeftButtons(container = g1)
    leftpane =gpanedgroup(g1)

    g2 = ggroup(horizontal=FALSE, expand=FALSE) # expand
    glabel("Step 3: Configure plots", cont=g2)
    gg = gframe(horizontal=FALSE, expand=FALSE, cont=g2) # expand
    .$makeRightBody(container = gg)
    .$makeRightButtons(container = g2)
    rightpane =gpanedgroup(g2)
    
    part.pg = gpanedgroup(leftpane, rightpane)

    add(.$window, part.pg)
  },

  makeLeftBody = function(., container)
  {
   
      tbl <- glayout(cont=container)
      size(tbl) <- c(getSubWidth()*0.5, getSubHeight()*0.8)

      cline <- 0
      for (i in 1: length(.$data.name))
      {
          cline <- cline + 1
          tbl[cline, 1, anchor = c(-1,-1)] = "Data Name"
          tbl[cline, 2, anchor = c(-1,-1)] = .$data.name[i]
          cline <- cline + 1
          tbl[cline, 1, anchor = c(-1,-1)] = "x"

          tbl[cline, 2, anchor = c(-1,-1)] = gdroplist(items=colnames(getCurrentData(.$data.name[i])))
          cline <- cline + 1
          tbl[cline, 1, anchor = c(-1,-1)] = "y"
          tbl[cline, 2, anchor = c(-1,-1)] = gdroplist(items=colnames(getCurrentData(.$data.name[i])))
      }

  },

  makeRightBody = function(., container)
  {
    gt <- gtable(items=c(letters[1:24]), cont=container, expand=TRUE)
    size(gt) <- c(getSubWidth()*0.5, getSubHeight()*0.8)
    gt[] <- .$right.table
  },

  makeLeftButtons = function(., container)
  {
    ## add buttons help, cancel, ok (if xxxButtonHandler is not NULL)
    #gseparator(cont=container)
    bg = ggroup(cont=container)

     addSpace(bg, getSubWidth()/5, horizontal=TRUE)

    if(!is.null(.$mapButtonHandler))
      okButton = gbutton("Map data", cont=bg,
        action = list(self=., super=.super),
        handler = .$okButtonHandler)
  },
  
  makeRightButtons = function(., container)
  {
    ## add buttons help, cancel, ok (if xxxButtonHandler is not NULL)
    #gseparator(cont=container)
    bg = ggroup(cont=container)

     addSpace(bg, getSubWidth()/5, horizontal=TRUE)

    if(!is.null(.$plotButtonHandler))
      okButton = gbutton("Plot data", cont=bg,
        action = list(self=., super=.super),
        handler = .$okButtonHandler)
  },

  ## Notice, the signature includes the initial "."
  helpButtonHandler = NULL,             # make a handler if interested
  cancelButtonHandler = NULL,           # make non-NULL handler
  clearButtonHandler = NULL,           # make non-NULL handler
  drawButtonHandler = NULL,
  saveImageHandler = NULL,
  ggobiImageHandler = NULL,
  closeImageHandler = NULL,
  mapButtonHandler = function(.,h,...)
  {

  },
  plotButtonHandler = function(.,h,...)
  {

  },

  saveImageHandler = function(.,h,...)
  {
  },

  closeImageHandler = function(.,h,...)
  {
      close(ggobi_get())
  },

  cancelButtonHandler = function(.,h,...)
  {
      dispose(.$window)
      ## others?
  },

  ## gwindow stuff
  window = NULL,                      # top-level gwindow
  parent.win = NULL,

  ## properties
  message = "Basic widget",
  props = list(),                     # for storing properties of widgets
  ## for generic use
  widgetList =  list(),
  widgets = list(),

  saveList = list(),
  savewd = list(),
  
  data.name = c(),
  
  ## right panel table
  right.table = c()
)

TableGUI = proto(
  new = function(., message = "Table GUI",...)
  {
     .$proto(message=message,...)
  },

  ## method to check if window has been drawn or destroyed
  isVisible = function(.,win = .$window)
  {
    if(!is.null(win) && is(win,"guiWidget") && isExtant(win)) return(TRUE)
    return(FALSE)
  },

  show = function(.,...)
  {
    if(.$isVisible()) return()

    winH <- getSubHeight()
    winW <- getSubWidth()*0.45
    .$window <- gwindow(title=.$message, parent=parent.win, height= winH, width=winW,...)
    g = ggroup(horizontal=FALSE, cont=.$window, use.scrollwindow = TRUE)

    ## container for body --e xpand = TRUE
    g1 = ggroup(horizontal=FALSE, container=g, expand=TRUE) # expand, spacing=0

    gframe1 = gframe(text = .$frameMessage, markup = FALSE, pos = 0, cont=g1, horizontal=TRUE)
    if (is.null(.$table.data) || length(.$table.data) < 1)
    {
        ErrorMessage("You don't have proper data!")
        return(invisible(NULL))
    }
    
    .$table.widget = gtable(.$table.data, multiple = TRUE, sort.columns = 1:2,
                              expand=TRUE, cont=gframe1)
    size(.$table.widget) <- c(winW*0.5, winH*0.8)

    if (!is.null(.$checkbox.data)) .$checkbox.widget = gcheckboxgroup(.$checkbox.data, horizontal=TRUE, cont=g1)

    gbutton(text=.$buttonMessage, horizontal=FALSE, cont=g1, handler = .$okButtonHandler)

  },
  
  ## gwindow stuff
  window = NULL,                      # top-level gwindow
  parent.win = NULL,

  ## properties
  message = "Basic widget",
  frameMessage = "My frame",
  buttonMessage = "OK",
  
  table.data = NULL,
  table.widget = NULL,
  checkbox.data = NULL,
  checkbox.widget = NULL,
  okButtonHandler = NULL

)
