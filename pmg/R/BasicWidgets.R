## Some helper functions to use proto to create widgets
## use proto to create widgets

## THIS USE PRE 0.4 proto version (.super instead of super())

                                        #require(proto)
#require(gWidgets)
#options("guiToolkit"="RGtk2")

## a Trait: uppercase
## a method, prototype: start with lowercase (new, show, makeButtons, ...)

## A Trait for a BasicGUI (with window, cancel, ok buttons)
## show() creates a new window we have 5 parts:
## menubar: set menubarList to a list for gmenu
## toolbar: set toolbarList to a list for gtoolbar
## body: override with makeBody(.,container).
## If widgetlist is provided, then generates "generic widget" like widget
## buttons: override [ok|cancel|help]ButtonHandler. Set NULL to not have
##   or write makeButtons(.,container)
## statusbar: set statusbarText to get

BasicGUI = proto(
  new = function(., message = "Basic GUI",...) {
     .$proto(message=message,...)
  },
  ## method to check if window has been drawn or destroyed
  isVisible = function(.,win = .$window) {
    if(!is.null(win) && is(win,"guiWidget") && isExtant(win)) return(TRUE)
    return(FALSE)
  },
  show = function(.,...) {
    ## ... passed to gwindow
    ## check if window is already there
    if(.$isVisible()) return()    
    ## window withing pmg, write this way to give flexibility outside of pmg
    if(exists("pmgWC"))
      .$window <- pmgWC$new(title = .$message,...)
    else
      .$window <- gwindow(title = .$message, ...)
    g = ggroup(horizontal=FALSE, container=.$window, expand=TRUE)
    ## group for toolbar and menubar
    if(!is.null(.$menubarList) || !is.null(.$toolbarList)) {
      g1 = ggroup(horizontal=FALSE, container=g, expand=FALSE) 
      if(!is.null(.$menubarList)) .$menubar <- gmenu(.$menubarList, container=g1)
      if(!is.null(.$toolbarList)) .$toolbar <- gtoolbar(.$toolbarList, style="icons",container=g1)
      gseparator(container=g1)
    }
    ## container for body --e xpand = TRUE
    g1 = ggroup(horizontal=FALSE, container=g, expand=TRUE) # expand
    .$makeBody(container = g1)
    .$makeButtons(container = g)
    if(!is.null(.$statusbarText))
      .$statusbar <- gstatusbar(.$statusbarText, container=g)    
  },
  makeBody = function(., container) {
    glabel(.$message, container=container)
    if(length(.$widgetList) > 0) {
      tbl <- glayout(container=container)
      ctr = 1; 
      for(i in names(.$widgetList)) {
        tmp = .$widgetList[[i]]
        FUN = tmp[[1]]
        tmp[[1]] <- NULL
        tbl[ctr,1] = i
        tbl[ctr,2] <-
          (.$widgets[[i]] <- do.call(FUN, c(tmp, container = tbl)))
        ctr = ctr + 1
    }
      visible(tbl) <- TRUE
    }
  },
  makeButtons = function(., container) {
    ## add buttons help, cancel, ok (if xxxButtonHandler is not NULL)
    gseparator(container=container)
    bg = ggroup(container=container)
    if(!is.null(.$helpButtonHandler)) 
      helpButton = gbutton("help", container=bg,  
        action = list(self=., super=.super),
        handler = .$helpButtonHandler)
    addSpring(bg)
    ## for these we take advantage of the fact that when we call
    ## the handlers this way the "." gets passed in via the first argument
    if(!is.null(.$cancelButtonHandler))
      cancelButton = gbutton("cancel", container=bg,  
        action = list(self=., super=.super),
        handler = .$cancelButtonHandler)
    if(!is.null(.$clearButtonHandler))
      clearButton = gbutton("clear", container=bg,  
        action = list(self=., super=.super),
        handler = .$clearButtonHandler)
    if(!is.null(.$okButtonHandler)) 
      okButton = gbutton("ok", container=bg, 
        action = list(self=., super=.super),
        handler = .$okButtonHandler)
  },
  ## Notice, the signature includes the initial "."
  helpButtonHandler = NULL,             # make a handler if interested
  cancelButtonHandler = NULL,           # make non-NULL handler
  clearButtonHandler = NULL,           # make non-NULL handler
  okButtonHandler = function(.,h,...) {
    for(i in names(.$widgetList))  {
      ## store vals in props of super
#      .$.super$props[[i]] <- svalue(.$widgets[[i]]) # pre 0.4-0
     h$action$super$props[[i]] <- svalue(.$widgets[[i]])
    }
    dispose(.$window)
    },
  cancelButtonHandler = function(.,h,...) {
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
  ## properties
  message = "Basic widget",
  props = list(),                     # for storing properties of widgets
  ## for generic use
  widgetList =  list(),
  widgets = list()
  )

## Test it
##   BGTest = BasicGUI$new(message="Basic Widget Test",
##   widgetList = list(
##     edit = list(type="gedit",text="starting text"),
##     droplist = list(type = "gdroplist", items = letters),
##     slider = list(type = "gslider", value = 10),
##     radio = list(type="gradio", items = 1:3, horizontal=FALSE)
##  ))
## ## override handler so we don't set values in parent
## BGTest$okButtonHandler = function(.,handler,...) {
##   print(sapply(.$widgets,svalue)) ## or whatever else
##   dispose(.$window)
## }
## BGTest$show()  ## show the widget




## A Trait for a basic widget. To be embedded in a container
## Override the makeBody to change
BasicWidget = proto(
  new = function(., container=NULL, ...) {
    .$container = container
    ## setup widget
  },
  show = function(., ...) {
    ## show widget
    .$makeBody(container=.$container)
  },
  makeBody = function(.,container) {
    glabel("This space for rent", container=container,...)
  },
  getValue = function(.,...) {
    if(is.null(.$widget))
      return(NA)
    else if(inherits(.$widget,"proto"))
      return(.$widget$getValue(...))
    else
      return(svalue(.$widget,...))
  },
  setValues = function(.,...) {},
  widget = NULL
)



## Make some Traits for extending gtable:
## SelectItemsWithOrder: two table panes, order is clear
## SelectItemsWithSelectionOrder: one table, order by click order
## UpDownTable: widget to move items up and down a table
## orderedGtable (return with order clicked, more subtle form of
## Up and Down Table (give buttons to move up and down an element)



## A Trait for a widget that allows one to select one or more from a
## list with order. -- only vectors, not data frames
SelectItemsWithOrder = BasicWidget$proto()
SelectItemsWithOrder$new = function(., container=NULL, allItems, curItems=c(), allItemsLabel = "", curItemsLabel = "") {
  if(missing(allItems)) {
    warn("Need to call with allItems and optionally  curItems")
    return()
  }
  .$proto(container=container, allItems=allItems, curItems=curItems,
          allItemsLabel = allItemsLabel, curItemsLabel=curItemsLabel)
} 
SelectItemsWithOrder$makeBody = function(.,container,...) {
  g = ggroup(container = container)
  g1 = ggroup(horizontal=FALSE, container=g)
  glabel(.$allItemsLabel, container=g1)
  .$tbl1 = gtable(setdiff(.$allItems,.$curItems), container= g1, expand=TRUE)
  .$leftRightArrow = gimage("rarrow",dirname="stock", container=g)
  g1 = ggroup(horizontal=FALSE, container=g)
  glabel(.$curItemsLabel, container=g1)
  .$tbl2 = gtable(.$allItems, container=g1, expand=TRUE)
  .$tbl2[] <- .$curItems
  bg = ggroup(horizontal=FALSE, container=g)
  addSpace(bg,50)
  .$upArrow = gimage("uarrow", dirname="stock", container=bg)
  .$downArrow = gimage("darrow", dirname="stock", container=bg)

  ## assign widget
  .$widget <- .$tbl2

  ## add handlers
  addHandlerClicked(.$tbl1, handler = function(h,...) {
    svalue(.$leftRightArrow) <- "rarrow"
    .$leftRightArrowState = "right"
  })
  addHandlerClicked(.$tbl2, handler = function(h,...) {
    svalue(.$leftRightArrow) <- "larrow"
    .$leftRightArrowState = "left"
  })
  addHandlerClicked(.$leftRightArrow, handler = function(h,...) {
    from = .$tbl1
    to = .$tbl2
    if(.$leftRightArrowState == "left") {
      from = .$tbl2; to = .$tbl1
    }
    curSelected = svalue(from)
    if(length(curSelected) > 0) {
      from[] <- setdiff(from[],curSelected)
      toVals = to[]; toVals = toVals[!is.na(toVals)]
      to[] <- c(toVals,curSelected)
    }
  })
  addHandlerClicked(.$upArrow,  handler = function(h,...) {
    curItems = .$tbl2[]
    curSelected = svalue(.$tbl2)
    if(length(curSelected) > 0) {
      curInd = which(curSelected == curItems)
      if(curInd !=1) {
        a = curItems[curInd-1]
        .$tbl2[curInd-1] <- curSelected
        .$tbl2[curInd] <- a
        svalue(.$tbl2, index=TRUE) <- curInd - 1
      }
    }
  })
  addHandlerClicked(.$downArrow, handler = function(h,...) {
    curItems = .$tbl2[]; n<- length(curItems)
    curSelected = svalue(.$tbl2)
    if(length(curSelected) > 0) {
      curInd = which(curSelected == curItems)
      if(curInd !=n) {
        a = curItems[curInd+1]
        .$tbl2[curInd+1] <- curSelected
        .$tbl2[curInd] <- a
        svalue(.$tbl2, index=TRUE) <- curInd + 1
      }
    }
  })
}
SelectItemsWithOrder$getValue = function(.,...) {
  .$tbl2[]                              # override svalue
}

### TEST IT
## Use this to select contrasts
## allC = c('contr.helmert', 'contr.poly', 'contr.sum',
##      'contr.treatment')
## b =SelectItemsWithOrder$new(container=gwindow("test"), allItems=allC,
##   allItemsLabel = "Avail. contrasts",curItemsLabel="Selected contrasts")
## b$show()


##################################################
## A Trait for selecting from a gtable with order
## data.frames or vectors for items
## This is a more subtle ordering so that user barely notices
SelectItemsWithSelectionOrder = BasicWidget$proto()                    
SelectItemsWithSelectionOrder$new = function(.,container=NULL,items=c(),label="",chosencol=1,...) {
  .$proto(container=container, items=items, label=label, chosencol=1, value=c())
}
SelectItemsWithSelectionOrder$makeBody = function(.,container,...) {
  g = ggroup(horizontal=FALSE, container=container,...)
  glabel(.$label, container=g)
  .$widget = gtable(.$items, multiple=TRUE, chosencol=.$chosencol,
    container=g, expand=TRUE)
  addHandlerClicked(.$widget, function(h,...) {
    ## set .$value based on number set. curvalue of value
    curVals = svalue(.$widget, index=TRUE)
    if(length(curVals) == 1)
      .$value = curVals
    else if(length(curVals) > 1) {
      ## add missing to value
      .$value = c(.$value, setdiff(curVals, .$value))
    }
    ## call click handler
    .$clickedHandler(h,...)
  })
  addHandlerDoubleclick(.$widget, function(h,...) .$doubleClickHandler(h,...))
}
SelectItemsWithSelectionOrder$clickedHandler = function(.,h,...) print(.$getValue(drop=FALSE))
SelectItemsWithSelectionOrder$doubleClickHandler = function(.,h,...) {}
SelectItemsWithSelectionOrder$getValue = function(.,...) {
  chosencol = tag(.$widget,"chosencol")
  return(.$widget[.$value,chosencol,...])
}
SelectItemsWithSelectionOrder$setValues = function(.,values) .$widget[,]<-values

## ## test it
## testit = SelectItemsWithSelectionOrder$new(
##   container=gwindow("test SelectItemsWithSelectionOrder"),
##   items = mtcars, label="mtcars")
## testit$show()

###################################################
## UpDownTable. This works with data frame
UpDownTable = BasicWidget$proto()
UpDownTable$new = function(., container=NULL, items=c(),label="") {
  .$proto(container=container, items=items, label=label)
} 
UpDownTable$makeBody = function(.,container,...) {
##  g = ggroup(container= .$container)
  g = gframe(.$label, container=container, expand=TRUE)
  lg = ggroup(horizontal=FALSE, container=g, expand=TRUE)
##  glabel(.$label, container=lg)
  .$widget = gtable(.$items, container=lg, expand=TRUE)
  bg = ggroup(horizontal=FALSE, container=g)
  addSpace(bg,50)
  .$upArrow = gimage("uarrow", dirname="stock", container=bg)
  .$downArrow = gimage("darrow", dirname="stock", container=bg)

  ## add handlers
  addHandlerClicked(.$upArrow,  handler = function(h,...) {
    curItems = .$widget[,]
    curInd = svalue(.$widget, index=TRUE)
    curSelected = curItems[curInd,,drop=FALSE]
    if(!is.null(curInd)) {
      if(curInd !=1) {
        a = curItems[curInd-1,,drop=FALSE]
        .$widget[curInd-1,] <- curSelected
        .$widget[curInd,] <- a
        svalue(.$widget, index=TRUE) <- curInd - 1
      }
    }
  })
  addHandlerClicked(.$downArrow, handler = function(h,...) {
    curItems = .$widget[,]
    if(is.data.frame(curItems))
      n <- dim(curItems)[1]
    else
      n<- length(curItems)
    curInd = svalue(.$widget, index=TRUE)
    curSelected = curItems[curInd,,drop=FALSE]
    if(!is.null(curInd)) {
      if(curInd !=n) {
        a = curItems[curInd+1,,drop=FALSE]
        .$widget[curInd+1,] <- curSelected
        .$widget[curInd,] <- a
        svalue(.$widget, index=TRUE) <- curInd + 1
      }
    }
  })
  addHandlerClicked(.$widget, action=list(self=.,super=.super), handler = .$clickedHandler)
  addHandlerDoubleclick(.$widget, action=list(self=.,super=.super), handler = .$doubleClickHandler)
}
UpDownTable$getValue = function(.,...) {
  .$widget[,]                              # override svalue
}
UpDownTable$setValues = function(.,value,...) .$widget[,]<-value
UpDownTable$clickedHandler = function(.,h,...) {}
UpDownTable$doubleClickHandler = function(.,h,...) {}


## ### TEST IT
## ## Use this to select contrasts
## allC = c('contr.helmert', 'contr.poly', 'contr.sum',
##      'contr.treatment')
## b =UpDownTable$new(container=gwindow("test"), items = mtcars, label="test")
## b$clickedHandler = function(.,h,...) print(.$getValue())
## b$show()



### Test
##g = SelectItemsWithSelectionOrder$new(items = letters, container=gwindow("test"))
##g$show()

