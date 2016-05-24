## Store windows for gWidgets

## the windowcollector is a place to organize windows within pmg
## methods
## new -- returns a new window
## getWindow(ID) -- returns the window with ID
## delete(win,[ID]) -- delete that window
## register(win) -- register window. Done by new()
## show() -- show the table with the windows. Double clicking an entry raises the window. SHould have a way to delete the window.


winCollector = BasicGUI$new("message"="Open windows")
winCollector$makeBody = function(.,container) {
  g = ggroup(horizontal = FALSE, container=container, expand=TRUE)
  glabel("window list", container=g)
  .$tbl <- gtable(.$summary(), chosencol=2,container=g, expand=TRUE)
  ## add stuff, handlers buttons
  ## double click -- raise
  addHandlerDoubleclick(.$tbl, function(h,...) {
    ID = svalue(h$obj)
    if(length(ID) == 0) return(TRUE)
    w <- .$getWindow(ID)
    focus(w) <- TRUE
  })
}
winCollector$makeButtons = function(.,container) {
  bg = ggroup(container=container)
  addSpring(bg)
  gbutton("cancel",container=bg, handler = function(h,...) dispose(.$window))
  addSpace(bg,10)
  gbutton("Raise", container=bg, handler = function(h,...) {
    ID = svalue(.$tbl)
    if(length(ID) == 0) return(TRUE)    
    w = .$getWindow(ID)
    focus(w) <- TRUE
  })
  gbutton("Delete window", container=bg, handler= function(h,...) {
    ID = svalue(.$tbl)
    if(length(ID) == 0) return(TRUE)    
    w = .$getWindow(ID)
    dispose(w)
  })
}
winCollector$updateBody = function(.) {
  ## check that it is visible
  if(.$isVisible())
    .$tbl[,] <- .$summary()
}

winCollector$ctr = 0                    # for the ID
winCollector$list = list()              # stores windows
winCollector$register = function(.,win) {    # register a window
  if(!.$isVisible(win)) return(NA)
  
  if(!is.null(tag(win,"wcID"))) {
    ## already added
    return(NA)
  }

  .$ctr = .$ctr + 1
  ID = as.character(.$ctr )
  tag(win,"wcID") <- ID
  .$list[[ID]] <- win
  .$updateBody()
  return(ID)
}
## delete from list, dispose is separate
winCollector$delete = function(.,win,ID=NULL) { # delete window from list
  if(is.null(ID))
    ID = tag(win,"wcID")
  if(is.null(ID)) {
    cat("Window not among list\n")
    return(FALSE)
  }
  win = .$list[[ID]]
  .$list[[ID]] <- NULL

  if(.$isVisible(win)) dispose(win)

  .$updateBody()  
  return(TRUE)
}
winCollector$summary = function(.) { ## return df with window names
  if(length(.$list) > 0) {
    d = data.frame(title = sapply(.$list,svalue),
      ID = sapply(.$list,function(o) tag(o,"wcID")),
      stringsAsFactors=FALSE
      )
    d = d[sapply(d[,2],function(ID) .$isVisible(.$getWindow(ID))),]
  } else {
    d = data.frame(title = c(""), ID = "", stringsAsFactors=FALSE)
  }
  return(d)
}

## new window, added to 
winCollector$new = function(.,...) {
  w <- gwindow(...)
  ID <- .$register(w)
  addHandlerUnrealize(w,action=ID, handler = function(h,...) {
    .$delete(ID=h$action)
  })
  addhandlerdestroy(w, action=ID, handler = function(h,...) {
    .$delete(ID=h$action)
  })
    
  return(w)                             # return window
}
## get window From ID
winCollector$getWindow = function(.,ID) {
  if(length(ID) == 0 || is.na(ID) || is.null(ID) )
    return(NA)
  .$list[[ID]]
}

## closeAll
winCollector$closeAll = function(.) {
  ID = tag(.$window,"wcID")
  d = .$summary
  sapply(d$ID, function(i) .$delete(i))
}

## give a shorter name
pmgWC = winCollector
