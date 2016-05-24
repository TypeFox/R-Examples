# For all queries to do with tcl/tk the best resource is http://wiki.tcl.tk/

handshake <- function(fun, ...) {
# this is a vital function for R and tcltk communication these days as R now 
# gets out of sync with tcltk due to its speed, thus tcltk does not behave 
# properly
  tclvalue <- NULL
  starttime <- proc.time()
  tclvalue <- fun(...)
  while (is.null(tclvalue)) {
  # while ( (is.null(tclvalue)) && (as.integer(proc.time()-starttime)[3] < 2) ){
    Sys.sleep(0.05)
  }
  # if (as.integer(proc.time()-starttime)[3] >= 2) { 
  #   stop("The connection to tcltk has failed.") 
  # }
  invisible(tclvalue)
}

handshakereverse <- function(fun, ...) {
# this is a vital function for R and tcltk communication these days as R now 
# gets out of sync with tcltk due to its speed, thus tcltk does not behave 
# properly
# this is used where no value is returned from the function - ie NULL
  tclvalue <- 1
  starttime <- proc.time()
  tclvalue <- fun(...)
  while (!is.null(tclvalue)) {
  # while ( (!is.null(tclvalue)) && (as.integer(proc.time()-starttime)[3] < 2) ){
    Sys.sleep(0.05)
  }
  # if (as.integer(proc.time()-starttime)[3] >= 2) { 
  #   stop("The connection to tcltk has failed.") 
  # }
  invisible(tclvalue)
}

handshake(.Tcl, 'package require BWidget')
.rpenv <- new.env()

.onAttach <- function(library, pkg) {
# First function run on opening the package
  packageStartupMessage("Package `rpanel', version 1.1-3: type help(rpanel) for summary information")
  assign("counter", 0 , envir=.rpenv)
  assign("getpanel", TRUE, envir=.rpenv)  # get the panel from the environment at the start of a function
  assign("setparent", TRUE, envir=.rpenv) # set the panel.panelname to the deparse-substitute of the panel
#  assign("setpanel", FALSE, envir=.rpenv)  # immediately save the panel back to the environment
  assign("savepanel", FALSE, envir=.rpenv) # save the panel back to the environment at the end of the function
  invisible()
}

.nc <- function()
  assign("counter", .rpenv$counter+1, envir=.rpenv)

rp.setup <- function(getpanel=.rpenv$getpanel, setparent=.rpenv$setparent) {
      # rp.setup <- function(getpanel=.rpenv$getpanel, setparent=.rpenv$setparent,
      # setpanel=.rpenv$setpanel, savepanel=.rpenv$savepanel)
   assign("getpanel",  getpanel,  envir=.rpenv)
   assign("setparent", setparent, envir=.rpenv)
   # assign("setpanel", setpanel, envir=.rpenv)
   assign("savepanel", savepanel, envir=.rpenv)
}

rp.settings <- function() {
  print(paste("getpanel is ", .rpenv$getpanel))
  print(paste("setparent is ", .rpenv$setparent))
  # print(paste("setpanel is ", .rpenv$setpanel))
  print(paste("savepanel is ", .rpenv$savepanel))
}

rp.panelname <- function() {
   # this function is retained for backword compatability only
   paste('.rpanel', .nc(), sep="")     
}

rp.env <- function()
  .rpenv

rp <- function()
  .rpenv

w.assign <- function(...) {
  params <- list(...)
  for(i in 1:length(params)) { 
    if (!is.null(names(params[i]))) { 
      assign(names(params[i]), params[[i]], .rpenv) 
    } 
  }  
}

# Is this function ever used?  It doesn't look right because of deparse(substitute()).

rp.assign <- function(panel, ...) {
  panelname <- deparse(substitute(panel))
  params <- list(...)
  for(i in 1:length(params)) { 
    if (!is.null(names(params[i]))) {
      rp.var.put(panelname, names(params[i]), params[[i]]) 
    } 
  }  
}

w.setfont <- function(widget, font)
  if (!is.null(font)) handshake(tkconfigure, widget, font=font)

w.setforeground <- function(widget, foreground)
  if (!is.null(foreground)) handshake(tkconfigure, widget, foreground=foreground)

w.setbackground <- function(widget, background)
  if (!is.null(background)) handshake(tkconfigure, widget, background=background)

w.createcontainer <- function(parent, pos, background, title = NULL, expand = "false", 
                              tkrplottype = FALSE, oneonly = TRUE) {

  # parent can be window, grid, notebooktab
  # place all widgets into a container - a frame
  
  if (!is.null(title)) {
     if (oneonly)
        container <- handshake(tkwidget, parent$.handle, "labelframe", text = title,
                               borderwidth = 2, padx = 0, pady = 0)
     else
        container <- handshake(tkwidget, parent$.handle, "labelframe", text = title)
  }
  else {

    if (is.list(pos) && (length(pos) > 1) &&
                ((!is.null(pos$width)) || (!is.null(pos$height)))) { 
      if (is.null(pos$width))
        container <- handshake(tkframe, parent$.handle, height = pos$height)
      else if (is.null(pos$height))
        container <- handshake(tkframe, parent$.handle, width = pos$width)
      else
        container <- handshake(tkframe, parent$.handle, width = pos$width,
                                 height = pos$height)

      # 06/08/2012 why done here and below? Is this going to be a problem if used?
      # From http://wiki.tcl.tk/9908
      # The grid geometry manager normally computes how large a master must be to just
      # exactly meet the needs of its slaves, and it sets the requested width and height
      # of the master to these dimensions. This causes geometry information to propagate
      # up through a window hierarchy to a top-level window so that the entire sub-tree
      # sizes itself to fit the needs of the leaf windows. However, the grid propagate
      # command may be used to turn off propagation for one or more masters. If propagation
      # is disabled then grid will not set the requested width and height of the master
      # window. This may be useful if, for example, you wish for a master window to have
      # a fixed size that you specify.

      # handshake(tkgrid.propagate, container, 0)
      pos$width  <- NULL 
      pos$height <- NULL
      # set to null to avoid messing up the next part
    }
    else
      container <- handshake(tkframe, parent$.handle, padx = 0, pady = 0)
  }

  if (!is.null(background)) w.setbackground(container, background)

  # 29/11/13 ppx added to options
  ppx <- if (oneonly) 0 else 2
  ppy <- if (oneonly) 0 else 2
  if ((is.null(pos)) || (length(pos) == 0) )
    handshake(tkpack, container, expand = expand, fill = "both",  padx=ppx, pady=ppy)
  else if (length(pos) == 1)
     handshake(tkpack, container, expand = expand, side = pos[[1]], padx=ppx, pady=ppy)
  else if (is.numeric(pos))
     handshake(tkplace, container,
               x = as.integer(pos[1]), y = as.integer(pos[2]),
               w = as.integer(pos[3]), h = as.integer(pos[4]))
  else {
    # 06/08 what is this? propogate seems to be causing a problem
    #handshake(tkgrid.propagate, container, 0); 
    if (is.null(pos$sticky))     pos$sticky     <- "news"
    if (is.null(pos$rowspan))    pos$rowspan    <- 1
    if (is.null(pos$columnspan)) pos$columnspan <- 1
    handshake(tkgrid, container, row = pos$row, column = pos$column,
              sticky = pos$sticky,
              rowspan = pos$rowspan, columnspan = pos$columnspan,
              ipadx = 1 - as.numeric(tkrplottype),
              ipady = 1 - as.numeric(tkrplottype))

    # 06/08 note that default weights of 1 will not do!  Nor will 0 (ie do not resize)
    if (!is.null(pos$rweight))
      handshake(tkgrid.rowconfigure, parent$.handle, pos$row, weight = pos$rweight)
    if (!is.null(pos$cweight))
      handshake(tkgrid.columnconfigure, parent$.handle, pos$column,
                weight = pos$cweight)
  }

  invisible(container)
}

w.layout <- function(widget, container, scr=NULL) {
  # Note - the "in" must be in quotes due to the fact that this is a reserved word in R.
  if (!is.null(scr)) { 
    handshake(tkgrid, widget$.widget, row = 0, column = 0, sticky = "ew",
              "in" = container, rowspan = 1, columnspan = 1) 
    handshake(tkgrid, scr, row=0, column=1, sticky="nws",
              "in" = container, rowspan = 1, columnspan = 1) 
  }
  else { 
    handshake(tkgrid, widget$.widget, row = 0, column = 0, sticky = "ew",
              "in" = container, rowspan = 1, columnspan = 1) 
    # Below was added 21/08/2012 - suspect row configure and column configure 
    # needed to size button to use full size of grid cell.
    # So - sticky may be irrelevant in this w.layout function since this
    # is already handled in the container function.
    # See this address for an explanation
    # http://www.wellho.net/mouth/1335_Expanding-a-grid-Tcl-Tk.html
    if (widget$.type %in% c("button", "slider")) {
       handshake(tkgrid.rowconfigure,    container, 0, weight=1)
       handshake(tkgrid.columnconfigure, container, 0, weight=1)    
    } 
  }
}

w.createwidget <- function(parent, pos = NULL, background = NULL, title = NULL, 
                           expand = "false", tkrplottype = FALSE) {
  widget <- list()
  widget$.handle <- w.createcontainer(parent, pos, background, title, expand, tkrplottype)
  widget$.parent <- parent$.handle
  invisible(widget)
}

w.appearancewidget <- function(widget, font, foreground, background, scr=NULL) {
  w.setfont(widget$.widget, font)
  w.setforeground(widget$.widget, foreground)
  w.setbackground(widget$.widget, background)
  if ((widget$.type != 'menu') && (widget$.type != 'table')) {
  	# if (is.numeric(widget$pos)) {
  	#   tkplace(widget$.widget, x = widget$pos[1], y = widget$pos[2],
  	#                           w = widget$pos[3], h = widget$pos[4])
  	#   }
  	# else
      w.layout(widget, widget$.handle, scr) 
  }
}

rp.widget.dispose <- function(panel, name) {
  if (!exists(panel$panelname, .rpenv, inherits = FALSE))
    panelname <- deparse(substitute(panel))
  else 
    panelname <- panel$panelname

  widget <- eval(parse(text = paste(panelname, ".", name, sep = "")), envir = .rpenv)
  # widget <- rp.widget.get(panelname, name)
  
  if (widget$.type == "doublebutton") { 
    handshakereverse(tkdestroy, widget$.text$.widget)
    handshakereverse(tkdestroy, widget$.inc$.widget)
    handshakereverse(tkdestroy, widget$.dec$.widget)
    handshakereverse(tkdestroy, widget$.text$.handle) 
    handshakereverse(tkdestroy, widget$.inc$.handle) 
    handshakereverse(tkdestroy, widget$.dec$.handle) 
    if (widget$.showvalue) { 
      handshakereverse(tkdestroy, widget$.show$.widget) 
      handshakereverse(tkdestroy, widget$.show$.handle) 
    }
  }

  else if (widget$.type == "sliders") {
    for (i in 1:length(widget$.sl)) {
      handshakereverse(tkdestroy, widget$.sl[[i]]$.widget)
      handshakereverse(tkdestroy, widget$.sl[[i]]$.handle)
      if (widget$.showvalue == "widgets") {
        handshakereverse(tkdestroy, widget$.text[[i]]$.widget)
        handshakereverse(tkdestroy, widget$.text[[i]]$.handle)
        handshakereverse(tkdestroy, widget$.show[[i]]$.widget)
        handshakereverse(tkdestroy, widget$.show[[i]]$.handle)
      }
    } 
    handshakereverse(tkdestroy, widget$.handle)       
  }  

  else if ( (widget$.type == "checkgroup") || (widget$.type == "radiogroup") ) {
    for (i in 1:length(widget$.cb)) {
      handshakereverse(tkdestroy, widget$.cb[[i]]$.widget)
      handshakereverse(tkdestroy, widget$.cb[[i]]$.handle)
    }
    handshakereverse(tkdestroy, widget$.handle)       
  }  

  else if (widget$.type == "textentrys") {
    for (i in 1:length(widget$.label)) {
      handshakereverse(tkdestroy, widget$.label[[i]]$.widget)
      handshakereverse(tkdestroy, widget$.label[[i]]$.handle)
      handshakereverse(tkdestroy, widget$.entry[[i]]$.widget)
      handshakereverse(tkdestroy, widget$.entry[[i]]$.handle)
    }
    handshakereverse(tkdestroy, widget$.handle)       
  }  

  else if (widget$.type == "combobox") {
    handshakereverse(tkdestroy, widget$.label[[1]]$.widget)
    handshakereverse(tkdestroy, widget$.label[[1]]$.handle)
    handshakereverse(tkdestroy, widget$.combo[[1]]$.widget)
    handshakereverse(tkdestroy, widget$.combo[[1]]$.handle)
    handshakereverse(tkdestroy, widget$.handle)       
  }  

  else {
  	# print(ls(.rpenv))
  	# print(get("window1.b1", envir = .rpenv))
  	# tkdestroy.temp <- function(win) {
  	   # tcl("destroy", win)
       # ID <- .Tk.ID(win)
       # print(ID)
       # env <- get("parent", envir = win$env)$env
  	   # print(ls(env))
  	   # print("ID")
       # print(get(ID, envir = env, inherits = FALSE))
       # if (exists(ID, envir = env, inherits = FALSE)) 
           # rm(list = ID, envir = env)
    # }
    
    # tcl("destroy", widget$.widget)
    # tcl("destroy", widget$.handle)
    handshakereverse(tkdestroy, widget$.widget) 
    handshakereverse(tkdestroy, widget$.handle)
  }
  rm(list = paste(panelname, ".", name, sep = ""), envir = .rpenv)

}

rp.var.get <- function(panelname, name) {
   # note you must send the panelname, not the panel as deparse substitute
   # only works one layer in, not two!
   pnm <- if (is.null(panelname)) "" else paste(panelname, "$", sep = "")
   eval(parse(text = paste(pnm, name, sep = "")), envir = .rpenv)
}

rp.var.put.new <- function(panel, name, val) {
   assign(paste(panel$panelname, "$", name, sep = ""), val, envir = .rpenv)
}

rp.var.put <- function(panelname, name, val, labels = NULL) {
   # note you must send the panelname, not the panel as deparse substitute
   # only works one layer in, not two!
   pnm <- if (is.null(panelname)) "" else paste(panelname, "$", sep = "")
   if ((is.character(val)) && (length(val) == 1))
      eval(parse(text=paste(pnm, name, " <- '", val, "'", sep="")),
                  envir=.rpenv)
   else {
      if ((!is.character(val)) && (length(val) == 1))
         eval(parse(text=paste(pnm, name, " <- ", val, sep="")),
                   envir=.rpenv)
      else { 
         if ((is.null(names(val))) && (is.null(labels))) {  
            eval(parse(text=paste(pnm, name, " <- ", list(val), sep="")),
                   envir=.rpenv) 
         }
         else { 
            for (j in 1:length(val)) {
               eval(parse(text=paste(pnm, name, "[", j, "] <- ",
                     val[j], sep="")), envir=.rpenv) 
               if (is.null(labels))
                  eval(parse(text=paste("names(", pnm, name,
                     ")[", j, "] <- '", names(val)[j], "'", sep="")), envir=.rpenv) 
               else {
                  eval(parse(text=paste("names(", pnm, name,
                     ")[", j, "] <- '", labels[j], "'", sep="")), envir=.rpenv) 
               }     
            }
         }
      }
   }
}

rp.matrix.put <- function(panelname, name, val, ncol, nrow) {
  eval(parse(text=paste(panelname, "$", name, " <- matrix(ncol=", ncol, ",nrow=", nrow, ",data=",
             list(val), ")", sep="")), envir=.rpenv) 
}

rp.isnull <- function(panelname, name) { 
  eval(parse(text=paste("is.null(", panelname, "$", name, ")", sep="")), envir=.rpenv)
}

rp.widget.exists <- function(panelname, name) { 
  eval(parse(text=paste("exists('", panelname, ".", name, "')", sep="")), envir=.rpenv)
}

rp.widget.get <- function(panelname, name) {
  eval(parse(text = paste(panelname, ".", name, sep = "")), envir = .rpenv)
}

rp.widget.put <- function(panelname, name, val) {
  assign(paste(panelname, ".", name, sep=""), val, envir=.rpenv)
}

rp.control.get <- function(panelname, panel = NULL) {
  if (exists(panelname, envir=.rpenv))
    panel <- eval(parse(text=panelname), envir=.rpenv)
  else
    assign(panelname, panel, envir=.rpenv)
  panel
}

rp.control.put <- function(panelname, panel) {
  assign(panelname, panel, envir=.rpenv)    
}

rp.screenresolution <- function() {
  width <- as.numeric(strsplit(handshake(tclvalue, handshake(tkwm.maxsize, '.'))," ")[[1]])[1]
  height <- as.numeric(strsplit(handshake(tclvalue, handshake(tkwm.maxsize, '.'))," ")[[1]])[2]
  return(list(width=width, height=height))
}


if(getRversion() >= "2.15.1")
   utils::globalVariables(c("aniso.angle", "aniso.ratio", "band", "case", "case.showing",
                "country", "curve.showing", "data.range", "Display", "display", "display.checks",
                "display.options", "fac.showing", "firth.true", "fitted.showing", "fplot",
                "gulls.image", "gx", "gy", "Influence", "intercept", "key", "lambda", "location",
                "model11", "model12", "model01", "model02", "model13", "model14", "model03",
                "model04", "mururoa.true", "n", "ngrid", "npts", "Nugget", "pars", "phi",
                "plot1", "plot1a", "plot2", "points.only", "prob", "prop", "pSill",
                "random.alignment", "Range", "ranges", "savepanel", "slope", "sp.type", "ssize",
                "stan", "str.type", "stype", "theta", "trend.setting", "tsb", "tsw",
                "var1", "var2", "var3", "vgm.checks", "year.ind"))
