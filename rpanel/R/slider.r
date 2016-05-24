rp.slider <- function(panel, variable, from, to, action = I, 
                      labels = NULL, names = NULL, title = NULL,
                      log = rep(FALSE, length(from)),
                      showvalue = FALSE, showvaluewidth = 4, 
                      resolution = 0, initval = from, pos = NULL, 
                      horizontal = TRUE, foreground = NULL, background = NULL, font = NULL, 
                      parentname = deparse(substitute(panel)),
                      name = paste("slider", .nc(), sep=""), ...) {

  panelname <- panel$panelname
  varname   <- deparse(substitute(variable))
  
  if (is.null(labels)) {
  	if (length(from) == 1) 
  	  labels <- varname
  	else
  	  labels <- paste(varname, 1:length(from), sep = "")
  }

  if (!rp.isnull(panelname, varname)) {
    variable <- rp.var.get(panelname, varname)
    if (is.null(names)) {
       if (!is.null(names(variable)))
          names <- names(variable)
       else
          names <- labels
    }
  }
  else {
    if (is.null(names)) names   <- labels
  	variable <- initval
  }
  names(variable)     <- names
  rp.var.put(panelname, varname, variable)
  
  if (is.null(pos) & length(list(...)) > 0) pos <- list(...)

  f <- function(val) {
    valexisting <- rp.var.get(panelname, varname)
    names(val)  <- names(valexisting)
    rp.var.put(panelname, varname, val)
    panel <- rp.control.get(panelname) 
    panel <- action(panel)
    rp.control.put(panelname, panel)
  }

  if (rp.widget.exists(panelname, parentname))
     parent <- rp.widget.get(panelname, parentname)
  else 
     parent <- panel
  if (is.list(pos) && !is.null(pos$grid))
     parent <- rp.widget.get(panelname, pos$grid)

  widget <- w.slider(parent, initval = variable, from, to, action = f, 
                     labels, names, title, log,
                     showvalue, showvaluewidth, resolution, pos, horizontal, foreground, 
                     background, font)
  rp.widget.put(panelname, name, widget)

  if (.rpenv$savepanel) rp.control.put(panelname, panel) # put the panel back into the environment

  invisible(panelname)
}

rp.slider.change <- function(panel, name, value, i = 1, do = TRUE) {
  if (!exists(panel$panelname, .rpenv, inherits = FALSE))
    panelname = deparse(substitute(panel))
  else
    panelname = panel$panelname 
  w.slider.change(rp.widget.get(panelname, name), value, i, do)
}

w.slider <- function(parent, initval, from, to, action = I,
                     labels = deparse(substitute(initval)), names = labels, title = NULL,
                     log = rep(FALSE, length(from)),
                     showvalue = FALSE, showvaluewidth = 4, resolution = 0, pos = NULL,
                     horizontal = TRUE, foreground = NULL, background = NULL, font = NULL) {
                     	
  widget                 <- w.createwidget(parent, pos, background, title)
  widget$.type           <- "sliders"
  widget$.showvalue      <- showvalue
  widget$.showvaluewidth <- showvaluewidth
  widget$.log            <- log
  widget$.var            <- c()
  if (showvalue == "widgets") {
    widget$.text <- list()   
    widget$.show <- list()   
  }  
  widget$.sl <- list()

  f <- function(...) { 
    variable <- c()
    for (j in (1:length(from))) { 
      variable[j] <- as.numeric(handshake(tclvalue, widget$.var[[j]]))
      if (log[j]) variable[j] <- exp(variable[j])
      if (showvalue == "widgets")
        w.text.change(widget$.show[[j]], signif(variable[j], showvaluewidth))
      }
    names(variable) <- names
    if (length(from) == 1)
      action(as.numeric(variable))
    else
      action(variable)
  }

  widget$.f <- f

  orient <- if (horizontal) "horizontal" else "vertical"

  for (i in 1:length(from)) { 
    if (log[i]) {
      wto <- log(to[i])
      wfrom <- log(from[i])
      initv <- log(initval[i])
    }
    else {
      wto <- to[i]
      wfrom <- from[i];
      initv <- initval[i]
    }
    widget$.var[i] <- list(handshake(tclVar, signif(initv, showvaluewidth)))

    if (showvalue != "widgets") {

      # Over-ride! This stops the incorrect 'unlogged' value being shown.
      if ((showvalue == TRUE) && log[i]) showvalue <- FALSE

      if (horizontal == TRUE) { 
        pos=list(
          column=0,
          row=i-1,
          sticky="news",
          # 2 times is to give height for the prompt (above) and slider (below)
          height= 2 * as.integer(handshake(.Tcl, 'font metrics systemfont -linespace')),
          cweight=1,
          rweight=1
        ) 
      }
      else { 
        pos=list(
          column=i-1,
          row=0,
          sticky="news",
          width = 2 * as.integer(handshake(.Tcl, 'font metrics systemfont -linespace')),
          cweight=1,
          rweight=1
        )
      }

      sl <- w.createwidget(widget, pos=pos, background)
      sl$.type <- "slider"
      if (length(from) == 1)
         sl$.widget <- handshake(tkscale, widget$.handle, from = wfrom, to = wto,
                                 showvalue = showvalue, orient = orient,
                                 resolution = resolution, variable = widget$.var[[i]])
      else
         sl$.widget <- handshake(tkscale, widget$.handle, from = wfrom, to = wto,
                                 showvalue = showvalue, orient = orient, label = labels[i],
                                 resolution = resolution, variable = widget$.var[[i]])

      #      sl$.widget <- handshake(tkscale, widget$.handle, from=wfrom, to=wto, showvalue=showvalue, orient=orient, resolution=resolution, variable=widget$.var[[i]], length=200)
      handshake(tkbind, sl$.widget, "<B1-Motion>", f)
      w.appearancewidget(sl, font, foreground, background)
      widget$.sl[i] <- list(sl)
    }    

    else {
      if (horizontal) pos <- "left" else pos <- "top"
      text <- w.text(widget, title[i], NA, 
        pos=list(
          column=0, 
          row=i-1, 
          sticky="w", 
          cweight=1,
          height=as.integer(handshake(.Tcl, 'font metrics systemfont -linespace')),
          width = as.integer(handshake(.Tcl, paste('font measure systemfont "', title[i], '"', sep="") ))
        ), foreground, background, font)
      if (horizontal==TRUE) { 
        pos=list(
          column=0,
          row=i-1,
          sticky="news",
          height= as.integer(handshake(.Tcl, 'font metrics systemfont -linespace')),
          cweight=100          
        ) 
      }
      else { 
        pos=list(
          column=i-1,
          row=0,
          sticky="news",
          width=as.integer(handshake(.Tcl, 'font metrics systemfont -linespace')),
          cweight=100
        )
      }
#      if (horizontal==TRUE) { pos=list(column=1,row=i-1,sticky="ew",cweight=100) }
#      else { pos=list(column=i-1,row=1,sticky="ew",cweight=100) }
      sl <- w.createwidget(widget, pos=pos, background)
      sl$.type <- "slider"

     # 29/11/13 if added to allow for displaying name if needed
     if (length(from) == 1)
        sl$.widget <- handshake(tkscale, widget$.handle, from = wfrom, to = wto,
                                showvalue = FALSE, orient = orient,
                                resolution = resolution, variable = widget$.var[[i]])
     else
        sl$.widget <- handshake(tkscale, widget$.handle, from = wfrom, to = wto,
                                showvalue = FALSE, orient = orient, label = labels[i],
                                resolution = resolution, variable = widget$.var[[i]])

      sl$.widget <- handshake(tkscale, widget$.handle, from=wfrom, to=wto,
                              showvalue=FALSE, orient=orient, resolution=resolution,
                              variable=widget$.var[[i]])

      handshake(tkbind, sl$.widget, "<B1-Motion>", f)
      #  The statement below applies to the slider only

      w.appearancewidget(sl, font, foreground, background)
      show <- w.text(widget, signif(initval, showvaluewidth), NA,
                     pos=list(column=2, row=i-1, sticky="e", cweight=1),
                     foreground, background, font, width=showvaluewidth+1)
      widget$.text[i] <- list(text)
      widget$.sl[i]   <- list(sl)
      widget$.show[i] <- list(show)
    }  
  }

  invisible(widget)
}

w.slider.change <- function(widget, value, i=1, do=TRUE) {

  if (widget$.log[i])
  	tclvalue(widget$.var[[i]]) <- log(value)
  else
    tclvalue(widget$.var[[i]]) <- value

  if (widget$.showvalue=="widgets") { w.text.change(widget$.show[[i]], signif(value, widget$.showvaluewidth)) }

  if (do) { widget$.f(widget$.var) }

}
