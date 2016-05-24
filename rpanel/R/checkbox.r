rp.checkbox <- function(panel, variable, action = I, labels = NULL, names = NULL, 
                        title = NULL, initval = rep(FALSE, length(labels)), pos = NULL,
                        doaction = FALSE, foreground = NULL, background = NULL, 
                        font = NULL, parentname = deparse(substitute(panel)),
                        name = paste("checkbox", .nc(), sep=""), ...) {
                        	
  if (!exists(panel$panelname, .rpenv, inherits = FALSE))
     panelname <- deparse(substitute(panel))
  else
     panelname <- panel$panelname 

  varname <- deparse(substitute(variable))
  if (is.null(labels)) labels <- varname
 
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
    if (length(initval) == 0) initval <- FALSE 
    if (is.null(names))       names   <- labels
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
  
  widget <- w.checkbox(parent, action = f, labels = labels, names = names, title = title,
                       initval = variable, pos = pos,
                       foreground = foreground, background = background, font = font)
  rp.widget.put(panelname, name, widget)  
  # rp.widget.put(parentname, varname, widget)  
  # rp.widget.put(parentname, name, widget)  
  if (doaction) f(initval)
  if (.rpenv$savepanel) rp.control.put(panelname, panel)

  invisible(panelname)
}

w.checkbox <- function(parent, action = I, labels, names, title, 
                       initval = rep(FALSE, length(labels)), 
                       pos = NULL, foreground = NULL, background = "white", font = NULL) {
                       	
  widget       <- w.createwidget(parent, pos, background, title)
  widget$.type <- "checkgroup" 
  widget$.var  <- c()
  widget$.cb   <- list()

  f <- function(...) { 
     variable <- c()
     for (j in (1:length(labels)))
        variable[j] <- !(handshake(tclvalue, widget$.var[[j]]) == '0')
     names(variable) <- names
     action(variable) # .rpenv, variable)
  }

  for (i in (1:length(labels))) {
    if (initval[i] == TRUE)
       widget$.var[i] <- list(handshake(tclVar, '1'))
    else
       widget$.var[i] <- list(handshake(tclVar, '0'))
    cb <- w.createwidget(widget, pos = list(column = 0, row = i - 1, sticky = "news",
                  width = as.integer(handshake(.Tcl, 
                     paste('font measure systemfont "', labels[[i]], '"', sep = "") )),
                  # 06/08/2012 NOTE - why do we need to add this offset?
                  # Note that 1 does not work
                  height = 2 + as.integer(handshake(.Tcl,
                             'font metrics systemfont -linespace'))
                  ), background)
    cb$.type   <- "checkbutton"
    cb$.widget <- handshake(tkcheckbutton, widget$.handle, command = f, text = labels[[i]],
                            variable = widget$.var[[i]])
    w.appearancewidget(cb, font, foreground, background)    
    widget$.cb[i] <- list(cb)
  } 

  invisible(widget)
}

# The checkbox change function would be useful but isn't working yet.

#rp.checkbox.change <- function(panel, name, varname, value, i = 1, action = I, do = TRUE) {
#   if (is.na(charmatch("window", panel$panelname))) # if the panelname is not set then
#      panelname <- deparse(substitute(panel)) # the panel name should be the panel deparse subst'ed
#      # panel <- rp.control.get(panelname, panel) # now get the panel
#      # panel$panelname = panelname # now set the panelname properly
#      # assign(panelname, panel, envir=.rpenv) # now send back the panel
#   else 
#      panelname <- panel$panelname 
#      # panel <- rp.control.get(panelname, panel) # now get the panel
#   widget <- rp.widget.get(panelname, name)
#   w.checkbox.change(widget, value, i, action = action, do = do)
#   if (do) {
#         valexisting <- rp.var.get(panelname, varname)
#         val         <- valexisting
#         val[i]      <- value
#         rp.var.put(panelname, varname, val)
#         panel <- rp.control.get(panelname)
#         panel <- action(panel)
#         rp.control.put(panelname, panel)
#   }
#   invisible()
#}
#
#w.checkbox.change <- function(widget, value, i = 1, action = I, do = TRUE) {
#   tclvalue(widget$.var[[i]]) <- as.numeric(value)
#}
