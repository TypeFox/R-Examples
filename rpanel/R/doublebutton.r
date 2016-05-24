w.doublebutton <- function(parent, step, title, action = I, parentvarname, initval, range = c(NA, NA), 
                           log = FALSE, showvalue = FALSE, showvaluewidth = 4,
                           repeatinterval = 100, repeatdelay = 100, pos = "left", 
                           foreground = NULL, background = NULL, font = NULL) {
  widget <- list()
  widget$.type <- "doublebutton"
  widget$.showvalue <- showvalue
  # varname <- paste("dbvar", .nc(), sep = "")
  varname <- parentvarname
  assign(varname, initval, .rpenv)
  
  fchange <- function(op) { # this returns the desired function
    function(...) {
      val <- assign(varname, 
                eval(parse(text = paste(varname, op, 
                           as.character(step))), .rpenv), .rpenv)
      if (!is.na(range[1])) val <- max(range[1], val)
      if (!is.na(range[2])) val <- min(range[2], val)
      assign(varname, val, .rpenv)
      if (showvalue)
         w.text.change(widget$.show, signif(val, showvaluewidth))
      action(val)
    }
  }

  finc <- fchange(if (log) "*" else "+")
  fdec <- fchange(if (log) "/" else "-")  

  widget$.cont <- w.createwidget(parent = parent, pos = pos, background = background)
#  widget$.text <- w.text(parent=widget$.cont, text=title, pos=list(row=0, column=0, cweight=10), foreground=foreground, background=background, font=font)
#  widget$.dec <- w.button(parent=widget$.cont, title="-", action=fdec, repeatdelay=repeatdelay, repeatinterval=repeatinterval, pos=list(row=0, column=1, cweight=1), foreground=foreground, background=background, font=font)
#  if (showvalue) { widget$.show <- w.text(parent=widget$.cont, text=signif(initval, showvaluewidth), pos=list(row=0, column=2, cweight=1), foreground=foreground, background=background, font=font, width=showvaluewidth+1) }
#  widget$.inc <- w.button(parent=widget$.cont, title="+", action=finc, repeatdelay=repeatdelay, repeatinterval=repeatinterval, pos=list(row=0, column=3, cweight=1), foreground=foreground, background=background, font=font)

  widget$.dec <- w.button(parent = widget$.cont, title = "-", action = fdec,
                          repeatdelay = repeatdelay, repeatinterval = repeatinterval,
                          pos = "left", foreground = foreground, background = background, 
                          font=font)
  if (showvalue) 
     widget$.show <- w.text(parent = widget$.cont, text = signif(initval, showvaluewidth),
                            pos = "left", foreground = foreground, background = background, 
                            font = font, width = showvaluewidth + 1)
  widget$.inc <- w.button(parent = widget$.cont, title = "+", action = finc,
                          repeatdelay = repeatdelay, repeatinterval = repeatinterval,
                          pos = "left", foreground = foreground, background = background,
                          font = font)
  widget$.text <- w.text(parent = widget$.cont, text = title, pos = "left",
                         foreground = foreground, background = background, font = font)

  invisible(widget)
}

rp.doublebutton <- function(panel, variable, step, title = deparse(substitute(variable)), 
                            action = I, initval = NULL, range = c(NA, NA), 
                            log = FALSE, showvalue = FALSE, showvaluewidth = 4, 
                            repeatinterval = 100, repeatdelay = 100, pos = NULL, 
                            foreground = NULL, background = NULL, font = NULL, 
                            parentname = deparse(substitute(panel)),
                            name=paste("doublebutton", .nc(), sep = ""), ...) {
# The following line appears to be redundant, but is not. There is a problem with deparse substitute and rtcltk where the wrong value is taken into tcltk. Printing or assigning the title gets round this problem.
  title <- paste(title)

  if (!exists(panel$panelname, .rpenv, inherits = FALSE)) { # if the panelname is not set then 
     panelname = deparse(substitute(panel)) # the panel name should be the panel deparse subst'ed
# 13/03/2012 these lines are not commented out in previous version
#    panel <- rp.control.get(panelname, panel) # now get the panel
#    panel$panelname = panelname # now set the panelname properly
#    assign(panelname, panel, envir=.rpenv) # now send back the panel
  } 
  else  { 
    panelname = panel$panelname 
# 13/03/2012 these lines are not commented out in previous version
#    panel <- rp.control.get(panelname, panel) # now get the panel
  }
  
  varname <- deparse(substitute(variable))
  if (!rp.isnull(panelname, varname))
     variable <- rp.var.get(panelname, varname)
  else {
     variable <- initval
     rp.var.put(panelname, varname, variable)
  }
  
  if (is.null(pos) & length(list(...)) > 0) pos <- list(...)
  f <- function(val) {
    if (!is.na(range[1])) val <- max(range[1], val)
    if (!is.na(range[2])) val <- min(range[2], val)
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
  
  widget <- w.doublebutton(parent, step, title, action = f, paste(panelname, "$", varname, sep=""), initval = variable, 
                           range, log, showvalue, showvaluewidth, repeatinterval,
                           repeatdelay, pos, foreground, background, font) 

  rp.widget.put(panelname, name, widget)
  if (.rpenv$savepanel) rp.control.put(panelname, panel) # put the panel back into the environment
  invisible(panelname)
}
