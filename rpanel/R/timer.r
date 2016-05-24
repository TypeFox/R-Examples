w.timer <- function(microseconds, action, where) {
  dothis <- function() {
    action()
    w.timer(microseconds, action, where)
  }

  if (eval(parse(text=where), .rpenv))
     timer <- handshake(tcl, "after", microseconds, dothis)
  else
     timer <- NULL

  invisible(timer)
}

rp.timer <- function(panel, microseconds, action, where) {

  if (!exists(panel$panelname, .rpenv, inherits = FALSE))
    panelname <- deparse(substitute(panel))
  else 
    panelname <- panel$panelname 

  dothis <- function() {
    panel <- rp.control.get(panelname)
    panel <- action(panel)
    rp.control.put(panelname, panel)
    # **************************************************
    #   need to return value back into various places  *
    # **************************************************
    if (where(panel))
       timer <- handshake(tcl, "after", microseconds, dothis)
    else
       timer <- NULL
  }

  if (where(panel))
     timer <- handshake(tcl, "after", microseconds, dothis)
  else
     timer <- NULL

  invisible(timer)
}

rp.do <- function(panel, action, x = NA, y = NA) {
  if (!exists(panel$panelname, .rpenv, inherits = FALSE))
    panelname <- deparse(substitute(panel))
  else 
    panelname <- panel$panelname 
  panel <- rp.control.get(panelname)
  if (is.na(x) & is.na(y))
     panel <- action(panel)
  else
     panel <- action(panel, x, y)
  rp.control.put(panelname, panel)
  invisible(panelname)
}

rp.block <- function(panel) {
  if (!exists(panel$panelname, .rpenv, inherits = FALSE))
    panelname <- deparse(substitute(panel))
  else 
    panelname <- panel$panelname

  open <- eval(parse(text=paste("!is.null(",panelname, ")", sep="")), envir=.rpenv)
  while(open) {
    Sys.sleep(0.01)
    # The commented code was added in 1.1-2 but seems to cause a problem
    # open <- exists(panelname)
    # if (open) {
      open <- exists(panelname, envir= .rpenv)
      if (open)
        open <- eval(parse(text=paste("!is.null(",panelname, ")", sep="")), envir=.rpenv)
    # }
  }

  invisible()
}
