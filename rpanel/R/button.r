rp.button <- function(panel, action, title=deparse(substitute(action)), 
                      repeatdelay = 0, repeatinterval = 0, quitbutton = FALSE, 
                      pos = NULL, foreground = NULL, background = NULL, font = NULL,
                      parentname = deparse(substitute(panel)), 
                      name = paste("button", .nc(), sep = ""), ...) {

   if (!exists(panel$panelname, .rpenv, inherits = FALSE))
     panelname <- deparse(substitute(panel))
   else
     panelname <- panel$panelname
  
   if (is.null(pos) && (length(list(...)) > 0)) pos <- list(...)

   f <- function() {
     panel <- rp.control.get(panelname)
     panel <- action(panel)
     rp.control.put(panelname, panel)
     if (quitbutton) {
       rp.control.dispose(panel)
       # if (exists(paste(panelname,"$.handle", sep = ""), envir = .rpenv)) 
       #    eval(parse(text=paste("try(tkdestroy(", panelname, "$.handle))",
       #    sep="")), envir=.rpenv)
     }
   }

   if (rp.widget.exists(panelname, parentname)) 
      parent <- rp.widget.get(panelname, parentname) 
   else 
      parent <- panel
   if (is.list(pos) && (!is.null(pos$grid)))
      parent <- rp.widget.get(panelname, pos$grid)

   widget <- w.button(parent, f, title, repeatdelay, repeatinterval, pos,
                      foreground, background, font)
   rp.widget.put(panelname, name, widget)

   if (.rpenv$savepanel) rp.control.put(panelname, panel)
   invisible(panelname)
}

w.button <- function(parent, action = I, title = deparse(substitute(action)), 
                     repeatdelay = 0, repeatinterval = 0, pos = NULL, 
                     foreground = NULL, background = NULL, font = NULL) {
   widget         <- w.createwidget(parent, pos, background)
   widget$.type   <- "button"
   f              <- function(...) action()
   widget$.widget <- handshake(tkbutton, parent$.handle, text=title, command=f, 
                               repeatdelay=repeatdelay, repeatinterval=repeatinterval)
   w.appearancewidget(widget, font, foreground, background)
   invisible(widget)
}
