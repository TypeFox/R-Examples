rp.notebook <- function(panel, tabs, tabnames = tabs, width = 600, height = 400,
                        pos = NULL, foreground = NULL, 
                        background = "lightgray", font = NULL,
                        parentname = deparse(substitute(panel)),
                        name = paste("notebook", .nc(), sep=""), ...)  {
  if (is.na(charmatch("window", panel$panelname)))
    panelname <- deparse(substitute(panel))
  else 
    panelname <- panel$panelname 

  if (rp.widget.exists(panelname, parentname))
    parent <- rp.widget.get(panelname, parentname)
  else
    parent <- panel
  
  if (is.null(pos) && (length(list(...)) > 0)) pos <- list(...)
  if (is.list(pos) && !is.null(pos$grid)) parent <- rp.widget.get(panelname, pos$grid)
  
  widget <- w.notebook(parent, width, height, pos, foreground, background, font)
  for (i in (1:length(tabs))) { 
    tab <- w.notebook.add(widget, tabs[i])
    rp.widget.put(panelname, gsub(' ', '_', tabnames[i]), tab)  
  }
  rp.widget.put(panelname, name, widget)  
  if (.rpenv$savepanel) rp.control.put(panelname, panel)
  
  invisible(panelname)
}

rp.notebook.raise <- function(panel, parentname, label) {
  if (is.na(charmatch("window", panel$panelname))) 
    panelname <- deparse(substitute(panel))
  else 
    panelname <- panel$panelname 
  w.notebook.raise(rp.widget.get(panelname, parentname), label)
}

w.notebook <- function(parent, width=NULL, height=NULL, pos=NULL, foreground=NULL, background="lightgray", font=NULL) {
  widget <- w.createwidget(parent, pos, background)
  widget$.type = "notebook"  
  handshake(.Tcl, 'package require BWidget')
  widget$.widget <- handshake(tkwidget, parent$.handle, "NoteBook")
  if ( (!is.null(width)) && (!is.null(height)) )
     handshake(tkconfigure, widget$.widget, width=width, height=height,
                  background=background)
  w.appearancewidget(widget, font=font, foreground=foreground, background=background)
  invisible(widget)
}

w.notebook.add <- function(parent, label) {
  tabpage <- list()
  page <- handshake(tkinsert, parent$.widget, "end", gsub(" ", "_", label), "-text", label)
  pagewin <- handshake(.Tk.newwin, page)
  tabpage$.handle <- handshake(tkframe, pagewin)
  handshake(tkgrid, tabpage$.handle, sticky="ew")
  invisible(tabpage)
}

w.notebook.raise <- function(parent, label)
  invisible(handshake(tcl, parent$.widget, "raise", gsub(" ", "_", label)))

