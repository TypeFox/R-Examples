w.grid <- function(parent, pos=NULL, background=NULL)
{
  widget <- w.createwidget(parent, pos, background, expand="true")
  widget$.type = "grid"
  invisible(widget)
}

rp.grid <- function(panel, name=paste("grid", .nc(), sep=""), pos=NULL, background=NULL, 
  parentname=deparse(substitute(panel)), ...) 
{

  if (!exists(panel$panelname, .rpenv, inherits = FALSE)) # if the panelname is not set then
  { 
    panelname = deparse(substitute(panel)) # the panel name should be the panel deparse subst'ed
# 13/03/2012 these lines are not commented out in previous version
#    panel <- rp.control.get(panelname, panel) # now get the panel
#    panel$panelname = panelname # now set the panelname properly
#    assign(panelname, panel, envir=.rpenv) # now send back the panel
  } 
  else 
  { 
    panelname = panel$panelname 
# 13/03/2012 these lines are not commented out in previous version
#    panel <- rp.control.get(panelname, panel) # now get the panel
  }
  
  if (rp.widget.exists(panelname, parentname)) { parent <- rp.widget.get(panelname, parentname) }
  else { parent <- panel }
  if (is.null(pos)) { if (length(list(...))>0) { pos <- list(...) } }
  if (is.list(pos)) { if (!is.null(pos$grid)) { parent <- rp.widget.get(panelname, pos$grid) } }
  
# 06/08/2012 duff width and height added
  widget <- w.grid(parent, pos, background)
  
  rp.widget.put(panelname, name, widget)  
  if (.rpenv$savepanel) { rp.control.put(panelname, panel) } # put the panel back into the environment
  invisible(panelname)
}
