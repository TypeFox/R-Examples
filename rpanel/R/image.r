w.image <- function(parent, filename, action=NA, mousedrag=NA, mouseup=NA, pos=NULL)
{
 widget <- w.createwidget(parent, pos, NULL, expand="true")
 widget$.type = "image"

 image <- handshake(tkimage.create, 'photo', file=filename)
 widget$.widget <- handshake(tkcanvas, parent$.handle, width=handshake(tcl, 'image', 'width', image), height=handshake(tcl, 'image', 'height', image))
 imageincanvas <- handshake(tkcreate, widget$.widget, 'image', 0, 0, image=image, anchor = 'nw')

 fdown <- function(x, y) { action(widget$.widget, as.numeric(x), as.numeric(y)) }
 fdrag <- function(x, y) { mousedrag(widget$.widget, as.numeric(x), as.numeric(y)) }
 fup <- function(x, y) { mouseup(widget$.widget, as.numeric(x), as.numeric(y)) }
 #returning the control should make it easier to act on the click - this first parameter could be extended to all widgets

 if (is.function(action)) { handshake(tkbind, widget$.widget, "<Button-1>", fdown) }
 if (is.function(mousedrag)) { handshake(tkbind, widget$.widget, "<B1-Motion>", fdrag) }
 if (is.function(mouseup)) { handshake(tkbind, widget$.widget, "<ButtonRelease-1>", fup) }

 w.appearancewidget(widget, NULL, NULL, NULL)
 invisible(widget)
}

w.line <- function(canvas, x1, y1, x2, y2, color="black", width=2, id="rpline") 
{
 rpline <- handshake(tkcreate, canvas, 'line' , x1, y1, x2, y2, fill=color, width=width)
 handshake(tkaddtag, canvas, 'rpline', 'withtag', rpline)
# 13/03/2012 this next line was missing in the previous version
 handshake(tkaddtag, canvas, id, 'withtag', rpline)
 invisible(rpline)
}

w.deleteline <- function(line) 
{
 handshake(tkdelete, line)
}

w.deletelinebyid <- function(image, id) 
{  # this method of deletion could have been used for all widget deletions
 handshake(tkdelete, image$.widget, 'withtag', id)
}

w.clearlines <- function(image) 
{
 handshake(tkdelete, image$.widget, 'withtag', 'rpline')
}

rp.image<-function(panel, filename, imagename, action=NA, mousedrag=NA, mouseup=NA, pos=NULL,
  parentname=deparse(substitute(panel)), ...) 
{
# formerly, but only in the new version, images were given numbered names as below. This fails for rp.gulls, thus the change.
#  if (is.null(deparse(substitute(imagename))))
#  {
#    imagename = paste("image", .nc(), sep="")
#  }
#  else
#  {
    imagename = deparse(substitute(imagename))
#  }

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
  
  if (is.function(action)) { fa <- function(w, x, y) { 
    panel <- rp.control.get(panelname); 
    panel <- action(panel, x ,y); 
    rp.control.put(panelname, panel)
  } }
  else { fa <- NA }  
  if (is.function(mousedrag)) { fd <- function(w, x, y) { 
    panel <- rp.control.get(panelname); 
    panel <- mousedrag(panel, x, y); rp.control.put(panelname, panel) } }
  else { fd <- NA }
  if (is.function(mouseup)) { fu <- function(w, x, y) { 
    panel <- rp.control.get(panelname); 
    panel <- mouseup(panel, x ,y); rp.control.put(panelname, panel) } }
  else { fu <- NA }

  if (rp.widget.exists(panelname, parentname)) { parent <- rp.widget.get(panelname, parentname) }
  else { parent <- panel }
  
  if (is.null(pos)) { if (length(list(...))>0) { pos <- list(...) } }
  if (is.list(pos)) { if (!is.null(pos$grid)) { parent <- rp.widget.get(panelname, pos$grid) } }
  
  widget <- w.image(parent, filename, action=fa, mousedrag=fd, mouseup=fu, pos)
  rp.widget.put(panelname, imagename, widget)
  if (.rpenv$savepanel) { rp.control.put(panelname, panel) } # put the panel back into the environment
#  hm, panelname or panel
  invisible(panelname)
}

rp.line <- function(panel, imagename, x1, y1, x2, y2, color="black", width=2, id="rpline") 
{
  if (!is.null(panel$panelname)) 
  { 
    panelname <- panel$panelname 
  } else { 
    panelname <- deparse(substitute(panel)) 
  } # this may need repeated elsewhere along with the adding this name logic
  imagename <- deparse(substitute(imagename))
  canvas = rp.widget.get(panelname, imagename)$.widget
  w.line(canvas, x1, y1, x2, y2, color, width, id)
}

rp.deleteline <- function(panel, imagename, id) 
{
  imagename <- deparse(substitute(imagename))
  if (!is.null(panel$panelname)) { panelname <- panel$panelname } else { panelname <- deparse(substitute(panel)) } # this may need repeated elsewhere along with the adding this name logic
  w.deletelinebyid(rp.widget.get(panelname, imagename), id)
}

rp.clearlines <- function(panel, imagename) 
{
  imagename <- deparse(substitute(imagename))
  if (!is.null(panel$panelname)) { panelname <- panel$panelname } else { panelname <- deparse(substitute(panel)) } # this may need repeated elsewhere along with the adding this name logic
  w.clearlines(rp.widget.get(panelname, imagename))
}

