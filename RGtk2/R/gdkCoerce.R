as.GdkAtom <-
function(x)
{
  if (is.integer(x))
    x <- as.numeric(x)
  else if (!inherits(x, "GdkAtom") && !is.numeric(x))
    x <- as.character(x)
  x
}

# either 'pixel' or ('red', 'green', 'blue') must exist (may be combined) 
as.GdkColor <-
function(x)
{
	if (is.character(x))
		return(gdkColorParse(x)$color)
	
  if (length(x) == 1) # only one field, must be 'pixel'
    fields <- "pixel"
  else { # otherwise, must have 'rgb' and possibly 'pixel'
    fields <- c("red", "green", "blue")
    if (length(x) > 3)
      fields <- c("pixel", fields)
  }
	
	x <- as.struct(x, "GdkColor", fields)
    
	if (length(x) > 3)
		x[[1]] <- as.numeric(x[[1]])
	else x[[1]] <- as.integer(x[[1]])
    
	x[-1] <- sapply(x[-1], as.integer) 

    return(x)
}

as.GdkRectangle <-
function(x)
{
    x <- as.struct(x, "GdkRectangle", c("x", "y", "width", "height"))
    x[[1]] <- as.integer(x[[1]])
    x[[2]] <- as.integer(x[[2]])
    x[[3]] <- as.integer(x[[3]])
    x[[4]] <- as.integer(x[[4]])

    return(x)
}

as.GdkTrapezoid <-
function(x)
{
  x <- as.struct("GdkTrapezoid", c("y1", "x11", "x21", "y2", "x12", "x22"))
  x[[1]] <- as.numeric(x[[1]])
  x[[2]] <- as.numeric(x[[2]])
  x[[3]] <- as.numeric(x[[3]])
  x[[4]] <- as.numeric(x[[4]])
  x[[5]] <- as.numeric(x[[5]])
  x[[6]] <- as.numeric(x[[6]])
  
  return(x)
}

as.GdkSpan <-
function(x)
{
    x <- as.struct(x, "GdkSpan", c("x", "y", "width"))
    x[[1]] <- as.integer(x[[1]])
    x[[2]] <- as.integer(x[[2]])
    x[[3]] <- as.integer(x[[3]])
    
    return(x)
}

as.GdkRgbCmap <-
function(x)
{
	x <- as.numeric(x)
	class(x) <- "GdkRgbCmap"
	x
}
as.GdkKeymapKey <-
function(x)
{
	x <- as.struct(x, "GdkKeymapKey", c("keycode", "group", "level"))
	x[[1]] <- as.numeric(x[[1]])
	x[[2]] <- as.integer(x[[2]])
	x[[3]] <- as.integer(x[[3]])
	return(x)
}
as.GdkGCValues <-
function(x)
{
	x <- as.struct(x, "GdkGCValues", c("foreground", "background", "font", "function", "fill", "tile", "stipple", 
		"clip.mask", "subwindow.mode", "ts.x.origin", "ts.y.origin", "clip.x.origin", "clip.y.origin", 
		"graphics.exposures", "line.width", "line.style", "cap.style", "join.style"))
	
	if (!is.null(x[[1]])) x[[1]] <- as.GdkColor(x[[1]])
	if (!is.null(x[[2]])) x[[2]] <- as.GdkColor(x[[2]])
	if (!is.null(x[[3]])) x[[3]] <- checkPtrType(x[[3]], "GdkFont")
	if (!is.null(x[[4]])) x[[4]] <- as.function(x[[4]])
	# GdkFill	    fill;
	if (!is.null(x[[6]])) x[[6]] <- checkPtrType(x[[6]], "GdkPixmap")
	if (!is.null(x[[7]])) x[[7]] <- checkPtrType(x[[7]], "GdkPixmap")
	if (!is.null(x[[8]])) x[[8]] <- checkPtrType(x[[8]], "GdkPixmap")
	# GdkSubwindowMode  subwindow_mode;
	if (!is.null(x[[10]])) x[[10]] <- as.integer(x[[10]])
	if (!is.null(x[[11]])) x[[11]] <- as.integer(x[[11]])
	if (!is.null(x[[12]])) x[[12]] <- as.integer(x[[12]])
	if (!is.null(x[[13]])) x[[13]] <- as.integer(x[[13]])
	if (!is.null(x[[14]])) x[[14]] <- as.integer(x[[14]])
	if (!is.null(x[[15]])) x[[15]] <- as.integer(x[[15]])
	# GdkLineStyle	    line_style;
	# GdkCapStyle	    cap_style;
	# GdkJoinStyle	    join_style;
	return(x)
}
as.GdkGeometry <-
function(x)
{
  x <- as.struct(x, "GdkGeometry", c("min.width", "min.height", "max.width", "max.height",
  	"base.width", "base.height", "width.inc", "height.inc", "min.aspect", "max.aspect", "win.gravity"))
  
	if (!is.null(x[[1]])) {
		x[[1]] <-  as.integer(x[[1]])
		x[[2]] <-  as.integer(x[[2]])
	}
	if (!is.null(x[[3]])) {
		x[[3]] <- as.integer(x[[3]])
		x[[4]] <- as.integer(x[[4]])
	}
	if (!is.null(x[[4]])) {
		 x[[5]] <- as.integer(x[[5]])
		 x[[6]] <- as.integer(x[[6]])
	}
	if (!is.null(x[[7]])) {
		x[[7]] <- as.integer(x[[7]])
		x[[8]] <- as.integer(x[[8]])
	}
	if (!is.null(x[[9]])) {
		x[[9]] <- as.numeric(x[[9]])
		x[[10]] <- as.numeric(x[[10]])
	}
  
  # GdkGravity win_gravity;
  
  return(x)
}
as.GdkWindowAttr <-
function(x)
{
	x <- as.struct(x, "GdkWindowAttr", c("title", "event.mask", "x", "y", "width", "height", "wclass",
		"visual", "colormap", "window.type", "cursor", "wmclass.name", "wmclass.class", "override.redirect"))
	if (!is.null(x[[1]])) x[[1]] <- as.character(x[[1]])
	x[[2]] <- as.integer(x[[2]])
	if (!is.null(x[[3]]))
    x[[3]] <- as.integer(x[[3]])
	if (!is.null(x[[4]]))
    x[[4]] <- as.integer(x[[4]])
	x[[5]] <- as.integer(x[[5]])
	x[[6]] <- as.integer(x[[6]])
	# wclass
	if (!is.null(x[[8]])) x[[8]] <- checkPtrType(x[[8]], "GdkVisual")
	if (!is.null(x[[9]])) x[[9]] <- checkPtrType(x[[9]], "GdkColormap")
  # window.type
	if (!is.null(x[[11]])) x[[11]] <- checkPtrType(x[[11]], "GdkCursor")
	if (!is.null(x[[12]])) {
		x[[12]] <- as.character(x[[12]])
		x[[13]] <- as.character(x[[13]])
	}
	if (!is.null(x[[14]])) x[[14]] <- as.logical(x[[14]])
	
	return(x)
}

as.GdkNativeWindow <-
function(x)
{
	class(x) <- "GdkNativeWindow"
	x
}

as.GdkSegment <-
function(x)
{
	x <- as.struct(x, "GdkSegment", c("x1", "y1", "x2", "y2"))
	x[[1]] <- as.integer(x[[1]])
	x[[2]] <- as.integer(x[[2]])
	x[[3]] <- as.integer(x[[3]])
	x[[4]] <- as.integer(x[[4]])
	
	return(x)
}

as.GdkPoint <-
function(x)
{
	x <- as.struct(x, "GdkPoint", c("x", "y"))
	x[[1]] <- as.integer(x[[1]])
	x[[2]] <- as.integer(x[[2]])
	return(x)
}
