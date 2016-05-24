if(!exists('.virtuals')) .virtuals <- new.env()
assign("GdkDisplay", c("get_display_name", "get_n_screens", "get_screen", "get_default_screen", "closed"), .virtuals)
assign("GdkDisplayManager", c("display_opened"), .virtuals)
assign("GdkDrawable", c("create_gc", "draw_rectangle", "draw_arc", "draw_polygon", "draw_text", "draw_text_wc", "draw_drawable", "draw_points", "draw_segments", "draw_lines", "draw_glyphs", "draw_image", "get_depth", "get_size", "set_colormap", "get_colormap", "get_visual", "get_screen", "get_image", "get_clip_region", "get_visible_region", "get_composite_drawable", "draw_pixbuf", "draw_glyphs_transformed", "draw_trapezoids", "ref_cairo_surface"), .virtuals)
assign("GdkGC", c("get_values", "set_values", "set_dashes"), .virtuals)
assign("GdkKeymap", c("direction_changed", "keys_changed"), .virtuals)
assign("GdkScreen", c("size_changed", "composited_changed"), .virtuals)
assign("GdkPixbufAnimation", c("is_static_image", "get_static_image", "get_size", "get_iter"), .virtuals)
assign("GdkPixbufAnimationIter", c("get_delay_time", "get_pixbuf", "on_currently_loading_frame", "advance"), .virtuals)
assign("GdkPixbufLoader", c("size_prepared", "area_prepared", "area_updated", "closed"), .virtuals)


gdkDisplayClassGetDisplayName <-
function(object.class, object)
{
  checkPtrType(object.class, "GdkDisplayClass")
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_class_get_display_name", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gdkDisplayClassGetNScreens <-
function(object.class, object)
{
  checkPtrType(object.class, "GdkDisplayClass")
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_class_get_n_screens", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gdkDisplayClassGetScreen <-
function(object.class, object, screen.num)
{
  checkPtrType(object.class, "GdkDisplayClass")
  checkPtrType(object, "GdkDisplay")
  screen.num <- as.integer(screen.num)

  w <- .RGtkCall("S_gdk_display_class_get_screen", object.class, object, screen.num, PACKAGE = "RGtk2")

  return(w)
}

gdkDisplayClassGetDefaultScreen <-
function(object.class, object)
{
  checkPtrType(object.class, "GdkDisplayClass")
  checkPtrType(object, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_class_get_default_screen", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gdkDisplayClassClosed <-
function(object.class, object, is.error)
{
  checkPtrType(object.class, "GdkDisplayClass")
  checkPtrType(object, "GdkDisplay")
  is.error <- as.logical(is.error)

  w <- .RGtkCall("S_gdk_display_class_closed", object.class, object, is.error, PACKAGE = "RGtk2")

  return(invisible(w))
}

gdkDisplayManagerClassDisplayOpened <-
function(object.class, object, display)
{
  checkPtrType(object.class, "GdkDisplayManagerClass")
  checkPtrType(object, "GdkDisplayManager")
  checkPtrType(display, "GdkDisplay")

  w <- .RGtkCall("S_gdk_display_manager_class_display_opened", object.class, object, display, PACKAGE = "RGtk2")

  return(invisible(w))
}

gdkDrawableClassCreateGc <-
function(object.class, object, values)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")
  values <- as.GdkGCValues(values)

  w <- .RGtkCall("S_gdk_drawable_class_create_gc", object.class, object, values, PACKAGE = "RGtk2")

  return(w)
}

gdkDrawableClassDrawRectangle <-
function(object.class, object, gc, filled, x, y, width, height)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  filled <- as.logical(filled)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gdk_drawable_class_draw_rectangle", object.class, object, gc, filled, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
}

gdkDrawableClassDrawArc <-
function(object.class, object, gc, filled, x, y, width, height, angle1, angle2)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  filled <- as.logical(filled)
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)
  angle1 <- as.integer(angle1)
  angle2 <- as.integer(angle2)

  w <- .RGtkCall("S_gdk_drawable_class_draw_arc", object.class, object, gc, filled, x, y, width, height, angle1, angle2, PACKAGE = "RGtk2")

  return(invisible(w))
}

gdkDrawableClassDrawPolygon <-
function(object.class, object, gc, filled, points, npoints)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  filled <- as.logical(filled)
  points <- as.GdkPoint(points)
  npoints <- as.integer(npoints)

  w <- .RGtkCall("S_gdk_drawable_class_draw_polygon", object.class, object, gc, filled, points, npoints, PACKAGE = "RGtk2")

  return(invisible(w))
}

gdkDrawableClassDrawText <-
function(object.class, object, font, gc, x, y, text, text.length)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")
  checkPtrType(font, "GdkFont")
  checkPtrType(gc, "GdkGC")
  x <- as.integer(x)
  y <- as.integer(y)
  text <- as.character(text)
  text.length <- as.integer(text.length)

  w <- .RGtkCall("S_gdk_drawable_class_draw_text", object.class, object, font, gc, x, y, text, text.length, PACKAGE = "RGtk2")

  return(invisible(w))
}

gdkDrawableClassDrawTextWc <-
function(object.class, object, font, gc, x, text)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")
  checkPtrType(font, "GdkFont")
  checkPtrType(gc, "GdkGC")
  x <- as.integer(x)
  text <- as.list(as.numeric(text))

  w <- .RGtkCall("S_gdk_drawable_class_draw_text_wc", object.class, object, font, gc, x, text, PACKAGE = "RGtk2")

  return(w)
}

gdkDrawableClassDrawDrawable <-
function(object.class, object, gc, src, xsrc, ysrc, xdest, ydest, width, height)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  checkPtrType(src, "GdkDrawable")
  xsrc <- as.integer(xsrc)
  ysrc <- as.integer(ysrc)
  xdest <- as.integer(xdest)
  ydest <- as.integer(ydest)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gdk_drawable_class_draw_drawable", object.class, object, gc, src, xsrc, ysrc, xdest, ydest, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
}

gdkDrawableClassDrawPoints <-
function(object.class, object, gc, points)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  points <- lapply(points, function(x) { x <- as.GdkPoint(x); x })

  w <- .RGtkCall("S_gdk_drawable_class_draw_points", object.class, object, gc, points, PACKAGE = "RGtk2")

  return(invisible(w))
}

gdkDrawableClassDrawSegments <-
function(object.class, object, gc, segs)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  segs <- lapply(segs, function(x) { x <- as.GdkSegment(x); x })

  w <- .RGtkCall("S_gdk_drawable_class_draw_segments", object.class, object, gc, segs, PACKAGE = "RGtk2")

  return(invisible(w))
}

gdkDrawableClassDrawLines <-
function(object.class, object, gc, points)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  points <- lapply(points, function(x) { x <- as.GdkPoint(x); x })

  w <- .RGtkCall("S_gdk_drawable_class_draw_lines", object.class, object, gc, points, PACKAGE = "RGtk2")

  return(invisible(w))
}

gdkDrawableClassDrawGlyphs <-
function(object.class, object, gc, font, x, y, glyphs)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  checkPtrType(font, "PangoFont")
  x <- as.integer(x)
  y <- as.integer(y)
  checkPtrType(glyphs, "PangoGlyphString")

  w <- .RGtkCall("S_gdk_drawable_class_draw_glyphs", object.class, object, gc, font, x, y, glyphs, PACKAGE = "RGtk2")

  return(invisible(w))
}

gdkDrawableClassDrawImage <-
function(object.class, object, gc, image, xsrc, ysrc, xdest, ydest, width, height)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  checkPtrType(image, "GdkImage")
  xsrc <- as.integer(xsrc)
  ysrc <- as.integer(ysrc)
  xdest <- as.integer(xdest)
  ydest <- as.integer(ydest)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gdk_drawable_class_draw_image", object.class, object, gc, image, xsrc, ysrc, xdest, ydest, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
}

gdkDrawableClassGetDepth <-
function(object.class, object)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")

  w <- .RGtkCall("S_gdk_drawable_class_get_depth", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gdkDrawableClassGetSize <-
function(object.class, object)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")

  w <- .RGtkCall("S_gdk_drawable_class_get_size", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gdkDrawableClassSetColormap <-
function(object.class, object, cmap)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")
  checkPtrType(cmap, "GdkColormap")

  w <- .RGtkCall("S_gdk_drawable_class_set_colormap", object.class, object, cmap, PACKAGE = "RGtk2")

  return(invisible(w))
}

gdkDrawableClassGetColormap <-
function(object.class, object)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")

  w <- .RGtkCall("S_gdk_drawable_class_get_colormap", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gdkDrawableClassGetVisual <-
function(object.class, object)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")

  w <- .RGtkCall("S_gdk_drawable_class_get_visual", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gdkDrawableClassGetScreen <-
function(object.class, object)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")

  w <- .RGtkCall("S_gdk_drawable_class_get_screen", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gdkDrawableClassGetImage <-
function(object.class, object, x, y, width, height)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gdk_drawable_class_get_image", object.class, object, x, y, width, height, PACKAGE = "RGtk2")

  return(w)
}

gdkDrawableClassGetClipRegion <-
function(object.class, object)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")

  w <- .RGtkCall("S_gdk_drawable_class_get_clip_region", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gdkDrawableClassGetVisibleRegion <-
function(object.class, object)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")

  w <- .RGtkCall("S_gdk_drawable_class_get_visible_region", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gdkDrawableClassGetCompositeDrawable <-
function(object.class, object, x, y, width, height)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gdk_drawable_class_get_composite_drawable", object.class, object, x, y, width, height, PACKAGE = "RGtk2")

  return(w)
}

gdkDrawableClassDrawPixbuf <-
function(object.class, object, gc, pixbuf, src.x, src.y, dest.x, dest.y, width, height, dither, x.dither, y.dither)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  checkPtrType(pixbuf, "GdkPixbuf")
  src.x <- as.integer(src.x)
  src.y <- as.integer(src.y)
  dest.x <- as.integer(dest.x)
  dest.y <- as.integer(dest.y)
  width <- as.integer(width)
  height <- as.integer(height)
  
  x.dither <- as.integer(x.dither)
  y.dither <- as.integer(y.dither)

  w <- .RGtkCall("S_gdk_drawable_class_draw_pixbuf", object.class, object, gc, pixbuf, src.x, src.y, dest.x, dest.y, width, height, dither, x.dither, y.dither, PACKAGE = "RGtk2")

  return(invisible(w))
}

gdkDrawableClassDrawGlyphsTransformed <-
function(object.class, object, gc, matrix, font, x, y, glyphs)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  checkPtrType(matrix, "PangoMatrix")
  checkPtrType(font, "PangoFont")
  x <- as.integer(x)
  y <- as.integer(y)
  checkPtrType(glyphs, "PangoGlyphString")

  w <- .RGtkCall("S_gdk_drawable_class_draw_glyphs_transformed", object.class, object, gc, matrix, font, x, y, glyphs, PACKAGE = "RGtk2")

  return(invisible(w))
}

gdkDrawableClassDrawTrapezoids <-
function(object.class, object, gc, trapezoids)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")
  checkPtrType(gc, "GdkGC")
  trapezoids <- lapply(trapezoids, function(x) { x <- as.GdkTrapezoid(x); x })

  w <- .RGtkCall("S_gdk_drawable_class_draw_trapezoids", object.class, object, gc, trapezoids, PACKAGE = "RGtk2")

  return(invisible(w))
}

gdkDrawableClassRefCairoSurface <-
function(object.class, object)
{
  checkPtrType(object.class, "GdkDrawableClass")
  checkPtrType(object, "GdkDrawable")

  w <- .RGtkCall("S_gdk_drawable_class_ref_cairo_surface", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gdkGCclassGetValues <-
function(object.class, object)
{
  checkPtrType(object.class, "GdkGCClass")
  checkPtrType(object, "GdkGC")

  w <- .RGtkCall("S_gdk_gcclass_get_values", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gdkGCclassSetValues <-
function(object.class, object, values)
{
  checkPtrType(object.class, "GdkGCClass")
  checkPtrType(object, "GdkGC")
  values <- as.GdkGCValues(values)

  w <- .RGtkCall("S_gdk_gcclass_set_values", object.class, object, values, PACKAGE = "RGtk2")

  return(invisible(w))
}

gdkGCclassSetDashes <-
function(object.class, object, dash.list)
{
  checkPtrType(object.class, "GdkGCClass")
  checkPtrType(object, "GdkGC")
  dash.list <- as.list(as.raw(dash.list))

  w <- .RGtkCall("S_gdk_gcclass_set_dashes", object.class, object, dash.list, PACKAGE = "RGtk2")

  return(w)
}

gdkKeymapClassDirectionChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GdkKeymapClass")
  checkPtrType(object, "GdkKeymap")

  w <- .RGtkCall("S_gdk_keymap_class_direction_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gdkKeymapClassKeysChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GdkKeymapClass")
  checkPtrType(object, "GdkKeymap")

  w <- .RGtkCall("S_gdk_keymap_class_keys_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gdkScreenClassSizeChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GdkScreenClass")
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_class_size_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gdkPixbufAnimationClassIsStaticImage <-
function(object.class, object)
{
  checkPtrType(object.class, "GdkPixbufAnimationClass")
  checkPtrType(object, "GdkPixbufAnimation")

  w <- .RGtkCall("S_gdk_pixbuf_animation_class_is_static_image", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gdkPixbufAnimationClassGetStaticImage <-
function(object.class, object)
{
  checkPtrType(object.class, "GdkPixbufAnimationClass")
  checkPtrType(object, "GdkPixbufAnimation")

  w <- .RGtkCall("S_gdk_pixbuf_animation_class_get_static_image", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gdkPixbufAnimationClassGetSize <-
function(object.class, object)
{
  checkPtrType(object.class, "GdkPixbufAnimationClass")
  checkPtrType(object, "GdkPixbufAnimation")

  w <- .RGtkCall("S_gdk_pixbuf_animation_class_get_size", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gdkPixbufAnimationClassGetIter <-
function(object.class, object, start.time)
{
  checkPtrType(object.class, "GdkPixbufAnimationClass")
  checkPtrType(object, "GdkPixbufAnimation")
  start.time <- as.GTimeVal(start.time)

  w <- .RGtkCall("S_gdk_pixbuf_animation_class_get_iter", object.class, object, start.time, PACKAGE = "RGtk2")

  return(w)
}

gdkPixbufAnimationIterClassGetDelayTime <-
function(object.class, object)
{
  checkPtrType(object.class, "GdkPixbufAnimationIterClass")
  checkPtrType(object, "GdkPixbufAnimationIter")

  w <- .RGtkCall("S_gdk_pixbuf_animation_iter_class_get_delay_time", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gdkPixbufAnimationIterClassGetPixbuf <-
function(object.class, object)
{
  checkPtrType(object.class, "GdkPixbufAnimationIterClass")
  checkPtrType(object, "GdkPixbufAnimationIter")

  w <- .RGtkCall("S_gdk_pixbuf_animation_iter_class_get_pixbuf", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gdkPixbufAnimationIterClassOnCurrentlyLoadingFrame <-
function(object.class, object)
{
  checkPtrType(object.class, "GdkPixbufAnimationIterClass")
  checkPtrType(object, "GdkPixbufAnimationIter")

  w <- .RGtkCall("S_gdk_pixbuf_animation_iter_class_on_currently_loading_frame", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

gdkPixbufAnimationIterClassAdvance <-
function(object.class, object, current.time)
{
  checkPtrType(object.class, "GdkPixbufAnimationIterClass")
  checkPtrType(object, "GdkPixbufAnimationIter")
  current.time <- as.GTimeVal(current.time)

  w <- .RGtkCall("S_gdk_pixbuf_animation_iter_class_advance", object.class, object, current.time, PACKAGE = "RGtk2")

  return(w)
}

gdkPixbufLoaderClassSizePrepared <-
function(object.class, object, width, height)
{
  checkPtrType(object.class, "GdkPixbufLoaderClass")
  checkPtrType(object, "GdkPixbufLoader")
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gdk_pixbuf_loader_class_size_prepared", object.class, object, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
}

gdkPixbufLoaderClassAreaPrepared <-
function(object.class, object)
{
  checkPtrType(object.class, "GdkPixbufLoaderClass")
  checkPtrType(object, "GdkPixbufLoader")

  w <- .RGtkCall("S_gdk_pixbuf_loader_class_area_prepared", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gdkPixbufLoaderClassAreaUpdated <-
function(object.class, object, x, y, width, height)
{
  checkPtrType(object.class, "GdkPixbufLoaderClass")
  checkPtrType(object, "GdkPixbufLoader")
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_gdk_pixbuf_loader_class_area_updated", object.class, object, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
}

gdkPixbufLoaderClassClosed <-
function(object.class, object)
{
  checkPtrType(object.class, "GdkPixbufLoaderClass")
  checkPtrType(object, "GdkPixbufLoader")

  w <- .RGtkCall("S_gdk_pixbuf_loader_class_closed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

gdkScreenClassCompositedChanged <-
function(object.class, object)
{
  checkPtrType(object.class, "GdkScreenClass")
  checkPtrType(object, "GdkScreen")

  w <- .RGtkCall("S_gdk_screen_class_composited_changed", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}