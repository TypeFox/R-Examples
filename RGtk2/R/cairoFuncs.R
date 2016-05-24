
cairoVersion <-
function()
{
  

  w <- .RGtkCall("S_cairo_version", PACKAGE = "RGtk2")

  return(w)
} 


cairoVersionString <-
function()
{
  

  w <- .RGtkCall("S_cairo_version_string", PACKAGE = "RGtk2")

  return(w)
} 


cairoCreate <-
function(target)
{
  checkPtrType(target, "CairoSurface")

  w <- .RGtkCall("S_cairo_create", target, PACKAGE = "RGtk2")

  return(w)
} 


cairoSave <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_save", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoRestore <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_restore", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoSetOperator <-
function(cr, op)
{
  checkPtrType(cr, "Cairo")
  

  w <- .RGtkCall("S_cairo_set_operator", cr, op, PACKAGE = "RGtk2")

  return(w)
} 


cairoSetSource <-
function(cr, source)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(source, "CairoPattern")

  w <- .RGtkCall("S_cairo_set_source", cr, source, PACKAGE = "RGtk2")

  return(w)
} 


cairoSetSourceRgb <-
function(cr, red, green, blue)
{
  checkPtrType(cr, "Cairo")
  red <- as.numeric(red)
  green <- as.numeric(green)
  blue <- as.numeric(blue)

  w <- .RGtkCall("S_cairo_set_source_rgb", cr, red, green, blue, PACKAGE = "RGtk2")

  return(w)
} 


cairoSetSourceRgba <-
function(cr, red, green, blue, alpha)
{
  checkPtrType(cr, "Cairo")
  red <- as.numeric(red)
  green <- as.numeric(green)
  blue <- as.numeric(blue)
  alpha <- as.numeric(alpha)

  w <- .RGtkCall("S_cairo_set_source_rgba", cr, red, green, blue, alpha, PACKAGE = "RGtk2")

  return(w)
} 


cairoSetSourceSurface <-
function(cr, surface, x, y)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(surface, "CairoSurface")
  x <- as.numeric(x)
  y <- as.numeric(y)

  w <- .RGtkCall("S_cairo_set_source_surface", cr, surface, x, y, PACKAGE = "RGtk2")

  return(w)
} 


cairoSetTolerance <-
function(cr, tolerance)
{
  checkPtrType(cr, "Cairo")
  tolerance <- as.numeric(tolerance)

  w <- .RGtkCall("S_cairo_set_tolerance", cr, tolerance, PACKAGE = "RGtk2")

  return(w)
} 


cairoSetAntialias <-
function(cr, antialias)
{
  checkPtrType(cr, "Cairo")
  

  w <- .RGtkCall("S_cairo_set_antialias", cr, antialias, PACKAGE = "RGtk2")

  return(w)
} 


cairoSetFillRule <-
function(cr, fill.rule)
{
  checkPtrType(cr, "Cairo")
  

  w <- .RGtkCall("S_cairo_set_fill_rule", cr, fill.rule, PACKAGE = "RGtk2")

  return(w)
} 


cairoSetLineWidth <-
function(cr, width)
{
  checkPtrType(cr, "Cairo")
  width <- as.numeric(width)

  w <- .RGtkCall("S_cairo_set_line_width", cr, width, PACKAGE = "RGtk2")

  return(w)
} 


cairoSetLineCap <-
function(cr, line.cap)
{
  checkPtrType(cr, "Cairo")
  

  w <- .RGtkCall("S_cairo_set_line_cap", cr, line.cap, PACKAGE = "RGtk2")

  return(w)
} 


cairoSetLineJoin <-
function(cr, line.join)
{
  checkPtrType(cr, "Cairo")
  

  w <- .RGtkCall("S_cairo_set_line_join", cr, line.join, PACKAGE = "RGtk2")

  return(w)
} 


cairoSetDash <-
function(cr, dashes, offset)
{
  checkPtrType(cr, "Cairo")
  dashes <- as.list(as.numeric(dashes))
  offset <- as.numeric(offset)

  w <- .RGtkCall("S_cairo_set_dash", cr, dashes, offset, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoSetMiterLimit <-
function(cr, limit)
{
  checkPtrType(cr, "Cairo")
  limit <- as.numeric(limit)

  w <- .RGtkCall("S_cairo_set_miter_limit", cr, limit, PACKAGE = "RGtk2")

  return(w)
} 


cairoTranslate <-
function(cr, tx, ty)
{
  checkPtrType(cr, "Cairo")
  tx <- as.numeric(tx)
  ty <- as.numeric(ty)

  w <- .RGtkCall("S_cairo_translate", cr, tx, ty, PACKAGE = "RGtk2")

  return(w)
} 


cairoScale <-
function(cr, sx, sy)
{
  checkPtrType(cr, "Cairo")
  sx <- as.numeric(sx)
  sy <- as.numeric(sy)

  w <- .RGtkCall("S_cairo_scale", cr, sx, sy, PACKAGE = "RGtk2")

  return(w)
} 


cairoRotate <-
function(cr, angle)
{
  checkPtrType(cr, "Cairo")
  angle <- as.numeric(angle)

  w <- .RGtkCall("S_cairo_rotate", cr, angle, PACKAGE = "RGtk2")

  return(w)
} 


cairoTransform <-
function(cr, matrix)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(matrix, "CairoMatrix")

  w <- .RGtkCall("S_cairo_transform", cr, matrix, PACKAGE = "RGtk2")

  return(w)
} 


cairoSetMatrix <-
function(cr, matrix)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(matrix, "CairoMatrix")

  w <- .RGtkCall("S_cairo_set_matrix", cr, matrix, PACKAGE = "RGtk2")

  return(w)
} 


cairoIdentityMatrix <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_identity_matrix", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoUserToDevice <-
function(cr, x, y)
{
  checkPtrType(cr, "Cairo")
  x <- as.list(as.numeric(x))
  y <- as.list(as.numeric(y))

  w <- .RGtkCall("S_cairo_user_to_device", cr, x, y, PACKAGE = "RGtk2")

  return(w)
} 


cairoUserToDeviceDistance <-
function(cr, dx, dy)
{
  checkPtrType(cr, "Cairo")
  dx <- as.list(as.numeric(dx))
  dy <- as.list(as.numeric(dy))

  w <- .RGtkCall("S_cairo_user_to_device_distance", cr, dx, dy, PACKAGE = "RGtk2")

  return(w)
} 


cairoDeviceToUser <-
function(cr, x, y)
{
  checkPtrType(cr, "Cairo")
  x <- as.list(as.numeric(x))
  y <- as.list(as.numeric(y))

  w <- .RGtkCall("S_cairo_device_to_user", cr, x, y, PACKAGE = "RGtk2")

  return(w)
} 


cairoDeviceToUserDistance <-
function(cr, dx, dy)
{
  checkPtrType(cr, "Cairo")
  dx <- as.list(as.numeric(dx))
  dy <- as.list(as.numeric(dy))

  w <- .RGtkCall("S_cairo_device_to_user_distance", cr, dx, dy, PACKAGE = "RGtk2")

  return(w)
} 


cairoNewPath <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_new_path", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoMoveTo <-
function(cr, x, y)
{
  checkPtrType(cr, "Cairo")
  x <- as.numeric(x)
  y <- as.numeric(y)

  w <- .RGtkCall("S_cairo_move_to", cr, x, y, PACKAGE = "RGtk2")

  return(w)
} 


cairoLineTo <-
function(cr, x, y)
{
  checkPtrType(cr, "Cairo")
  x <- as.numeric(x)
  y <- as.numeric(y)

  w <- .RGtkCall("S_cairo_line_to", cr, x, y, PACKAGE = "RGtk2")

  return(w)
} 


cairoCurveTo <-
function(cr, x1, y1, x2, y2, x3, y3)
{
  checkPtrType(cr, "Cairo")
  x1 <- as.numeric(x1)
  y1 <- as.numeric(y1)
  x2 <- as.numeric(x2)
  y2 <- as.numeric(y2)
  x3 <- as.numeric(x3)
  y3 <- as.numeric(y3)

  w <- .RGtkCall("S_cairo_curve_to", cr, x1, y1, x2, y2, x3, y3, PACKAGE = "RGtk2")

  return(w)
} 


cairoArc <-
function(cr, xc, yc, radius, angle1, angle2)
{
  checkPtrType(cr, "Cairo")
  xc <- as.numeric(xc)
  yc <- as.numeric(yc)
  radius <- as.numeric(radius)
  angle1 <- as.numeric(angle1)
  angle2 <- as.numeric(angle2)

  w <- .RGtkCall("S_cairo_arc", cr, xc, yc, radius, angle1, angle2, PACKAGE = "RGtk2")

  return(w)
} 


cairoArcNegative <-
function(cr, xc, yc, radius, angle1, angle2)
{
  checkPtrType(cr, "Cairo")
  xc <- as.numeric(xc)
  yc <- as.numeric(yc)
  radius <- as.numeric(radius)
  angle1 <- as.numeric(angle1)
  angle2 <- as.numeric(angle2)

  w <- .RGtkCall("S_cairo_arc_negative", cr, xc, yc, radius, angle1, angle2, PACKAGE = "RGtk2")

  return(w)
} 


cairoRelMoveTo <-
function(cr, dx, dy)
{
  checkPtrType(cr, "Cairo")
  dx <- as.numeric(dx)
  dy <- as.numeric(dy)

  w <- .RGtkCall("S_cairo_rel_move_to", cr, dx, dy, PACKAGE = "RGtk2")

  return(w)
} 


cairoRelLineTo <-
function(cr, dx, dy)
{
  checkPtrType(cr, "Cairo")
  dx <- as.numeric(dx)
  dy <- as.numeric(dy)

  w <- .RGtkCall("S_cairo_rel_line_to", cr, dx, dy, PACKAGE = "RGtk2")

  return(w)
} 


cairoRelCurveTo <-
function(cr, dx1, dy1, dx2, dy2, dx3, dy3)
{
  checkPtrType(cr, "Cairo")
  dx1 <- as.numeric(dx1)
  dy1 <- as.numeric(dy1)
  dx2 <- as.numeric(dx2)
  dy2 <- as.numeric(dy2)
  dx3 <- as.numeric(dx3)
  dy3 <- as.numeric(dy3)

  w <- .RGtkCall("S_cairo_rel_curve_to", cr, dx1, dy1, dx2, dy2, dx3, dy3, PACKAGE = "RGtk2")

  return(w)
} 


cairoRectangle <-
function(cr, x, y, width, height)
{
  checkPtrType(cr, "Cairo")
  x <- as.numeric(x)
  y <- as.numeric(y)
  width <- as.numeric(width)
  height <- as.numeric(height)

  w <- .RGtkCall("S_cairo_rectangle", cr, x, y, width, height, PACKAGE = "RGtk2")

  return(w)
} 


cairoClosePath <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_close_path", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoPaint <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_paint", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoPaintWithAlpha <-
function(cr, alpha)
{
  checkPtrType(cr, "Cairo")
  alpha <- as.numeric(alpha)

  w <- .RGtkCall("S_cairo_paint_with_alpha", cr, alpha, PACKAGE = "RGtk2")

  return(w)
} 


cairoMask <-
function(cr, pattern)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(pattern, "CairoPattern")

  w <- .RGtkCall("S_cairo_mask", cr, pattern, PACKAGE = "RGtk2")

  return(w)
} 


cairoMaskSurface <-
function(cr, surface, surface.x, surface.y)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(surface, "CairoSurface")
  surface.x <- as.numeric(surface.x)
  surface.y <- as.numeric(surface.y)

  w <- .RGtkCall("S_cairo_mask_surface", cr, surface, surface.x, surface.y, PACKAGE = "RGtk2")

  return(w)
} 


cairoStroke <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_stroke", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoStrokePreserve <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_stroke_preserve", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoFill <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_fill", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoFillPreserve <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_fill_preserve", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoCopyPage <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_copy_page", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoShowPage <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_show_page", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoInStroke <-
function(cr, x, y)
{
  checkPtrType(cr, "Cairo")
  x <- as.numeric(x)
  y <- as.numeric(y)

  w <- .RGtkCall("S_cairo_in_stroke", cr, x, y, PACKAGE = "RGtk2")

  return(w)
} 


cairoInFill <-
function(cr, x, y)
{
  checkPtrType(cr, "Cairo")
  x <- as.numeric(x)
  y <- as.numeric(y)

  w <- .RGtkCall("S_cairo_in_fill", cr, x, y, PACKAGE = "RGtk2")

  return(w)
} 


cairoStrokeExtents <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_stroke_extents", cr, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoFillExtents <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_fill_extents", cr, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoResetClip <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_reset_clip", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoClip <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_clip", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoClipPreserve <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_clip_preserve", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoSelectFontFace <-
function(cr, family, slant, weight)
{
  checkPtrType(cr, "Cairo")
  family <- as.character(family)
  
  

  w <- .RGtkCall("S_cairo_select_font_face", cr, family, slant, weight, PACKAGE = "RGtk2")

  return(w)
} 


cairoSetFontSize <-
function(cr, size)
{
  checkPtrType(cr, "Cairo")
  size <- as.numeric(size)

  w <- .RGtkCall("S_cairo_set_font_size", cr, size, PACKAGE = "RGtk2")

  return(w)
} 


cairoSetFontMatrix <-
function(cr, matrix)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(matrix, "CairoMatrix")

  w <- .RGtkCall("S_cairo_set_font_matrix", cr, matrix, PACKAGE = "RGtk2")

  return(w)
} 


cairoGetFontMatrix <-
function(cr, matrix)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(matrix, "CairoMatrix")

  w <- .RGtkCall("S_cairo_get_font_matrix", cr, matrix, PACKAGE = "RGtk2")

  return(w)
} 


cairoShowText <-
function(cr, utf8)
{
  checkPtrType(cr, "Cairo")
  utf8 <- as.character(utf8)

  w <- .RGtkCall("S_cairo_show_text", cr, utf8, PACKAGE = "RGtk2")

  return(w)
} 


cairoShowGlyphs <-
function(cr, glyphs, num.glyphs)
{
  checkPtrType(cr, "Cairo")
  glyphs <- as.CairoGlyph(glyphs)
  num.glyphs <- as.integer(num.glyphs)

  w <- .RGtkCall("S_cairo_show_glyphs", cr, glyphs, num.glyphs, PACKAGE = "RGtk2")

  return(w)
} 


cairoGetFontFace <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_get_font_face", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoFontExtents <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_font_extents", cr, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoSetFontFace <-
function(cr, font.face)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(font.face, "CairoFontFace")

  w <- .RGtkCall("S_cairo_set_font_face", cr, font.face, PACKAGE = "RGtk2")

  return(w)
} 


cairoTextExtents <-
function(cr, utf8)
{
  checkPtrType(cr, "Cairo")
  utf8 <- as.character(utf8)

  w <- .RGtkCall("S_cairo_text_extents", cr, utf8, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoGlyphExtents <-
function(cr, glyphs)
{
  checkPtrType(cr, "Cairo")
  glyphs <- lapply(glyphs, function(x) { x <- as.CairoGlyph(x); x })

  w <- .RGtkCall("S_cairo_glyph_extents", cr, glyphs, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoTextPath <-
function(cr, utf8)
{
  checkPtrType(cr, "Cairo")
  utf8 <- as.character(utf8)

  w <- .RGtkCall("S_cairo_text_path", cr, utf8, PACKAGE = "RGtk2")

  return(w)
} 


cairoGlyphPath <-
function(cr, glyphs)
{
  checkPtrType(cr, "Cairo")
  glyphs <- lapply(glyphs, function(x) { x <- as.CairoGlyph(x); x })

  w <- .RGtkCall("S_cairo_glyph_path", cr, glyphs, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoSetFontOptions <-
function(cr, options)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(options, "CairoFontOptions")

  w <- .RGtkCall("S_cairo_set_font_options", cr, options, PACKAGE = "RGtk2")

  return(w)
} 


cairoGetFontOptions <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_get_font_options", cr, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoFontFaceGetUserData <-
function(font.face, key)
{
  checkPtrType(font.face, "CairoFontFace")
  checkPtrType(key, "CairoUserDataKey")

  w <- .RGtkCall("S_cairo_font_face_get_user_data", font.face, key, PACKAGE = "RGtk2")

  return(w)
} 


cairoFontFaceSetUserData <-
function(font.face, key, user.data)
{
  checkPtrType(font.face, "CairoFontFace")
  checkPtrType(key, "CairoUserDataKey")
  

  w <- .RGtkCall("S_cairo_font_face_set_user_data", font.face, key, user.data, PACKAGE = "RGtk2")

  return(w)
} 


cairoFontFaceStatus <-
function(font.face)
{
  checkPtrType(font.face, "CairoFontFace")

  w <- .RGtkCall("S_cairo_font_face_status", font.face, PACKAGE = "RGtk2")

  return(w)
} 


cairoScaledFontCreate <-
function(font.face, font.matrix, ctm, option)
{
  checkPtrType(font.face, "CairoFontFace")
  checkPtrType(font.matrix, "CairoMatrix")
  checkPtrType(ctm, "CairoMatrix")
  checkPtrType(option, "CairoFontOptions")

  w <- .RGtkCall("S_cairo_scaled_font_create", font.face, font.matrix, ctm, option, PACKAGE = "RGtk2")

  return(w)
} 


cairoScaledFontExtents <-
function(scaled.font)
{
  checkPtrType(scaled.font, "CairoScaledFont")

  w <- .RGtkCall("S_cairo_scaled_font_extents", scaled.font, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoScaledFontGlyphExtents <-
function(scaled.font, glyphs, num.glyphs)
{
  checkPtrType(scaled.font, "CairoScaledFont")
  glyphs <- as.CairoGlyph(glyphs)
  num.glyphs <- as.integer(num.glyphs)

  w <- .RGtkCall("S_cairo_scaled_font_glyph_extents", scaled.font, glyphs, num.glyphs, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoScaledFontStatus <-
function(scaled.font)
{
  checkPtrType(scaled.font, "CairoScaledFont")

  w <- .RGtkCall("S_cairo_scaled_font_status", scaled.font, PACKAGE = "RGtk2")

  return(w)
} 


cairoGetOperator <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_get_operator", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoGetSource <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_get_source", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoGetTolerance <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_get_tolerance", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoGetAntialias <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_get_antialias", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoGetCurrentPoint <-
function(cr, x, y)
{
  checkPtrType(cr, "Cairo")
  x <- as.list(as.numeric(x))
  y <- as.list(as.numeric(y))

  w <- .RGtkCall("S_cairo_get_current_point", cr, x, y, PACKAGE = "RGtk2")

  return(w)
} 


cairoGetFillRule <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_get_fill_rule", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoGetLineWidth <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_get_line_width", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoGetLineCap <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_get_line_cap", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoGetLineJoin <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_get_line_join", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoGetMiterLimit <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_get_miter_limit", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoGetMatrix <-
function(cr, matrix)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(matrix, "CairoMatrix")

  w <- .RGtkCall("S_cairo_get_matrix", cr, matrix, PACKAGE = "RGtk2")

  return(w)
} 


cairoGetTarget <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_get_target", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoCopyPath <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_copy_path", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoCopyPathFlat <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_copy_path_flat", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoAppendPath <-
function(cr, path)
{
  checkPtrType(cr, "Cairo")
  path <- as.CairoPath(path)

  w <- .RGtkCall("S_cairo_append_path", cr, path, PACKAGE = "RGtk2")

  return(w)
} 


cairoStatus <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_status", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoStatusToString <-
function(status)
{
  

  w <- .RGtkCall("S_cairo_status_to_string", status, PACKAGE = "RGtk2")

  return(w)
} 


cairoSurfaceCreateSimilar <-
function(other, content, width, height)
{
  checkPtrType(other, "CairoSurface")
  
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_cairo_surface_create_similar", other, content, width, height, PACKAGE = "RGtk2")

  return(w)
} 


cairoSurfaceDestroy <-
function(surface)
{
  checkPtrType(surface, "CairoSurface")

  w <- .RGtkCall("S_cairo_surface_destroy", surface, PACKAGE = "RGtk2")

  return(w)
} 


cairoSurfaceFinish <-
function(surface)
{
  checkPtrType(surface, "CairoSurface")

  w <- .RGtkCall("S_cairo_surface_finish", surface, PACKAGE = "RGtk2")

  return(w)
} 


cairoSurfaceWriteToPng <-
function(surface, filename)
{
  checkPtrType(surface, "CairoSurface")
  filename <- as.character(filename)

  w <- .RGtkCall("S_cairo_surface_write_to_png", surface, filename, PACKAGE = "RGtk2")

  return(w)
} 


cairoSurfaceGetUserData <-
function(surface, key)
{
  checkPtrType(surface, "CairoSurface")
  checkPtrType(key, "CairoUserDataKey")

  w <- .RGtkCall("S_cairo_surface_get_user_data", surface, key, PACKAGE = "RGtk2")

  return(w)
} 


cairoSurfaceSetUserData <-
function(surface, key, user.data)
{
  checkPtrType(surface, "CairoSurface")
  checkPtrType(key, "CairoUserDataKey")
  

  w <- .RGtkCall("S_cairo_surface_set_user_data", surface, key, user.data, PACKAGE = "RGtk2")

  return(w)
} 


cairoSurfaceFlush <-
function(surface)
{
  checkPtrType(surface, "CairoSurface")

  w <- .RGtkCall("S_cairo_surface_flush", surface, PACKAGE = "RGtk2")

  return(w)
} 


cairoSurfaceMarkDirty <-
function(surface)
{
  checkPtrType(surface, "CairoSurface")

  w <- .RGtkCall("S_cairo_surface_mark_dirty", surface, PACKAGE = "RGtk2")

  return(w)
} 


cairoSurfaceMarkDirtyRectangle <-
function(surface, x, y, width, height)
{
  checkPtrType(surface, "CairoSurface")
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_cairo_surface_mark_dirty_rectangle", surface, x, y, width, height, PACKAGE = "RGtk2")

  return(w)
} 


cairoSurfaceSetDeviceOffset <-
function(surface, x.offset, y.offset)
{
  checkPtrType(surface, "CairoSurface")
  x.offset <- as.numeric(x.offset)
  y.offset <- as.numeric(y.offset)

  w <- .RGtkCall("S_cairo_surface_set_device_offset", surface, x.offset, y.offset, PACKAGE = "RGtk2")

  return(w)
} 


cairoSurfaceGetFontOptions <-
function(surface)
{
  checkPtrType(surface, "CairoSurface")

  w <- .RGtkCall("S_cairo_surface_get_font_options", surface, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoSurfaceStatus <-
function(surface)
{
  checkPtrType(surface, "CairoSurface")

  w <- .RGtkCall("S_cairo_surface_status", surface, PACKAGE = "RGtk2")

  return(w)
} 


cairoImageSurfaceCreate <-
function(format, width, height)
{
  
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_cairo_image_surface_create", format, width, height, PACKAGE = "RGtk2")

  return(w)
} 


cairoImageSurfaceCreateForData <-
function(data, format, width, height, stride)
{
  data <- as.list(as.raw(data))
  
  width <- as.integer(width)
  height <- as.integer(height)
  stride <- as.integer(stride)

  w <- .RGtkCall("S_cairo_image_surface_create_for_data", data, format, width, height, stride, PACKAGE = "RGtk2")

  return(w)
} 


cairoImageSurfaceGetWidth <-
function(surface)
{
  checkPtrType(surface, "CairoSurface")

  w <- .RGtkCall("S_cairo_image_surface_get_width", surface, PACKAGE = "RGtk2")

  return(w)
} 


cairoImageSurfaceGetHeight <-
function(surface)
{
  checkPtrType(surface, "CairoSurface")

  w <- .RGtkCall("S_cairo_image_surface_get_height", surface, PACKAGE = "RGtk2")

  return(w)
} 


cairoImageSurfaceCreateFromPng <-
function(filename)
{
  filename <- as.character(filename)

  w <- .RGtkCall("S_cairo_image_surface_create_from_png", filename, PACKAGE = "RGtk2")

  return(w)
} 


cairoPatternCreateRgb <-
function(red, green, blue)
{
  red <- as.numeric(red)
  green <- as.numeric(green)
  blue <- as.numeric(blue)

  w <- .RGtkCall("S_cairo_pattern_create_rgb", red, green, blue, PACKAGE = "RGtk2")

  return(w)
} 


cairoPatternCreateRgba <-
function(red, green, blue, alpha)
{
  red <- as.numeric(red)
  green <- as.numeric(green)
  blue <- as.numeric(blue)
  alpha <- as.numeric(alpha)

  w <- .RGtkCall("S_cairo_pattern_create_rgba", red, green, blue, alpha, PACKAGE = "RGtk2")

  return(w)
} 


cairoPatternCreateForSurface <-
function(surface)
{
  checkPtrType(surface, "CairoSurface")

  w <- .RGtkCall("S_cairo_pattern_create_for_surface", surface, PACKAGE = "RGtk2")

  return(w)
} 


cairoPatternCreateLinear <-
function(x0, y0, x1, y1)
{
  x0 <- as.numeric(x0)
  y0 <- as.numeric(y0)
  x1 <- as.numeric(x1)
  y1 <- as.numeric(y1)

  w <- .RGtkCall("S_cairo_pattern_create_linear", x0, y0, x1, y1, PACKAGE = "RGtk2")

  return(w)
} 


cairoPatternCreateRadial <-
function(cx0, cy0, radius0, cx1, cy1, radius1)
{
  cx0 <- as.numeric(cx0)
  cy0 <- as.numeric(cy0)
  radius0 <- as.numeric(radius0)
  cx1 <- as.numeric(cx1)
  cy1 <- as.numeric(cy1)
  radius1 <- as.numeric(radius1)

  w <- .RGtkCall("S_cairo_pattern_create_radial", cx0, cy0, radius0, cx1, cy1, radius1, PACKAGE = "RGtk2")

  return(w)
} 


cairoPatternStatus <-
function(pattern)
{
  checkPtrType(pattern, "CairoPattern")

  w <- .RGtkCall("S_cairo_pattern_status", pattern, PACKAGE = "RGtk2")

  return(w)
} 


cairoPatternAddColorStopRgb <-
function(pattern, offset, red, green, blue)
{
  checkPtrType(pattern, "CairoPattern")
  offset <- as.numeric(offset)
  red <- as.numeric(red)
  green <- as.numeric(green)
  blue <- as.numeric(blue)

  w <- .RGtkCall("S_cairo_pattern_add_color_stop_rgb", pattern, offset, red, green, blue, PACKAGE = "RGtk2")

  return(w)
} 


cairoPatternAddColorStopRgba <-
function(pattern, offset, red, green, blue, alpha)
{
  checkPtrType(pattern, "CairoPattern")
  offset <- as.numeric(offset)
  red <- as.numeric(red)
  green <- as.numeric(green)
  blue <- as.numeric(blue)
  alpha <- as.numeric(alpha)

  w <- .RGtkCall("S_cairo_pattern_add_color_stop_rgba", pattern, offset, red, green, blue, alpha, PACKAGE = "RGtk2")

  return(w)
} 


cairoPatternSetMatrix <-
function(pattern, matrix)
{
  checkPtrType(pattern, "CairoPattern")
  checkPtrType(matrix, "CairoMatrix")

  w <- .RGtkCall("S_cairo_pattern_set_matrix", pattern, matrix, PACKAGE = "RGtk2")

  return(w)
} 


cairoPatternGetMatrix <-
function(pattern, matrix)
{
  checkPtrType(pattern, "CairoPattern")
  checkPtrType(matrix, "CairoMatrix")

  w <- .RGtkCall("S_cairo_pattern_get_matrix", pattern, matrix, PACKAGE = "RGtk2")

  return(w)
} 


cairoPatternSetExtend <-
function(pattern, extend)
{
  checkPtrType(pattern, "CairoPattern")
  

  w <- .RGtkCall("S_cairo_pattern_set_extend", pattern, extend, PACKAGE = "RGtk2")

  return(w)
} 


cairoPatternGetExtend <-
function(pattern)
{
  checkPtrType(pattern, "CairoPattern")

  w <- .RGtkCall("S_cairo_pattern_get_extend", pattern, PACKAGE = "RGtk2")

  return(w)
} 


cairoPatternSetFilter <-
function(pattern, filter)
{
  checkPtrType(pattern, "CairoPattern")
  

  w <- .RGtkCall("S_cairo_pattern_set_filter", pattern, filter, PACKAGE = "RGtk2")

  return(w)
} 


cairoPatternGetFilter <-
function(pattern)
{
  checkPtrType(pattern, "CairoPattern")

  w <- .RGtkCall("S_cairo_pattern_get_filter", pattern, PACKAGE = "RGtk2")

  return(w)
} 


cairoMatrixInit <-
function(xx, yx, xy, yy, x0, y0)
{
  xx <- as.numeric(xx)
  yx <- as.numeric(yx)
  xy <- as.numeric(xy)
  yy <- as.numeric(yy)
  x0 <- as.numeric(x0)
  y0 <- as.numeric(y0)

  w <- .RGtkCall("S_cairo_matrix_init", xx, yx, xy, yy, x0, y0, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoMatrixInitIdentity <-
function()
{
  

  w <- .RGtkCall("S_cairo_matrix_init_identity", PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoMatrixInitTranslate <-
function(tx, ty)
{
  tx <- as.numeric(tx)
  ty <- as.numeric(ty)

  w <- .RGtkCall("S_cairo_matrix_init_translate", tx, ty, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoMatrixInitScale <-
function(sx, sy)
{
  sx <- as.numeric(sx)
  sy <- as.numeric(sy)

  w <- .RGtkCall("S_cairo_matrix_init_scale", sx, sy, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoMatrixInitRotate <-
function(radians)
{
  radians <- as.numeric(radians)

  w <- .RGtkCall("S_cairo_matrix_init_rotate", radians, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoMatrixTranslate <-
function(matrix, tx, ty)
{
  checkPtrType(matrix, "CairoMatrix")
  tx <- as.numeric(tx)
  ty <- as.numeric(ty)

  w <- .RGtkCall("S_cairo_matrix_translate", matrix, tx, ty, PACKAGE = "RGtk2")

  return(w)
} 


cairoMatrixScale <-
function(matrix, sx, sy)
{
  checkPtrType(matrix, "CairoMatrix")
  sx <- as.numeric(sx)
  sy <- as.numeric(sy)

  w <- .RGtkCall("S_cairo_matrix_scale", matrix, sx, sy, PACKAGE = "RGtk2")

  return(w)
} 


cairoMatrixRotate <-
function(matrix, radians)
{
  checkPtrType(matrix, "CairoMatrix")
  radians <- as.numeric(radians)

  w <- .RGtkCall("S_cairo_matrix_rotate", matrix, radians, PACKAGE = "RGtk2")

  return(w)
} 


cairoMatrixInvert <-
function(matrix)
{
  checkPtrType(matrix, "CairoMatrix")

  w <- .RGtkCall("S_cairo_matrix_invert", matrix, PACKAGE = "RGtk2")

  return(w)
} 


cairoMatrixMultiply <-
function(result, a, b)
{
  checkPtrType(result, "CairoMatrix")
  checkPtrType(a, "CairoMatrix")
  checkPtrType(b, "CairoMatrix")

  w <- .RGtkCall("S_cairo_matrix_multiply", result, a, b, PACKAGE = "RGtk2")

  return(w)
} 


cairoMatrixTransformDistance <-
function(matrix, dx, dy)
{
  checkPtrType(matrix, "CairoMatrix")
  dx <- as.list(as.numeric(dx))
  dy <- as.list(as.numeric(dy))

  w <- .RGtkCall("S_cairo_matrix_transform_distance", matrix, dx, dy, PACKAGE = "RGtk2")

  return(w)
} 


cairoMatrixTransformPoint <-
function(matrix, x, y)
{
  checkPtrType(matrix, "CairoMatrix")
  x <- as.list(as.numeric(x))
  y <- as.list(as.numeric(y))

  w <- .RGtkCall("S_cairo_matrix_transform_point", matrix, x, y, PACKAGE = "RGtk2")

  return(w)
} 


cairoFontOptionsCreate <-
function()
{
  

  w <- .RGtkCall("S_cairo_font_options_create", PACKAGE = "RGtk2")

  return(w)
} 


cairoFontOptionsCopy <-
function(original)
{
  checkPtrType(original, "CairoFontOptions")

  w <- .RGtkCall("S_cairo_font_options_copy", original, PACKAGE = "RGtk2")

  return(w)
} 


cairoFontOptionsStatus <-
function(options)
{
  checkPtrType(options, "CairoFontOptions")

  w <- .RGtkCall("S_cairo_font_options_status", options, PACKAGE = "RGtk2")

  return(w)
} 


cairoFontOptionsMerge <-
function(options, other)
{
  checkPtrType(options, "CairoFontOptions")
  checkPtrType(other, "CairoFontOptions")

  w <- .RGtkCall("S_cairo_font_options_merge", options, other, PACKAGE = "RGtk2")

  return(w)
} 


cairoFontOptionsEqual <-
function(options, other)
{
  checkPtrType(options, "CairoFontOptions")
  checkPtrType(other, "CairoFontOptions")

  w <- .RGtkCall("S_cairo_font_options_equal", options, other, PACKAGE = "RGtk2")

  return(w)
} 


cairoFontOptionsSetAntialias <-
function(options, antialias)
{
  checkPtrType(options, "CairoFontOptions")
  

  w <- .RGtkCall("S_cairo_font_options_set_antialias", options, antialias, PACKAGE = "RGtk2")

  return(w)
} 


cairoFontOptionsGetAntialias <-
function(options)
{
  checkPtrType(options, "CairoFontOptions")

  w <- .RGtkCall("S_cairo_font_options_get_antialias", options, PACKAGE = "RGtk2")

  return(w)
} 


cairoFontOptionsSetSubpixelOrder <-
function(options, subpixel.order)
{
  checkPtrType(options, "CairoFontOptions")
  

  w <- .RGtkCall("S_cairo_font_options_set_subpixel_order", options, subpixel.order, PACKAGE = "RGtk2")

  return(w)
} 


cairoFontOptionsGetSubpixelOrder <-
function(options)
{
  checkPtrType(options, "CairoFontOptions")

  w <- .RGtkCall("S_cairo_font_options_get_subpixel_order", options, PACKAGE = "RGtk2")

  return(w)
} 


cairoFontOptionsSetHintStyle <-
function(options, hint.style)
{
  checkPtrType(options, "CairoFontOptions")
  

  w <- .RGtkCall("S_cairo_font_options_set_hint_style", options, hint.style, PACKAGE = "RGtk2")

  return(w)
} 


cairoFontOptionsGetHintStyle <-
function(options)
{
  checkPtrType(options, "CairoFontOptions")

  w <- .RGtkCall("S_cairo_font_options_get_hint_style", options, PACKAGE = "RGtk2")

  return(w)
} 


cairoFontOptionsSetHintMetrics <-
function(options, hint.metrics)
{
  checkPtrType(options, "CairoFontOptions")
  

  w <- .RGtkCall("S_cairo_font_options_set_hint_metrics", options, hint.metrics, PACKAGE = "RGtk2")

  return(w)
} 


cairoFontOptionsGetHintMetrics <-
function(options)
{
  checkPtrType(options, "CairoFontOptions")

  w <- .RGtkCall("S_cairo_font_options_get_hint_metrics", options, PACKAGE = "RGtk2")

  return(w)
} 


cairoPushGroup <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_push_group", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoPushGroupWithContent <-
function(cr, content)
{
  checkPtrType(cr, "Cairo")
  

  w <- .RGtkCall("S_cairo_push_group_with_content", cr, content, PACKAGE = "RGtk2")

  return(w)
} 


cairoPopGroup <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_pop_group", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoPopGroupToSource <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_pop_group_to_source", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoGetGroupTarget <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_get_group_target", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoNewSubPath <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_new_sub_path", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoSetScaledFont <-
function(cr, scaled.font)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(scaled.font, "CairoScaledFont")

  w <- .RGtkCall("S_cairo_set_scaled_font", cr, scaled.font, PACKAGE = "RGtk2")

  return(w)
} 


cairoScaledFontGetFontFace <-
function(scaled.font)
{
  checkPtrType(scaled.font, "CairoScaledFont")

  w <- .RGtkCall("S_cairo_scaled_font_get_font_face", scaled.font, PACKAGE = "RGtk2")

  return(w)
} 


cairoScaledFontGetFontMatrix <-
function(scaled.font)
{
  checkPtrType(scaled.font, "CairoScaledFont")

  w <- .RGtkCall("S_cairo_scaled_font_get_font_matrix", scaled.font, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoScaledFontGetCtm <-
function(scaled.font)
{
  checkPtrType(scaled.font, "CairoScaledFont")

  w <- .RGtkCall("S_cairo_scaled_font_get_ctm", scaled.font, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoScaledFontGetFontOptions <-
function(scaled.font)
{
  checkPtrType(scaled.font, "CairoScaledFont")

  w <- .RGtkCall("S_cairo_scaled_font_get_font_options", scaled.font, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoScaledFontTextExtents <-
function(scaled.font, utf8)
{
  checkPtrType(scaled.font, "CairoScaledFont")
  utf8 <- as.character(utf8)

  w <- .RGtkCall("S_cairo_scaled_font_text_extents", scaled.font, utf8, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoScaledFontGetType <-
function(scaled.font)
{
  checkPtrType(scaled.font, "CairoScaledFont")

  w <- .RGtkCall("S_cairo_scaled_font_get_type", scaled.font, PACKAGE = "RGtk2")

  return(w)
} 


cairoFontFaceGetType <-
function(font.face)
{
  checkPtrType(font.face, "CairoFontFace")

  w <- .RGtkCall("S_cairo_font_face_get_type", font.face, PACKAGE = "RGtk2")

  return(w)
} 


cairoSurfaceGetType <-
function(surface)
{
  checkPtrType(surface, "CairoSurface")

  w <- .RGtkCall("S_cairo_surface_get_type", surface, PACKAGE = "RGtk2")

  return(w)
} 


cairoSurfaceGetDeviceOffset <-
function(surface)
{
  checkPtrType(surface, "CairoSurface")

  w <- .RGtkCall("S_cairo_surface_get_device_offset", surface, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoSurfaceSetFallbackResolution <-
function(surface, x.pixels.per.inch, y.pixels.per.inch)
{
  checkPtrType(surface, "CairoSurface")
  x.pixels.per.inch <- as.numeric(x.pixels.per.inch)
  y.pixels.per.inch <- as.numeric(y.pixels.per.inch)

  w <- .RGtkCall("S_cairo_surface_set_fallback_resolution", surface, x.pixels.per.inch, y.pixels.per.inch, PACKAGE = "RGtk2")

  return(w)
} 


cairoSurfaceGetContent <-
function(surface)
{
  checkPtrType(surface, "CairoSurface")

  w <- .RGtkCall("S_cairo_surface_get_content", surface, PACKAGE = "RGtk2")

  return(w)
} 


cairoImageSurfaceGetFormat <-
function(surface)
{
  checkPtrType(surface, "CairoSurface")

  w <- .RGtkCall("S_cairo_image_surface_get_format", surface, PACKAGE = "RGtk2")

  return(w)
} 


cairoImageSurfaceGetStride <-
function(surface)
{
  checkPtrType(surface, "CairoSurface")

  w <- .RGtkCall("S_cairo_image_surface_get_stride", surface, PACKAGE = "RGtk2")

  return(w)
} 


cairoImageSurfaceGetData <-
function(surface)
{
  checkPtrType(surface, "CairoSurface")

  w <- .RGtkCall("S_cairo_image_surface_get_data", surface, PACKAGE = "RGtk2")

  return(w)
} 


cairoPatternGetType <-
function(pattern)
{
  checkPtrType(pattern, "CairoPattern")

  w <- .RGtkCall("S_cairo_pattern_get_type", pattern, PACKAGE = "RGtk2")

  return(w)
} 


cairoPdfSurfaceCreate <-
function(filename, width.in.points, height.in.points)
{
  filename <- as.character(filename)
  width.in.points <- as.numeric(width.in.points)
  height.in.points <- as.numeric(height.in.points)

  w <- .RGtkCall("S_cairo_pdf_surface_create", filename, width.in.points, height.in.points, PACKAGE = "RGtk2")

  return(w)
} 


cairoPdfSurfaceCreateForStream <-
function(write.func, closure = NULL, width.in.points, height.in.points)
{
  write.func <- as.function(write.func)
  
  width.in.points <- as.numeric(width.in.points)
  height.in.points <- as.numeric(height.in.points)

  w <- .RGtkCall("S_cairo_pdf_surface_create_for_stream", write.func, closure, width.in.points, height.in.points, PACKAGE = "RGtk2")

  return(w)
} 


cairoPdfSurfaceSetSize <-
function(surface, width.in.points, height.in.points)
{
  checkPtrType(surface, "CairoSurface")
  width.in.points <- as.numeric(width.in.points)
  height.in.points <- as.numeric(height.in.points)

  w <- .RGtkCall("S_cairo_pdf_surface_set_size", surface, width.in.points, height.in.points, PACKAGE = "RGtk2")

  return(w)
} 


cairoPsSurfaceCreate <-
function(filename, width.in.points, height.in.points)
{
  filename <- as.character(filename)
  width.in.points <- as.numeric(width.in.points)
  height.in.points <- as.numeric(height.in.points)

  w <- .RGtkCall("S_cairo_ps_surface_create", filename, width.in.points, height.in.points, PACKAGE = "RGtk2")

  return(w)
} 


cairoPsSurfaceCreateForStream <-
function(write.func, closure = NULL, width.in.points, height.in.points)
{
  write.func <- as.function(write.func)
  
  width.in.points <- as.numeric(width.in.points)
  height.in.points <- as.numeric(height.in.points)

  w <- .RGtkCall("S_cairo_ps_surface_create_for_stream", write.func, closure, width.in.points, height.in.points, PACKAGE = "RGtk2")

  return(w)
} 


cairoPsSurfaceSetSize <-
function(surface, width.in.points, height.in.points)
{
  checkPtrType(surface, "CairoSurface")
  width.in.points <- as.numeric(width.in.points)
  height.in.points <- as.numeric(height.in.points)

  w <- .RGtkCall("S_cairo_ps_surface_set_size", surface, width.in.points, height.in.points, PACKAGE = "RGtk2")

  return(w)
} 


cairoPsSurfaceDscComment <-
function(surface, comment)
{
  checkPtrType(surface, "CairoSurface")
  comment <- as.character(comment)

  w <- .RGtkCall("S_cairo_ps_surface_dsc_comment", surface, comment, PACKAGE = "RGtk2")

  return(w)
} 


cairoPsSurfaceDscBeginSetup <-
function(surface)
{
  checkPtrType(surface, "CairoSurface")

  w <- .RGtkCall("S_cairo_ps_surface_dsc_begin_setup", surface, PACKAGE = "RGtk2")

  return(w)
} 


cairoPsSurfaceDscBeginPageSetup <-
function(surface)
{
  checkPtrType(surface, "CairoSurface")

  w <- .RGtkCall("S_cairo_ps_surface_dsc_begin_page_setup", surface, PACKAGE = "RGtk2")

  return(w)
} 


cairoSvgSurfaceCreate <-
function(filename, width.in.points, height.in.points)
{
  filename <- as.character(filename)
  width.in.points <- as.numeric(width.in.points)
  height.in.points <- as.numeric(height.in.points)

  w <- .RGtkCall("S_cairo_svg_surface_create", filename, width.in.points, height.in.points, PACKAGE = "RGtk2")

  return(w)
} 


cairoSvgSurfaceCreateForStream <-
function(write.func, closure = NULL, width.in.points, height.in.points)
{
  write.func <- as.function(write.func)
  
  width.in.points <- as.numeric(width.in.points)
  height.in.points <- as.numeric(height.in.points)

  w <- .RGtkCall("S_cairo_svg_surface_create_for_stream", write.func, closure, width.in.points, height.in.points, PACKAGE = "RGtk2")

  return(w)
} 


cairoSvgSurfaceRestrictToVersion <-
function(surface, version)
{
  checkPtrType(surface, "CairoSurface")
  

  w <- .RGtkCall("S_cairo_svg_surface_restrict_to_version", surface, version, PACKAGE = "RGtk2")

  return(w)
} 


cairoSvgGetVersions <-
function()
{
  

  w <- .RGtkCall("S_cairo_svg_get_versions", PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoSvgVersionToString <-
function(version)
{
  

  w <- .RGtkCall("S_cairo_svg_version_to_string", version, PACKAGE = "RGtk2")

  return(w)
} 


cairoClipExtents <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_clip_extents", cr, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoCopyClipRectangleList <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_copy_clip_rectangle_list", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoGetDashCount <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_get_dash_count", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoGetDash <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_get_dash", cr, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoPatternGetRgba <-
function(pattern)
{
  checkPtrType(pattern, "CairoPattern")

  w <- .RGtkCall("S_cairo_pattern_get_rgba", pattern, PACKAGE = "RGtk2")

  return(w)
} 


cairoPatternGetSurface <-
function(pattern)
{
  checkPtrType(pattern, "CairoPattern")

  w <- .RGtkCall("S_cairo_pattern_get_surface", pattern, PACKAGE = "RGtk2")

  return(w)
} 


cairoPatternGetColorStopRgba <-
function(pattern, index)
{
  checkPtrType(pattern, "CairoPattern")
  index <- as.integer(index)

  w <- .RGtkCall("S_cairo_pattern_get_color_stop_rgba", pattern, index, PACKAGE = "RGtk2")

  return(w)
} 


cairoPatternGetColorStopCount <-
function(pattern)
{
  checkPtrType(pattern, "CairoPattern")

  w <- .RGtkCall("S_cairo_pattern_get_color_stop_count", pattern, PACKAGE = "RGtk2")

  return(w)
} 


cairoPatternGetLinearPoints <-
function(pattern)
{
  checkPtrType(pattern, "CairoPattern")

  w <- .RGtkCall("S_cairo_pattern_get_linear_points", pattern, PACKAGE = "RGtk2")

  return(w)
} 


cairoPatternGetRadialCircles <-
function(pattern)
{
  checkPtrType(pattern, "CairoPattern")

  w <- .RGtkCall("S_cairo_pattern_get_radial_circles", pattern, PACKAGE = "RGtk2")

  return(w)
} 


cairoGetScaledFont <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_get_scaled_font", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoSetUserData <-
function(cr, key, user.data)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(key, "CairoUserDataKey")
  

  w <- .RGtkCall("S_cairo_set_user_data", cr, key, user.data, PACKAGE = "RGtk2")

  return(w)
} 


cairoGetUserData <-
function(cr, key)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(key, "CairoUserDataKey")

  w <- .RGtkCall("S_cairo_get_user_data", cr, key, PACKAGE = "RGtk2")

  return(w)
} 


cairoScaledFontGetUserData <-
function(scaled.font, key)
{
  checkPtrType(scaled.font, "CairoScaledFont")
  checkPtrType(key, "CairoUserDataKey")

  w <- .RGtkCall("S_cairo_scaled_font_get_user_data", scaled.font, key, PACKAGE = "RGtk2")

  return(w)
} 


cairoScaledFontSetUserData <-
function(scaled.font, key, user.data)
{
  checkPtrType(scaled.font, "CairoScaledFont")
  checkPtrType(key, "CairoUserDataKey")
  

  w <- .RGtkCall("S_cairo_scaled_font_set_user_data", scaled.font, key, user.data, PACKAGE = "RGtk2")

  return(w)
} 


cairoPatternGetUserData <-
function(pattern, key)
{
  checkPtrType(pattern, "CairoPattern")
  checkPtrType(key, "CairoUserDataKey")

  w <- .RGtkCall("S_cairo_pattern_get_user_data", pattern, key, PACKAGE = "RGtk2")

  return(w)
} 


cairoPatternSetUserData <-
function(pattern, key, user.data)
{
  checkPtrType(pattern, "CairoPattern")
  checkPtrType(key, "CairoUserDataKey")
  

  w <- .RGtkCall("S_cairo_pattern_set_user_data", pattern, key, user.data, PACKAGE = "RGtk2")

  return(w)
} 


cairoFormatStrideForWidth <-
function(format, width)
{
  
  width <- as.integer(width)

  w <- .RGtkCall("S_cairo_format_stride_for_width", format, width, PACKAGE = "RGtk2")

  return(w)
} 


cairoHasCurrentPoint <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_has_current_point", cr, PACKAGE = "RGtk2")

  return(w)
} 


cairoPathExtents <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_cairo_path_extents", cr, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoSurfaceCopyPage <-
function(surface)
{
  checkPtrType(surface, "CairoSurface")

  w <- .RGtkCall("S_cairo_surface_copy_page", surface, PACKAGE = "RGtk2")

  return(w)
} 


cairoSurfaceShowPage <-
function(surface)
{
  checkPtrType(surface, "CairoSurface")

  w <- .RGtkCall("S_cairo_surface_show_page", surface, PACKAGE = "RGtk2")

  return(w)
} 


cairoPsSurfaceRestrictToLevel <-
function(surface, level)
{
  checkPtrType(surface, "CairoSurface")
  

  w <- .RGtkCall("S_cairo_ps_surface_restrict_to_level", surface, level, PACKAGE = "RGtk2")

  return(w)
} 


cairoPsGetLevels <-
function()
{
  

  w <- .RGtkCall("S_cairo_ps_get_levels", PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoPsLevelToString <-
function(level)
{
  

  w <- .RGtkCall("S_cairo_ps_level_to_string", level, PACKAGE = "RGtk2")

  return(w)
} 


cairoPsSurfaceSetEps <-
function(surface, eps)
{
  checkPtrType(surface, "CairoSurface")
  eps <- as.logical(eps)

  w <- .RGtkCall("S_cairo_ps_surface_set_eps", surface, eps, PACKAGE = "RGtk2")

  return(w)
} 


cairoPsSurfaceGetEps <-
function(surface)
{
  checkPtrType(surface, "CairoSurface")

  w <- .RGtkCall("S_cairo_ps_surface_get_eps", surface, PACKAGE = "RGtk2")

  return(w)
} 


cairoToyFontFaceCreate <-
function(family, slant, weight)
{
  family <- as.character(family)
  
  

  w <- .RGtkCall("S_cairo_toy_font_face_create", family, slant, weight, PACKAGE = "RGtk2")

  return(w)
} 


cairoToyFontFaceGetFamily <-
function(font.face)
{
  checkPtrType(font.face, "CairoFontFace")

  w <- .RGtkCall("S_cairo_toy_font_face_get_family", font.face, PACKAGE = "RGtk2")

  return(w)
} 


cairoToyFontFaceGetSlant <-
function(font.face)
{
  checkPtrType(font.face, "CairoFontFace")

  w <- .RGtkCall("S_cairo_toy_font_face_get_slant", font.face, PACKAGE = "RGtk2")

  return(w)
} 


cairoToyFontFaceGetWeight <-
function(font.face)
{
  checkPtrType(font.face, "CairoFontFace")

  w <- .RGtkCall("S_cairo_toy_font_face_get_weight", font.face, PACKAGE = "RGtk2")

  return(w)
} 


cairoSurfaceGetFallbackResolution <-
function(surface)
{
  checkPtrType(surface, "CairoSurface")

  w <- .RGtkCall("S_cairo_surface_get_fallback_resolution", surface, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoSurfaceHasShowTextGlyphs <-
function(surface)
{
  checkPtrType(surface, "CairoSurface")

  w <- .RGtkCall("S_cairo_surface_has_show_text_glyphs", surface, PACKAGE = "RGtk2")

  return(w)
} 


cairoShowTextGlyphs <-
function(cr, utf8, glyphs, clusters, cluster.flags)
{
  checkPtrType(cr, "Cairo")
  utf8 <- as.character(utf8)
  glyphs <- lapply(glyphs, function(x) { x <- as.CairoGlyph(x); x })
  clusters <- lapply(clusters, function(x) { x <- as.CairoTextCluster(x); x })
  

  w <- .RGtkCall("S_cairo_show_text_glyphs", cr, utf8, glyphs, clusters, cluster.flags, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoScaledFontGetScaleMatrix <-
function(scaled.font)
{
  checkPtrType(scaled.font, "CairoScaledFont")

  w <- .RGtkCall("S_cairo_scaled_font_get_scale_matrix", scaled.font, PACKAGE = "RGtk2")

  return(invisible(w))
} 


cairoScaledFontTextToGlyphs <-
function(scaled.font, x, y, utf8, utf8.len = -1)
{
  checkPtrType(scaled.font, "CairoScaledFont")
  x <- as.numeric(x)
  y <- as.numeric(y)
  utf8 <- as.character(utf8)
  utf8.len <- as.integer(utf8.len)

  w <- .RGtkCall("S_cairo_scaled_font_text_to_glyphs", scaled.font, x, y, utf8, utf8.len, PACKAGE = "RGtk2")

  return(w)
} 


cairoUserFontFaceCreate <-
function()
{
  

  w <- .RGtkCall("S_cairo_user_font_face_create", PACKAGE = "RGtk2")

  return(w)
} 


cairoUserFontFaceSetInitFunc <-
function(font.face, init.func)
{
  checkPtrType(font.face, "CairoFontFace")
  init.func <- as.function(init.func)

  w <- .RGtkCall("S_cairo_user_font_face_set_init_func", font.face, init.func, PACKAGE = "RGtk2")

  return(w)
} 


cairoUserFontFaceSetRenderGlyphFunc <-
function(font.face, render.glyph.func)
{
  checkPtrType(font.face, "CairoFontFace")
  render.glyph.func <- as.function(render.glyph.func)

  w <- .RGtkCall("S_cairo_user_font_face_set_render_glyph_func", font.face, render.glyph.func, PACKAGE = "RGtk2")

  return(w)
} 


cairoUserFontFaceSetUnicodeToGlyphFunc <-
function(font.face, unicode.to.glyph.func)
{
  checkPtrType(font.face, "CairoFontFace")
  unicode.to.glyph.func <- as.function(unicode.to.glyph.func)

  w <- .RGtkCall("S_cairo_user_font_face_set_unicode_to_glyph_func", font.face, unicode.to.glyph.func, PACKAGE = "RGtk2")

  return(w)
} 


cairoUserFontFaceSetTextToGlyphsFunc <-
function(font.face, text.to.glyphs.func)
{
  checkPtrType(font.face, "CairoFontFace")
  text.to.glyphs.func <- as.function(text.to.glyphs.func)

  w <- .RGtkCall("S_cairo_user_font_face_set_text_to_glyphs_func", font.face, text.to.glyphs.func, PACKAGE = "RGtk2")

  return(w)
} 


cairoUserFontFaceGetInitFunc <-
function(font.face)
{
  checkPtrType(font.face, "CairoFontFace")

  w <- .RGtkCall("S_cairo_user_font_face_get_init_func", font.face, PACKAGE = "RGtk2")

  return(w)
} 


cairoUserFontFaceGetRenderGlyphFunc <-
function(font.face)
{
  checkPtrType(font.face, "CairoFontFace")

  w <- .RGtkCall("S_cairo_user_font_face_get_render_glyph_func", font.face, PACKAGE = "RGtk2")

  return(w)
} 


cairoUserFontFaceGetUnicodeToGlyphFunc <-
function(font.face)
{
  checkPtrType(font.face, "CairoFontFace")

  w <- .RGtkCall("S_cairo_user_font_face_get_unicode_to_glyph_func", font.face, PACKAGE = "RGtk2")

  return(w)
} 


cairoUserFontFaceGetTextToGlyphsFunc <-
function(font.face)
{
  checkPtrType(font.face, "CairoFontFace")

  w <- .RGtkCall("S_cairo_user_font_face_get_text_to_glyphs_func", font.face, PACKAGE = "RGtk2")

  return(w)
} 

