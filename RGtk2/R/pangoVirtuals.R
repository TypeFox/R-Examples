if(!exists('.virtuals')) .virtuals <- new.env()
assign("PangoAttr", c("copy", "destroy", "equal"), .virtuals)
assign("PangoFontFamily", c("list_faces", "get_name", "is_monospace"), .virtuals)
assign("PangoFontFace", c("get_face_name", "describe", "list_sizes"), .virtuals)
assign("PangoFont", c("describe", "get_coverage", "get_glyph_extents", "get_metrics", "get_font_map"), .virtuals)
assign("PangoFontMap", c("load_font", "list_families", "load_fontset"), .virtuals)
assign("PangoFontset", c("get_font", "get_metrics", "get_language", "foreach"), .virtuals)
assign("PangoRenderer", c("draw_glyphs", "draw_rectangle", "draw_error_underline", "draw_shape", "draw_trapezoid", "draw_glyph", "part_changed", "begin", "end", "prepare_run", "draw_glyph_item"), .virtuals)


pangoAttrClassCopy <-
function(object.class, object)
{
  checkPtrType(object.class, "PangoAttrClass")
  checkPtrType(object, "PangoAttr")

  w <- .RGtkCall("S_pango_attr_class_copy", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

pangoAttrClassDestroy <-
function(object.class, object)
{
  checkPtrType(object.class, "PangoAttrClass")
  checkPtrType(object, "PangoAttr")

  w <- .RGtkCall("S_pango_attr_class_destroy", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

pangoAttrClassEqual <-
function(object.class, object, attr2)
{
  checkPtrType(object.class, "PangoAttrClass")
  checkPtrType(object, "PangoAttr")
  checkPtrType(attr2, "PangoAttribute")

  w <- .RGtkCall("S_pango_attr_class_equal", object.class, object, attr2, PACKAGE = "RGtk2")

  return(w)
}

pangoFontFamilyClassListFaces <-
function(object.class, object)
{
  checkPtrType(object.class, "PangoFontFamilyClass")
  checkPtrType(object, "PangoFontFamily")

  w <- .RGtkCall("S_pango_font_family_class_list_faces", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

pangoFontFamilyClassGetName <-
function(object.class, object)
{
  checkPtrType(object.class, "PangoFontFamilyClass")
  checkPtrType(object, "PangoFontFamily")

  w <- .RGtkCall("S_pango_font_family_class_get_name", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

pangoFontFamilyClassIsMonospace <-
function(object.class, object)
{
  checkPtrType(object.class, "PangoFontFamilyClass")
  checkPtrType(object, "PangoFontFamily")

  w <- .RGtkCall("S_pango_font_family_class_is_monospace", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

pangoFontFaceClassGetFaceName <-
function(object.class, object)
{
  checkPtrType(object.class, "PangoFontFaceClass")
  checkPtrType(object, "PangoFontFace")

  w <- .RGtkCall("S_pango_font_face_class_get_face_name", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

pangoFontFaceClassDescribe <-
function(object.class, object)
{
  checkPtrType(object.class, "PangoFontFaceClass")
  checkPtrType(object, "PangoFontFace")

  w <- .RGtkCall("S_pango_font_face_class_describe", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

pangoFontFaceClassListSizes <-
function(object.class, object)
{
  checkPtrType(object.class, "PangoFontFaceClass")
  checkPtrType(object, "PangoFontFace")

  w <- .RGtkCall("S_pango_font_face_class_list_sizes", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

pangoFontClassDescribe <-
function(object.class, object)
{
  checkPtrType(object.class, "PangoFontClass")
  checkPtrType(object, "PangoFont")

  w <- .RGtkCall("S_pango_font_class_describe", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

pangoFontClassGetCoverage <-
function(object.class, object, lang)
{
  checkPtrType(object.class, "PangoFontClass")
  checkPtrType(object, "PangoFont")
  checkPtrType(lang, "PangoLanguage")

  w <- .RGtkCall("S_pango_font_class_get_coverage", object.class, object, lang, PACKAGE = "RGtk2")

  return(w)
}

pangoFontClassGetGlyphExtents <-
function(object.class, object, glyph)
{
  checkPtrType(object.class, "PangoFontClass")
  checkPtrType(object, "PangoFont")
  glyph <- as.numeric(glyph)

  w <- .RGtkCall("S_pango_font_class_get_glyph_extents", object.class, object, glyph, PACKAGE = "RGtk2")

  return(w)
}

pangoFontClassGetMetrics <-
function(object.class, object, language)
{
  checkPtrType(object.class, "PangoFontClass")
  checkPtrType(object, "PangoFont")
  checkPtrType(language, "PangoLanguage")

  w <- .RGtkCall("S_pango_font_class_get_metrics", object.class, object, language, PACKAGE = "RGtk2")

  return(w)
}

pangoFontClassGetFontMap <-
function(object.class, object)
{
  checkPtrType(object.class, "PangoFontClass")
  checkPtrType(object, "PangoFont")

  w <- .RGtkCall("S_pango_font_class_get_font_map", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

pangoFontMapClassLoadFont <-
function(object.class, object, context, desc)
{
  checkPtrType(object.class, "PangoFontMapClass")
  checkPtrType(object, "PangoFontMap")
  checkPtrType(context, "PangoContext")
  checkPtrType(desc, "PangoFontDescription")

  w <- .RGtkCall("S_pango_font_map_class_load_font", object.class, object, context, desc, PACKAGE = "RGtk2")

  return(w)
}

pangoFontMapClassListFamilies <-
function(object.class, object)
{
  checkPtrType(object.class, "PangoFontMapClass")
  checkPtrType(object, "PangoFontMap")

  w <- .RGtkCall("S_pango_font_map_class_list_families", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

pangoFontMapClassLoadFontset <-
function(object.class, object, context, desc, language)
{
  checkPtrType(object.class, "PangoFontMapClass")
  checkPtrType(object, "PangoFontMap")
  checkPtrType(context, "PangoContext")
  checkPtrType(desc, "PangoFontDescription")
  checkPtrType(language, "PangoLanguage")

  w <- .RGtkCall("S_pango_font_map_class_load_fontset", object.class, object, context, desc, language, PACKAGE = "RGtk2")

  return(w)
}

pangoFontsetClassGetFont <-
function(object.class, object, wc)
{
  checkPtrType(object.class, "PangoFontsetClass")
  checkPtrType(object, "PangoFontset")
  wc <- as.numeric(wc)

  w <- .RGtkCall("S_pango_fontset_class_get_font", object.class, object, wc, PACKAGE = "RGtk2")

  return(w)
}

pangoFontsetClassGetMetrics <-
function(object.class, object)
{
  checkPtrType(object.class, "PangoFontsetClass")
  checkPtrType(object, "PangoFontset")

  w <- .RGtkCall("S_pango_fontset_class_get_metrics", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

pangoFontsetClassGetLanguage <-
function(object.class, object)
{
  checkPtrType(object.class, "PangoFontsetClass")
  checkPtrType(object, "PangoFontset")

  w <- .RGtkCall("S_pango_fontset_class_get_language", object.class, object, PACKAGE = "RGtk2")

  return(w)
}

pangoFontsetClassForeach <-
function(object.class, object, func, data)
{
  checkPtrType(object.class, "PangoFontsetClass")
  checkPtrType(object, "PangoFontset")
  func <- as.function(func)
  

  w <- .RGtkCall("S_pango_fontset_class_foreach", object.class, object, func, data, PACKAGE = "RGtk2")

  return(invisible(w))
}

pangoRendererClassDrawGlyphs <-
function(object.class, object, font, glyphs, x, y)
{
  checkPtrType(object.class, "PangoRendererClass")
  checkPtrType(object, "PangoRenderer")
  checkPtrType(font, "PangoFont")
  checkPtrType(glyphs, "PangoGlyphString")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_pango_renderer_class_draw_glyphs", object.class, object, font, glyphs, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
}

pangoRendererClassDrawRectangle <-
function(object.class, object, part, x, y, width, height)
{
  checkPtrType(object.class, "PangoRendererClass")
  checkPtrType(object, "PangoRenderer")
  
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_pango_renderer_class_draw_rectangle", object.class, object, part, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
}

pangoRendererClassDrawErrorUnderline <-
function(object.class, object, x, y, width, height)
{
  checkPtrType(object.class, "PangoRendererClass")
  checkPtrType(object, "PangoRenderer")
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_pango_renderer_class_draw_error_underline", object.class, object, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
}

pangoRendererClassDrawShape <-
function(object.class, object, attr, x, y)
{
  checkPtrType(object.class, "PangoRendererClass")
  checkPtrType(object, "PangoRenderer")
  checkPtrType(attr, "PangoAttrShape")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_pango_renderer_class_draw_shape", object.class, object, attr, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
}

pangoRendererClassDrawTrapezoid <-
function(object.class, object, part, y1., x11, x21, y2, x12, x22)
{
  checkPtrType(object.class, "PangoRendererClass")
  checkPtrType(object, "PangoRenderer")
  
  y1. <- as.numeric(y1.)
  x11 <- as.numeric(x11)
  x21 <- as.numeric(x21)
  y2 <- as.numeric(y2)
  x12 <- as.numeric(x12)
  x22 <- as.numeric(x22)

  w <- .RGtkCall("S_pango_renderer_class_draw_trapezoid", object.class, object, part, y1., x11, x21, y2, x12, x22, PACKAGE = "RGtk2")

  return(invisible(w))
}

pangoRendererClassDrawGlyph <-
function(object.class, object, font, glyph, x, y)
{
  checkPtrType(object.class, "PangoRendererClass")
  checkPtrType(object, "PangoRenderer")
  checkPtrType(font, "PangoFont")
  glyph <- as.numeric(glyph)
  x <- as.numeric(x)
  y <- as.numeric(y)

  w <- .RGtkCall("S_pango_renderer_class_draw_glyph", object.class, object, font, glyph, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
}

pangoRendererClassPartChanged <-
function(object.class, object, part)
{
  checkPtrType(object.class, "PangoRendererClass")
  checkPtrType(object, "PangoRenderer")
  

  w <- .RGtkCall("S_pango_renderer_class_part_changed", object.class, object, part, PACKAGE = "RGtk2")

  return(invisible(w))
}

pangoRendererClassBegin <-
function(object.class, object)
{
  checkPtrType(object.class, "PangoRendererClass")
  checkPtrType(object, "PangoRenderer")

  w <- .RGtkCall("S_pango_renderer_class_begin", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

pangoRendererClassEnd <-
function(object.class, object)
{
  checkPtrType(object.class, "PangoRendererClass")
  checkPtrType(object, "PangoRenderer")

  w <- .RGtkCall("S_pango_renderer_class_end", object.class, object, PACKAGE = "RGtk2")

  return(invisible(w))
}

pangoRendererClassPrepareRun <-
function(object.class, object, run)
{
  checkPtrType(object.class, "PangoRendererClass")
  checkPtrType(object, "PangoRenderer")
  checkPtrType(run, "PangoGlyphItem")

  w <- .RGtkCall("S_pango_renderer_class_prepare_run", object.class, object, run, PACKAGE = "RGtk2")

  return(invisible(w))
}

pangoRendererClassDrawGlyphItem <-
function(object.class, object, text, glyph.item, x, y)
{
  checkPtrType(object.class, "PangoRendererClass")
  checkPtrType(object, "PangoRenderer")
  text <- as.character(text)
  checkPtrType(glyph.item, "PangoGlyphItem")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_pango_renderer_class_draw_glyph_item", object.class, object, text, glyph.item, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
}