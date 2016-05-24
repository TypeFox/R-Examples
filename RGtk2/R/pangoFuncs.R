
pangoColorGetType <-
function()
{
  

  w <- .RGtkCall("S_pango_color_get_type", PACKAGE = "RGtk2")

  return(w)
} 


pangoColorCopy <-
function(object)
{
  checkPtrType(object, "PangoColor")

  w <- .RGtkCall("S_pango_color_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoColorFree <-
function(object)
{
  checkPtrType(object, "PangoColor")

  w <- .RGtkCall("S_pango_color_free", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoColorParse <-
function(spec)
{
  spec <- as.character(spec)

  w <- .RGtkCall("S_pango_color_parse", spec, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrTypeRegister <-
function(name)
{
  name <- as.character(name)

  w <- .RGtkCall("S_pango_attr_type_register", name, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttributeCopy <-
function(object)
{
  checkPtrType(object, "PangoAttribute")

  w <- .RGtkCall("S_pango_attribute_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttributeEqual <-
function(object, attr2)
{
  checkPtrType(object, "PangoAttribute")
  checkPtrType(attr2, "PangoAttribute")

  w <- .RGtkCall("S_pango_attribute_equal", object, attr2, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrLanguageNew <-
function(language)
{
  checkPtrType(language, "PangoLanguage")

  w <- .RGtkCall("S_pango_attr_language_new", language, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrFamilyNew <-
function(family)
{
  family <- as.character(family)

  w <- .RGtkCall("S_pango_attr_family_new", family, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrForegroundNew <-
function(red, green, blue)
{
  red <- as.integer(red)
  green <- as.integer(green)
  blue <- as.integer(blue)

  w <- .RGtkCall("S_pango_attr_foreground_new", red, green, blue, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrBackgroundNew <-
function(red, green, blue)
{
  red <- as.integer(red)
  green <- as.integer(green)
  blue <- as.integer(blue)

  w <- .RGtkCall("S_pango_attr_background_new", red, green, blue, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrStrikethroughColorNew <-
function(red, green, blue)
{
  red <- as.integer(red)
  green <- as.integer(green)
  blue <- as.integer(blue)

  w <- .RGtkCall("S_pango_attr_strikethrough_color_new", red, green, blue, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrUnderlineColorNew <-
function(red, green, blue)
{
  red <- as.integer(red)
  green <- as.integer(green)
  blue <- as.integer(blue)

  w <- .RGtkCall("S_pango_attr_underline_color_new", red, green, blue, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrSizeNew <-
function(size)
{
  size <- as.integer(size)

  w <- .RGtkCall("S_pango_attr_size_new", size, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrSizeNewAbsolute <-
function(size)
{
  size <- as.integer(size)

  w <- .RGtkCall("S_pango_attr_size_new_absolute", size, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrStyleNew <-
function(style)
{
  

  w <- .RGtkCall("S_pango_attr_style_new", style, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrWeightNew <-
function(weight)
{
  

  w <- .RGtkCall("S_pango_attr_weight_new", weight, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrVariantNew <-
function(variant)
{
  

  w <- .RGtkCall("S_pango_attr_variant_new", variant, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrStretchNew <-
function(stretch)
{
  

  w <- .RGtkCall("S_pango_attr_stretch_new", stretch, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrFontDescNew <-
function(desc)
{
  checkPtrType(desc, "PangoFontDescription")

  w <- .RGtkCall("S_pango_attr_font_desc_new", desc, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrUnderlineNew <-
function(underline)
{
  

  w <- .RGtkCall("S_pango_attr_underline_new", underline, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrStrikethroughNew <-
function(strikethrough)
{
  strikethrough <- as.logical(strikethrough)

  w <- .RGtkCall("S_pango_attr_strikethrough_new", strikethrough, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrRiseNew <-
function(rise)
{
  rise <- as.integer(rise)

  w <- .RGtkCall("S_pango_attr_rise_new", rise, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrShapeNew <-
function(ink.rect, logical.rect)
{
  ink.rect <- as.PangoRectangle(ink.rect)
  logical.rect <- as.PangoRectangle(logical.rect)

  w <- .RGtkCall("S_pango_attr_shape_new", ink.rect, logical.rect, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrShapeNewWithData <-
function(ink.rect, logical.rect, data)
{
  ink.rect <- as.PangoRectangle(ink.rect)
  logical.rect <- as.PangoRectangle(logical.rect)
  

  w <- .RGtkCall("S_pango_attr_shape_new_with_data", ink.rect, logical.rect, data, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrLetterSpacingNew <-
function(letter.spacing)
{
  letter.spacing <- as.integer(letter.spacing)

  w <- .RGtkCall("S_pango_attr_letter_spacing_new", letter.spacing, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrScaleNew <-
function(scale.factor)
{
  scale.factor <- as.numeric(scale.factor)

  w <- .RGtkCall("S_pango_attr_scale_new", scale.factor, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrFallbackNew <-
function(fallback)
{
  fallback <- as.logical(fallback)

  w <- .RGtkCall("S_pango_attr_fallback_new", fallback, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrListGetType <-
function()
{
  

  w <- .RGtkCall("S_pango_attr_list_get_type", PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrListNew <-
function()
{
  

  w <- .RGtkCall("S_pango_attr_list_new", PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrListCopy <-
function(object)
{
  checkPtrType(object, "PangoAttrList")

  w <- .RGtkCall("S_pango_attr_list_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrListInsert <-
function(object, attr)
{
  checkPtrType(object, "PangoAttrList")
  checkPtrType(attr, "PangoAttribute")

  w <- .RGtkCall("S_pango_attr_list_insert", object, attr, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoAttrListInsertBefore <-
function(object, attr)
{
  checkPtrType(object, "PangoAttrList")
  checkPtrType(attr, "PangoAttribute")

  w <- .RGtkCall("S_pango_attr_list_insert_before", object, attr, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoAttrListChange <-
function(object, attr)
{
  checkPtrType(object, "PangoAttrList")
  checkPtrType(attr, "PangoAttribute")

  w <- .RGtkCall("S_pango_attr_list_change", object, attr, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoAttrListSplice <-
function(object, other, pos, len)
{
  checkPtrType(object, "PangoAttrList")
  checkPtrType(other, "PangoAttrList")
  pos <- as.integer(pos)
  len <- as.integer(len)

  w <- .RGtkCall("S_pango_attr_list_splice", object, other, pos, len, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoAttrListGetIterator <-
function(object)
{
  checkPtrType(object, "PangoAttrList")

  w <- .RGtkCall("S_pango_attr_list_get_iterator", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrListFilter <-
function(object, func, data)
{
  checkPtrType(object, "PangoAttrList")
  func <- as.function(func)
  

  w <- .RGtkCall("S_pango_attr_list_filter", object, func, data, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrIteratorRange <-
function(object)
{
  checkPtrType(object, "PangoAttrIterator")

  w <- .RGtkCall("S_pango_attr_iterator_range", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoAttrIteratorNext <-
function(object)
{
  checkPtrType(object, "PangoAttrIterator")

  w <- .RGtkCall("S_pango_attr_iterator_next", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrIteratorCopy <-
function(object)
{
  checkPtrType(object, "PangoAttrIterator")

  w <- .RGtkCall("S_pango_attr_iterator_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrIteratorGet <-
function(object, type)
{
  checkPtrType(object, "PangoAttrIterator")
  

  w <- .RGtkCall("S_pango_attr_iterator_get", object, type, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrIteratorGetFont <-
function(object)
{
  checkPtrType(object, "PangoAttrIterator")

  w <- .RGtkCall("S_pango_attr_iterator_get_font", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoAttrIteratorGetAttrs <-
function(object)
{
  checkPtrType(object, "PangoAttrIterator")

  w <- .RGtkCall("S_pango_attr_iterator_get_attrs", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoParseMarkup <-
function(markup.text, length = -1, accel.marker = 0, .errwarn = TRUE)
{
  markup.text <- as.character(markup.text)
  length <- as.integer(length)
  accel.marker <- as.numeric(accel.marker)

  w <- .RGtkCall("S_pango_parse_markup", markup.text, length, accel.marker, PACKAGE = "RGtk2")

  w <- handleError(w, .errwarn)

  return(w)
} 


pangoFindParagraphBoundary <-
function(text, length = -1)
{
  text <- as.character(text)
  length <- as.integer(length)

  w <- .RGtkCall("S_pango_find_paragraph_boundary", text, length, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoCairoFontMapGetType <-
function()
{
  

  w <- .RGtkCall("S_pango_cairo_font_map_get_type", PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoFontMapNew <-
function()
{
  

  w <- .RGtkCall("S_pango_cairo_font_map_new", PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoFontMapGetDefault <-
function()
{
  

  w <- .RGtkCall("S_pango_cairo_font_map_get_default", PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoFontMapSetResolution <-
function(object, dpi)
{
  checkPtrType(object, "PangoCairoFontMap")
  dpi <- as.numeric(dpi)

  w <- .RGtkCall("S_pango_cairo_font_map_set_resolution", object, dpi, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoCairoFontMapGetResolution <-
function(object)
{
  checkPtrType(object, "PangoCairoFontMap")

  w <- .RGtkCall("S_pango_cairo_font_map_get_resolution", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoFontMapCreateContext <-
function(object)
{
  checkPtrType(object, "PangoCairoFontMap")

  w <- .RGtkCall("S_pango_cairo_font_map_create_context", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoUpdateContext <-
function(cr, context)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(context, "PangoContext")

  w <- .RGtkCall("S_pango_cairo_update_context", cr, context, PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoContextSetFontOptions <-
function(context, options)
{
  checkPtrType(context, "PangoContext")
  checkPtrType(options, "CairoFontOptions")

  w <- .RGtkCall("S_pango_cairo_context_set_font_options", context, options, PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoContextGetFontOptions <-
function(context)
{
  checkPtrType(context, "PangoContext")

  w <- .RGtkCall("S_pango_cairo_context_get_font_options", context, PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoContextSetResolution <-
function(context, dpi)
{
  checkPtrType(context, "PangoContext")
  dpi <- as.numeric(dpi)

  w <- .RGtkCall("S_pango_cairo_context_set_resolution", context, dpi, PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoContextGetResolution <-
function(context)
{
  checkPtrType(context, "PangoContext")

  w <- .RGtkCall("S_pango_cairo_context_get_resolution", context, PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoCreateLayout <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_pango_cairo_create_layout", cr, PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoUpdateLayout <-
function(cr, layout)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(layout, "PangoLayout")

  w <- .RGtkCall("S_pango_cairo_update_layout", cr, layout, PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoShowGlyphString <-
function(cr, font, glyphs)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(font, "PangoFont")
  checkPtrType(glyphs, "PangoGlyphString")

  w <- .RGtkCall("S_pango_cairo_show_glyph_string", cr, font, glyphs, PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoShowLayoutLine <-
function(cr, line)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(line, "PangoLayoutLine")

  w <- .RGtkCall("S_pango_cairo_show_layout_line", cr, line, PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoShowLayout <-
function(cr, layout)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(layout, "PangoLayout")

  w <- .RGtkCall("S_pango_cairo_show_layout", cr, layout, PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoGlyphStringPath <-
function(cr, font, glyphs)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(font, "PangoFont")
  checkPtrType(glyphs, "PangoGlyphString")

  w <- .RGtkCall("S_pango_cairo_glyph_string_path", cr, font, glyphs, PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoLayoutLinePath <-
function(cr, line)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(line, "PangoLayoutLine")

  w <- .RGtkCall("S_pango_cairo_layout_line_path", cr, line, PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoLayoutPath <-
function(cr, layout)
{
  checkPtrType(cr, "Cairo")
  checkPtrType(layout, "PangoLayout")

  w <- .RGtkCall("S_pango_cairo_layout_path", cr, layout, PACKAGE = "RGtk2")

  return(w)
} 


pangoContextSetFontMap <-
function(object, font.map)
{
  checkPtrType(object, "PangoContext")
  checkPtrType(font.map, "PangoFontMap")

  w <- .RGtkCall("S_pango_context_set_font_map", object, font.map, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoContextGetFontMap <-
function(object)
{
  checkPtrType(object, "PangoContext")

  w <- .RGtkCall("S_pango_context_get_font_map", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoContextListFamilies <-
function(object)
{
  checkPtrType(object, "PangoContext")

  w <- .RGtkCall("S_pango_context_list_families", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoGetMirrorChar <-
function(ch)
{
  if(getOption("depwarn"))
    .Deprecated("nothing", "RGtk2")

  ch <- as.numeric(ch)

  w <- .RGtkCall("S_pango_get_mirror_char", ch, PACKAGE = "RGtk2")

  return(w)
} 


pangoUnicharDirection <-
function(ch)
{
  ch <- as.numeric(ch)

  w <- .RGtkCall("S_pango_unichar_direction", ch, PACKAGE = "RGtk2")

  return(w)
} 


pangoFindBaseDir <-
function(text, length = -1)
{
  text <- as.character(text)
  length <- as.integer(length)

  w <- .RGtkCall("S_pango_find_base_dir", text, length, PACKAGE = "RGtk2")

  return(w)
} 


pangoContextLoadFont <-
function(object, desc)
{
  checkPtrType(object, "PangoContext")
  checkPtrType(desc, "PangoFontDescription")

  w <- .RGtkCall("S_pango_context_load_font", object, desc, PACKAGE = "RGtk2")

  return(w)
} 


pangoContextLoadFontset <-
function(object, desc, language)
{
  checkPtrType(object, "PangoContext")
  checkPtrType(desc, "PangoFontDescription")
  checkPtrType(language, "PangoLanguage")

  w <- .RGtkCall("S_pango_context_load_fontset", object, desc, language, PACKAGE = "RGtk2")

  return(w)
} 


pangoContextSetMatrix <-
function(object, matrix)
{
  checkPtrType(object, "PangoContext")
  checkPtrType(matrix, "PangoMatrix")

  w <- .RGtkCall("S_pango_context_set_matrix", object, matrix, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoContextGetMatrix <-
function(object)
{
  checkPtrType(object, "PangoContext")

  w <- .RGtkCall("S_pango_context_get_matrix", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoContextGetMetrics <-
function(object, desc, language = NULL)
{
  checkPtrType(object, "PangoContext")
  checkPtrType(desc, "PangoFontDescription")
  if (!is.null( language )) checkPtrType(language, "PangoLanguage")

  w <- .RGtkCall("S_pango_context_get_metrics", object, desc, language, PACKAGE = "RGtk2")

  return(w)
} 


pangoContextSetFontDescription <-
function(object, desc)
{
  checkPtrType(object, "PangoContext")
  checkPtrType(desc, "PangoFontDescription")

  w <- .RGtkCall("S_pango_context_set_font_description", object, desc, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoContextGetFontDescription <-
function(object)
{
  checkPtrType(object, "PangoContext")

  w <- .RGtkCall("S_pango_context_get_font_description", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoContextGetLanguage <-
function(object)
{
  checkPtrType(object, "PangoContext")

  w <- .RGtkCall("S_pango_context_get_language", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoContextSetLanguage <-
function(object, language)
{
  checkPtrType(object, "PangoContext")
  checkPtrType(language, "PangoLanguage")

  w <- .RGtkCall("S_pango_context_set_language", object, language, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoContextSetBaseDir <-
function(object, direction)
{
  checkPtrType(object, "PangoContext")
  

  w <- .RGtkCall("S_pango_context_set_base_dir", object, direction, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoContextGetBaseDir <-
function(object)
{
  checkPtrType(object, "PangoContext")

  w <- .RGtkCall("S_pango_context_get_base_dir", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoItemize <-
function(context, text, start.index, length, attrs, cached.iter = NULL)
{
  checkPtrType(context, "PangoContext")
  text <- as.character(text)
  start.index <- as.integer(start.index)
  length <- as.integer(length)
  checkPtrType(attrs, "PangoAttrList")
  if (!is.null( cached.iter )) checkPtrType(cached.iter, "PangoAttrIterator")

  w <- .RGtkCall("S_pango_itemize", context, text, start.index, length, attrs, cached.iter, PACKAGE = "RGtk2")

  return(w)
} 


pangoItemizeWithBaseDir <-
function(context, base.dir, text, start.index, length, attrs, cached.iter = NULL)
{
  checkPtrType(context, "PangoContext")
  
  text <- as.character(text)
  start.index <- as.integer(start.index)
  length <- as.integer(length)
  checkPtrType(attrs, "PangoAttrList")
  if (!is.null( cached.iter )) checkPtrType(cached.iter, "PangoAttrIterator")

  w <- .RGtkCall("S_pango_itemize_with_base_dir", context, base.dir, text, start.index, length, attrs, cached.iter, PACKAGE = "RGtk2")

  return(w)
} 


pangoCoverageNew <-
function()
{
  

  w <- .RGtkCall("S_pango_coverage_new", PACKAGE = "RGtk2")

  return(w)
} 


pangoCoverageCopy <-
function(object)
{
  checkPtrType(object, "PangoCoverage")

  w <- .RGtkCall("S_pango_coverage_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoCoverageGet <-
function(object, index)
{
  checkPtrType(object, "PangoCoverage")
  index <- as.integer(index)

  w <- .RGtkCall("S_pango_coverage_get", object, index, PACKAGE = "RGtk2")

  return(w)
} 


pangoCoverageSet <-
function(object, index, level)
{
  checkPtrType(object, "PangoCoverage")
  index <- as.integer(index)
  

  w <- .RGtkCall("S_pango_coverage_set", object, index, level, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoCoverageMax <-
function(object, other)
{
  checkPtrType(object, "PangoCoverage")
  checkPtrType(other, "PangoCoverage")

  w <- .RGtkCall("S_pango_coverage_max", object, other, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoCoverageToBytes <-
function(object)
{
  checkPtrType(object, "PangoCoverage")

  w <- .RGtkCall("S_pango_coverage_to_bytes", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoCoverageFromBytes <-
function(bytes)
{
  bytes <- as.list(as.raw(bytes))

  w <- .RGtkCall("S_pango_coverage_from_bytes", bytes, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontDescriptionNew <-
function()
{
  

  w <- .RGtkCall("S_pango_font_description_new", PACKAGE = "RGtk2")

  return(w)
} 


pangoFontDescriptionCopy <-
function(object)
{
  checkPtrType(object, "PangoFontDescription")

  w <- .RGtkCall("S_pango_font_description_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontDescriptionCopyStatic <-
function(object)
{
  checkPtrType(object, "PangoFontDescription")

  w <- .RGtkCall("S_pango_font_description_copy_static", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontDescriptionHash <-
function(object)
{
  checkPtrType(object, "PangoFontDescription")

  w <- .RGtkCall("S_pango_font_description_hash", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontDescriptionEqual <-
function(object, desc2)
{
  checkPtrType(object, "PangoFontDescription")
  checkPtrType(desc2, "PangoFontDescription")

  w <- .RGtkCall("S_pango_font_description_equal", object, desc2, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontDescriptionSetFamily <-
function(object, family)
{
  checkPtrType(object, "PangoFontDescription")
  family <- as.character(family)

  w <- .RGtkCall("S_pango_font_description_set_family", object, family, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoFontDescriptionSetFamilyStatic <-
function(object, family)
{
  checkPtrType(object, "PangoFontDescription")
  family <- as.character(family)

  w <- .RGtkCall("S_pango_font_description_set_family_static", object, family, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoFontDescriptionGetFamily <-
function(object)
{
  checkPtrType(object, "PangoFontDescription")

  w <- .RGtkCall("S_pango_font_description_get_family", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontDescriptionSetStyle <-
function(object, style)
{
  checkPtrType(object, "PangoFontDescription")
  

  w <- .RGtkCall("S_pango_font_description_set_style", object, style, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoFontDescriptionGetStyle <-
function(object)
{
  checkPtrType(object, "PangoFontDescription")

  w <- .RGtkCall("S_pango_font_description_get_style", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontDescriptionSetVariant <-
function(object, variant)
{
  checkPtrType(object, "PangoFontDescription")
  

  w <- .RGtkCall("S_pango_font_description_set_variant", object, variant, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoFontDescriptionGetVariant <-
function(object)
{
  checkPtrType(object, "PangoFontDescription")

  w <- .RGtkCall("S_pango_font_description_get_variant", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontDescriptionSetWeight <-
function(object, weight)
{
  checkPtrType(object, "PangoFontDescription")
  

  w <- .RGtkCall("S_pango_font_description_set_weight", object, weight, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoFontDescriptionGetWeight <-
function(object)
{
  checkPtrType(object, "PangoFontDescription")

  w <- .RGtkCall("S_pango_font_description_get_weight", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontDescriptionSetStretch <-
function(object, stretch)
{
  checkPtrType(object, "PangoFontDescription")
  

  w <- .RGtkCall("S_pango_font_description_set_stretch", object, stretch, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoFontDescriptionGetStretch <-
function(object)
{
  checkPtrType(object, "PangoFontDescription")

  w <- .RGtkCall("S_pango_font_description_get_stretch", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontDescriptionSetAbsoluteSize <-
function(object, size)
{
  checkPtrType(object, "PangoFontDescription")
  size <- as.numeric(size)

  w <- .RGtkCall("S_pango_font_description_set_absolute_size", object, size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoFontDescriptionGetSizeIsAbsolute <-
function(object)
{
  checkPtrType(object, "PangoFontDescription")

  w <- .RGtkCall("S_pango_font_description_get_size_is_absolute", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontDescriptionSetSize <-
function(object, size)
{
  checkPtrType(object, "PangoFontDescription")
  size <- as.integer(size)

  w <- .RGtkCall("S_pango_font_description_set_size", object, size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoFontDescriptionGetSize <-
function(object)
{
  checkPtrType(object, "PangoFontDescription")

  w <- .RGtkCall("S_pango_font_description_get_size", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontDescriptionGetSetFields <-
function(object)
{
  checkPtrType(object, "PangoFontDescription")

  w <- .RGtkCall("S_pango_font_description_get_set_fields", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontDescriptionUnsetFields <-
function(object, to.unset)
{
  checkPtrType(object, "PangoFontDescription")
  

  w <- .RGtkCall("S_pango_font_description_unset_fields", object, to.unset, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoFontDescriptionMerge <-
function(object, desc.to.merge, replace.existing)
{
  checkPtrType(object, "PangoFontDescription")
  checkPtrType(desc.to.merge, "PangoFontDescription")
  replace.existing <- as.logical(replace.existing)

  w <- .RGtkCall("S_pango_font_description_merge", object, desc.to.merge, replace.existing, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoFontDescriptionBetterMatch <-
function(object, old.match = NULL, new.match)
{
  checkPtrType(object, "PangoFontDescription")
  if (!is.null( old.match )) checkPtrType(old.match, "PangoFontDescription")
  checkPtrType(new.match, "PangoFontDescription")

  w <- .RGtkCall("S_pango_font_description_better_match", object, old.match, new.match, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontDescriptionFromString <-
function(str)
{
  str <- as.character(str)

  w <- .RGtkCall("S_pango_font_description_from_string", str, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontDescriptionToString <-
function(object)
{
  checkPtrType(object, "PangoFontDescription")

  w <- .RGtkCall("S_pango_font_description_to_string", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontDescriptionToFilename <-
function(object)
{
  checkPtrType(object, "PangoFontDescription")

  w <- .RGtkCall("S_pango_font_description_to_filename", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontMetricsGetType <-
function()
{
  

  w <- .RGtkCall("S_pango_font_metrics_get_type", PACKAGE = "RGtk2")

  return(w)
} 


pangoFontMetricsGetAscent <-
function(object)
{
  checkPtrType(object, "PangoFontMetrics")

  w <- .RGtkCall("S_pango_font_metrics_get_ascent", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontMetricsGetDescent <-
function(object)
{
  checkPtrType(object, "PangoFontMetrics")

  w <- .RGtkCall("S_pango_font_metrics_get_descent", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontMetricsGetApproximateCharWidth <-
function(object)
{
  checkPtrType(object, "PangoFontMetrics")

  w <- .RGtkCall("S_pango_font_metrics_get_approximate_char_width", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontMetricsGetApproximateDigitWidth <-
function(object)
{
  checkPtrType(object, "PangoFontMetrics")

  w <- .RGtkCall("S_pango_font_metrics_get_approximate_digit_width", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontMetricsGetStrikethroughPosition <-
function(object)
{
  checkPtrType(object, "PangoFontMetrics")

  w <- .RGtkCall("S_pango_font_metrics_get_strikethrough_position", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontMetricsGetStrikethroughThickness <-
function(object)
{
  checkPtrType(object, "PangoFontMetrics")

  w <- .RGtkCall("S_pango_font_metrics_get_strikethrough_thickness", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontMetricsGetUnderlinePosition <-
function(object)
{
  checkPtrType(object, "PangoFontMetrics")

  w <- .RGtkCall("S_pango_font_metrics_get_underline_position", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontMetricsGetUnderlineThickness <-
function(object)
{
  checkPtrType(object, "PangoFontMetrics")

  w <- .RGtkCall("S_pango_font_metrics_get_underline_thickness", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontFamilyGetType <-
function()
{
  

  w <- .RGtkCall("S_pango_font_family_get_type", PACKAGE = "RGtk2")

  return(w)
} 


pangoFontFamilyListFaces <-
function(object)
{
  checkPtrType(object, "PangoFontFamily")

  w <- .RGtkCall("S_pango_font_family_list_faces", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoFontFamilyGetName <-
function(object)
{
  checkPtrType(object, "PangoFontFamily")

  w <- .RGtkCall("S_pango_font_family_get_name", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontFamilyIsMonospace <-
function(object)
{
  checkPtrType(object, "PangoFontFamily")

  w <- .RGtkCall("S_pango_font_family_is_monospace", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontFaceGetType <-
function()
{
  

  w <- .RGtkCall("S_pango_font_face_get_type", PACKAGE = "RGtk2")

  return(w)
} 


pangoFontFaceDescribe <-
function(object)
{
  checkPtrType(object, "PangoFontFace")

  w <- .RGtkCall("S_pango_font_face_describe", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontFaceGetFaceName <-
function(object)
{
  checkPtrType(object, "PangoFontFace")

  w <- .RGtkCall("S_pango_font_face_get_face_name", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontFaceListSizes <-
function(object)
{
  checkPtrType(object, "PangoFontFace")

  w <- .RGtkCall("S_pango_font_face_list_sizes", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoFontGetType <-
function()
{
  

  w <- .RGtkCall("S_pango_font_get_type", PACKAGE = "RGtk2")

  return(w)
} 


pangoFontDescribe <-
function(object)
{
  checkPtrType(object, "PangoFont")

  w <- .RGtkCall("S_pango_font_describe", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontGetCoverage <-
function(object, language)
{
  checkPtrType(object, "PangoFont")
  checkPtrType(language, "PangoLanguage")

  w <- .RGtkCall("S_pango_font_get_coverage", object, language, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontGetMetrics <-
function(object, language = NULL)
{
  checkPtrType(object, "PangoFont")
  if (!is.null( language )) checkPtrType(language, "PangoLanguage")

  w <- .RGtkCall("S_pango_font_get_metrics", object, language, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontGetGlyphExtents <-
function(object, glyph)
{
  checkPtrType(object, "PangoFont")
  glyph <- as.numeric(glyph)

  w <- .RGtkCall("S_pango_font_get_glyph_extents", object, glyph, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoFontGetFontMap <-
function(object)
{
  checkPtrType(object, "PangoFont")

  w <- .RGtkCall("S_pango_font_get_font_map", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontMapLoadFont <-
function(object, context, desc)
{
  checkPtrType(object, "PangoFontMap")
  checkPtrType(context, "PangoContext")
  checkPtrType(desc, "PangoFontDescription")

  w <- .RGtkCall("S_pango_font_map_load_font", object, context, desc, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontMapLoadFontset <-
function(object, context, desc, language)
{
  checkPtrType(object, "PangoFontMap")
  checkPtrType(context, "PangoContext")
  checkPtrType(desc, "PangoFontDescription")
  checkPtrType(language, "PangoLanguage")

  w <- .RGtkCall("S_pango_font_map_load_fontset", object, context, desc, language, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontMapListFamilies <-
function(object)
{
  checkPtrType(object, "PangoFontMap")

  w <- .RGtkCall("S_pango_font_map_list_families", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoFontsetGetFont <-
function(object, wc)
{
  checkPtrType(object, "PangoFontset")
  wc <- as.numeric(wc)

  w <- .RGtkCall("S_pango_fontset_get_font", object, wc, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontsetGetMetrics <-
function(object)
{
  checkPtrType(object, "PangoFontset")

  w <- .RGtkCall("S_pango_fontset_get_metrics", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontsetForeach <-
function(object, func, data)
{
  checkPtrType(object, "PangoFontset")
  func <- as.function(func)
  

  w <- .RGtkCall("S_pango_fontset_foreach", object, func, data, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoGlyphStringNew <-
function()
{
  

  w <- .RGtkCall("S_pango_glyph_string_new", PACKAGE = "RGtk2")

  return(w)
} 


pangoGlyphStringSetSize <-
function(object, new.len)
{
  checkPtrType(object, "PangoGlyphString")
  new.len <- as.integer(new.len)

  w <- .RGtkCall("S_pango_glyph_string_set_size", object, new.len, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoGlyphStringGetType <-
function()
{
  

  w <- .RGtkCall("S_pango_glyph_string_get_type", PACKAGE = "RGtk2")

  return(w)
} 


pangoGlyphStringCopy <-
function(object)
{
  checkPtrType(object, "PangoGlyphString")

  w <- .RGtkCall("S_pango_glyph_string_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoGlyphStringExtents <-
function(object, font)
{
  checkPtrType(object, "PangoGlyphString")
  checkPtrType(font, "PangoFont")

  w <- .RGtkCall("S_pango_glyph_string_extents", object, font, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoGlyphStringExtentsRange <-
function(object, start, end, font)
{
  checkPtrType(object, "PangoGlyphString")
  start <- as.integer(start)
  end <- as.integer(end)
  checkPtrType(font, "PangoFont")

  w <- .RGtkCall("S_pango_glyph_string_extents_range", object, start, end, font, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoGlyphItemSplit <-
function(orig, text, split.index)
{
  checkPtrType(orig, "PangoGlyphItem")
  text <- as.character(text)
  split.index <- as.integer(split.index)

  w <- .RGtkCall("S_pango_glyph_item_split", orig, text, split.index, PACKAGE = "RGtk2")

  return(w)
} 


pangoGlyphItemApplyAttrs <-
function(glyph.item, text, list)
{
  checkPtrType(glyph.item, "PangoGlyphItem")
  text <- as.character(text)
  checkPtrType(list, "PangoAttrList")

  w <- .RGtkCall("S_pango_glyph_item_apply_attrs", glyph.item, text, list, PACKAGE = "RGtk2")

  return(w)
} 


pangoGlyphItemLetterSpace <-
function(glyph.item, text, log.attrs)
{
  checkPtrType(glyph.item, "PangoGlyphItem")
  text <- as.character(text)
  log.attrs <- lapply(log.attrs, function(x) { checkPtrType(x, "PangoLogAttr"); x })

  w <- .RGtkCall("S_pango_glyph_item_letter_space", glyph.item, text, log.attrs, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoMatrixTranslate <-
function(object, tx, ty)
{
  checkPtrType(object, "PangoMatrix")
  tx <- as.numeric(tx)
  ty <- as.numeric(ty)

  w <- .RGtkCall("S_pango_matrix_translate", object, tx, ty, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoMatrixScale <-
function(object, scale.x, scale.y)
{
  checkPtrType(object, "PangoMatrix")
  scale.x <- as.numeric(scale.x)
  scale.y <- as.numeric(scale.y)

  w <- .RGtkCall("S_pango_matrix_scale", object, scale.x, scale.y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoMatrixRotate <-
function(object, degrees)
{
  checkPtrType(object, "PangoMatrix")
  degrees <- as.numeric(degrees)

  w <- .RGtkCall("S_pango_matrix_rotate", object, degrees, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoMatrixConcat <-
function(object, new.matrix)
{
  checkPtrType(object, "PangoMatrix")
  checkPtrType(new.matrix, "PangoMatrix")

  w <- .RGtkCall("S_pango_matrix_concat", object, new.matrix, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoMatrixCopy <-
function(object)
{
  checkPtrType(object, "PangoMatrix")

  w <- .RGtkCall("S_pango_matrix_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoItemCopy <-
function(item)
{
  checkPtrType(item, "PangoItem")

  w <- .RGtkCall("S_pango_item_copy", item, PACKAGE = "RGtk2")

  return(w)
} 


pangoItemNew <-
function()
{
  

  w <- .RGtkCall("S_pango_item_new", PACKAGE = "RGtk2")

  return(w)
} 


pangoItemSplit <-
function(orig, split.index, split.offset)
{
  checkPtrType(orig, "PangoItem")
  split.index <- as.integer(split.index)
  split.offset <- as.integer(split.offset)

  w <- .RGtkCall("S_pango_item_split", orig, split.index, split.offset, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutGetType <-
function()
{
  

  w <- .RGtkCall("S_pango_layout_get_type", PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutNew <-
function(context)
{
  checkPtrType(context, "PangoContext")

  w <- .RGtkCall("S_pango_layout_new", context, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutCopy <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutGetContext <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_context", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutSetAttributes <-
function(object, attrs)
{
  checkPtrType(object, "PangoLayout")
  checkPtrType(attrs, "PangoAttrList")

  w <- .RGtkCall("S_pango_layout_set_attributes", object, attrs, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutGetAttributes <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_attributes", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutSetText <-
function(object, text, length = -1)
{
  checkPtrType(object, "PangoLayout")
  text <- as.character(text)
  length <- as.integer(length)

  w <- .RGtkCall("S_pango_layout_set_text", object, text, length, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutGetText <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_text", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutSetMarkup <-
function(object, markup, length = -1)
{
  checkPtrType(object, "PangoLayout")
  markup <- as.character(markup)
  length <- as.integer(length)

  w <- .RGtkCall("S_pango_layout_set_markup", object, markup, length, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutSetFontDescription <-
function(object, desc = NULL)
{
  checkPtrType(object, "PangoLayout")
  if (!is.null( desc )) checkPtrType(desc, "PangoFontDescription")

  w <- .RGtkCall("S_pango_layout_set_font_description", object, desc, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutGetFontDescription <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_font_description", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutSetWidth <-
function(object, width)
{
  checkPtrType(object, "PangoLayout")
  width <- as.integer(width)

  w <- .RGtkCall("S_pango_layout_set_width", object, width, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutGetWidth <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_width", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutSetWrap <-
function(object, wrap)
{
  checkPtrType(object, "PangoLayout")
  

  w <- .RGtkCall("S_pango_layout_set_wrap", object, wrap, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutGetWrap <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_wrap", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutSetIndent <-
function(object, indent)
{
  checkPtrType(object, "PangoLayout")
  indent <- as.integer(indent)

  w <- .RGtkCall("S_pango_layout_set_indent", object, indent, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutGetIndent <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_indent", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutSetSpacing <-
function(object, spacing)
{
  checkPtrType(object, "PangoLayout")
  spacing <- as.integer(spacing)

  w <- .RGtkCall("S_pango_layout_set_spacing", object, spacing, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutGetSpacing <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_spacing", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutSetJustify <-
function(object, justify)
{
  checkPtrType(object, "PangoLayout")
  justify <- as.logical(justify)

  w <- .RGtkCall("S_pango_layout_set_justify", object, justify, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutGetJustify <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_justify", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutSetAutoDir <-
function(object, auto.dir)
{
  checkPtrType(object, "PangoLayout")
  auto.dir <- as.logical(auto.dir)

  w <- .RGtkCall("S_pango_layout_set_auto_dir", object, auto.dir, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutGetAutoDir <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_auto_dir", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutSetAlignment <-
function(object, alignment)
{
  checkPtrType(object, "PangoLayout")
  

  w <- .RGtkCall("S_pango_layout_set_alignment", object, alignment, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutGetAlignment <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_alignment", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutSetTabs <-
function(object, tabs = NULL)
{
  checkPtrType(object, "PangoLayout")
  if (!is.null( tabs )) checkPtrType(tabs, "PangoTabArray")

  w <- .RGtkCall("S_pango_layout_set_tabs", object, tabs, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutGetTabs <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_tabs", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutSetSingleParagraphMode <-
function(object, setting)
{
  checkPtrType(object, "PangoLayout")
  setting <- as.logical(setting)

  w <- .RGtkCall("S_pango_layout_set_single_paragraph_mode", object, setting, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutGetSingleParagraphMode <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_single_paragraph_mode", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutSetEllipsize <-
function(object, ellipsize)
{
  checkPtrType(object, "PangoLayout")
  

  w <- .RGtkCall("S_pango_layout_set_ellipsize", object, ellipsize, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutGetEllipsize <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_ellipsize", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutContextChanged <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_context_changed", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutGetLogAttrs <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_log_attrs", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutIndexToPos <-
function(object, index, pos)
{
  checkPtrType(object, "PangoLayout")
  index <- as.integer(index)
  pos <- as.PangoRectangle(pos)

  w <- .RGtkCall("S_pango_layout_index_to_pos", object, index, pos, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutGetCursorPos <-
function(object, index)
{
  checkPtrType(object, "PangoLayout")
  index <- as.integer(index)

  w <- .RGtkCall("S_pango_layout_get_cursor_pos", object, index, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutMoveCursorVisually <-
function(object, strong, old.index, old.trailing, direction)
{
  checkPtrType(object, "PangoLayout")
  strong <- as.logical(strong)
  old.index <- as.integer(old.index)
  old.trailing <- as.integer(old.trailing)
  direction <- as.integer(direction)

  w <- .RGtkCall("S_pango_layout_move_cursor_visually", object, strong, old.index, old.trailing, direction, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutXyToIndex <-
function(object, x, y)
{
  checkPtrType(object, "PangoLayout")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_pango_layout_xy_to_index", object, x, y, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutGetExtents <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_extents", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutGetPixelExtents <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_pixel_extents", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutGetSize <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_size", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutGetPixelSize <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_pixel_size", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutGetLineCount <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_line_count", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutGetLine <-
function(object, line)
{
  checkPtrType(object, "PangoLayout")
  line <- as.integer(line)

  w <- .RGtkCall("S_pango_layout_get_line", object, line, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutGetLines <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_lines", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutLineXToIndex <-
function(object, x.pos)
{
  checkPtrType(object, "PangoLayoutLine")
  x.pos <- as.integer(x.pos)

  w <- .RGtkCall("S_pango_layout_line_x_to_index", object, x.pos, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutLineIndexToX <-
function(object, index, trailing)
{
  checkPtrType(object, "PangoLayoutLine")
  index <- as.integer(index)
  trailing <- as.logical(trailing)

  w <- .RGtkCall("S_pango_layout_line_index_to_x", object, index, trailing, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutLineGetXRanges <-
function(object, start.index, end.index)
{
  checkPtrType(object, "PangoLayoutLine")
  start.index <- as.integer(start.index)
  end.index <- as.integer(end.index)

  w <- .RGtkCall("S_pango_layout_line_get_x_ranges", object, start.index, end.index, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutLineGetExtents <-
function(object)
{
  checkPtrType(object, "PangoLayoutLine")

  w <- .RGtkCall("S_pango_layout_line_get_extents", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutLineGetPixelExtents <-
function(object)
{
  checkPtrType(object, "PangoLayoutLine")

  w <- .RGtkCall("S_pango_layout_line_get_pixel_extents", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutIterGetType <-
function()
{
  

  w <- .RGtkCall("S_pango_layout_iter_get_type", PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutGetIter <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_iter", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutIterGetIndex <-
function(object)
{
  checkPtrType(object, "PangoLayoutIter")

  w <- .RGtkCall("S_pango_layout_iter_get_index", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutIterGetRun <-
function(object)
{
  checkPtrType(object, "PangoLayoutIter")

  w <- .RGtkCall("S_pango_layout_iter_get_run", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutIterGetLine <-
function(object)
{
  checkPtrType(object, "PangoLayoutIter")

  w <- .RGtkCall("S_pango_layout_iter_get_line", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutIterAtLastLine <-
function(object)
{
  checkPtrType(object, "PangoLayoutIter")

  w <- .RGtkCall("S_pango_layout_iter_at_last_line", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutIterNextChar <-
function(object)
{
  checkPtrType(object, "PangoLayoutIter")

  w <- .RGtkCall("S_pango_layout_iter_next_char", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutIterNextCluster <-
function(object)
{
  checkPtrType(object, "PangoLayoutIter")

  w <- .RGtkCall("S_pango_layout_iter_next_cluster", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutIterNextRun <-
function(object)
{
  checkPtrType(object, "PangoLayoutIter")

  w <- .RGtkCall("S_pango_layout_iter_next_run", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutIterNextLine <-
function(object)
{
  checkPtrType(object, "PangoLayoutIter")

  w <- .RGtkCall("S_pango_layout_iter_next_line", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutIterGetCharExtents <-
function(object)
{
  checkPtrType(object, "PangoLayoutIter")

  w <- .RGtkCall("S_pango_layout_iter_get_char_extents", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutIterGetClusterExtents <-
function(object)
{
  checkPtrType(object, "PangoLayoutIter")

  w <- .RGtkCall("S_pango_layout_iter_get_cluster_extents", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutIterGetRunExtents <-
function(object)
{
  checkPtrType(object, "PangoLayoutIter")

  w <- .RGtkCall("S_pango_layout_iter_get_run_extents", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutIterGetLineExtents <-
function(object)
{
  checkPtrType(object, "PangoLayoutIter")

  w <- .RGtkCall("S_pango_layout_iter_get_line_extents", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutIterGetLineYrange <-
function(object)
{
  checkPtrType(object, "PangoLayoutIter")

  w <- .RGtkCall("S_pango_layout_iter_get_line_yrange", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutIterGetLayoutExtents <-
function(object)
{
  checkPtrType(object, "PangoLayoutIter")

  w <- .RGtkCall("S_pango_layout_iter_get_layout_extents", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoLayoutIterGetBaseline <-
function(object)
{
  checkPtrType(object, "PangoLayoutIter")

  w <- .RGtkCall("S_pango_layout_iter_get_baseline", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoRendererGetType <-
function()
{
  

  w <- .RGtkCall("S_pango_renderer_get_type", PACKAGE = "RGtk2")

  return(w)
} 


pangoRendererDrawLayout <-
function(object, layout, x, y)
{
  checkPtrType(object, "PangoRenderer")
  checkPtrType(layout, "PangoLayout")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_pango_renderer_draw_layout", object, layout, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoRendererDrawLayoutLine <-
function(object, line, x, y)
{
  checkPtrType(object, "PangoRenderer")
  checkPtrType(line, "PangoLayoutLine")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_pango_renderer_draw_layout_line", object, line, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoRendererDrawGlyphs <-
function(object, font, glyphs, x, y)
{
  checkPtrType(object, "PangoRenderer")
  checkPtrType(font, "PangoFont")
  checkPtrType(glyphs, "PangoGlyphString")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_pango_renderer_draw_glyphs", object, font, glyphs, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoRendererDrawRectangle <-
function(object, part, x, y, width, height)
{
  checkPtrType(object, "PangoRenderer")
  
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_pango_renderer_draw_rectangle", object, part, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoRendererDrawErrorUnderline <-
function(object, x, y, width, height)
{
  checkPtrType(object, "PangoRenderer")
  x <- as.integer(x)
  y <- as.integer(y)
  width <- as.integer(width)
  height <- as.integer(height)

  w <- .RGtkCall("S_pango_renderer_draw_error_underline", object, x, y, width, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoRendererDrawTrapezoid <-
function(object, part, y1., x11, x21, y2, x12, x22)
{
  checkPtrType(object, "PangoRenderer")
  
  y1. <- as.numeric(y1.)
  x11 <- as.numeric(x11)
  x21 <- as.numeric(x21)
  y2 <- as.numeric(y2)
  x12 <- as.numeric(x12)
  x22 <- as.numeric(x22)

  w <- .RGtkCall("S_pango_renderer_draw_trapezoid", object, part, y1., x11, x21, y2, x12, x22, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoRendererDrawGlyph <-
function(object, font, glyph, x, y)
{
  checkPtrType(object, "PangoRenderer")
  checkPtrType(font, "PangoFont")
  glyph <- as.numeric(glyph)
  x <- as.numeric(x)
  y <- as.numeric(y)

  w <- .RGtkCall("S_pango_renderer_draw_glyph", object, font, glyph, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoRendererActivate <-
function(object)
{
  checkPtrType(object, "PangoRenderer")

  w <- .RGtkCall("S_pango_renderer_activate", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoRendererDeactivate <-
function(object)
{
  checkPtrType(object, "PangoRenderer")

  w <- .RGtkCall("S_pango_renderer_deactivate", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoRendererPartChanged <-
function(object, part)
{
  checkPtrType(object, "PangoRenderer")
  

  w <- .RGtkCall("S_pango_renderer_part_changed", object, part, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoRendererSetColor <-
function(object, part, color)
{
  checkPtrType(object, "PangoRenderer")
  
  checkPtrType(color, "PangoColor")

  w <- .RGtkCall("S_pango_renderer_set_color", object, part, color, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoRendererGetColor <-
function(object, part)
{
  checkPtrType(object, "PangoRenderer")
  

  w <- .RGtkCall("S_pango_renderer_get_color", object, part, PACKAGE = "RGtk2")

  return(w)
} 


pangoRendererSetMatrix <-
function(object, matrix)
{
  checkPtrType(object, "PangoRenderer")
  checkPtrType(matrix, "PangoMatrix")

  w <- .RGtkCall("S_pango_renderer_set_matrix", object, matrix, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoRendererGetMatrix <-
function(object)
{
  checkPtrType(object, "PangoRenderer")

  w <- .RGtkCall("S_pango_renderer_get_matrix", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoTabArrayNew <-
function(initial.size, positions.in.pixels)
{
  initial.size <- as.integer(initial.size)
  positions.in.pixels <- as.logical(positions.in.pixels)

  w <- .RGtkCall("S_pango_tab_array_new", initial.size, positions.in.pixels, PACKAGE = "RGtk2")

  return(w)
} 


pangoTabArrayGetType <-
function()
{
  

  w <- .RGtkCall("S_pango_tab_array_get_type", PACKAGE = "RGtk2")

  return(w)
} 


pangoTabArrayCopy <-
function(object)
{
  checkPtrType(object, "PangoTabArray")

  w <- .RGtkCall("S_pango_tab_array_copy", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoTabArrayGetSize <-
function(object)
{
  checkPtrType(object, "PangoTabArray")

  w <- .RGtkCall("S_pango_tab_array_get_size", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoTabArrayResize <-
function(object, new.size)
{
  checkPtrType(object, "PangoTabArray")
  new.size <- as.integer(new.size)

  w <- .RGtkCall("S_pango_tab_array_resize", object, new.size, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoTabArraySetTab <-
function(object, tab.index, alignment, location)
{
  checkPtrType(object, "PangoTabArray")
  tab.index <- as.integer(tab.index)
  
  location <- as.integer(location)

  w <- .RGtkCall("S_pango_tab_array_set_tab", object, tab.index, alignment, location, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoTabArrayGetTab <-
function(object, tab.index)
{
  checkPtrType(object, "PangoTabArray")
  tab.index <- as.integer(tab.index)

  w <- .RGtkCall("S_pango_tab_array_get_tab", object, tab.index, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoTabArrayGetTabs <-
function(object)
{
  checkPtrType(object, "PangoTabArray")

  w <- .RGtkCall("S_pango_tab_array_get_tabs", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoTabArrayGetPositionsInPixels <-
function(object)
{
  checkPtrType(object, "PangoTabArray")

  w <- .RGtkCall("S_pango_tab_array_get_positions_in_pixels", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLanguageFromString <-
function(language)
{
  language <- as.character(language)

  w <- .RGtkCall("S_pango_language_from_string", language, PACKAGE = "RGtk2")

  return(w)
} 


pangoLanguageMatches <-
function(object, range.list)
{
  checkPtrType(object, "PangoLanguage")
  range.list <- as.character(range.list)

  w <- .RGtkCall("S_pango_language_matches", object, range.list, PACKAGE = "RGtk2")

  return(w)
} 


pangoLanguageToString <-
function(object)
{
  checkPtrType(object, "PangoLanguage")

  w <- .RGtkCall("S_pango_language_to_string", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoPixels <-
function(size)
{
  size <- as.integer(size)

  w <- .RGtkCall("S_PANGO_PIXELS", size, PACKAGE = "RGtk2")

  return(w)
} 


pangoAscent <-
function(rect)
{
  rect <- as.PangoRectangle(rect)

  w <- .RGtkCall("S_PANGO_ASCENT", rect, PACKAGE = "RGtk2")

  return(w)
} 


pangoDescent <-
function(rect)
{
  rect <- as.PangoRectangle(rect)

  w <- .RGtkCall("S_PANGO_DESCENT", rect, PACKAGE = "RGtk2")

  return(w)
} 


pangoLbearing <-
function(rect)
{
  rect <- as.PangoRectangle(rect)

  w <- .RGtkCall("S_PANGO_LBEARING", rect, PACKAGE = "RGtk2")

  return(w)
} 


pangoRbearing <-
function(rect)
{
  rect <- as.PangoRectangle(rect)

  w <- .RGtkCall("S_PANGO_RBEARING", rect, PACKAGE = "RGtk2")

  return(w)
} 


pangoScriptForUnichar <-
function(ch)
{
  ch <- as.numeric(ch)

  w <- .RGtkCall("S_pango_script_for_unichar", ch, PACKAGE = "RGtk2")

  return(w)
} 


pangoScriptIterNew <-
function(text, length)
{
  text <- as.character(text)
  length <- as.integer(length)

  w <- .RGtkCall("S_pango_script_iter_new", text, length, PACKAGE = "RGtk2")

  return(w)
} 


pangoScriptIterGetRange <-
function(object)
{
  checkPtrType(object, "PangoScriptIter")

  w <- .RGtkCall("S_pango_script_iter_get_range", object, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoScriptIterNext <-
function(object)
{
  checkPtrType(object, "PangoScriptIter")

  w <- .RGtkCall("S_pango_script_iter_next", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoScriptGetSampleLanguage <-
function(script)
{
  

  w <- .RGtkCall("S_pango_script_get_sample_language", script, PACKAGE = "RGtk2")

  return(w)
} 


pangoLanguageIncludesScript <-
function(object, script)
{
  checkPtrType(object, "PangoLanguage")
  

  w <- .RGtkCall("S_pango_language_includes_script", object, script, PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoShowErrorUnderline <-
function(cr, x, y, width, height)
{
  checkPtrType(cr, "Cairo")
  x <- as.numeric(x)
  y <- as.numeric(y)
  width <- as.numeric(width)
  height <- as.numeric(height)

  w <- .RGtkCall("S_pango_cairo_show_error_underline", cr, x, y, width, height, PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoErrorUnderlinePath <-
function(cr, x, y, width, height)
{
  checkPtrType(cr, "Cairo")
  x <- as.numeric(x)
  y <- as.numeric(y)
  width <- as.numeric(width)
  height <- as.numeric(height)

  w <- .RGtkCall("S_pango_cairo_error_underline_path", cr, x, y, width, height, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontDescribeWithAbsoluteSize <-
function(object)
{
  checkPtrType(object, "PangoFont")

  w <- .RGtkCall("S_pango_font_describe_with_absolute_size", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoGlyphStringGetWidth <-
function(object)
{
  checkPtrType(object, "PangoGlyphString")

  w <- .RGtkCall("S_pango_glyph_string_get_width", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoMatrixGetFontScaleFactor <-
function(object)
{
  checkPtrType(object, "PangoMatrix")

  w <- .RGtkCall("S_pango_matrix_get_font_scale_factor", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutIndexToLineX <-
function(object, index., trailing)
{
  checkPtrType(object, "PangoLayout")
  index. <- as.integer(index.)
  trailing <- as.logical(trailing)

  w <- .RGtkCall("S_pango_layout_index_to_line_x", object, index., trailing, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoGravityToRotation <-
function(base.gravity)
{
  

  w <- .RGtkCall("S_pango_gravity_to_rotation", base.gravity, PACKAGE = "RGtk2")

  return(w)
} 


pangoGravityGetForMatrix <-
function(matrix)
{
  checkPtrType(matrix, "PangoMatrix")

  w <- .RGtkCall("S_pango_gravity_get_for_matrix", matrix, PACKAGE = "RGtk2")

  return(w)
} 


pangoGravityGetForScript <-
function(script, base.gravity, hint)
{
  
  
  

  w <- .RGtkCall("S_pango_gravity_get_for_script", script, base.gravity, hint, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrGravityNew <-
function(gravity)
{
  

  w <- .RGtkCall("S_pango_attr_gravity_new", gravity, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrGravityHintNew <-
function(hint)
{
  

  w <- .RGtkCall("S_pango_attr_gravity_hint_new", hint, PACKAGE = "RGtk2")

  return(w)
} 


pangoContextSetBaseGravity <-
function(object, gravity)
{
  checkPtrType(object, "PangoContext")
  

  w <- .RGtkCall("S_pango_context_set_base_gravity", object, gravity, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoContextGetBaseGravity <-
function(object)
{
  checkPtrType(object, "PangoContext")

  w <- .RGtkCall("S_pango_context_get_base_gravity", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoContextGetGravity <-
function(object)
{
  checkPtrType(object, "PangoContext")

  w <- .RGtkCall("S_pango_context_get_gravity", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoContextSetGravityHint <-
function(object, hint)
{
  checkPtrType(object, "PangoContext")
  

  w <- .RGtkCall("S_pango_context_set_gravity_hint", object, hint, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoContextGetGravityHint <-
function(object)
{
  checkPtrType(object, "PangoContext")

  w <- .RGtkCall("S_pango_context_get_gravity_hint", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontDescriptionSetGravity <-
function(object, gravity)
{
  checkPtrType(object, "PangoFontDescription")
  

  w <- .RGtkCall("S_pango_font_description_set_gravity", object, gravity, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoFontDescriptionGetGravity <-
function(object)
{
  checkPtrType(object, "PangoFontDescription")

  w <- .RGtkCall("S_pango_font_description_get_gravity", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutGetLineReadonly <-
function(object, line)
{
  checkPtrType(object, "PangoLayout")
  line <- as.integer(line)

  w <- .RGtkCall("S_pango_layout_get_line_readonly", object, line, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutGetLinesReadonly <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_lines_readonly", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutIterGetLineReadonly <-
function(object)
{
  checkPtrType(object, "PangoLayoutIter")

  w <- .RGtkCall("S_pango_layout_iter_get_line_readonly", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutIterGetRunReadonly <-
function(object)
{
  checkPtrType(object, "PangoLayoutIter")

  w <- .RGtkCall("S_pango_layout_iter_get_run_readonly", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoColorToString <-
function(object)
{
  checkPtrType(object, "PangoColor")

  w <- .RGtkCall("S_pango_color_to_string", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoMatrixTransformPoint <-
function(object, x, y)
{
  checkPtrType(object, "PangoMatrix")
  x <- as.list(as.numeric(x))
  y <- as.list(as.numeric(y))

  w <- .RGtkCall("S_pango_matrix_transform_point", object, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoMatrixTransformDistance <-
function(object, dx, dy)
{
  checkPtrType(object, "PangoMatrix")
  dx <- as.list(as.numeric(dx))
  dy <- as.list(as.numeric(dy))

  w <- .RGtkCall("S_pango_matrix_transform_distance", object, dx, dy, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoMatrixTransformRectangle <-
function(object, rect)
{
  checkPtrType(object, "PangoMatrix")
  rect <- as.PangoRectangle(rect)

  w <- .RGtkCall("S_pango_matrix_transform_rectangle", object, rect, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoMatrixTransformPixelRectangle <-
function(object, rect)
{
  checkPtrType(object, "PangoMatrix")
  rect <- as.PangoRectangle(rect)

  w <- .RGtkCall("S_pango_matrix_transform_pixel_rectangle", object, rect, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoUnitsFromDouble <-
function(d)
{
  d <- as.numeric(d)

  w <- .RGtkCall("S_pango_units_from_double", d, PACKAGE = "RGtk2")

  return(w)
} 


pangoUnitsToDouble <-
function(i)
{
  i <- as.integer(i)

  w <- .RGtkCall("S_pango_units_to_double", i, PACKAGE = "RGtk2")

  return(w)
} 


pangoExtentsToPixels <-
function(inclusive, nearest)
{
  inclusive <- as.PangoRectangle(inclusive)
  nearest <- as.PangoRectangle(nearest)

  w <- .RGtkCall("S_pango_extents_to_pixels", inclusive, nearest, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutIsWrapped <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_is_wrapped", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutIsEllipsized <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_is_ellipsized", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutGetUnknownGlyphsCount <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_unknown_glyphs_count", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoVersion <-
function()
{
  

  w <- .RGtkCall("S_pango_version", PACKAGE = "RGtk2")

  return(w)
} 


pangoVersionString <-
function()
{
  

  w <- .RGtkCall("S_pango_version_string", PACKAGE = "RGtk2")

  return(w)
} 


pangoVersionCheck <-
function(required.major, required.minor, required.micro)
{
  required.major <- as.integer(required.major)
  required.minor <- as.integer(required.minor)
  required.micro <- as.integer(required.micro)

  w <- .RGtkCall("S_pango_version_check", required.major, required.minor, required.micro, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutGetHeight <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_height", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutSetHeight <-
function(object, height)
{
  checkPtrType(object, "PangoLayout")
  height <- as.integer(height)

  w <- .RGtkCall("S_pango_layout_set_height", object, height, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoAttributeInit <-
function(attr, klass)
{
  checkPtrType(attr, "PangoAttribute")
  checkPtrType(klass, "PangoAttrClass")

  w <- .RGtkCall("S_pango_attribute_init", attr, klass, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutIterGetLayout <-
function(iter)
{
  checkPtrType(iter, "PangoLayoutIter")

  w <- .RGtkCall("S_pango_layout_iter_get_layout", iter, PACKAGE = "RGtk2")

  return(w)
} 


pangoRendererGetLayout <-
function(renderer)
{
  checkPtrType(renderer, "PangoRenderer")

  w <- .RGtkCall("S_pango_renderer_get_layout", renderer, PACKAGE = "RGtk2")

  return(w)
} 


pangoRendererGetLayoutLine <-
function(renderer)
{
  checkPtrType(renderer, "PangoRenderer")

  w <- .RGtkCall("S_pango_renderer_get_layout_line", renderer, PACKAGE = "RGtk2")

  return(w)
} 


pangoFontFaceIsSynthesized <-
function(object)
{
  checkPtrType(object, "PangoFontFace")

  w <- .RGtkCall("S_pango_font_face_is_synthesized", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoFontGetScaledFont <-
function(object)
{
  checkPtrType(object, "PangoCairoFont")

  w <- .RGtkCall("S_pango_cairo_font_get_scaled_font", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoFontMapNewForFontType <-
function(fonttype)
{
  

  w <- .RGtkCall("S_pango_cairo_font_map_new_for_font_type", fonttype, PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoFontMapGetFontType <-
function(object)
{
  checkPtrType(object, "PangoCairoFontMap")

  w <- .RGtkCall("S_pango_cairo_font_map_get_font_type", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoContextSetShapeRenderer <-
function(object, func, data)
{
  checkPtrType(object, "PangoContext")
  func <- as.function(func)
  

  w <- .RGtkCall("S_pango_cairo_context_set_shape_renderer", object, func, data, PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoContextGetShapeRenderer <-
function(object)
{
  checkPtrType(object, "PangoContext")

  w <- .RGtkCall("S_pango_cairo_context_get_shape_renderer", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLanguageGetDefault <-
function()
{
  

  w <- .RGtkCall("S_pango_language_get_default", PACKAGE = "RGtk2")

  return(w)
} 


pangoLanguageGetSampleString <-
function(object)
{
  checkPtrType(object, "PangoLanguage")

  w <- .RGtkCall("S_pango_language_get_sample_string", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoBidiTypeForUnichar <-
function(ch)
{
  ch <- as.numeric(ch)

  w <- .RGtkCall("S_pango_bidi_type_for_unichar", ch, PACKAGE = "RGtk2")

  return(w)
} 


pangoAttrTypeGetName <-
function(type)
{
  

  w <- .RGtkCall("S_pango_attr_type_get_name", type, PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoCreateContext <-
function(cr)
{
  checkPtrType(cr, "Cairo")

  w <- .RGtkCall("S_pango_cairo_create_context", cr, PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoFontMapSetDefault <-
function(fontmap)
{
  checkPtrType(fontmap, "PangoCairoFontMap")

  w <- .RGtkCall("S_pango_cairo_font_map_set_default", fontmap, PACKAGE = "RGtk2")

  return(w)
} 


pangoCairoShowGlyphItem <-
function(cr, text, glyph.item)
{
  checkPtrType(cr, "Cairo")
  text <- as.character(text)
  checkPtrType(glyph.item, "PangoGlyphItem")

  w <- .RGtkCall("S_pango_cairo_show_glyph_item", cr, text, glyph.item, PACKAGE = "RGtk2")

  return(w)
} 


pangoRendererDrawGlyphItem <-
function(object, text, glyph.item, x, y)
{
  checkPtrType(object, "PangoRenderer")
  text <- as.character(text)
  checkPtrType(glyph.item, "PangoGlyphItem")
  x <- as.integer(x)
  y <- as.integer(y)

  w <- .RGtkCall("S_pango_renderer_draw_glyph_item", object, text, glyph.item, x, y, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoFontMapCreateContext <-
function(object)
{
  checkPtrType(object, "PangoFontMap")

  w <- .RGtkCall("S_pango_font_map_create_context", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoGlyphItemIterInitStart <-
function(object, glyph.item, text)
{
  checkPtrType(object, "PangoGlyphItemIter")
  checkPtrType(glyph.item, "PangoGlyphItem")
  text <- as.character(text)

  w <- .RGtkCall("S_pango_glyph_item_iter_init_start", object, glyph.item, text, PACKAGE = "RGtk2")

  return(w)
} 


pangoGlyphItemIterInitEnd <-
function(object, glyph.item, text)
{
  checkPtrType(object, "PangoGlyphItemIter")
  checkPtrType(glyph.item, "PangoGlyphItem")
  text <- as.character(text)

  w <- .RGtkCall("S_pango_glyph_item_iter_init_end", object, glyph.item, text, PACKAGE = "RGtk2")

  return(w)
} 


pangoGlyphItemIterNextCluster <-
function(object)
{
  checkPtrType(object, "PangoGlyphItemIter")

  w <- .RGtkCall("S_pango_glyph_item_iter_next_cluster", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoGlyphItemIterPrevCluster <-
function(object)
{
  checkPtrType(object, "PangoGlyphItemIter")

  w <- .RGtkCall("S_pango_glyph_item_iter_prev_cluster", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLanguageGetScripts <-
function(object)
{
  checkPtrType(object, "PangoLanguage")

  w <- .RGtkCall("S_pango_language_get_scripts", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoLayoutGetBaseline <-
function(object)
{
  checkPtrType(object, "PangoLayout")

  w <- .RGtkCall("S_pango_layout_get_baseline", object, PACKAGE = "RGtk2")

  return(w)
} 


pangoGlyphItemGetLogicalWidths <-
function(glyph.item, text)
{
  checkPtrType(glyph.item, "PangoGlyphItem")
  text <- as.character(text)

  w <- .RGtkCall("S_pango_glyph_item_get_logical_widths", glyph.item, text, PACKAGE = "RGtk2")

  return(invisible(w))
} 


pangoGravityGetForScriptAndWidth <-
function(script, wide, base.gravity, hint)
{
  
  wide <- as.logical(wide)
  
  

  w <- .RGtkCall("S_pango_gravity_get_for_script_and_width", script, wide, base.gravity, hint, PACKAGE = "RGtk2")

  return(w)
} 

