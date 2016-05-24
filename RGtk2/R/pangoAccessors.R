pangoColorGetRed <-
function(obj)
{
  checkPtrType(obj, 'PangoColor')
  v <- .Call('S_PangoColorGetRed', obj, PACKAGE = "RGtk2")
  v
} 
pangoColorGetGreen <-
function(obj)
{
  checkPtrType(obj, 'PangoColor')
  v <- .Call('S_PangoColorGetGreen', obj, PACKAGE = "RGtk2")
  v
} 
pangoColorGetBlue <-
function(obj)
{
  checkPtrType(obj, 'PangoColor')
  v <- .Call('S_PangoColorGetBlue', obj, PACKAGE = "RGtk2")
  v
} 
pangoGlyphStringGetNumGlyphs <-
function(obj)
{
  checkPtrType(obj, 'PangoGlyphString')
  v <- .Call('S_PangoGlyphStringGetNumGlyphs', obj, PACKAGE = "RGtk2")
  v
} 
pangoGlyphStringGetGlyphs <-
function(obj)
{
  checkPtrType(obj, 'PangoGlyphString')
  v <- .Call('S_PangoGlyphStringGetGlyphs', obj, PACKAGE = "RGtk2")
  v
} 
pangoGlyphStringGetLogClusters <-
function(obj)
{
  checkPtrType(obj, 'PangoGlyphString')
  v <- .Call('S_PangoGlyphStringGetLogClusters', obj, PACKAGE = "RGtk2")
  v
} 
pangoItemGetOffset <-
function(obj)
{
  checkPtrType(obj, 'PangoItem')
  v <- .Call('S_PangoItemGetOffset', obj, PACKAGE = "RGtk2")
  v
} 
pangoItemGetLength <-
function(obj)
{
  checkPtrType(obj, 'PangoItem')
  v <- .Call('S_PangoItemGetLength', obj, PACKAGE = "RGtk2")
  v
} 
pangoItemGetNumChars <-
function(obj)
{
  checkPtrType(obj, 'PangoItem')
  v <- .Call('S_PangoItemGetNumChars', obj, PACKAGE = "RGtk2")
  v
} 
pangoItemGetAnalysis <-
function(obj)
{
  checkPtrType(obj, 'PangoItem')
  v <- .Call('S_PangoItemGetAnalysis', obj, PACKAGE = "RGtk2")
  v
} 
pangoLayoutLineGetLayout <-
function(obj)
{
  checkPtrType(obj, 'PangoLayoutLine')
  v <- .Call('S_PangoLayoutLineGetLayout', obj, PACKAGE = "RGtk2")
  v
} 
pangoLayoutLineGetStartIndex <-
function(obj)
{
  checkPtrType(obj, 'PangoLayoutLine')
  v <- .Call('S_PangoLayoutLineGetStartIndex', obj, PACKAGE = "RGtk2")
  v
} 
pangoLayoutLineGetLength <-
function(obj)
{
  checkPtrType(obj, 'PangoLayoutLine')
  v <- .Call('S_PangoLayoutLineGetLength', obj, PACKAGE = "RGtk2")
  v
} 
pangoLayoutLineGetRuns <-
function(obj)
{
  checkPtrType(obj, 'PangoLayoutLine')
  v <- .Call('S_PangoLayoutLineGetRuns', obj, PACKAGE = "RGtk2")
  v
} 
pangoLayoutLineGetIsParagraphStart <-
function(obj)
{
  checkPtrType(obj, 'PangoLayoutLine')
  v <- .Call('S_PangoLayoutLineGetIsParagraphStart', obj, PACKAGE = "RGtk2")
  v
} 
pangoLayoutLineGetResolvedDir <-
function(obj)
{
  checkPtrType(obj, 'PangoLayoutLine')
  v <- .Call('S_PangoLayoutLineGetResolvedDir', obj, PACKAGE = "RGtk2")
  v
} 
pangoAnalysisGetFont <-
function(obj)
{
  checkPtrType(obj, 'PangoAnalysis')
  v <- .Call('S_PangoAnalysisGetFont', obj, PACKAGE = "RGtk2")
  v
} 
pangoAnalysisGetLevel <-
function(obj)
{
  checkPtrType(obj, 'PangoAnalysis')
  v <- .Call('S_PangoAnalysisGetLevel', obj, PACKAGE = "RGtk2")
  v
} 
pangoAnalysisGetLanguage <-
function(obj)
{
  checkPtrType(obj, 'PangoAnalysis')
  v <- .Call('S_PangoAnalysisGetLanguage', obj, PACKAGE = "RGtk2")
  v
} 
pangoAnalysisGetExtraAttrs <-
function(obj)
{
  checkPtrType(obj, 'PangoAnalysis')
  v <- .Call('S_PangoAnalysisGetExtraAttrs', obj, PACKAGE = "RGtk2")
  v
} 
pangoLogAttrGetIsLineBreak <-
function(obj)
{
  checkPtrType(obj, 'PangoLogAttr')
  v <- .Call('S_PangoLogAttrGetIsLineBreak', obj, PACKAGE = "RGtk2")
  v
} 
pangoLogAttrGetIsMandatoryBreak <-
function(obj)
{
  checkPtrType(obj, 'PangoLogAttr')
  v <- .Call('S_PangoLogAttrGetIsMandatoryBreak', obj, PACKAGE = "RGtk2")
  v
} 
pangoLogAttrGetIsCharBreak <-
function(obj)
{
  checkPtrType(obj, 'PangoLogAttr')
  v <- .Call('S_PangoLogAttrGetIsCharBreak', obj, PACKAGE = "RGtk2")
  v
} 
pangoLogAttrGetIsWhite <-
function(obj)
{
  checkPtrType(obj, 'PangoLogAttr')
  v <- .Call('S_PangoLogAttrGetIsWhite', obj, PACKAGE = "RGtk2")
  v
} 
pangoLogAttrGetIsCursorPosition <-
function(obj)
{
  checkPtrType(obj, 'PangoLogAttr')
  v <- .Call('S_PangoLogAttrGetIsCursorPosition', obj, PACKAGE = "RGtk2")
  v
} 
pangoLogAttrGetIsWordStart <-
function(obj)
{
  checkPtrType(obj, 'PangoLogAttr')
  v <- .Call('S_PangoLogAttrGetIsWordStart', obj, PACKAGE = "RGtk2")
  v
} 
pangoLogAttrGetIsWordEnd <-
function(obj)
{
  checkPtrType(obj, 'PangoLogAttr')
  v <- .Call('S_PangoLogAttrGetIsWordEnd', obj, PACKAGE = "RGtk2")
  v
} 
pangoLogAttrGetIsSentenceBoundary <-
function(obj)
{
  checkPtrType(obj, 'PangoLogAttr')
  v <- .Call('S_PangoLogAttrGetIsSentenceBoundary', obj, PACKAGE = "RGtk2")
  v
} 
pangoLogAttrGetIsSentenceStart <-
function(obj)
{
  checkPtrType(obj, 'PangoLogAttr')
  v <- .Call('S_PangoLogAttrGetIsSentenceStart', obj, PACKAGE = "RGtk2")
  v
} 
pangoLogAttrGetIsSentenceEnd <-
function(obj)
{
  checkPtrType(obj, 'PangoLogAttr')
  v <- .Call('S_PangoLogAttrGetIsSentenceEnd', obj, PACKAGE = "RGtk2")
  v
} 
pangoLogAttrGetBackspaceDeletesCharacter <-
function(obj)
{
  checkPtrType(obj, 'PangoLogAttr')
  v <- .Call('S_PangoLogAttrGetBackspaceDeletesCharacter', obj, PACKAGE = "RGtk2")
  v
} 
pangoAttributeGetKlass <-
function(obj)
{
  checkPtrType(obj, 'PangoAttribute')
  v <- .Call('S_PangoAttributeGetKlass', obj, PACKAGE = "RGtk2")
  v
} 
pangoAttributeGetStartIndex <-
function(obj)
{
  checkPtrType(obj, 'PangoAttribute')
  v <- .Call('S_PangoAttributeGetStartIndex', obj, PACKAGE = "RGtk2")
  v
} 
pangoAttributeGetEndIndex <-
function(obj)
{
  checkPtrType(obj, 'PangoAttribute')
  v <- .Call('S_PangoAttributeGetEndIndex', obj, PACKAGE = "RGtk2")
  v
} 
pangoAttrClassGetType <-
function(obj)
{
  checkPtrType(obj, 'PangoAttrClass')
  v <- .Call('S_PangoAttrClassGetType', obj, PACKAGE = "RGtk2")
  v
} 
pangoGlyphItemGetItem <-
function(obj)
{
  checkPtrType(obj, 'PangoGlyphItem')
  v <- .Call('S_PangoGlyphItemGetItem', obj, PACKAGE = "RGtk2")
  v
} 
pangoGlyphItemGetGlyphs <-
function(obj)
{
  checkPtrType(obj, 'PangoGlyphItem')
  v <- .Call('S_PangoGlyphItemGetGlyphs', obj, PACKAGE = "RGtk2")
  v
} 
pangoGlyphInfoGetGlyph <-
function(obj)
{
  checkPtrType(obj, 'PangoGlyphInfo')
  v <- .Call('S_PangoGlyphInfoGetGlyph', obj, PACKAGE = "RGtk2")
  v
} 
pangoGlyphInfoGetGeometry <-
function(obj)
{
  checkPtrType(obj, 'PangoGlyphInfo')
  v <- .Call('S_PangoGlyphInfoGetGeometry', obj, PACKAGE = "RGtk2")
  v
} 
pangoGlyphInfoGetAttr <-
function(obj)
{
  checkPtrType(obj, 'PangoGlyphInfo')
  v <- .Call('S_PangoGlyphInfoGetAttr', obj, PACKAGE = "RGtk2")
  v
} 
pangoGlyphGeometryGetWidth <-
function(obj)
{
  checkPtrType(obj, 'PangoGlyphGeometry')
  v <- .Call('S_PangoGlyphGeometryGetWidth', obj, PACKAGE = "RGtk2")
  v
} 
pangoGlyphGeometryGetXOffset <-
function(obj)
{
  checkPtrType(obj, 'PangoGlyphGeometry')
  v <- .Call('S_PangoGlyphGeometryGetXOffset', obj, PACKAGE = "RGtk2")
  v
} 
pangoGlyphGeometryGetYOffset <-
function(obj)
{
  checkPtrType(obj, 'PangoGlyphGeometry')
  v <- .Call('S_PangoGlyphGeometryGetYOffset', obj, PACKAGE = "RGtk2")
  v
} 
pangoGlyphVisAttrGetIsClusterStart <-
function(obj)
{
  checkPtrType(obj, 'PangoGlyphVisAttr')
  v <- .Call('S_PangoGlyphVisAttrGetIsClusterStart', obj, PACKAGE = "RGtk2")
  v
} 
pangoAttrStringGetValue <-
function(obj)
{
  checkPtrType(obj, 'PangoAttrString')
  v <- .Call('S_PangoAttrStringGetValue', obj, PACKAGE = "RGtk2")
  v
} 
pangoAttrLanguageGetValue <-
function(obj)
{
  checkPtrType(obj, 'PangoAttrLanguage')
  v <- .Call('S_PangoAttrLanguageGetValue', obj, PACKAGE = "RGtk2")
  v
} 
pangoAttrColorGetColor <-
function(obj)
{
  checkPtrType(obj, 'PangoAttrColor')
  v <- .Call('S_PangoAttrColorGetColor', obj, PACKAGE = "RGtk2")
  v
} 
pangoAttrIntGetValue <-
function(obj)
{
  checkPtrType(obj, 'PangoAttrInt')
  v <- .Call('S_PangoAttrIntGetValue', obj, PACKAGE = "RGtk2")
  v
} 
pangoAttrFloatGetValue <-
function(obj)
{
  checkPtrType(obj, 'PangoAttrFloat')
  v <- .Call('S_PangoAttrFloatGetValue', obj, PACKAGE = "RGtk2")
  v
} 
pangoAttrFontDescGetDesc <-
function(obj)
{
  checkPtrType(obj, 'PangoAttrFontDesc')
  v <- .Call('S_PangoAttrFontDescGetDesc', obj, PACKAGE = "RGtk2")
  v
} 
pangoAttrShapeGetInkRect <-
function(obj)
{
  checkPtrType(obj, 'PangoAttrShape')
  v <- .Call('S_PangoAttrShapeGetInkRect', obj, PACKAGE = "RGtk2")
  v
} 
pangoAttrShapeGetLogicalRect <-
function(obj)
{
  checkPtrType(obj, 'PangoAttrShape')
  v <- .Call('S_PangoAttrShapeGetLogicalRect', obj, PACKAGE = "RGtk2")
  v
} 
pangoAttrSizeGetSize <-
function(obj)
{
  checkPtrType(obj, 'PangoAttrSize')
  v <- .Call('S_PangoAttrSizeGetSize', obj, PACKAGE = "RGtk2")
  v
} 
pangoAttrSizeGetAbsolute <-
function(obj)
{
  checkPtrType(obj, 'PangoAttrSize')
  v <- .Call('S_PangoAttrSizeGetAbsolute', obj, PACKAGE = "RGtk2")
  v
} 
