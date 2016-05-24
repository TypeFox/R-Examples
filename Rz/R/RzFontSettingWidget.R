fontSettingWidget <- 
setRefClass("RzFontSettingWidget",
  fields = c("title", "fontName", "fontFamily", "fontBox", "showSize", "showStyle"),
  methods = list(
    initialize  = function(...) {
      initFields(...)
      fontBox <<- gtkHBoxNew(spacing=5)
      font.label <- gtkLabelNew(title)
      font.button <- gtkFontButtonNew()
      font.button$setFontName(fontName)
      font.button$setUseFont(TRUE)
      font.button$setTitle(title)
      font.button$setShowSize(showSize)
      font.button$setUseSize(showSize)
#      font.button$setShowStyle(showStyle)
      fontBox$packStart(font.label, expand=FALSE)
      fontBox$packStart(font.button)
      fontBox$showAll()
      gSignalConnect(font.button, "font-set", .self$onSelectFont)
    },

    onSelectFont = function(button){
      fontName   <<- button$getFontName()
      fontFamily <<- pangoFontDescriptionFromString(fontName)$getFamily()
    }
))
fontSettingWidget$accessors(c("fontName", "fontFamily", "fontBox"))
