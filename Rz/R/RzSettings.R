settings <- 
setRefClass("RzSettings",
  fields = c("RzPath", "theme",
             "globalFont", "variableViewFont", "monospaceFont", "monospaceFontFamily",
             "useEmbededDevice", "embededDeviceOn", "popupOff",
             "psFont", "pdfFont", "autosave"),
  methods = list(
    load = function(){
      path <- file.path(rzConfPath(), "Rz.conf")
      RzPath <<- path.package("Rz")
      embededDeviceOn <<- NULL
      # Does a setting file exists?
      if(file.exists(path)) {
        con <- file(path)
        open(con)
        settings <- dget(con)
        close(con)
        # platform dependent settings
        if(grepl("linux", R.Version()$os)){
          theme            <<- ifelse(is.null(settings$theme)   , "Default", settings$theme)
        } else {
          theme            <<- ifelse(is.null(settings$theme)   , "kde42-oxygen", settings$theme)
        }
        globalFont       <<- ifelse(is.null(settings$globalFont)   , "sans 10"     , settings$globalFont)
        variableViewFont <<- ifelse(is.null(settings$variableViewFont)   , "sans 10"     , settings$variableViewFont)

        # settings        
        monospaceFont    <<- ifelse(is.null(settings$monospaceFont), "monospace 10", settings$monospaceFont)
        monospaceFontFamily <<- pangoFontDescriptionFromString(monospaceFont)$getFamily()
        psFont           <<- ifelse(is.null(settings$psFont)       , "sans"        , settings$psFont)
        pdfFont          <<- ifelse(is.null(settings$pdfFont)      , "sans"        , settings$pdfFont)
        useEmbededDevice <<- ifelse(is.null(settings$useEmbededDevice), FALSE, settings$useEmbededDevice)
        popupOff         <<- ifelse(is.null(settings$popupOff),         FALSE, settings$popupOff)
        autosave         <<- ifelse(is.null(settings$autosave),         TRUE , settings$autosave)
      } else {
        # initialize settings
        if(grepl("linux", R.Version()$os)){
          theme      <<- "Default"
        } else {
          theme      <<- "kde42-oxygen"          
        }
        globalFont <<- "sans 10"
        variableViewFont <<- "sans 10"
        
        monospaceFont <<- "monospace 10"
        monospaceFontFamily <<- pangoFontDescriptionFromString(monospaceFont)$getFamily()
        psFont <<- "sans"
        pdfFont <<- "sans"
        useEmbededDevice <<- FALSE
        popupOff         <<- FALSE
        autosave         <<- TRUE
      }
      theme.path <- system.file("themes", theme, "gtk-2.0", "gtkrc", package="Rz")
      gtkRcParse(theme.path) 
    },
    
    runDialog = function(win){
      dialog <- gtkDialogNewWithButtons(gettext("Settings"), win,
                                        c("modal", "destroy-with-parent"), 
                                        "gtk-ok", GtkResponseType["accept"], 
                                        "gtk-cancel", GtkResponseType["reject"],
                                        show=FALSE)
      
      themes.label <- gtkLabelNew(gettext("Theme (requires restart R)"))
      themesCombo  <- gtkComboBoxNewText()
      themesCombo$getCells()[[1]]$setAlignment(0.5, 0.5)
      themes  <- sapply(list.dirs(system.file("themes", package="Rz"), recursive=FALSE), basename)
      
      for(i in themes) themesCombo$appendText(i)
      themesCombo$setActive(which(theme==themes) - 1)
      themes.hbox <- gtkHBoxNew(spacing=5)
      themes.hbox$packStart(themes.label, expand=FALSE)
      themes.hbox$packStart(themesCombo)
            
      checkButtonUseEmbededDevice <- gtkCheckButtonNewWithLabel(gettext("Use embeded graphics divice (requires cairoDevice package)"))
      checkButtonUseEmbededDevice$setActive(useEmbededDevice)
      checkButtonPopupOff <- gtkCheckButtonNewWithLabel(gettext("Don't Popup Summary"))
      checkButtonPopupOff$setActive(popupOff)
      checkButtonAutosave <- gtkCheckButtonNewWithLabel(gettext("Automatically save the session"))
      checkButtonAutosave$setActive(autosave)
      
      general.tab <- gtkVBoxNew()
      general.tab["border-width"] <- 2
      general.tab$packStart(themes.hbox, expand=FALSE)
      general.tab$packStart(checkButtonUseEmbededDevice, expand=FALSE)
      general.tab$packStart(checkButtonPopupOff, expand=FALSE)
      general.tab$packStart(checkButtonAutosave, expand=FALSE)
      
      rzFontSettingWidget1 <- new("RzFontSettingWidget", title = gettext("Global Font"), fontName = globalFont, showSize = TRUE, showStyle=TRUE)
      rzFontSettingWidget4 <- new("RzFontSettingWidget", title = gettext("Variable View Font"), fontName = variableViewFont, showSize = TRUE, showStyle=TRUE)
      rzFontSettingWidget2 <- new("RzFontSettingWidget", title = gettext("Monospace Font"), fontName = monospaceFont, showSize = TRUE, showStyle=TRUE)
      
      pdf.font.label <- gtkLabelNew(gettext("PDF Font"))
      pdfFontCombo <- gtkComboBoxNewText()
      pdfFontCombo$getCells()[[1]]$setAlignment(0.5, 0.5)
      for(i in names(pdfFonts())) pdfFontCombo$appendText(i)
      pdfFontCombo$setActive(which(pdfFont==names(pdfFonts())) - 1)
      
      ps.font.label <- gtkLabelNew(gettext("PostScript Font"))
      psFontCombo <- gtkComboBoxNewText()
      psFontCombo$getCells()[[1]]$setAlignment(0.5, 0.5)
      for(i in names(postscriptFonts())) psFontCombo$appendText(i)
      psFontCombo$setActive(which(psFont==names(postscriptFonts())) - 1)
      
      pdffont.hbox <- gtkHBoxNew(spacing=5)
      pdffont.hbox$packStart(pdf.font.label, expand=FALSE)
      pdffont.hbox$packStart(pdfFontCombo)
      psfont.hbox <- gtkHBoxNew(spacing=5)
      psfont.hbox$packStart(ps.font.label, expand=FALSE)
      psfont.hbox$packStart(psFontCombo)
      
      font.tab <- gtkVBoxNew(spacing=2)
      font.tab["border-width"] <- 2
      font.tab$packStart(rzFontSettingWidget1$getFontBox(), fill=FALSE, expand=FALSE)
      font.tab$packStart(rzFontSettingWidget4$getFontBox(), fill=FALSE, expand=FALSE)
      font.tab$packStart(rzFontSettingWidget2$getFontBox(), fill=FALSE, expand=FALSE)
      font.tab$packStart(pdffont.hbox, fill=FALSE, expand=FALSE)
      font.tab$packStart(psfont.hbox, fill=FALSE, expand=FALSE)
            
      note <- gtkNotebookNew()
      note$appendPage(general.tab, gtkLabelNew(gettext("General")))
      note$appendPage(font.tab, gtkLabelNew(gettext("Font")))
      dialog[["vbox"]]$packStart(note)
      
      onResponse <- function(dialog, response.id){
        if(response.id == GtkResponseType["accept"]) {
          theme            <<- localize(themesCombo$getActiveText())
          globalFont       <<- localize(rzFontSettingWidget1$getFontName())
          variableViewFont <<- localize(rzFontSettingWidget4$getFontName())
          monospaceFont    <<- localize(rzFontSettingWidget2$getFontName())
          monospaceFontFamily <<- pangoFontDescriptionFromString(monospaceFont)$getFamily()
          psFont           <<- localize(psFontCombo$getActiveText())
          pdfFont          <<- localize(pdfFontCombo$getActiveText())
          useEmbededDevice <<- checkButtonUseEmbededDevice$getActive()
          popupOff         <<- checkButtonPopupOff$getActive()
          autosave         <<- checkButtonAutosave$getActive()
          settings <- gtkSettingsGetDefault()
          settings$setStringProperty("gtk-font-name", rzSettings$getGlobalFont(), NULL)
          
          if (checkConfDir()) {
            path <- file.path(rzConfPath(), "Rz.conf")
            
            con <- file(path, open="w")
            dput(list(
              theme            = theme,
              globalFont       = globalFont,
              variableViewFont = variableViewFont,
              monospaceFont    = monospaceFont,
              psFont           = psFont,
              pdfFont          = pdfFont,
              useEmbededDevice = useEmbededDevice,
              popupOff         = popupOff,
              autosave         = autosave
            ),
                 file=con, control=NULL)
            close(con)
            if (autosave) {
              saveSession()
            } else {
              unlink(file.path(rzConfPath(), "session.rzs"))
            }
          }
          dialog$hide()
        } else {
          dialog$hide()
        }
      }
      gSignalConnect(dialog, "response", onResponse)
      
      dialog$run()
    }
  )
)
settings$accessors(c("RzPath", "globalFont", "variableViewFont", "monospaceFont", "monospaceFontFamily",
                     "useEmbededDevice", "embededDeviceOn", "popupOff", "autosave",
                     "psFont", "pdfFont"))

