menu <- 
setRefClass("RzMenu",
  fields = c("uimanager", "menu.bar", "a.file", "a.open",
             "a.import", "a.importsav", "a.importpor", "a.importdta", "a.importt",
             "a.save", "a.importg", "a.quit", "a.menusettings", "a.settings",
             "a.help", "a.tutorial", "a.load.sample"),
  methods = list(
    initialize            = function(...) {
      initFields(...)
      path <- file.path(path.package("Rz"), "ui", "menu.ui")
      
      a.file      <<- gtkActionNew("MenuFile", gettext("_File"))
      a.open      <<- gtkActionNew("Open", gettext("_Open Data..."))
      a.import    <<-  gtkActionNew("MenuImport", gettext("_Import Data"))
      a.importsav <<- gtkActionNew("ImportSav", gettext("From SPSS System File..."))
      a.importpor <<- gtkActionNew("ImportPor", gettext("From SPSS Portable File..."))
      a.importdta <<- gtkActionNew("ImportDta", gettext("From Stata Data File..."))
      a.importt   <<- gtkActionNew("ImportText", gettext("From Text File..."))
      a.importg   <<- gtkActionNew("ImportFromGlobalEnv", gettext("From Grobal Environment..."))
      a.save      <<- gtkActionNew("Save As", gettext("Save _As..."))
      a.quit      <<- gtkActionNew("Close", gettext("_Close"))
      
      a.menusettings <<- gtkActionNew("MenuSettings", gettext("_Preferences"))
      a.settings  <<- gtkActionNew("Settings", gettext("_Preferences..."), gettext("Preferences..."))
      
      a.help      <<- gtkActionNew("MenuHelp", gettext("_Help"))
      a.tutorial  <<- gtkActionNew("Tutorial", gettext("Tutorial on Web"), gettext("Tutorial on Web"))
      a.load.sample <<- gtkActionNew("LoadSample", gettext("Load Sample Dataset"), gettext("Load Sample Dataset"))
      
      a.open$setAccelPath("<Rz-Menu>/File/Open")
      a.save$setAccelPath("<Rz-Menu>/File/Save As")
      a.importg$setAccelPath("<Rz-Menu>/File/Import")
      a.quit$setAccelPath("<Rz-Menu>/File/Close")
      
      action.group  <- gtkActionGroupNew()
      action.group$setTranslationDomain("pkg-RGtk2")
      
      action.group$addAction(a.file)
      action.group$addAction(a.open)
      action.group$addAction(a.import)
      action.group$addAction(a.importsav)
      action.group$addAction(a.importpor)
      action.group$addAction(a.importdta)
      action.group$addAction(a.importt)
      action.group$addAction(a.importg)
      action.group$addAction(a.save)
      action.group$addAction(a.quit)
      
      action.group$addAction(a.menusettings)
      action.group$addAction(a.settings)
      
      action.group$addAction(a.help)
      action.group$addAction(a.tutorial)
      action.group$addAction(a.load.sample)
      
      uimanager     <<- gtkUIManagerNew()
      uimanager$insertActionGroup(action.group, 0)
      uimanager$addUiFromFile(path)
      menu.bar     <<- uimanager$getWidget("/MenuBar")
    }
))
menu$accessors(c("uimanager", "menu.bar", "a.file", "a.open",
                 "a.import", "a.importsav", "a.importpor", "a.importdta", "a.importt",
                 "a.save", "a.importg", "a.quit", "a.menusettings", "a.settings",
                 "a.help", "a.tutorial", "a.load.sample"))
