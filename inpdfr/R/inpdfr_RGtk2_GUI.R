#' RGtk2 GUI function: switch on the "Processing..." message.
#'
#' This function is provided so that you can easily see its content. It is not
#'   intended to be used, prefer \code{loadGUI()} to load the \code{RGtk2} GUI.
#' @import RGtk2
#' @export
switchOnDialogWait <- function(){
  dialogX <- gtkMessageDialog(parent = NULL, flags = 0, type = "info", 
    buttons = "none", "Processing...", show = FALSE)
  # dialogX$setDeletable(FALSE) # not working on Ubuntu LTS 14.04
  dialogX$setModal(TRUE)
  hbox <- gtkHBox()
  spinner <- gtkSpinner()
  spinner$start()
  hbox$packStart(spinner)
  vbox <- dialogX$getContentArea()
  vbox$packStart(hbox)
  dialogX$show()

  return(dialogX)
}

#' RGtk2 GUI function: switch off the "Processing..." message.
#'
#' This function is provided so that you can easily see its content. It is not
#'   intended to be used, prefer \code{loadGUI()} to load the \code{RGtk2} GUI.
#' @param dialogX The dialog window to be closed.
#' @import RGtk2
#' @export
switchOffDialogWait <- function(dialogX){
  dialogX$destroy()
}

#' RGtk2 GUI function: ask confirmation to quit if topright button is used to quit.
#'
#' This function is provided so that you can easily see its content. It is not
#'   intended to be used, prefer \code{loadGUI()} to load the \code{RGtk2} GUI.
#' @param myobject The parent window to be closed.
#' @import RGtk2
#' @export
askQuit <- function(myobject) {
  gSignalConnect(myobject, "delete-event", function(event, ...) {
    dialog <- gtkMessageDialog(parent = myobject, flags = 0, type = "question", 
      buttons = "yes-no", "Are you sure you want to quit?")
    out <- dialog$run()
    dialog$destroy()
    out != GtkResponseType["yes"]
  })
}

#' RGtk2 GUI function: check data validity in entries.
#'
#' This function is provided so that you can easily see its content. It is not
#'   intended to be used, prefer \code{loadGUI()} to load the \code{RGtk2} GUI.
#' @param validatedEntry Entry to be checked.
#' @param myRegEx Regular expression to test the \code{validatedEntry}.
#' @import RGtk2
#' @export
checkEntry <- function(validatedEntry,myRegEx){
  gSignalConnect(validatedEntry, "changed", function(entry) {
    mytext <- entry$getText()
    if (nzchar(gsub(myRegEx, "", mytext))) { #  "^([a-zA-Z])([a-zA-Z0-9]*)"
      entry$setIconFromStock("primary", "gtk-no")
      # entry$setIconTooltipText("primary","Only letters and numbers are allowed")
    } else {
      entry$setIconFromStock("primary", "gtk-yes")
      # entry$setIconTooltipText("primary", NULL)
    }
  })
}

#' RGtk2 GUI function: open a window in order to choose the working directory.
#'
#' This function is provided so that you can easily see its content. It is not
#'   intended to be used, prefer \code{loadGUI()} to load the \code{RGtk2} GUI.
#' @param widget Widget to open.
#' @param window A window to contain the widget.
#' @return The path to the user-defeined working directory.
#' @import RGtk2
#' @export

open_cb <- function(widget, window) {
  dialog <- gtkFileChooserDialog("working directory", window, "select-folder","gtk-cancel", 
    GtkResponseType["cancel"], "gtk-open", GtkResponseType["accept"])
  if (dialog$run() == GtkResponseType["accept"]) {
    currentObjectName <- dialog$getFilename()
    message(dialog$getFilename())
  }
  dialog$destroy()
  return(currentObjectName)
}

#' RGtk2 GUI function: open a window in order to choose an R data file (rda, RDA, RData).
#'
#' This function is provided so that you can easily see its content. It is not
#'   intended to be used, prefer \code{loadGUI()} to load the \code{RGtk2} GUI.
#' @param widget Widget to open.
#' @param window A window to contain the widget.
#' @return The path to the user-defeined file.
#' @import RGtk2
#' @export
open_cbFile <- function(widget, window) {
  dialog <- gtkFileChooserDialog("File", window, "open","gtk-cancel", 
    GtkResponseType["cancel"], "gtk-open", GtkResponseType["accept"])
  fileFilter <- gtkFileFilter()
  fileFilter$setName("R files")
  fileFilter$addPattern("*.rda")
  fileFilter$addPattern("*.RDA")
  fileFilter$addPattern("*.RData")
  dialog$addFilter(fileFilter)
  if (dialog$run() == GtkResponseType["accept"]) {
    currentObjectNameFile <- dialog$getFilename()
    message(dialog$getFilename())
  }
  dialog$destroy()
  return(currentObjectNameFile)
}

#' RGtk2 GUI function: main window.
#'
#' This function is provided so that you can easily see its content. It is not
#'   intended to be used, prefer \code{loadGUI()} to load the \code{RGtk2} GUI.
#' @param main_window Main window where the menu is created.
#' @import RGtk2
#' @export
makeMenuMainWindow <- function(main_window){
  new_cb <- function(widget, window){
    ### not used
  }

  open_cb <- function(widget, window) {
    dialog <- gtkFileChooserDialog(
      "working directory",
      window,
      "select-folder",
      "gtk-cancel",
      GtkResponseType["cancel"],
      "gtk-open",
      GtkResponseType["accept"]
    )
    if (dialog$run() == GtkResponseType["accept"]) {
      currentObjectName <- dialog$getFilename()
      message(dialog$getFilename())
    }
    dialog$destroy()
    return(currentObjectName)
  }

  quit_cb <- function(widget, window) window$destroy() ### close the GUI

  activate_email <- function(about, link, data){
    utils::browseURL(paste("mailto:", link, sep = ""))
  }

  activate_url <- function(about, link, data){
    utils::browseURL(link)
  }

  about_cb <- function(widget, window) { ### credits
    gtkAboutDialogSetEmailHook(activate_email)
    gtkAboutDialogSetUrlHook(activate_url)
    gtkShowAboutDialog(window,
                       program_name = "inpdfr GUI",
                       version = "version beta 0.49",
                       copyright = "Copyright (C) IRD / Fran\u00e7ois Rebaudo",
                       license = "GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007 (https://www.r-project.org/Licenses/GPL-3)",
                       website = "",
                       comments = "",
                       authors = c("Fran\u00e7ois Rebaudo", "<francois.rebaudo@ird.fr>"),
                       documenters = "Fran\u00e7ois Rebaudo"#,
                       #translator_credits = "Fran\u00e7ois Rebaudo"#,
                       #logo = myIMAGE
    )
  }

  actions <- list( ### actions associated with the menu items and tool items
    list("FileMenu", NULL, "File"),
    list("About", NULL, "About"),
    list("About_GUI", "gtk-about", "About", "<control>T", "About...", about_cb),
    list("Quit", "gtk-quit", "Quit", "<control>Q","Quit", quit_cb)
  )
  action_group <- gtkActionGroup("spreadsheetActions")
  action_group$addActions(actions, main_window)
  ui_manager <- gtkUIManager()
  ui_manager$insertActionGroup(action_group, 0)

  merge <- ui_manager$newMergeId()
  ui_manager$addUi(merge.id = merge, path = "/", name = "menubar",action = NULL, type = "menubar", top = FALSE)
  ui_manager$addUi(merge, "/menubar", "file", "FileMenu", "menu", FALSE)
  ui_manager$addUi(merge, "/menubar/file", "quit", "Quit", "menuitem", FALSE)
  ui_manager$addUi(merge, "/menubar", "about", "About", "menu", FALSE)
  ui_manager$addUi(merge, "/menubar/about", "about", "About_GUI", "menuitem", FALSE)
  menubar <- ui_manager$getWidget("/menubar")
  main_window$addAccelGroup(ui_manager$getAccelGroup()) ### keyboard shortcuts in the main window.
  vbox <- gtkVBox(homogeneous = FALSE, spacing = 0)
  vbox$packStart(menubar, expand = FALSE, fill = FALSE, padding = 0)
  return(vbox)
}

#' RGtk2 GUI function: dynamic content of main window.
#'
#' This function is provided so that you can easily see its content. It is not
#'   intended to be used, prefer \code{loadGUI()} to load the \code{RGtk2} GUI.
#' @param main_window Main window where the content is created.
#' @import RGtk2
#' @export
makeMainWindowsContent <- function(main_window){

  mergedD <- NULL # makeMainWindowsContent : <anonymous>: no visible binding for global variable 'mergedD'

  myPadding <- 1

  vbox <- gtkVBox(homogeneous = FALSE, spacing = 0)

  hboxGEN <- gtkHBoxNew(homogeneous = FALSE, spacing = 0)
  vboxSUB <- gtkVBox(homogeneous = FALSE, spacing = 0)
  hboxWD <- gtkHBoxNew(homogeneous = FALSE, spacing = 0)
  labelWD <- gtkLabelNew("Select folder: ") ; gtkMiscSetAlignment(labelWD, xalign = 1, yalign = 0.5)
  buttonWD <- gtkButton(stock.id = "gtk-find")
  entryWD <- gtkEntryNew()
  hboxWD$packStart(labelWD, expand = FALSE, fill = FALSE, padding = myPadding)
  hboxWD$packStart(buttonWD, expand = FALSE, fill = FALSE, padding = myPadding)
  hboxWD$packStart(entryWD, expand = TRUE, fill = TRUE, padding = myPadding)
  hboxSW <- gtkHBoxNew(homogeneous = FALSE, spacing = 0)
  labelStopWords <- gtkLabelNew("Stop words: ")
  gtkMiscSetAlignment(labelStopWords, xalign = 1, yalign = 0.5)
  choicesStopWords<-rGtkDataFrame(c("none","French","English","Spanish","user-defined"))
  comboBoxStopWords <- gtkComboBox(choicesStopWords)
  crtMAP <- gtkCellRendererText()
  comboBoxStopWords$packStart(crtMAP)
  comboBoxStopWords$addAttribute(crtMAP, "text", 0)
  gtkComboBoxSetActive(comboBoxStopWords,0)
  labelTrunc <- gtkLabelNew("Max. num. words: ")
  gtkMiscSetAlignment(labelTrunc, xalign = 1, yalign = 0.5)
  entryTrunc <- gtkEntryNew()
  hboxSW$packStart(labelStopWords, expand = FALSE, fill = FALSE, padding = myPadding)
  hboxSW$packStart(comboBoxStopWords, expand = FALSE, fill = FALSE, padding = myPadding)
  hboxSW$packStart(labelTrunc, expand = FALSE, fill = FALSE, padding = myPadding)
  hboxSW$packStart(entryTrunc, expand = TRUE, fill = TRUE, padding = myPadding)
  vboxSUB$packStart(hboxWD, expand = TRUE, fill = TRUE, padding = myPadding)
  vboxSUB$packStart(hboxSW, expand = TRUE, fill = TRUE, padding = myPadding)
  hboxGEN$packStart(vboxSUB, expand = TRUE, fill = TRUE, padding = myPadding)
  buttonUpload <- gtkButton(stock.id="gtk-ok")
  hboxGEN$packStart(buttonUpload, expand = FALSE, fill = FALSE, padding = myPadding)

  vbox$packStart(hboxGEN, expand = FALSE, fill = FALSE, padding = myPadding)
  #############################
  ### TAB DESCRIBE
  blocDESC <- {
    vboxDesOptions <- gtkVBoxNew(homogeneous = FALSE, spacing = 0)
    hboxDesMostFreq <- gtkHBoxNew(homogeneous = TRUE, spacing = 0)
    labelDesMostFreq <- gtkLabelNew("X most frequent words: ")
    gtkMiscSetAlignment(labelDesMostFreq, xalign = 1, yalign = 0.5)
    entryDesMostFreq <- gtkEntryNew()
    buttonDesMostFreq <- gtkButton(stock.id = "gtk-ok")
    hboxDesMostFreq$packStart(labelDesMostFreq, expand = FALSE, fill = TRUE, padding = 1)
    hboxDesMostFreq$packStart(entryDesMostFreq, expand = TRUE, fill = TRUE, padding = 1)
    hboxDesMostFreq$packStart(buttonDesMostFreq, expand = TRUE, fill = TRUE, padding = 1)
    hboxDesXFreq <- gtkHBoxNew(homogeneous = TRUE, spacing = 0)
    labelDesXFreq <- gtkLabelNew("Words with occurrences > X: ")
    gtkMiscSetAlignment(labelDesXFreq, xalign = 1, yalign = 0.5)
    entryDesXFreq <- gtkEntryNew()
    buttonDesXFreq <- gtkButton(stock.id = "gtk-ok")
    hboxDesXFreq$packStart(labelDesXFreq, expand = FALSE, fill = TRUE, padding = 1)
    hboxDesXFreq$packStart(entryDesXFreq, expand = TRUE, fill = TRUE, padding = 1)
    hboxDesXFreq$packStart(buttonDesXFreq, expand = TRUE, fill = TRUE, padding = 1)
    hboxHisto <- gtkHBoxNew(homogeneous = TRUE, spacing = 0)
    labelHisto <- gtkLabelNew("Histogram (number of words): ")
    gtkMiscSetAlignment(labelHisto, xalign = 1, yalign = 0.5)
    buttonOK_Histo <- gtkButton(stock.id = "gtk-ok")
    hboxHisto$packStart(labelHisto, expand = FALSE, fill = TRUE, padding = 1)
    hboxHisto$packStart(buttonOK_Histo, expand = TRUE, fill = TRUE, padding = 1)
    hboxBarplot <- gtkHBoxNew(homogeneous = TRUE, spacing = 0)
    labelBarplot <- gtkLabelNew("Barplot (number of words): ")
    gtkMiscSetAlignment(labelBarplot, xalign = 1, yalign = 0.5)
    buttonOK_Barplot <- gtkButton(stock.id = "gtk-ok")
    hboxBarplot$packStart(labelBarplot, expand = FALSE, fill = TRUE, padding = 1)
    hboxBarplot$packStart(buttonOK_Barplot, expand = TRUE, fill = TRUE, padding = 1)
    hboxOccurr <- gtkHBoxNew(homogeneous = TRUE, spacing = 0)
    labelOccurr <- gtkLabelNew("Words occurrence: \n(may take some time)")
    gtkMiscSetAlignment(labelOccurr, xalign = 1, yalign = 0.5)
    buttonOK_Occurr <- gtkButton(stock.id = "gtk-ok")
    hboxOccurr$packStart(labelOccurr, expand = FALSE, fill = TRUE, padding = 1)
    hboxOccurr$packStart(buttonOK_Occurr, expand = TRUE, fill = TRUE, padding = 1)
    vboxDesOptions$add(hboxHisto)
    vboxDesOptions$add(hboxBarplot)
    vboxDesOptions$add(hboxOccurr)
    vboxDesOptions$add(hboxDesMostFreq)
    vboxDesOptions$add(hboxDesXFreq)
  }
  #############################
  ### TAB WORDCLOUD
  blocWORD <- {
    # plot
    hboxPlot <- gtkHBoxNew(homogeneous = FALSE, spacing = 0)
    checkboxFile <- gtkCheckButton(label = "Plot one wc per document")
    checkboxFile$active <- TRUE
    checkboxSet <- gtkCheckButton(label = "Plot global wc")
    checkboxSet$active <- TRUE
    hboxPlot$packStart(checkboxFile, expand = FALSE, fill = FALSE, padding = 1)
    hboxPlot$packStart(checkboxSet, expand = FALSE, fill = FALSE, padding = 1)
    # Min. word frequency
    hboxMinFreq <- gtkHBoxNew(homogeneous = FALSE, spacing = 0)
    labelMinFreq <- gtkLabelNew("Min. word occurrence:")
    gtkMiscSetAlignment(labelMinFreq, xalign = 1, yalign = 0.5)
    choicesMinFreq <- rGtkDataFrame(c(1:5, 10, 15, 20, 50, 100, 500))
    comboBoxMinFreq <- gtkComboBox(choicesMinFreq)
    crtMAP <- gtkCellRendererText()
    comboBoxMinFreq$packStart(crtMAP)
    comboBoxMinFreq$addAttribute(crtMAP, "text", 0)
    gtkComboBoxSetActive(comboBoxMinFreq, 0)
    hboxMinFreq$packStart(labelMinFreq, expand = FALSE, fill = FALSE, padding = 1)
    hboxMinFreq$packStart(comboBoxMinFreq, expand = TRUE, fill = TRUE, padding = 1)
    # Max. number of words
    hboxMaxWords <- gtkHBoxNew(homogeneous = FALSE, spacing = 0)
    labelMaxWords <- gtkLabelNew("Max. number of words:")
    gtkMiscSetAlignment(labelMaxWords, xalign = 1, yalign = 0.5)
    choicesMaxWords <- rGtkDataFrame(c("Inf", seq(from = 10, to = 100, by = 10),
      seq(from = 200, to = 500, by = 100), 1000))
    comboBoxMaxWords <- gtkComboBox(choicesMaxWords)
    crtMAP <- gtkCellRendererText()
    comboBoxMaxWords$packStart(crtMAP)
    comboBoxMaxWords$addAttribute(crtMAP, "text", 0)
    gtkComboBoxSetActive(comboBoxMaxWords, 0)
    hboxMaxWords$packStart(labelMaxWords, expand = FALSE, fill = FALSE, padding = 1)
    hboxMaxWords$packStart(comboBoxMaxWords, expand = TRUE, fill = TRUE, padding = 1)
    # Random order
    hboxRandOrder <- gtkHBoxNew(homogeneous = FALSE, spacing = 0) #TRUE/FALSE
    labelRandOrder <- gtkLabelNew("Random order:")
    gtkMiscSetAlignment(labelRandOrder, xalign = 1, yalign = 0.5)
    choicesRandOrder <- rGtkDataFrame(c("FALSE", "TRUE"))
    comboBoxRandOrder <- gtkComboBox(choicesRandOrder)
    crtMAP <- gtkCellRendererText()
    comboBoxRandOrder$packStart(crtMAP)
    comboBoxRandOrder$addAttribute(crtMAP, "text", 0)
    gtkComboBoxSetActive(comboBoxRandOrder, 0)
    hboxRandOrder$packStart(labelRandOrder, expand = FALSE, fill = FALSE, padding = 1)
    hboxRandOrder$packStart(comboBoxRandOrder, expand = TRUE, fill = TRUE, padding = 1)
    # colors
    hboxCol <- gtkHBoxNew(homogeneous = FALSE, spacing = 0) #palette
    labelCol <- gtkLabelNew("Colors:")
    gtkMiscSetAlignment(labelCol, xalign = 1, yalign = 0.5)
    choicesCol <- rGtkDataFrame(seq(from = 3, to = 8, by = 1))
    choicesPal <- rGtkDataFrame(c("Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", "Oranges",
      "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr",
      "YlOrRd", "BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", "RdYlBu", "RdYlGn", "Spectral",
      "Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3"))
    comboBoxCol <- gtkComboBox(choicesCol)
    comboBoxPal <- gtkComboBox(choicesPal)
    crtMAP <- gtkCellRendererText()
    comboBoxCol$packStart(crtMAP)
    comboBoxPal$packStart(crtMAP)
    comboBoxCol$addAttribute(crtMAP, "text", 0)
    comboBoxPal$addAttribute(crtMAP, "text", 0)
    gtkComboBoxSetActive(comboBoxCol, 0)
    gtkComboBoxSetActive(comboBoxPal, 0)
    hboxCol$packStart(labelCol, expand = FALSE, fill = FALSE, padding = 1)
    hboxCol$packStart(comboBoxCol, expand = TRUE, fill = TRUE, padding = 1)
    hboxCol$packStart(comboBoxPal, expand = TRUE, fill = TRUE, padding = 1)
    # output file format
    hboxOutFile <- gtkHBoxNew(homogeneous = FALSE, spacing = 0) #png;pdf;jpg;...
    labelOutFile <- gtkLabelNew("Output file format:")
    gtkMiscSetAlignment(labelOutFile, xalign = 1, yalign = 0.5)
    choicesOutFile <- rGtkDataFrame(c("png", "eps", "pdf", "svg", "tiff", "jpeg", "bmp"))
    comboBoxOutFile <- gtkComboBox(choicesOutFile)
    crtMAP <- gtkCellRendererText()
    comboBoxOutFile$packStart(crtMAP)
    comboBoxOutFile$addAttribute(crtMAP, "text", 0)
    gtkComboBoxSetActive(comboBoxOutFile, 0)
    hboxOutFile$packStart(labelOutFile, expand = FALSE, fill = FALSE, padding = 1)
    hboxOutFile$packStart(comboBoxOutFile, expand = TRUE, fill = TRUE, padding = 1)
    buttonOK_WC <- gtkButton(stock.id = "gtk-ok")

    vboxWcOptions <- gtkVBox(homogeneous = FALSE, spacing = 0)
    vboxWcOptions$add(hboxPlot)
    vboxWcOptions$add(hboxMinFreq)
    vboxWcOptions$add(hboxMaxWords)
    vboxWcOptions$add(hboxRandOrder)
    vboxWcOptions$add(hboxCol)
    vboxWcOptions$add(hboxOutFile)
    vboxWcOptions$add(buttonOK_WC)
  }
  #############################
  ### TAB H-CLUST
  blocHCLU <- {
    # normalize num words
    hboxHclustcheckBox <- gtkHBoxNew(homogeneous = FALSE, spacing = 0)
    checkboxStd1 <- gtkCheckButton(label = "Normalize number of words")
    checkboxStd1$active <- TRUE
    hboxHclustcheckBox$packStart(checkboxStd1, expand = TRUE, fill = TRUE, padding = 1)
    # Method
    hboxHclustMethod <- gtkHBoxNew(homogeneous = TRUE, spacing = 0)
    labelHclustMethod <- gtkLabelNew("Method:") ; gtkMiscSetAlignment(labelHclustMethod, xalign = 1, yalign = 0.5)
    choicesHclustMethod <- rGtkDataFrame(c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"))
    comboBoxHclustMethod <- gtkComboBox(choicesHclustMethod)
    crtMAP <- gtkCellRendererText()
    comboBoxHclustMethod$packStart(crtMAP)
    comboBoxHclustMethod$addAttribute(crtMAP, "text", 0)
    gtkComboBoxSetActive(comboBoxHclustMethod, 0)
    hboxHclustMethod$packStart(labelHclustMethod, expand = FALSE, fill = TRUE, padding = 1)
    hboxHclustMethod$packStart(comboBoxHclustMethod, expand = TRUE, fill = TRUE, padding = 1)
    # Grouping
    hboxHclustGrouping <- gtkHBoxNew(homogeneous = TRUE, spacing = 0)
    labelHclustGrouping <- gtkLabelNew("Grouping:") ; gtkMiscSetAlignment(labelHclustGrouping, xalign = 1, yalign = 0.5)
    choicesHclustGrouping <- rGtkDataFrame(c("TRUE", "FALSE"))
    comboBoxHclustGrouping <- gtkComboBox(choicesHclustGrouping)
    crtMAP <- gtkCellRendererText()
    comboBoxHclustGrouping$packStart(crtMAP)
    comboBoxHclustGrouping$addAttribute(crtMAP, "text", 0)
    gtkComboBoxSetActive(comboBoxHclustGrouping, 0)
    hboxHclustGrouping$packStart(labelHclustGrouping, expand = FALSE, fill = TRUE, padding = 1)
    hboxHclustGrouping$packStart(comboBoxHclustGrouping, expand = TRUE, fill = TRUE, padding = 1)
    # NB groups
    hboxHclustNbGroup <- gtkHBoxNew(homogeneous = TRUE, spacing = 0)
    labelHclustNbGroup <- gtkLabelNew("Nb groups:") ; gtkMiscSetAlignment(labelHclustNbGroup, xalign = 1, yalign = 0.5)
    choicesHclustNbGroup <- rGtkDataFrame(1:10)
    comboBoxHclustNbGroup <- gtkComboBox(choicesHclustNbGroup)
    crtMAP <- gtkCellRendererText()
    comboBoxHclustNbGroup$packStart(crtMAP)
    comboBoxHclustNbGroup$addAttribute(crtMAP, "text", 0)
    gtkComboBoxSetActive(comboBoxHclustNbGroup, 0)
    hboxHclustNbGroup$packStart(labelHclustNbGroup, expand = FALSE, fill = TRUE, padding = 1)
    hboxHclustNbGroup$packStart(comboBoxHclustNbGroup, expand = TRUE, fill = TRUE, padding = 1)
    # button validate
    buttonOK_Hclust <- gtkButton(stock.id = "gtk-ok")

    vboxHclustOptions <- gtkVBox(homogeneous = FALSE, spacing = 0)
    vboxHclustOptions$add(hboxHclustcheckBox)
    vboxHclustOptions$add(hboxHclustMethod)
    vboxHclustOptions$add(hboxHclustGrouping)
    vboxHclustOptions$add(hboxHclustNbGroup)
    vboxHclustOptions$add(buttonOK_Hclust)
  }
  #############################
  ### TAB K-MEAN CLUST
  blocKCLU <- {
    # normalize num words
    hboxHclustcheckBox2 <- gtkHBoxNew(homogeneous = FALSE, spacing = 0)
    checkboxStd2 <- gtkCheckButton(label = "Normalize number of words")
    checkboxStd2$active <- TRUE
    hboxHclustcheckBox2$packStart(checkboxStd2, expand = TRUE, fill = TRUE, padding = 1)
    # centers
    hboxKmeanCenter <- gtkHBoxNew(homogeneous = TRUE, spacing = 0)
    labelKmeanCenter <- gtkLabelNew("Nb clusters:") ; gtkMiscSetAlignment(labelKmeanCenter, xalign = 1, yalign = 0.5)
    choicesKmeanCenter <- rGtkDataFrame(1:10)
    comboBoxKmeanCenter <- gtkComboBox(choicesKmeanCenter)
    crtMAP <- gtkCellRendererText()
    comboBoxKmeanCenter$packStart(crtMAP)
    comboBoxKmeanCenter$addAttribute(crtMAP, "text", 0)
    gtkComboBoxSetActive(comboBoxKmeanCenter, 0)
    hboxKmeanCenter$packStart(labelKmeanCenter, expand = FALSE, fill = TRUE, padding = 1)
    hboxKmeanCenter$packStart(comboBoxKmeanCenter, expand = TRUE, fill = TRUE, padding = 1)
    # iter.max
    hboxKmeanIterMax <- gtkHBoxNew(homogeneous = TRUE, spacing = 0)
    labelKmeanIterMax <- gtkLabelNew("Nb iterations:") ; gtkMiscSetAlignment(labelKmeanIterMax, xalign = 1, yalign = 0.5)
    choicesKmeanIterMax <- rGtkDataFrame(c(1, 5, 10, 20, 50, 100))
    comboBoxKmeanIterMax <- gtkComboBox(choicesKmeanIterMax)
    crtMAP <- gtkCellRendererText()
    comboBoxKmeanIterMax$packStart(crtMAP)
    comboBoxKmeanIterMax$addAttribute(crtMAP, "text", 0)
    gtkComboBoxSetActive(comboBoxKmeanIterMax, 0)
    hboxKmeanIterMax$packStart(labelKmeanIterMax, expand = FALSE, fill = TRUE, padding = 1)
    hboxKmeanIterMax$packStart(comboBoxKmeanIterMax, expand = TRUE, fill = TRUE, padding = 1)
    # algorithm
    hboxKmeanAlgo <- gtkHBoxNew(homogeneous = TRUE, spacing = 0)
    labelKmeanAlgo <- gtkLabelNew("Algorithm:") ; gtkMiscSetAlignment(labelKmeanAlgo, xalign = 1, yalign = 0.5)
    choicesKmeanAlgo <- rGtkDataFrame(c("Hartigan-Wong", "Lloyd", "Forgy", "MacQueen"))
    comboBoxKmeanAlgo <- gtkComboBox(choicesKmeanAlgo)
    crtMAP <- gtkCellRendererText()
    comboBoxKmeanAlgo$packStart(crtMAP)
    comboBoxKmeanAlgo$addAttribute(crtMAP, "text", 0)
    gtkComboBoxSetActive(comboBoxKmeanAlgo, 0)
    hboxKmeanAlgo$packStart(labelKmeanAlgo, expand = FALSE, fill = TRUE, padding = 1)
    hboxKmeanAlgo$packStart(comboBoxKmeanAlgo, expand = TRUE, fill = TRUE, padding = 1)
    # button validate
    buttonOK_KmeanClust <- gtkButton(stock.id = "gtk-ok")

    vboxKmeanClustOptions <- gtkVBox(homogeneous = FALSE, spacing = 0)
    vboxKmeanClustOptions$add(hboxHclustcheckBox2)
    vboxKmeanClustOptions$add(hboxKmeanCenter)
    vboxKmeanClustOptions$add(hboxKmeanIterMax)
    vboxKmeanClustOptions$add(hboxKmeanAlgo)
    vboxKmeanClustOptions$add(buttonOK_KmeanClust)
  }
  #############################
  ### TAB Corresp Analysis and Correl Analysis
  blocCORR <- {
    # normalize num words
    frameCorrespAna <- gtkFrameNew("Correspondance Analysis")
    vboxFrameCA <- gtkVBoxNew(homogeneous = FALSE, spacing = 0)
    hboxCAcheckBox3 <- gtkHBoxNew(homogeneous = FALSE, spacing = 0)
    checkboxStd3 <- gtkCheckButton(label = "Normalize number of words")
    checkboxStd3$active <- TRUE
    hboxCAcheckBox3$packStart(checkboxStd3, expand = TRUE, fill = TRUE, padding = 1)
    # button validate
    buttonOK_CA <- gtkButton(stock.id = "gtk-ok")
    vboxFrameCA$add(hboxCAcheckBox3)
    vboxFrameCA$packStart(buttonOK_CA, expand = TRUE, fill = TRUE, padding = 1)
    frameCorrespAna$add(vboxFrameCA)

    frameCorrelAna <- gtkFrameNew("Correlation Analysis")
    vboxFrameCorrelA <- gtkVBoxNew(homogeneous = FALSE, spacing = 0)
    hboxFrameCorrelA <- gtkHBoxNew(homogeneous = FALSE, spacing = 0)
    labelCorrelA <- gtkLabelNew("Nb words:") ; gtkMiscSetAlignment(labelCorrelA, xalign = 1, yalign = 0.5)
    choicesCorrelA <- rGtkDataFrame(c(1:10, 20, 30, 40, 50, 100))
    comboBoxCorrelA <- gtkComboBox(choicesCorrelA)
    crtMAP <- gtkCellRendererText()
    comboBoxCorrelA$packStart(crtMAP)
    comboBoxCorrelA$addAttribute(crtMAP, "text", 0)
    gtkComboBoxSetActive(comboBoxCorrelA, 0)
    hboxFrameCorrelA$packStart(labelCorrelA, expand = FALSE, fill = FALSE, padding = 1)
    hboxFrameCorrelA$packStart(comboBoxCorrelA, expand = TRUE, fill = TRUE, padding = 1)
    # button validate
    buttonOK_CorrelA <- gtkButton(stock.id = "gtk-ok")
    vboxFrameCorrelA$packStart(hboxFrameCorrelA, expand = FALSE, fill = FALSE, padding = 1)
    vboxFrameCorrelA$packStart(buttonOK_CorrelA, expand = TRUE, fill = TRUE, padding = 1)
    frameCorrelAna$add(vboxFrameCorrelA)

    vboxCAOptions <- gtkVBox(homogeneous = FALSE, spacing = 0)
    vboxCAOptions$packStart(frameCorrespAna, expand = TRUE, fill = TRUE, padding = 1)
    vboxCAOptions$packStart(frameCorrelAna, expand = TRUE, fill = TRUE, padding = 1)
  }
  #############################
  ### TAB MetacomMetacom
  blocMETA <- {
    hboxMetacomGen <- gtkHBoxNew(homogeneous = FALSE, spacing = 0)
    labelMetacom <- gtkLabelNew("Limit the incidence matrix to the x most used words: ")
    gtkMiscSetAlignment(labelMetacom, xalign = 1, yalign = 0.5)
    entryMetacom <- gtkEntryNew()
    entryMetacom$setText("Inf")
    hboxMetacomGen$packStart(labelMetacom, expand = FALSE, fill = FALSE, padding =1 )
    hboxMetacomGen$packStart(entryMetacom, expand = TRUE, fill = TRUE, padding =1 )

    hboxMetacomDiv <- gtkHBoxNew(homogeneous = FALSE, spacing = 0)
    labelMetacomDiv <- gtkLabelNew("Functions to calculate alpha, beta and gamma diversity of communities: ")
    gtkMiscSetAlignment(labelMetacomDiv, xalign = 1, yalign = 0.5)
    buttonOK_MetacomDiv <- gtkButton(stock.id = "gtk-ok")
    hboxMetacomDiv$packStart(labelMetacomDiv, expand = FALSE, fill = FALSE, padding = 1)
    hboxMetacomDiv$packStart(buttonOK_MetacomDiv, expand = TRUE, fill = TRUE, padding = 1)

    hboxMetacomStr <- gtkHBoxNew(homogeneous = FALSE, spacing = 0)
    labelMetacomStr <- gtkLabelNew("Functions for the analysis of the elements of metacommunity structure: \n(this may take some hours on large datasets...)")
    gtkMiscSetAlignment(labelMetacomStr, xalign = 1, yalign = 0.5)
    buttonOK_MetacomStr <- gtkButton(stock.id = "gtk-ok")
    hboxMetacomStr$packStart(labelMetacomStr, expand = FALSE, fill = FALSE, padding = 1)
    hboxMetacomStr$packStart(buttonOK_MetacomStr, expand = TRUE, fill = TRUE, padding = 1)

    vboxMetacomOptions <- gtkVBox(homogeneous = FALSE, spacing = 0)
    vboxMetacomOptions$add(hboxMetacomDiv)
    vboxMetacomOptions$add(hboxMetacomGen)
    vboxMetacomOptions$add(hboxMetacomStr)
  }
  #############################
  ### TAB Save
  blocSAVE <- {
    ## save
    hboxSaveIM <- gtkHBoxNew(homogeneous = FALSE, spacing = 0)
    labelSaveIM <- gtkLabelNew("Save incidence-matrix as: ")
    gtkMiscSetAlignment(labelSaveIM, xalign = 1, yalign = 0.5)
    entrySaveIM <- gtkEntryNew()
    buttonOK_SaveIM <- gtkButton(stock.id = "gtk-ok")
    hboxSaveIM$packStart(labelSaveIM, expand = FALSE, fill = FALSE, padding = 1)
    hboxSaveIM$packStart(entrySaveIM, expand = TRUE, fill = TRUE, padding = 1)
    hboxSaveIM$packStart(buttonOK_SaveIM, expand = FALSE, fill = FALSE, padding = 1)
    ## compare label
    hboxCompLabel <- gtkHBoxNew(homogeneous = FALSE, spacing = 0)
    labelCompLabel <- gtkLabelNew("Compare previously saved incidence-matrices: ")
    gtkMiscSetAlignment(labelCompLabel, xalign = 1, yalign = 0.5)
    hboxCompLabel$packStart(labelCompLabel, expand = FALSE, fill = FALSE, padding = 1)
    ## compare merge?
    hboxCompLabelMerge <- gtkHBoxNew(homogeneous = FALSE, spacing = 0)
    checkboxCompLabelMerge <- gtkCheckButton(label = "Merge documents per group")
    checkboxCompLabelMerge$active <- TRUE
    hboxCompLabelMerge$packStart(checkboxCompLabelMerge, expand = FALSE, fill = FALSE, padding = 1)
    ## compare function number of groups
    vboxCompGroup <- gtkVBox(homogeneous = FALSE, spacing = 0)
    for(k in 1:5){ # loop
      assign(x = paste0("hboxCompGroup", k), value = gtkHBoxNew(homogeneous = FALSE, spacing = 0))
      assign(x = paste0("labelCompGroup", k), value = gtkLabelNew(paste0("Group ", k, ": ")))
      gtkMiscSetAlignment(get(paste0("labelCompGroup", k)), xalign = 1, yalign = 0.5)
      assign(x = paste0("entryCompGroup", k), value = gtkEntryNew())
      assign(x = paste0("buttonFind_CompGroup", k), value = gtkButton(stock.id = "gtk-find"))
      get(paste0("hboxCompGroup", k))$packStart(get(paste0("labelCompGroup", k)), expand = FALSE, fill = FALSE, padding = 1)
      get(paste0("hboxCompGroup", k))$packStart(get(paste0("entryCompGroup", k)), expand = TRUE, fill = TRUE, padding = 1)
      get(paste0("hboxCompGroup", k))$packStart(get(paste0("buttonFind_CompGroup", k)), expand = FALSE, fill = FALSE, padding = 1)
      vboxCompGroup$add(get(paste0("hboxCompGroup", k)))
    }# end loop
    ## compare validation
    hboxCompValid <- gtkHBoxNew(homogeneous = FALSE, spacing = 0)
    buttonOK_CompValid <- gtkButton(stock.id = "gtk-ok")
    hboxCompValid$packStart(buttonOK_CompValid, expand = TRUE, fill = TRUE, padding = 1)
    vboxCompGroup$add(hboxCompValid)

    vboxSaveOptions <- gtkVBox(homogeneous = FALSE, spacing = 0)
    vboxSaveOptions$add(hboxSaveIM)
    vboxSaveOptions$add(hboxCompLabel)
    vboxSaveOptions$add(hboxCompLabelMerge)
    vboxSaveOptions$add(vboxCompGroup)
  }
  #############################

  notebook <- gtkNotebook()
  notebook$setTabPos("top")
  notebook$setHomogeneousTabs(TRUE)
  vbox$packStart(notebook, expand = FALSE, fill = TRUE, padding = 0)

  gSignalConnect(buttonUpload, "clicked", function(widget) {
    mergedD <- NULL

    dialogX <- switchOnDialogWait()

    if (length(entryWD$getText())>0){
      mywd <- entryWD$getText()
      listFilesExt <- getListFiles(mywd)
      quitSpaceFromChars(vectxt = listFilesExt$pdf)
      listFilesExt <- getListFiles(mywd)
      wordFreqPDF <- getPDF(myPDFs = listFilesExt$pdf)
      wordFreqTXT <- getTXT(myTXTs = listFilesExt$txt)
      wordFreq <- append(wordFreqPDF, wordFreqTXT)
      print("Message from GUI: texts extracted from files.")

      ### [1] exclude words from exclusion lists
      language <- choicesStopWords[gtkComboBoxGetActive(comboBoxStopWords)+1]
      if(is.factor(language)){language <- as.character(language)}
      wordFreq <- excludeStopWords(wordF = wordFreq, lang = language)
      print("Message from GUI: stop words excluded according to the specified list.")

      ### [2] select a number of words for analysis (Inf if !numeric)
      maxWords <- as.numeric(entryTrunc$getText())
      wordFreq <- truncNumWords(maxWords = maxWords, wordF = wordFreq)
      print("Message from GUI: words truncated to the specified limit.")

      mergedD <- mergeWordFreq(wordF = wordFreq)
      print("Message from GUI: documents merged into a data.frame (mergedD).")
    }
    mergedD<<-mergedD

    subDir <- "RESULTS"
    dir.create(file.path(getwd(), subDir), showWarnings = FALSE) # create subdirectory for output files
    print("Message from GUI: sub-directory RESULTS created.")

    switchOffDialogWait(dialogX) # close "Processing..." window

    for(i in 1:gtkNotebookGetNPages(notebook)){ # suppr tabs if already loaded
      gtkNotebookRemovePage(notebook, -1)
    }
    notebook$appendPage(vboxWcOptions, gtkLabel("WordCloud"))
    if(ncol(mergedD)>2){ # create tabs according to input
      notebook$appendPage(vboxDesOptions, gtkLabel("Describe"))
      notebook$appendPage(vboxHclustOptions, gtkLabel("Hcluster"))
      notebook$appendPage(vboxKmeanClustOptions, gtkLabel("KmeanClust"))
      notebook$appendPage(vboxCAOptions, gtkLabel("CorAnalysis"))
      notebook$appendPage(vboxMetacomOptions, gtkLabel("MetaCom"))
      notebook$appendPage(vboxSaveOptions, gtkLabel("Compare"))
    }
  })

  gSignalConnect(buttonWD, "clicked", function(widget) { # choose working dir.
    currentObjectName <- open_cb(window = main_window)
    entryWD$setText(currentObjectName)
  })

  gSignalConnect(buttonOK_WC, "clicked", function(widget, wordFreq = mergedD) {
    doAllWc <- checkboxFile$getActive()
    doOneWc <- checkboxSet$getActive()
    myMinFreq <- choicesMinFreq[gtkComboBoxGetActive(comboBoxMinFreq)+1]
    myMaxWords <- choicesMaxWords[gtkComboBoxGetActive(comboBoxMaxWords)+1]
    myRandOrder <- choicesRandOrder[gtkComboBoxGetActive(comboBoxRandOrder)+1]
    myCol <- choicesCol[gtkComboBoxGetActive(comboBoxCol)+1]
    myPal <- choicesPal[gtkComboBoxGetActive(comboBoxPal)+1]
    myOutFile <- toString(choicesOutFile[gtkComboBoxGetActive(comboBoxOutFile)+1])
    if(is.factor(myMaxWords)){myMaxWords <- as.character(myMaxWords)}
    if(is.factor(myRandOrder)){myRandOrder <- as.character(myRandOrder)}
    if(is.factor(myPal)){myPal <- as.character(myPal)}
    if(is.factor(myOutFile)){myOutFile <- as.character(myOutFile)}

    dialogX <- switchOnDialogWait() #!#

    makeWordcloud(wordF = wordFreq, wcFormat = "png", wcminFreq = myMinFreq, 
      wcmaxWords = myMaxWords, wcRandOrder = myRandOrder, 
      wcCol = RColorBrewer::brewer.pal(myCol, myPal),
      getPlot = c(doAllWc, doOneWc), formatType = myOutFile)

    switchOffDialogWait(dialogX) # close "Processing..." window

    print("Message from GUI: word cloud done.")
  })

  gSignalConnect(buttonOK_Hclust, "clicked", function(widget, mydbWords = mergedD) {
    myMethod <- choicesHclustMethod[gtkComboBoxGetActive(comboBoxHclustMethod)+1]
    if (is.factor(myMethod)){myMethod <- as.character(myMethod)}
    gp <- choicesHclustGrouping[gtkComboBoxGetActive(comboBoxHclustGrouping)+1]
    if (is.factor(gp)){
      gp <- as.character(gp)
      gp <- as.logical(gp)
    }
    nbGp <- choicesHclustNbGroup[gtkComboBoxGetActive(comboBoxHclustNbGroup)+1]
    if (is.factor(nbGp)){nbGp <- as.character(nbGp)}
    standard <- checkboxStd1$getActive()
    # print(standard)
    if (standard==TRUE){
      for (i in 2:length(mydbWords[1,])){
        mydbWords[,i] <- round(mydbWords[,i]/sum(mydbWords[,i])*1000, digits = 0)
      }
    }
    doCluster(wordF = mydbWords, myMethod = myMethod, gp = gp, nbGp = nbGp)
    print("Message from GUI: Hierarchical cluster analysis done.")
  })

  gSignalConnect(buttonOK_KmeanClust, "clicked", function(widget, mydbWords = mergedD) {
    nbClust <- choicesKmeanCenter[gtkComboBoxGetActive(comboBoxKmeanCenter)+1]
    if (is.factor(nbClust)){nbClust <- as.character(nbClust)}
    nbIter <- choicesKmeanIterMax[gtkComboBoxGetActive(comboBoxKmeanIterMax)+1]
    if (is.factor(nbIter)){nbIter <- as.character(nbIter)}
    algo <- choicesKmeanAlgo[gtkComboBoxGetActive(comboBoxKmeanAlgo)+1]
    if (is.factor(algo)){algo <- as.character(algo)}
    standard <- checkboxStd2$getActive()
    # print(standard)
    if (standard==TRUE){
      for (i in 2:length(mydbWords[1,])){
        mydbWords[,i] <- round(mydbWords[,i]/sum(mydbWords[,i])*1000, digits = 0)
      }
    }
    doKmeansClust(wordF = mydbWords, nbClust = nbClust, nbIter = nbIter, algo = algo)
    print("Message from GUI: k-means clustering analysis done.")

  })

  gSignalConnect(buttonOK_CA, "clicked", function(widget, mydbWords = mergedD) {
    standard <- checkboxStd3$getActive()
    # print(standard)
    if (standard==TRUE){
      for (i in 2:length(mydbWords[1,])){
        mydbWords[,i] <- round(mydbWords[,i]/sum(mydbWords[,i])*1000, digits = 0)
      }
    }
    doCA(wordF = mydbWords)
    print("Message from GUI: Simple correspondence analysis done.")

  })

  gSignalConnect(buttonOK_CorrelA, "clicked", function(widget, mydbWords = mergedD) {
    nbWordsCA <- choicesCorrelA[gtkComboBoxGetActive(comboBoxCorrelA)+1]
    if (is.factor(nbWordsCA)){nbWordsCA <- as.character(nbWordsCA)}
    try(getMostFreqWordCor(wordF = mydbWords, numWords = nbWordsCA))
    print("Message from GUI: correlation analysis done.")

  })

  gSignalConnect(buttonOK_MetacomDiv, "clicked", function(widget, mydbWords = mergedD) {

    dialogX <- switchOnDialogWait()

    doMetacomEntropart(wordF = mydbWords)

    switchOffDialogWait(dialogX) # close "Processing..." window

    cat("\n")
    print("Message from GUI: metacommunity analysis with 'entropart' package done.")

  })

  gSignalConnect(buttonOK_MetacomStr, "clicked", function(widget, mydbWords = mergedD) {

    limit <- entryMetacom$getText()
    tryCatch(limit <- as.numeric(limit), warning = function(e){print(paste0(limit, " is not a valid entry"))})

    dialogX <- switchOnDialogWait() #!#

    try(doMetacomMetacom(wordF = mydbWords, limit = limit), silent = TRUE)

    switchOffDialogWait(dialogX) # close "Processing..." window

    print("Message from GUI: metacommunity analysis with 'metacom' package done.")

  })

  gSignalConnect(buttonOK_Histo, "clicked", function(widget, mydbWords = mergedD) {
    getSummaryStatsHISTO(wordF = mydbWords)
    print("Message from GUI: histogram done.")
  })
  gSignalConnect(buttonOK_Barplot, "clicked", function(widget, mydbWords = mergedD) {
    getSummaryStatsBARPLOT(wordF = mydbWords)
    print("Message from GUI: bar plot done.")
  })
  gSignalConnect(buttonOK_Occurr, "clicked", function(widget, mydbWords = mergedD) {
    getSummaryStatsOCCUR(wordF = mydbWords)
    print("Message from GUI: occurrences plot done.")
  })

  gSignalConnect(buttonDesMostFreq, "clicked", function(widget, mydbWords = mergedD) {
    numWords <- entryDesMostFreq$getText()
    tryCatch(numWords <- as.numeric(numWords), warning = function(e){print(paste0("Warning from GUI: ", numWords, " is not a valid entry"))})
    sink(paste0('RESULTS/MostFreqWords_', numWords, '.txt'))
      cat('\n#######################\n### X=', numWords, '       ###\n#######################\n')
      try(print(getMostFreqWord(wordF = mergedD, numWords = numWords)), silent = TRUE)
    sink()
  })
  gSignalConnect(buttonDesXFreq, "clicked", function(widget, mydbWords = mergedD) {
    occuWords <- entryDesXFreq$getText()
    tryCatch(occuWords <- as.numeric(occuWords), warning = function(e){print(paste0("Warning from GUI: ", occuWords, " is not a valid entry"))})
    sink(paste0('RESULTS/OccurrenceWords_', occuWords, '.txt'))
      cat('\n#######################\n### X>=', occuWords, '       ###\n#######################\n')
      try(print(getXFreqWord(wordF = mergedD, occuWords = occuWords)), silent = TRUE)
    sink()
  })

  gSignalConnect(buttonOK_SaveIM, "clicked", function(widget, mydbWords = mergedD) {
    nameMergedD <- entrySaveIM$getText()
    assign(x = entrySaveIM$getText(), value = mergedD)
    do.call(save, list(nameMergedD, file = paste0('RESULTS/', paste(nameMergedD, "RData", sep = "."))))
    print(paste0("Message from GUI: word occurrence data.frame saved successfully as ", nameMergedD, ".RData."))
  })

  gSignalConnect(get('buttonFind_CompGroup1'), "clicked", function(widget, mydbWords = mergedD) { # buttonFind_CompGroup1
    currentObjectNameFile <- open_cbFile(window = main_window)
    currentObjectNameFile <- chartr("\\", "/", currentObjectNameFile)
    entryCompGroup1 <- get('entryCompGroup1')
    entryCompGroup1$setText(currentObjectNameFile)
  })
  gSignalConnect(get('buttonFind_CompGroup2'), "clicked", function(widget, mydbWords = mergedD) { # buttonFind_CompGroup2
    currentObjectNameFile <- open_cbFile(window = main_window)
    currentObjectNameFile <- chartr("\\", "/", currentObjectNameFile)
    entryCompGroup2 <- get('entryCompGroup2')
    entryCompGroup2$setText(currentObjectNameFile)
  })
  gSignalConnect(get('buttonFind_CompGroup3'), "clicked", function(widget, mydbWords = mergedD) { # buttonFind_CompGroup3
    currentObjectNameFile <- open_cbFile(window = main_window)
    currentObjectNameFile <- chartr("\\", "/", currentObjectNameFile)
    entryCompGroup3 <- get('entryCompGroup3')
    entryCompGroup3$setText(currentObjectNameFile)
  })
  gSignalConnect(get('buttonFind_CompGroup4'), "clicked", function(widget, mydbWords = mergedD) { # buttonFind_CompGroup4
    currentObjectNameFile <- open_cbFile(window = main_window)
    currentObjectNameFile <- chartr("\\", "/", currentObjectNameFile)
    entryCompGroup4 <- get('entryCompGroup4')
    entryCompGroup4$setText(currentObjectNameFile)
  })
  gSignalConnect(get('buttonFind_CompGroup5'), "clicked", function(widget, mydbWords = mergedD) { # buttonFind_CompGroup5
    currentObjectNameFile <- open_cbFile(window = main_window)
    currentObjectNameFile <- chartr("\\", "/", currentObjectNameFile)
    entryCompGroup5 <- get('entryCompGroup5')
    entryCompGroup5$setText(currentObjectNameFile)
  })

  gSignalConnect(buttonOK_CompValid, "clicked", function(widget, mydbWords = mergedD) {

    entryCompGroup1 <- get('entryCompGroup1')
    entryCompGroup2 <- get('entryCompGroup2')
    entryCompGroup3 <- get('entryCompGroup3')
    entryCompGroup4 <- get('entryCompGroup4')
    entryCompGroup5 <- get('entryCompGroup5')

    t1 <- utils::tail(strsplit(entryCompGroup1$getText(), split = "\\/")[[1]], 1)
    t2 <- utils::tail(strsplit(entryCompGroup2$getText(), split = "\\/")[[1]], 1)
    t3 <- utils::tail(strsplit(entryCompGroup3$getText(), split = "\\/")[[1]], 1)
    t4 <- utils::tail(strsplit(entryCompGroup4$getText(), split = "\\/")[[1]], 1)
    t5 <- utils::tail(strsplit(entryCompGroup5$getText(), split = "\\/")[[1]], 1)

    try(load(entryCompGroup1$getText()), silent = TRUE)
    try(load(entryCompGroup2$getText()), silent = TRUE)
    try(load(entryCompGroup3$getText()), silent = TRUE)
    try(load(entryCompGroup4$getText()), silent = TRUE)
    try(load(entryCompGroup5$getText()), silent = TRUE)

    g1 <- g2 <- g3 <- g4 <- g5 <- NULL
    try(g1 <- get(strsplit(utils::tail(strsplit(entryCompGroup1$getText(), split = "\\/")[[1]], 1), split = "\\.")[[1]][1]), silent = TRUE)
    try(g2 <- get(strsplit(utils::tail(strsplit(entryCompGroup2$getText(), split = "\\/")[[1]], 1), split = "\\.")[[1]][1]), silent = TRUE)
    try(g3 <- get(strsplit(utils::tail(strsplit(entryCompGroup3$getText(), split = "\\/")[[1]], 1), split = "\\.")[[1]][1]), silent = TRUE)
    try(g4 <- get(strsplit(utils::tail(strsplit(entryCompGroup4$getText(), split = "\\/")[[1]], 1), split = "\\.")[[1]][1]), silent = TRUE)
    try(g5 <- get(strsplit(utils::tail(strsplit(entryCompGroup5$getText(), split = "\\/")[[1]], 1), split = "\\.")[[1]][1]), silent = TRUE)

    listMergedD <- list(g1, g2, g3, g4, g5)
    names(listMergedD) <- c(t1, t2, t3, t4, t5)
    listMergedD <- Filter(Negate(is.null), listMergedD)

    if(length(listMergedD)>1){

      fileNames <- NULL

      std = checkboxCompLabelMerge$getActive()
      if(std==TRUE){
        listMergedD <- lapply(listMergedD, function(i){
          data.frame(word = as.character(i[,1]), freq = apply(i[,2:ncol(i)], MARGIN = 1, FUN = sum))
        })
        fileNames <- names(listMergedD)
      }else{
        fileNames <- as.vector(sapply(listMergedD, function(z){names(z)[2:ncol(z)]}))
      }

      words <- unique(as.vector(unique(unlist(sapply(listMergedD, function(i){unique(i[,1])})))))

      mydbWords <- data.frame(word = words)
      for(k in 1:length(listMergedD)){
        mydbWords <- try(merge(mydbWords, listMergedD[[k]], by.x = 1, by.y = 1, all = TRUE, suffixes = c(paste0(".x", k), paste0(".y", k))), silent = TRUE)
      }
      mydbWords[is.na(mydbWords)] <- 0
      colnames(mydbWords) <- c("word", fileNames)

      ### stem over the incidence-matrix
      dd <- NULL
      # Stem words
      dd$stem <- SnowballC::wordStem(as.character(mydbWords[,1]), language = "english") ###
      dd$word <- as.character(mydbWords[,1])
      dd$freq <- mydbWords[,2:ncol(mydbWords)]

      agg_freq <- sapply(1:ncol(dd$freq), function(i){stats::aggregate(freq[,i] ~ stem, data = dd, sum)}) ### a optimiser
      agg_freq <- cbind(data.frame(unlist(agg_freq[1,1])), data.frame(agg_freq[2,]))
      agg_word <- stats::aggregate(word ~ stem, data = dd, function(x) x[1]) ###
      dd <- cbind(word = agg_word[,2], agg_freq[,2:ncol(agg_freq)]) ###
      names(dd) <- names(mydbWords)
      mydbWords <- dd

      mydbWords <- mydbWords[order(apply(mydbWords[,2:ncol(mydbWords)], MARGIN = 1, FUN = sum), decreasing = T), ] # order matrix by incidence

      mergedD<<-mydbWords
      print("Message from GUI: groups of documents merged into a single data.frame.")
    }
  })

  return(vbox)
}

#' RGtk2 GUI function: Load the Graphical user Interface
#'
#' Load the Graphical user Interface in order to use \code{inpdfr} package through
#'   a user-friendly interface.
#' @details \code{inpdfr} package uses \code{RGtk2} package for its GUI. Non-linux
#'   users may need to download additional files such as the "gtk-file" icon, or
#'   the "hicolor" theme, which can be found by downloading GTK+ from http://www.gtk.org/.
#'   They are not needed for the GUI to work as intended, but you may get a "GTK-WARNING"
#'   when using \code{loadGUI()}. Feel free to ignore this warning. The RGtk2 GUI is
#'   not needed to access all funcionalities of \code{inpdfr} package. Some options
#'   are only available through the command line interface.
#' @examples
#' \dontrun{
#' loadGUI()
#' }
#' @export
loadGUI <- function(){
  main_window <- gtkWindow("toplevel", show = FALSE)
  main_window["title"] <- "GUI for inpdfr package."
  main_window$setDefaultSize(300, 250)
  ICONimage <- gdkPixbuf(filename = imagefile("gtk-logo-rgb.gif"))[[1]]
  main_window$set(icon = ICONimage)
  askQuit(main_window)
  vbox <- makeMenuMainWindow(main_window)
  vboxContent <- suppressWarnings(makeMainWindowsContent(main_window))
  vbox$add(vboxContent)
  main_window$add(vbox)
  main_window$show()
}
