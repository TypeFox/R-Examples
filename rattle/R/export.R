# Gnome R Data Miner: GNOME interface to R for Data Mining
#
# Time-stamp: <2015-05-17 08:56:28 gjw>
#
# Implement functionality associated with the Export button and Menu.
#
# Copyright (c) 2009-2013 Togaware Pty Ltd
#
# This files is part of Rattle.
#
# Rattle is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# Rattle is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Rattle. If not, see <http://www.gnu.org/licenses/>.

on_export_button_clicked <- function(action, window)
{
  # Wrap the actual call with a "try" so that the watch cursor turns
  # off even on error.
  
  setStatusBar()
  set.cursor("watch")
  on.exit(set.cursor())
  dispatchExportButton()

}

dispatchExportButton <- function()
{
  # Check which tab of notebook and dispatch to appropriate execute action

  ct <- getCurrentPageLabel(crv$NOTEBOOK)

  if (ct == crv$NOTEBOOK.CLUSTER.NAME)
  {  
    exportClusterTab()
  }
  else if (ct == crv$NOTEBOOK.MODEL.NAME)
  {
    exportModelTab()
  }
  else if (ct == crv$NOTEBOOK.ASSOCIATE.NAME)
  {
    exportAssociateTab()
  }
  else if (ct == crv$NOTEBOOK.DATA.NAME ||
           ct == crv$NOTEBOOK.TRANSFORM.NAME)
  {
    # For any of the DATA, SELECT, or TRANSFORM tabs, the logical
    # thing to EXPORT is the dataset.
    
    exportDataTab()
  }
##   else if (ct == crv$NOTEBOOK.EVALUATE.NAME)
##   {
##     exportEvaluateTab()
##   }
  else  if (ct == crv$NOTEBOOK.LOG.NAME) 
  {
    exportLogTab()
  }
  else if (ct == crv$NOTEBOOK.EVALUATE.NAME &&
           theWidget("evaluate_score_radiobutton")$getActive())
    # 091123 Users expect the Export to save the scores, so make this
    # the same as clicking the Execute button.
    executeEvaluateTab()
  else
    {
      # 100424 This is required for MS/Windows and Japanese and
      # sprintf for some reason - presumably a bug.
      if (isJapanese()) Encoding(ct) <- "unknown"
      infoDialog(sprintf(Rtxt("No export functionality is available for the",
                              "%s tab. Nothing done."), ct))
    }
}

## This is handled by the Cairo device save button now. Might want to
## interpret Export differently for these now.

## exportExploreTab <- function()
## {
##   if (theWidget("explore_distr_radiobutton")$getActive())
##     exportPlot("dist")
##   else if (theWidget("explore_correlation_radiobutton")$getActive())
##     exportPlot("corr")
##   else if (theWidget("explore_correlation_hier_radiobutton")$getActive())
##     exportPlot("hiercorr")
##   else if (theWidget("prcomp_radiobutton")$getActive())
##     exportPlot("prcomp")
##   else
##     infoDialog("No export functionality is available for the",
##                "selected option.")
## }

## ########################################################################

## exportEvaluateTab <- function()
## {
##   if (theWidget("risk_radiobutton")$getActive())
##     exportPlot("risk")
##   else if (theWidget("lift_radiobutton")$getActive())
##     exportPlot("lift")
##   else if (theWidget("roc_radiobutton")$getActive())
##     exportPlot("roc")
##   else if (theWidget("precision_radiobutton")$getActive())
##     exportPlot("precision")
##   else if (theWidget("sensitivity_radiobutton")$getActive())
##     exportPlot("sensitivity")
##   else
##     infoDialog("No export functionality from the Evaluate tab for",
##                "the selected option is yet available.")
## }

## exportPlot <- function(type="plot", devices=NULL)
## {
##   if (is.null(dev.list()))
##   {
##     warnDialog("There are currently no active graphics devices.",
##                "So there is nothing to export!",
##                "Please Execute (F2) to obtain a plot to export.")
##     return()
##   }

##   # Obtain a filename to save to. Ideally, this would also prompt for
##   # the device to export, and the fontsize, etc.

##   dialog <- gtkFileChooserDialog("Export Graphics (pdf, png, jpg)",
##                                  NULL, "save",
##                                  "gtk-cancel", RGtk2::GtkResponseType["cancel"],
##                                  "gtk-save", RGtk2::GtkResponseType["accept"])

##   if(not.null(crs$dataname))
##     dialog$setCurrentName(paste(get.stem(crs$dataname),
##                                 "_", type, ".pdf", sep=""))

##   ff <- gtkFileFilterNew()
##   ff$setName("Graphics Files")
##   ff$addPattern("*.pdf")
##   ff$addPattern("*.png")
##   ff$addPattern("*.jpg")
##   dialog$addFilter(ff)

##   ff <- gtkFileFilterNew()
##   ff$setName("All Files")
##   ff$addPattern("*")
##   dialog$addFilter(ff)
  
##   if (dialog$run() == RGtk2::GtkResponseType["accept"])
##   {
##     save.name <- dialog$getFilename()
##     dialog$destroy()
##   }
##   else
##   {
##     dialog$destroy()
##     return()
##   }

##   if (get.extension(save.name) == "") save.name <- sprintf("%s.pdf", save.name)
    
##   if (file.exists(save.name))
##     if ( ! questionDialog("A Graphics file of the name", save.name,
##                                 "already exists. Do you want to overwrite",
##                                 "this file?"))
##       return()
  
##   cur <- dev.cur()
##   ext <- get.extension(save.name)
##   if (ext == "pdf")
##     dev.copy(pdf, file=save.name, width=7, height=7)
##   else if (ext == "png")
##     dev.copy(png, file=save.name, width=700, height=700)
##   else if (ext == "jpg")
##     dev.copy(jpeg, file=save.name, width=700, height=700)
##   dev.off()
##   dev.set(cur)
  
##   infoDialog(sprintf("R Graphics: Device %d (ACTIVE)", cur),
##              "has been exported to", save.name)
## }

getWidgetOrObject <- function(dialog, name)
{
  if (crv$useGtkBuilder)
    return(dialog$getObject(name))
  else
    return(dialog$getWidget(name))
}

getExportSaveName <- function(mtype)
{
  # 090117 Request a filename to save the model to and return the
  # filename as a string.

  # Require the pmml package for either exporting to PMML or C (which
  # goes via PMML).
  
  lib.cmd <- "library(pmml, quietly=TRUE)"
  if (! (exists("pmml") ||
         packageIsAvailable("pmml", sprintf(Rtxt("export a %s model"),
                                            commonName(mtype)))))
      return(NULL)
  appendLog(packageProvides("pmml", "pmml"), lib.cmd)
  # Load the package unless we already have a pmml defined (through source).
  if (! exists("pmml")) eval(parse(text=lib.cmd))

  # Obtain filename to write the PMML or C code to.  081218 We use the
  # glade generated filechooser rather than my original hand-coded
  # one. It is much simpler to handle the formatting. It has been
  # modified to offer a choice of Class/Score and PMML/Info to the C
  # file.

  if (crv$useGtkBuilder)
  {
    dialogGUI <- RGtk2::gtkBuilderNew()
    dialogGUI$setTranslationDomain("R-rattle")
  }
  
  result <- try(etc <- file.path(path.package(package="rattle")[1], "etc"),
                silent=TRUE)
  if (inherits(result, "try-error"))
    if (crv$useGtkBuilder)
        dialogGUI$addFromFile(crv$rattleUI)
    else
      dialogGUI <- gladeXMLNew("rattle.glade",
                               root="export_filechooserdialog")
  else
      if (crv$useGtkBuilder)
        dialogGUI$addFromFile(file.path(etc, crv$rattleUI))
      else
        dialogGUI <- gladeXMLNew(file.path(etc,"rattle.glade"),
                                 root="export_filechooserdialog")
  if (crv$useGtkBuilder)
  {
    dialogGUI$getObject("export_filechooserdialog")$show()
    dialogGUI$connectSignals()
  }

  if (! crv$export.to.c.available)
    getWidgetOrObject(dialogGUI, "export_filechooser_options_table")$hide()

  dialog <- getWidgetOrObject(dialogGUI, "export_filechooserdialog")

  if (crv$export.to.c.available) dialog$setTitle(Rtxt("Export C or PMML"))
  #dialog$setIcon(crv$icon)
  
  if(not.null(crs$dataname))
    dialog$setCurrentName(paste(get.stem(crs$dataname), "_", mtype,
                                ifelse(crv$export.to.c.available, ".c", ".xml"),
                                sep=""))

  if (crv$export.to.c.available)
  {
    ff <- RGtk2::gtkFileFilterNew()
    ff$setName(Rtxt("C Files"))
    ff$addPattern("*.c")
    dialog$addFilter(ff)
  }

  ff <- RGtk2::gtkFileFilterNew()
  ff$setName(Rtxt("PMML Files"))
  ff$addPattern("*.xml")
  dialog$addFilter(ff)

  ff <- RGtk2::gtkFileFilterNew()
  ff$setName(Rtxt("All Files"))
  ff$addPattern("*")
  dialog$addFilter(ff)

  if (crv$export.to.c.available)
  {
    if (mtype %in% c("glm"))
      # 090629 The default for multinomial is class but we don't allow
      # the user to choose probablity - was that always the case?
      if (multinomialTarget())
        getWidgetOrObject(dialogGUI,
                          "export_filechooser_class_radiobutton")$setActive(TRUE)
      else
        getWidgetOrObject(dialogGUI,
                          "export_filechooser_probabilities_radiobutton")$setActive(TRUE)

    # 081218 Add glm when implemented.
    
    if (!binomialTarget() || mtype %notin% c("rpart", "glm"))
    {
      getWidgetOrObject(dialogGUI,
                        "export_filechooser_target_label")$setSensitive(FALSE)

      getWidgetOrObject(dialogGUI,
                        "export_filechooser_class_radiobutton")$setSensitive(FALSE)

      getWidgetOrObject(dialogGUI,
                        "export_filechooser_probabilities_radiobutton")$
      setSensitive(FALSE)
    }

    if (mtype %in% c("survival"))
    {
      getWidgetOrObject(dialogGUI,
                        "export_filechooser_class_radiobutton")$setLabel(Rtxt("Time"))

      getWidgetOrObject(dialogGUI,
                        "export_filechooser_probabilities_radiobutton")$
      setLabel(Rtxt("Risk"))

      getWidgetOrObject(dialogGUI,
                        "export_filechooser_probabilities_radiobutton")$
      setActive(class(crs$survival) == "coxph")

      if (class(crs$survival) == "coxph")
      {
        getWidgetOrObject(dialogGUI,
                          "export_filechooser_probabilities_radiobutton")$
        setSensitive(TRUE)
      }
      else
      {
        getWidgetOrObject(dialogGUI,
                          "export_filechooser_class_radiobutton")$
        setSensitive(TRUE)
      }
    }
  }

  if (dialog$run() == RGtk2::GtkResponseType["accept"])
  {
    save.name <- dialog$getFilename()
    save.type <- dialog$getFilter()$getName()
    if (crv$export.to.c.available)
    {
      includePMML <- getWidgetOrObject(dialogGUI,
                                       "export_filechooser_pmml_checkbutton")$getActive()

      includeMetaData <- getWidgetOrObject(dialogGUI,
                                           "export_filechooser_metadata_checkbutton")$
      getActive()

      exportClass <- getWidgetOrObject(dialogGUI,
                                       "export_filechooser_class_radiobutton")$
      getActive()
    }
    dialog$destroy()
  }
  else
  {
    dialog$destroy()
    return(NULL)
  }

  ext <- tolower(get.extension(save.name))

  if (crv$export.to.c.available)
  {
    attr(save.name, "includePMML") <- includePMML
    attr(save.name, "includeMetaData") <- includeMetaData
    attr(save.name, "exportClass") <- exportClass
  }

  Encoding(save.name) <- "UTF-8"
  return(save.name)
}

generateExportPMMLtoC <- function(model.name, save.name, TV)
{
  # 110116 Introduce an encoding to the file saved. Then under
  # MS/Windows the file is saved as UTF-8 rather than SHIFT-JIS. On
  # Linux it is saved as UTF-8 anyhow. However, I note that we seem to
  # be working pretty hard to make everything UTF-8 on MS/Windows. It
  # seems that we are fighting against "nature" here - is there
  # something about encodings and R that we are missing. Maybe that is
  # the key as to why everything works okay on Linux (UTF-8) but we
  # battle with MS/Windows.
  #
  # From cran.r-project.org/doc/manuals/R-data.html: "The hard part is
  # to know what file encoding to use. For use on Windows, it is best
  # to use what Windows calls `Unicode' that is "UTF-16LE". (Even
  # then, Windows applications may expect a Byte Order Mark which the
  # implementation of iconv used by R may or may not add depending on
  # the platform.), Using UTF-8 is a good way to make portable files
  # that will not easily be confused with any other encoding, but even
  # Mac OS X applications (where UTF-8 is the system encoding) may not
  # recognize them, and Windows applications are most unlikely
  # to. Apparently Excel:mac 2004/8 expects .csv files in "macroman"
  # encoding (the encoding used in much earlier versions of Mac OS)."

  #110122 IBI suggest the C code be SJIS for now.
  
  export.cmd <- paste('con <- file("%s", open="w")', # 110122, encoding="UTF8")',
                      "\ncat(pmmltoc(%s%s%s,",
                      '\n            name="%s",',
                      '\n            includePMML=%s,',
                      '\n            includeMetaData=%s,',
                      '\n            exportClass=%s),',
                      '\n    file=con)',
                      '\nclose(con)', sep="")

  export.cmd <- sprintf(export.cmd,
                        fixWindowsSlash(save.name),
                        ifelse(isWindows() && isJapanese(),
                               paste("paste('<?xml version=\"1.0\"",
                                     "encoding=\"shift_jis\"?>\\n',",
                                     "\n                  "), ""),
                        "toString(eval(parse(text=pmml.cmd)))",
                        ifelse(isWindows() && isJapanese(), ")", ""),
                        model.name,
                        ifelse(attr(save.name, "includePMML"), "TRUE", "FALSE"),
                        ifelse(attr(save.name, "includeMetaData"),
                               sprintf('rattle:::getTextviewContent("%s")', TV),
                               '"\\"Not Included\\""'),
                        ifelse(attr(save.name, "exportClass"), "TRUE", "FALSE"))
  return(export.cmd)
}
