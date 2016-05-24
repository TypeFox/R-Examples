# if(grepl("darwin",R.Version()$os)) {
#   gettext  <- c
#   gettextf <- sprintf  
# } else {
  gettext  <- function(...) base::gettext(..., domain = "R-Rz")
  gettextf <- function(...) base::gettextf(..., domain = "R-Rz")
  stop     <- function(...) base::stop(..., call.=FALSE, domain = "R-Rz")  
# }

units   <- c("npc", "cm", "inches", "mm", "points", "picas", "bigpts",
             "dida", "cicero", "scaledpts", "lines", "char", "native", "snpc")

#fixTranslations <- function(w){
#  if ("GtkLabel" %in% class(w))
#    w$setLabel(gettext(w$getLabel()))
#  else if ("GtkNotebook" %in% class(w))
#    lapply(gtkChildren(w),
#           function(wc)
#             w$getTabLabel(wc)$setLabel(gettext(w$getTabLabelText(wc))))
#
#  if ("GtkContainer" %in% class(w))
#    lapply(gtkChildren(w), fixTranslations)
#  
#  return()
#}

fileCheck <- function(filename, parent){
  if (file.exists(filename)){
    dialog <- gtkMessageDialogNew(parent, "destroy-with-parent",
                                   GtkMessageType["question"],
                                   GtkButtonsType["ok-cancel"],
                                   paste(gettext("File Path: "), filename, "\n",
                                         gettext("This file already exists. Overwrite it?"), sep=""))
    response <- dialog$run()
    dialog$hide()
    if (response==GtkResponseType["ok"]){
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    return(TRUE)
  }
}

#script.buffer <-
#setRefClass("RzScriptBuffer",
#  fields = c("script", "buffer")
#  )

#sb.obj <- script.buffer$new(script="test", buffer=gtkTextBufferNew())

localize <- function(char) iconv(char, "UTF-8", localeToCharset()[1])


fontsRegisterScript <- function(fonts){
  fonts.win <- paste(sprintf("    \"%s\" = windowsFont(\"%s\")", fonts, fonts), collapse=",\n")
  fonts.win <- sprintf("  windowsFonts(\n%s\n  )", fonts.win)
  
  fonts.x11 <- paste(
    sprintf("    \"%s\" = X11Font(\"-*-%s-*-*-*-*-*-*-*-*-*-*-*-*\")", fonts, fonts),
    collapse=",\n")
  fonts.x11 <- sprintf("  X11Fonts(\n%s\n  )", fonts.x11)
  
  fonts.quartz <- paste(sprintf("    \"%s\" = quartzFont(rep(\"%s\", 4))", fonts, fonts),
                       collapse=",\n")
  fonts.quartz <- sprintf("  quartzFonts(\n%s\n  )", fonts.quartz)

  script <- paste(
    "if (grepl(\"mingw\", R.Version()$os)) {", fonts.win,
    "} else if (grepl(\"darwin\", R.Version()$os)) {", fonts.x11, fonts.quartz,"}",
    sep="\n")
  return(script)
}

gtkFileChooserDialogFilteredNew <- function(title, parent=NULL,
                                            action=GtkFileChooserAction["open"],
                                            file.type.list){
  dialog <- gtkFileChooserDialogNew(title=title, parent=parent,
                                    action=action,
                                    "gtk-open", GtkResponseType["accept"],
                                    "gtk-cancel", GtkResponseType["cancel"], 
                                    show=FALSE)
  for (i in seq_along(file.type.list)) {
    filter <- gtkFileFilterNew()
    filter$setName(file.type.list[[i]]$name)
    filter$addPattern(file.type.list[[i]]$pattern)
    dialog$addFilter(filter)
  }
  class(dialog) <- c("GtkFileChooserDialogFiltered", class(dialog))
  return(dialog)
}

gtkFileChooserDialogFilteredActivate <- function(obj){
  if (! exists("theme_rz", envir=.GlobalEnv)) {
    rzTools$sync("theme_rz", theme_grey)        
  }
  
  if (obj$run() == GtkResponseType["accept"]) {
    filename <- localize(obj$getFilename())
    filetype <- localize(obj$getFilter()$getName())
    return(list(filename=filename, filetype=filetype))
  } else {
    return(NULL)
  }
}

gtkFileChooserDialogFilteredRun <- function(obj) gtkDialogRun(obj)

Rz <- function(...){
  if (is.null(rzTools$getMain())) {
    rzTools$setMain(new("RzMain"))    
  } else {
    rzTools$getMain()$show()
  }
}

RzThemeEditor <- function(new=FALSE, ...) {
  if (new | is.null(rzTools$getThemeEditor())) {
    rzTools$setThemeEditor(new("RzPlotTheme"))    
  }
  rzTools$getThemeEditor()$show()
}

gtkInfoBarRzNew <- function(show=TRUE){
    obj <- gtkInfoBarNew(show=TRUE)
    class(obj) <- c("gtkInfoBarRz", class(obj))
    return(obj)
}

gtkInfoBarRzSetText <- function(obj, txt){
  label <- obj$getContentArea()$getChildren()[[1]]
  label$setText(txt)
}

write.spss <-
function (df, datafile, codefile, varlabels, varnames = NULL) {
    dfn <- lapply(df, function(x) if (is.factor(x)) as.numeric(x) else x)
    write.table(dfn, file = datafile, row.names = FALSE, col.names = FALSE, 
        sep = ",", quote = FALSE, na = "", eol = ",\n")
    if (is.null(varnames)) {
        varnames <- names(df)
    }
    varnames <- gsub("[^[:alnum:]_\\$@#]", "\\.", varnames)
    dl.varnames <- varnames
    if (any(chv <- sapply(df, is.character))) {
        lengths <- sapply(df[chv], function(v) max(nchar(v)))
        if (any(lengths > 255L)) 
            stop("Cannot handle character variables longer than 255")
        lengths <- paste("(A", lengths, ")", sep = "")
        star <- ifelse(c(FALSE, diff(which(chv) > 1L)), " *", 
            " ")
        dl.varnames[chv] <- paste(star, dl.varnames[chv], lengths)
    }
    cat("DATA LIST FILE=", foreign:::adQuote(datafile), " free (\",\")\n", 
        file = codefile)
    cat("/", dl.varnames, " .\n\n", file = codefile, append = TRUE)
    cat("VARIABLE LABELS\n", file = codefile, append = TRUE)
    cat(paste(varnames, foreign:::adQuote(varlabels), "\n"), ".\n", file = codefile, 
        append = TRUE)
    factors <- sapply(df, is.factor)
    if (any(factors)) {
        cat("\nVALUE LABELS\n", file = codefile, append = TRUE)
        for (v in which(factors)) {
            cat("/\n", file = codefile, append = TRUE)
            cat(varnames[v], " \n", file = codefile, append = TRUE)
            levs <- levels(df[[v]])
            cat(paste(seq_along(levs), foreign:::adQuote(levs), "\n", sep = " "), 
                file = codefile, append = TRUE)
        }
        cat(".\n", file = codefile, append = TRUE)
    }
    cat("\nEXECUTE.\n", file = codefile, append = TRUE)
}

write.stata <- function (df, datafile, codefile, varlabels) 
{
    write.table(df, file = datafile, row.names = FALSE, col.names = FALSE, 
        sep = ",", quote = FALSE, na = ".")
    nms <- names(df)
    varlabels <- paste("label variable ", nms, " \"",varlabels, "\"", sep="", collapse="\n")
    factors <- sapply(df, is.factor) | sapply(df, is.character)
    formats <- paste(nms, "fmt", sep = "_")
    nms <- ifelse(factors, paste(nms, formats, sep = ":"), nms)
    cat("infile", nms, " using ", datafile, ", automatic\n", varlabels,
        file = codefile)
}

check.class <- function(obj, class){
  result <- NULL
  if (class=="numeric") {
    result <- try(eval(parse(text=(sprintf("c(%s)", obj)))), silent=TRUE)
    if(class(result) != "numeric") return(NULL)
    
  } else if (class=="any") {
    result <- try(eval(parse(text=(sprintf("c(%s)", obj)))), silent=TRUE)
    if(class(result) != "numeric") return(obj)
    
  } else if (class=="formula") {
    result <- try(as.formula(obj), silent=TRUE)
    if(class(result) != "formula") return(NULL)
    
  } else if (class=="logical") {
    result <- try(as.logical(obj), silent=TRUE)
    if(class(result) != "logical" | is.na(result)) return(NULL)    
    
  } else {
    result <- obj
  }
  return(result)
}

calc.hist.breaks <- function(breaks, x){
  suppressWarnings(binwidth <- as.numeric(breaks))
  if(breaks[1]=="based on Sturges"){
    breaks <- nclass.Sturges(x)
  } else if (breaks[1]=="based on Freedman-Diaconis"){
    breaks <- nclass.FD(x)
  } else if (breaks[1]=="based on Scott"){
    breaks <- nclass.scott(x)
  } else {
    breaks <- NA
  }
  
  if(all((!is.na(breaks)))) {
    breaks <- pretty(range(as.numeric(x), na.rm=TRUE), n = breaks, min.n = 1)
  }
  
  if (all(!is.na(binwidth))) {
    return(sprintf("binwidth=%s", deparse(binwidth)))
  } else if (all(!is.na(breaks))) {
    breaks <- deparse(breaks)
    breaks <- paste(breaks, collapse="")
    return(sprintf("breaks=%s", breaks))
  }
}

gtkScrolledWindowWithViewportNew <- function(hadjustment = NULL, vadjustment = NULL, show = TRUE){
  viewport <- gtkViewportNew()
  viewport$setShadowType(GtkShadowType["none"])
  scrolledWindow <- gtkScrolledWindowNew(hadjustment, vadjustment, show)
  scrolledWindow$setShadowType(GtkShadowType["none"])
  scrolledWindow$setPolicy(GtkPolicyType["automatic"], GtkPolicyType["automatic"])
  scrolledWindow$setBorderWidth(3)
  scrolledWindow$add(viewport)
  class(scrolledWindow) <- c("GtkScrolledWindowWithViewport", class(scrolledWindow))
  return(scrolledWindow)
}

gtkScrolledWindowWithViewportAdd <- function(object, child){
  object$getChild()$add(child)
}

buildPlotOptionPage <- function(widget){
  main <- gtkScrolledWindowWithViewportNew()
  vbox <- gtkVBoxNew()
  vbox$packStart(widget, expand=FALSE)
  main$add(vbox)
  return(main)
}


gtkColorSelectionForPlotNew <- function(show = TRUE){
  obj <- gtkColorSelectionNew(show)
  obj$setHasPalette(TRUE)
  obj$setHasOpacityControl(FALSE)
  class(obj) <- c("GtkColorSelectionForPlot", class(obj))
  return(obj)
}

gtkColorSelectionForPlotGetColor <- function(object){
  color <- object$getCurrentColor()
  color <- rgb(color$color$red, color$color$green,color$color$blue, 0, maxColorValue=65535)
  color <- paste(strsplit(color, "")[[1]][1:7], collapse="")
  return(color)
}

gtkColorSelectionForPlotSetColor <- function(object, color){
  color <- gdkColorParse(color)$color
  color <- object$setCurrentColor(color)
}

gtkColorSelectionDialogForPlotNew <- function(show=TRUE, parent){
  colorSelection <- gtkColorSelectionForPlotNew()
  obj <- gtkDialogNewWithButtons(gettext("Select colour"), parent,
                                 "destroy-with-parent",
                                 "gtk-ok", GtkResponseType["ok"], 
                                 "gtk-cancel", GtkResponseType["cancel"],
                                 show=show)
  obj["resizable"] <- FALSE
  obj$getContentArea()$setBorderWidth(5)
  obj$getContentArea()$packStart(colorSelection, expand=TRUE, fill=TRUE)
  obj$setData("colorSelection", colorSelection)
  class(obj) <- c("GtkColorSelectionDialogForPlot", class(obj))  
  return(obj)
}

gtkColorSelectionDialogForPlotGetColor <- function(object){
  colorSelection <- object$getData("colorSelection")
  colorSelection$getColor()
}

gtkColorSelectionDialogForPlotSetColor <- function(object, color){
  colorSelection <- object$getData("colorSelection")
  colorSelection$setColor(color)
}

gtkColorSelectionWidgetNew <- function(homogeneous = NULL, spacing = NULL, show = TRUE, parent, colour = "#FFFFFF"){
  model.color <- rGtkDataFrameNew(data.frame(colors()))
  comboentry <- gtkComboBoxEntryNewWithModel(model.color, 0)
  comboentry["width-request"] <- 1
  button <- gtkButtonNew()
  button$setImage(gtkImageNewFromStock(GTK_STOCK_SELECT_COLOR, GtkIconSize["button"]))
  dialog <- gtkColorSelectionDialogForPlotNew(show=FALSE, parent=parent)
  hbox <- gtkHBoxNew(homogeneous, spacing, show)
  hbox$packStart(comboentry, expand=TRUE, fill=TRUE)
  hbox$packStart(button, expand=FALSE)

  gSignalConnect(dialog, "response", function(obj, response){
    if (response == GtkResponseType["ok"]) {
      comboentry$getChild()$setText(dialog$getColor())
    }
  })
  
  gSignalConnect(button, "clicked", function(...){
    dialog$setColor(comboentry$getActiveText())
    dialog$run()
    dialog$hide()
  })
  
  hbox$setData("entry", comboentry$getChild())
  class(hbox) <- c("GtkColorSelectionWidget", class(hbox))
  return(hbox)
}

gtkColorSelectionWidgetSetColor <- function(object, colour){
  entry <- object$getData("entry")
  entry$setText(colour)
}

gtkColorSelectionWidgetGetEntry <- function(object){
  entry <- object$getData("entry")
  return(entry)
}


gtkFontSelectionDialog2New <- function(fontname=NULL, show=FALSE, parent=NULL){
  fontSelection <- gtkFontSelectionNew()
  if (!is.null(fontname)) {
    fontSelection$setFontName(fontname)
  }
  fontFamily  <- fontSelection$getFamilyList()$getParent()
  fontPreview <- fontSelection$getPreviewEntry()$getParent()
  fontSelection$getChildren()[[1]]$remove(fontFamily)
  fontSelection$getChildren()[[2]]$remove(fontPreview)
  
  label1 <- gtkLabelNew(gettext("Family:"))
  hbox1  <- gtkHBoxNew()
  hbox1$packStart(label1, expand=FALSE, padding=2)
  label2 <- gtkLabelNew(gettext("Preview:"))
  hbox2  <- gtkHBoxNew()
  hbox2$packStart(label2, expand=FALSE, padding=2)
  vbox   <- gtkVBoxNew(spacing=5)
  vbox$packStart(hbox1      , expand=FALSE)
  vbox$packStart(fontFamily , expand=TRUE)
  vbox$packStart(hbox2      , expand=FALSE)
  vbox$packStart(fontPreview, expand=FALSE)
  
  win <- gtkWindowNew(show=FALSE)
  win$add(fontSelection) # do NOT remove!
  
  obj <- gtkDialogNewWithButtons(gettext("Select Font Family"), parent,
                                 "destroy-with-parent",
                                 "gtk-ok", GtkResponseType["ok"], 
                                 "gtk-cancel", GtkResponseType["cancel"],
                                 show=show)
  obj["width-request"] <- 300
  obj$getContentArea()$setBorderWidth(5)
  obj$getContentArea()$packStart(vbox, expand=TRUE, fill=TRUE)
  obj$setData("fontSelection", fontSelection)
  class(obj) <- c("GtkFontSelectionDialog2", class(obj))
  return(obj)
}

gtkFontSelectionDialog2GetFontName <- function(object){
  fontSelection <- object$getData("fontSelection")
  return(fontSelection$getFamily()$getName())
}
gtkFontSelectionDialog2SetFontName <- function(object, fontname){
  fontSelection <- object$getData("fontSelection")
  fontSelection$setFontName(sprintf("%s 10", fontname))
}

gtkFontSelectionButtonNew <- function(fontname=NULL, parent){
  font.dialog <- gtkFontSelectionDialog2New(fontname, parent=parent)
  fontname    <- font.dialog$getFontName()
  obj <- gtkButtonNewWithLabel(fontname)
  obj["width-request"] <- 1
  obj$getChild()$modifyFont(pangoFontDescriptionFromString(sprintf("%s 10", fontname)))
  
  gSignalConnect(font.dialog, "response", function(dialog, response){
    if (response == GtkResponseType["ok"]) {
      new.fontname <- font.dialog$getFontName()
      obj$setLabel(new.fontname)
      obj$getChild()$modifyFont(pangoFontDescriptionFromString(sprintf("%s 10", new.fontname)))
    }
  })
  
  gSignalConnect(obj, "clicked", function(...){
    font.dialog$setFontName(obj$getLabel())
    font.dialog$run()
    font.dialog$hide()
  })
  
  obj$setData("font.dialog", font.dialog)
  class(obj) <- c("GtkFontSelectionButton", class(obj))
  return(obj)
}

gtkFontSelectionButtonGetFontName <- function(object){
  font.dialog <- object$getData("font.dialog")
  return(font.dialog$getFontName())
}

gtkFontSelectionButtonSetFontName <- function(object, fontname){
  font.dialog <- object$getData("font.dialog")
  font.dialog$setFontName(fontname)  
}

gtkPlotUnitNew <- function(name, value=NULL, unit="cm", inherit.from=NULL){
  inherit <- NULL
  if (is.null(value)) {
    value <- 1
    inherit <- TRUE
  } else {
    inherit <- FALSE
  }
  adj     <- gtkAdjustmentNew(value, -100, 100, 0.1)
  button  <- gtkSpinButtonNew(adj, climb.rate=0.1, digits=3)
  button$setValue(value)
  button["width-request"] <- 1
  combo <- gtkComboBoxEntryNewText()
  combo["width-request"] <- 80
  for(i in units) combo$appendText(i)
  combo$setActive(which(units == unit) - 1)
  
  hbox <- gtkHBoxNew(spacing=2)
  hbox$packStart(button)
  hbox$packStart(combo, expand=FALSE)
  hbox2 <- gtkHBoxNew(spacing=2)
  hbox2$setBorderWidth(5)
  hbox2$packStart(hbox)
  
  tooltip.text <- character()
  if (is.null(inherit.from)) {
    tooltip.text <- gettext("Inherit from the current theme")
  } else {
    tooltip.text <- gettextf("Inherit from <span font_style='italic' font_weight='bold' font_size='large'>%s</span>", inherit.from)
  }

  check  <- gtkToggleButtonNew()
  check$setImage(gtkImageNewFromStock(GTK_STOCK_GO_DOWN, GtkIconSize["button"]))
  check["tooltip-markup"] <- tooltip.text
  hbox2$packStart(check, expand=FALSE)    
  
  gSignalConnect(check, "toggled", function(obj){
    inherit <- obj$getActive()
    hbox$setSensitive( !inherit )
    obj$setData("inherit", inherit)
  })
  
  check$setActive(inherit)
  
  vbox <- gtkVBoxNew()
  vbox$packStart(hbox2, expand=FALSE)
  obj <- gtkFrameNew(name)
  obj$setShadowType(GtkShadowType["etched-in"])
  obj$add(vbox)
  obj$showAll()
  
  obj$setData("check", check)
  obj$setData("button", button)
  obj$setData("combo", combo)
  
  gSignalConnect(button, "value-changed", function(...){
    script <- sprintf("%s = unit(%s, %s)", name, button$getValue(), deparse(combo$getActiveText()))
    obj$setData("script", script)    
  })
  gSignalConnect(combo, "changed", function(...){
    script <- sprintf("%s = unit(%s, %s)", name, button$getValue(), deparse(combo$getActiveText()))
    obj$setData("script", script)    
  })
  
  script <- sprintf("%s = unit(%s, %s)", name, button$getValue(), deparse(combo$getActiveText()))
  obj$setData("script", script)
  
  class(obj) <- c("GtkPlotUnit", "GktPlotThemeWidgets", class(obj))
  return(obj)
}

gtkPlotUnitSetValue <- function(object, value){
  check  <- object$getData("check")
  button <- object$getData("button")
  combo  <- object$getData("combo")
  if (is.null(value)) {
    check$setActive(TRUE)          
  } else {
    button$setValue(value)
    combo$setActive( which(units == attr(value, "unit")) - 1)
    check$setActive(FALSE)
  }
}

gtkPlotUnit2New <- function(name, value=NULL, unit="cm"){
  inherit <- NULL
  if (is.null(value)) {
    value <- c(1,1,1,1)
    inherit <- TRUE
  } else {
    inherit <- FALSE
  }
  units   <- c("npc", "cm", "inches", "mm", "points", "picas", "bigpts",
               "dida", "cicero", "scaledpts", "lines", "char", "native", "snpc")
  adj1  <- gtkAdjustmentNew(value[1], -100, 100, 0.1)
  adj2  <- gtkAdjustmentNew(value[2], -100, 100, 0.1)
  adj3  <- gtkAdjustmentNew(value[3], -100, 100, 0.1)
  adj4  <- gtkAdjustmentNew(value[4], -100, 100, 0.1)
  button1  <- gtkSpinButtonNew(adj1, climb.rate=0.1, digits=3)
  button1$setValue(value[1])
  button1["width-request"] <- 1
  button2  <- gtkSpinButtonNew(adj2, climb.rate=0.1, digits=3)
  button2$setValue(value[2])
  button2["width-request"] <- 1
  button3  <- gtkSpinButtonNew(adj3, climb.rate=0.1, digits=3)
  button3$setValue(value[3])
  button3["width-request"] <- 1
  button4  <- gtkSpinButtonNew(adj4, climb.rate=0.1, digits=3)
  button4$setValue(value[4])
  button4["width-request"] <- 1
  combo <- gtkComboBoxEntryNewText()
  combo["width-request"] <- 80
  for(i in units) combo$appendText(i)
  combo$setActive(which(units == unit) - 1)
  
  check  <- gtkToggleButtonNew()
  check$setImage(gtkImageNewFromStock(GTK_STOCK_GO_DOWN, GtkIconSize["button"]))
  check["tooltip-markup"] <- gettext("Inherit from the current theme")
  
  hbox <- gtkHBoxNew(spacing=2)
  hbox$packStart(combo)
  hbox$packStart(check, expand=FALSE)
  
  table <- gtkTableNew()
  table$setBorderWidth(5)
  table$setColSpacings(2)
  table$setRowSpacings(5)
  table$attach(gtkLabelNew(gettext("Top"))   , 0, 1, 0, 1, "shrink", "shrink")
  table$attach(gtkLabelNew(gettext("Right")) , 2, 3, 0, 1, "shrink", "shrink")
  table$attach(gtkLabelNew(gettext("Bottom")), 0, 1, 1, 2, "shrink", "shrink")
  table$attach(gtkLabelNew(gettext("Left"))  , 2, 3, 1, 2, "shrink", "shrink")
  table$attach(gtkLabelNew(gettext("Unit"))  , 0, 1, 2, 3, "shrink", "shrink")
  table$attach(button1, 1, 2, 0, 1, 5, "shrink")
  table$attach(button2, 3, 4, 0, 1, 5, "shrink")
  table$attach(button3, 1, 2, 1, 2, 5, "shrink")
  table$attach(button4, 3, 4, 1, 2, 5, "shrink")  
  table$attach(hbox,   1, 4, 2, 3, 5, "shrink")  

  gSignalConnect(check, "toggled", function(obj){
    inherit <- obj$getActive()
    button1$setSensitive( !inherit )
    button2$setSensitive( !inherit )
    button3$setSensitive( !inherit )
    button4$setSensitive( !inherit )
    combo$setSensitive( !inherit )
  })
  
  check$setActive(inherit)
  
  
  obj <- gtkFrameNew(name)
  obj$setShadowType(GtkShadowType["etched-in"])
  obj$add(table)
  obj$showAll()
  
  obj$setData("inherit", inherit)
  obj$setData("check", check)
  obj$setData("button1", button1)
  obj$setData("button2", button2)
  obj$setData("button3", button3)
  obj$setData("button4", button4)
  obj$setData("combo", combo)
  
  gSignalConnect(button1, "value-changed", function(...){
    values <- c(button1$getValue(), button2$getValue(), button3$getValue(), button4$getValue())
    script <- sprintf("%s = unit(%s, %s)", name, deparse(values), deparse(combo$getActiveText()))
    obj$setData("script", script)
  })
  gSignalConnect(button2, "value-changed", function(...){
    values <- c(button1$getValue(), button2$getValue(), button3$getValue(), button4$getValue())
    script <- sprintf("%s = unit(%s, %s)", name, deparse(values), deparse(combo$getActiveText()))
    obj$setData("script", script)
  })
  gSignalConnect(button3, "value-changed", function(...){
    values <- c(button1$getValue(), button2$getValue(), button3$getValue(), button4$getValue())
    script <- sprintf("%s = unit(%s, %s)", name, deparse(values), deparse(combo$getActiveText()))
    obj$setData("script", script)
  })
  gSignalConnect(button4, "value-changed", function(...){
    values <- c(button1$getValue(), button2$getValue(), button3$getValue(), button4$getValue())
    script <- sprintf("%s = unit(%s, %s)", name, deparse(values), deparse(combo$getActiveText()))
    obj$setData("script", script)
  })
  gSignalConnect(combo, "changed", function(...){
    values <- c(button1$getValue(), button2$getValue(), button3$getValue(), button4$getValue())
    script <- sprintf("%s = unit(%s, %s)", name, deparse(values), deparse(combo$getActiveText()))
    obj$setData("script", script)    
  })
  
  values <- c(button1$getValue(), button2$getValue(), button3$getValue(), button4$getValue())
  script <- sprintf("%s = unit(%s, %s)", name, deparse(values), deparse(combo$getActiveText()))
  obj$setData("script", script)
  
  class(obj) <- c("GtkPlotUnit2", "GktPlotThemeWidgets", class(obj))
  return(obj)
}

gtkPlotUnit2SetValue <- function(object, value){
  check  <- object$getData("check")
  button1 <- object$getData("button1")
  button2 <- object$getData("button2")
  button3 <- object$getData("button3")
  button4 <- object$getData("button4")
  combo  <- object$getData("combo")
  if (is.null(value)) {
    check$setActive(TRUE)          
  } else {
    button1$setValue(value[1])
    button2$setValue(value[2])
    button3$setValue(value[3])
    button4$setValue(value[4])
    combo$setActive( which(units == attr(value[1], "unit")) - 1)
    check$setActive(FALSE)
  }
}



gtkPlotAlignNew <- function(name, value=NULL){
  inherit <- NULL
  if (is.null(value)) {
    value <- 0
    inherit <- TRUE
  } else {
    inherit <- FALSE
  }
  
  adj    <- gtkAdjustmentNew(value, 0, 1, 0.01)
  button <- gtkSpinButtonNew(adj, climb.rate=0.01, digits=2)
  button$setValue(value)
  button["width-request"] <- 1

  
  check  <- gtkToggleButtonNew()
  check$setImage(gtkImageNewFromStock(GTK_STOCK_GO_DOWN, GtkIconSize["button"]))
  check["tooltip-markup"] <- gettext("Inherit from the current theme")
  
  hbox <- gtkHBoxNew(spacing=2)
  hbox$packStart(button)
  hbox$packStart(check, expand=FALSE)

  vbox <- gtkVBoxNew(spacing=2)
  vbox$packStart(hbox, expand=FALSE)
  vbox$setBorderWidth(5)
  
  gSignalConnect(check, "toggled", function(obj){
    inherit <- obj$getActive()
    button$setSensitive( !inherit )
  })
  
  check$setActive(inherit)
  
  obj <- gtkFrameNew(name)
  obj$setShadowType(GtkShadowType["etched-in"])
  obj$add(vbox)
  obj$showAll()
  
  obj$setData("inherit", inherit)
  obj$setData("check", check)
  obj$setData("button", button)
  
  
  gSignalConnect(button, "value-changed", function(...){
    script <- sprintf("%s = %s", name, button$getValue())
    obj$setData("script", script)    
  })
  
  script <- sprintf("%s = %s", name, button$getValue())
  obj$setData("script", script)
  
  class(obj) <- c("GtkPlotAlign", "GktPlotThemeWidgets", class(obj))
  return(obj)
}

gtkPlotAlignSetValue <- function(object, value){
  check  <- object$getData("check")
  button1 <- object$getData("button")
  if (is.null(value)) {
    check$setActive(TRUE)
  } else {
    button1$setValue(value)
    check$setActive(FALSE)
  }
}


gtkPlotPositionNew <- function(name, value=NULL){
  labels <- NULL
  inherit <- NULL
  
  combo <- gtkComboBoxNewText()
  combo$show()
  if (name == "legend.position") {
    labels <- c("top", "right", "bottom", "left", "specify")
    for(i in labels) combo$appendText(i)
  } else {
    labels <- c("center", "specify")
    for(i in labels) combo$appendText(i)    
  }
  if (is.null(value)) {
    combo$setActive(0)
    value <- c(0, 0)
    inherit <- TRUE
  } else {
    if (is.numeric(value)) {
      combo$setActive(which("specify" == labels) - 1)
    } else {
      combo$setActive(which(value == labels) - 1)      
      value <- c(0, 0)
    }
    inherit <- FALSE
  }
  
  adj1    <- gtkAdjustmentNew(value[1], 0, 1, 0.01)
  button1 <- gtkSpinButtonNew(adj1, climb.rate=0.01, digits=2)
  button1$setValue(value[1])
  button1["width-request"] <- 1
  adj2    <- gtkAdjustmentNew(value[2], 0, 1, 0.01)
  button2 <- gtkSpinButtonNew(adj2, climb.rate=0.01, digits=2)
  button2$setValue(value[2])
  button2["width-request"] <- 1
  
  table <- gtkTableNew()
  table$setBorderWidth(5)
  table$setColSpacings(2)
  table$setRowSpacings(5)
  table$attach(combo,            0, 4, 0, 1, 5, "shrink")
  table$attach(gtkLabelNew("X"), 0, 1, 1, 2, "shrink", "shrink")
  table$attach(button1,          1, 2, 1, 2, 5, "shrink")
  table$attach(gtkLabelNew("Y"), 2, 3, 1, 2, "shrink", "shrink")
  table$attach(button2,          3, 4, 1, 2, 5, "shrink")
  
  check  <- gtkToggleButtonNew()
  check$setImage(gtkImageNewFromStock(GTK_STOCK_GO_DOWN, GtkIconSize["button"]))
  check["tooltip-markup"] <- gettext("Inherit from the current theme")
  
  table$attach(check, 4, 5, 0, 1, "shrink", "shrink")  
  
  gSignalConnect(combo, "changed", function(obj){
    text <- localize(combo$getActiveText())
    if (text == "specify") {
      button1$setSensitive(TRUE)
      button2$setSensitive(TRUE)      
    } else {
      button1$setSensitive(FALSE)
      button2$setSensitive(FALSE)
    }
  })
  
  gSignalConnect(check, "toggled", function(obj){
    inherit <- obj$getActive()
    combo$setSensitive( !inherit )    
    button1$setSensitive( !inherit )
    button2$setSensitive( !inherit )
    
    text <- localize(combo$getActiveText())
    if (!inherit && text == "specify") {
      button1$setSensitive(TRUE)
      button2$setSensitive(TRUE)      
    } else {
      button1$setSensitive(FALSE)
      button2$setSensitive(FALSE)
    }
  })
  
  check$setActive(inherit)
  
  obj <- gtkFrameNew(name)
  obj$setShadowType(GtkShadowType["etched-in"])
  obj$add(table)
  obj$showAll()
  
  obj$setData("inherit", inherit)
  obj$setData("check", check)
  obj$setData("name", name)
  obj$setData("combo", combo)
  obj$setData("button1", button1)
  obj$setData("button2", button2)
  
  gSignalConnect(button1, "value-changed", function(...){
    text <- localize(combo$getActiveText())
    if (text == "specify") {
      values <- c(button1$getValue(), button2$getValue())
      script <- sprintf("%s = %s", name, deparse(values))
    } else {
      script <- sprintf("%s = %s", name, deparse(text))
    }
    obj$setData("script", script)    
  })
  gSignalConnect(button2, "value-changed", function(...){
    text <- localize(combo$getActiveText())
    if (text == "specify") {
      values <- c(button1$getValue(), button2$getValue())
      script <- sprintf("%s = %s", name, deparse(values))
    } else {
      script <- sprintf("%s = %s", name, deparse(text))
    }
    obj$setData("script", script)    
  })
  gSignalConnect(combo, "changed", function(...){
    text <- localize(combo$getActiveText())
    if (text == "specify") {
      values <- c(button1$getValue(), button2$getValue())
      script <- sprintf("%s = %s", name, deparse(values))
    } else {
      script <- sprintf("%s = %s", name, deparse(text))
    }
    obj$setData("script", script)    
  })
    
  text <- localize(combo$getActiveText())
  script <- NULL
  if (text == "specify") {
    values <- c(button1$getValue(), button2$getValue())
    script <- sprintf("%s = %s", name, deparse(values))
  } else {
    script <- sprintf("%s = %s", name, deparse(text))
  }
  obj$setData("script", script)    
  
  class(obj) <- c("GtkPlotPosition", "GktPlotThemeWidgets", class(obj))
  return(obj)
}

gtkPlotPositionSetValue <- function(object, value){
  check   <- object$getData("check")
  name    <- object$getData("name")
  combo   <- object$getData("combo")
  button1 <- object$getData("button1")
  button2 <- object$getData("button2")
  
  labels <- NULL
  if (name == "legend.position") {
    labels <- c("top", "right", "bottom", "left", "specify")
  } else {
    labels <- c("center", "specify")
  }
  
  if (is.null(value)) {
    check$setActive(TRUE)          
  } else if (is.character(value)){
    combo$setActive( which(labels == value) - 1)
    check$setActive(FALSE)
  } else {
    combo$setActive( which(labels == "specify") - 1)
    button1$setValue(value[1])
    button2$setValue(value[2])
    check$setActive(FALSE)
  }
}

gtkPlotDirectionNew <- function(name, direction=NULL){
  inherit <- NULL
  if (is.null(direction)) {
    direction <- "horizontal"
    inherit <- TRUE
  } else {
    inherit <- FALSE
  }
  combo <- gtkComboBoxEntryNewText()
  combo["width-request"] <- 80
  combo$appendText("horizontal")
  combo$appendText("vertical")
  combo$setActive(which(c("horizontal", "vertical") == direction) - 1)
  
  check  <- gtkToggleButtonNew()
  check$setImage(gtkImageNewFromStock(GTK_STOCK_GO_DOWN, GtkIconSize["button"]))
  check["tooltip-markup"] <- gettext("Inherit from the current theme")
  
  hbox <- gtkHBoxNew(spacing=2)
  hbox$packStart(combo)
  hbox$packStart(check, expand=FALSE)
  
  gSignalConnect(check, "toggled", function(obj){
    inherit <- obj$getActive()
    combo$setSensitive( !inherit )
  })
  
  check$setActive(inherit)
  
  vbox <- gtkVBoxNew()
  vbox$packStart(hbox, expand=FALSE)
  vbox$setBorderWidth(5)

  obj <- gtkFrameNew(name)
  obj$setShadowType(GtkShadowType["etched-in"])
  obj$add(vbox)
  obj$showAll()
  
  obj$setData("inherit", inherit)
  obj$setData("check", check)
  obj$setData("combo", combo)
  
  gSignalConnect(combo, "changed", function(...){
    text <- localize(combo$getActiveText())
    script <- sprintf("%s = %s", name, deparse(text))
    obj$setData("script", script)    
  })
  
  text <- localize(combo$getActiveText())
  script <- sprintf("%s = %s", name, deparse(text))
  obj$setData("script", script)    
  
  class(obj) <- c("GtkPlotDirection", "GktPlotThemeWidgets", class(obj))
  return(obj)
}

gtkPlotDirectionSetValue <- function(object, value){
  check  <- object$getData("check")
  combo  <- object$getData("combo")
  if (is.null(value)) {
    check$setActive(TRUE)          
  } else {
    combo$setActive(which(c("horizontal", "vertical") == value) - 1)
    check$setActive(FALSE)
  }
}


gktPlotThemeWidgetsReset <- function(object) {
  check <- object$getData("check")
  check$setActive(TRUE)
}

gktPlotThemeWidgetsGetScript <- function(object) {
  check <- object$getData("check")
  if (check$getActive()) {
    return(NULL)
  } else {
    return(object$getData("script"))
  }
}

gtkActionSetIconFromFile <- function(widget, path, filename) {
  image <- gFileIconNew(gFileNewForPath(file.path(path, filename)))
  widget$setGicon(image)
}

spinStart <- function() {
  spinner <- rzTools$getSpinner()
  if (is.null(spinner)) {
    return()
  } else {
    spinner$getToplevel()$setSensitive(FALSE)
    spinner$start()
    spinner$show()
  }
}
spinStop  <- function() {
  spinner <- rzTools$getSpinner()
  if (is.null(spinner)) {
    return()
  } else {
    spinner$hide()
    spinner$stop() 
    spinner$getToplevel()$setSensitive(TRUE)
  }
}

rzAddData <- function(data.set, name=NULL){
  name <- ifelse(is.null(name),
                 sprintf("%s [from Global Environment]", as.character(match.call()[2])), name)
  rzTools$addData(data.set, name)
}

rzReloadData <- function(data.set.name = NULL, ask = TRUE){
  rzTools$reloadData(data.set.name, ask = ask)
}

rzAddItem <- function(item, name = as.character(substitute(item)), data.set.name = NULL, description = name,
                      measurement = c("auto", "nominal", "ordinal", "interval", "ratio"),
                      overwrite = FALSE, ask = FALSE){
  rzTools$addItem(item = item, name = name, data.set.name, description = description,
                  measurement = measurement, overwrite = overwrite, ask = ask)
}


saveSession <- function(file=NULL){
  data.collection.obj$saveSession(file)
}

loadSession <- function(file=NULL){
  data.collection.obj$loadSession(file)
}

checkConfDir <- function(){
  confDir <- normalizePath(rzConfPath(), winslash="/", mustWork=FALSE)
  if (file.exists(confDir)) {
    return(TRUE)
  } else {
    response <- rzTools$runDialog(gettextf(
      'Create the directory below to save the config file and/or the session.\n\t%s\nIs it OK?\n\nOr you can change the directory by\n\toptions(RzConfPath="/path/to/dir")',
      confDir
      ), type="question")
    if(response == GtkResponseType["ok"]) {
      dir.create(confDir, recursive=TRUE)
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
}

rzConfPath <- function() getOption("RzConfPath", path.expand("~/Rz"))
