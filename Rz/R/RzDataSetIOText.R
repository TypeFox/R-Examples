datasetiotext <- 
  setRefClass("RzDataSetIoText",
              fields = c("widget", "win", "filepath", "data"),
              methods = list(
                initialize  = function(...) {
                  initFields(...)
                  data <<- NULL
                  if(class(win)[[1]] == "uninitializedField") {
                    win <<- NULL
                  }
                  encode.label <- gtkLabelNew(gettext("Character code"))
                  encode.sep1 <- gtkHSeparatorNew()
                  encode.lab.hbox <- gtkHBoxNew(spacing=4)
                  encode.lab.hbox$packStart(encode.label, expand=FALSE)
                  encode.lab.hbox$packStart(encode.sep1, expand=TRUE, fill=TRUE)
                  encode1 <- gtkComboBoxNewText()
                  for(i in iconvlist()) encode1$appendText(i)
                  index <- which(localeToCharset()[1]==iconvlist()) - 1
                  if(length(index)==0) index <- -1
                  encode1$setActive(index)
                  encode.hbox <- gtkHBoxNew(spacing=5)
                  encode.hbox$packStart(encode1, expand=FALSE)                  
                  
                  header.label <- gtkLabelNew(gettext("Header"))
                  header.sep1 <- gtkHSeparatorNew()
                  header.lab.hbox <- gtkHBoxNew(spacing=4)
                  header.lab.hbox$packStart(header.label, expand=FALSE)
                  header.lab.hbox$packStart(header.sep1, expand=TRUE, fill=TRUE)
                  header1 <- gtkRadioButtonNewWithLabel(NULL, "TRUE")
                  header1$setData("value", TRUE)
                  header2 <- gtkRadioButtonNewWithLabelFromWidget(header1, "FALSE")
                  header2$setData("value", FALSE)
                  header.hbox <- gtkHBoxNew(spacing=5)
                  header.hbox$packStart(header1, expand=FALSE)
                  header.hbox$packStart(header2, expand=FALSE)
                  
                  sepchar.label <- gtkLabelNew(gettext("Separator"))
                  sepchar.sep1 <- gtkHSeparatorNew()
                  sepchar.lab.hbox <- gtkHBoxNew(spacing=4)
                  sepchar.lab.hbox$packStart(sepchar.label, expand=FALSE)
                  sepchar.lab.hbox$packStart(sepchar.sep1, expand=TRUE, fill=TRUE)                  
                  sepchar1 <- gtkRadioButtonNewWithLabel(NULL, gettext("Comma"))
                  sepchar1$setData("value", ",")
                  sepchar2 <- gtkRadioButtonNewWithLabelFromWidget(sepchar1, gettext("Whitespace"))
                  sepchar2$setData("value", " ")
                  sepchar3 <- gtkRadioButtonNewWithLabelFromWidget(sepchar1, gettext("Tab"))
                  sepchar3$setData("value", "\t")
                  sepchar4 <- gtkRadioButtonNewWithLabelFromWidget(sepchar1, gettext("Semicolon"))
                  sepchar4$setData("value", ";")
                  sepchar5 <- gtkRadioButtonNewWithLabelFromWidget(sepchar1, gettext("Other:"))
                  sepchar5$setData("value", "")
                  sepchar.entry <- gtkEntryNew()
                  sepchar.entry$setWidthChars(4)
                  sepchar.hbox <- gtkHBoxNew(spacing=5)
                  sepchar.hbox$packStart(sepchar5, expand=FALSE)
                  sepchar.hbox$packStart(sepchar.entry, expand=FALSE)
                  
                  missing.label <- gtkLabelNew(gettext("Missing"))
                  missing.sep1 <- gtkHSeparatorNew()
                  missing.lab.hbox <- gtkHBoxNew(spacing=4)
                  missing.lab.hbox$packStart(missing.label, expand=FALSE)
                  missing.lab.hbox$packStart(missing.sep1, expand=TRUE, fill=TRUE)
                  missing1 <- gtkRadioButtonNewWithLabel(NULL, gettext("Whitespace"))
                  missing1$setData("value", " ")
                  missing2 <- gtkRadioButtonNewWithLabelFromWidget(missing1, gettext("Blank"))
                  missing2$setData("value", "")
                  missing3 <- gtkRadioButtonNewWithLabelFromWidget(missing1, gettext("Period"))
                  missing3$setData("value", ".")
                  missing4 <- gtkRadioButtonNewWithLabelFromWidget(missing1, gettext("\"NA\""))
                  missing4$setData("value", "NA")
                  missing5 <- gtkRadioButtonNewWithLabelFromWidget(missing1, gettext("Other:"))
                  missing5$setData("value", "")
                  missing.entry <- gtkEntryNew()
                  missing.entry$setWidthChars(4)
                  missing.hbox <- gtkHBoxNew(spacing=5)
                  missing.hbox$packStart(missing5, expand=FALSE)
                  missing.hbox$packStart(missing.entry, expand=FALSE)
                  
                  decimal.label <- gtkLabelNew(gettext("Decimal"))
                  decimal.sep1 <- gtkHSeparatorNew()
                  decimal.lab.hbox <- gtkHBoxNew(spacing=4)
                  decimal.lab.hbox$packStart(decimal.label, expand=FALSE)
                  decimal.lab.hbox$packStart(decimal.sep1, expand=TRUE, fill=TRUE)
                  decimal1 <- gtkRadioButtonNewWithLabel(NULL, gettext("Period"))
                  decimal1$setData("value", ".")
                  decimal2 <- gtkRadioButtonNewWithLabelFromWidget(decimal1, gettext("Comma"))
                  decimal2$setData("value", ",")
                  decimal.hbox <- gtkHBoxNew(spacing=5)
                  decimal.hbox$packStart(decimal1, expand=FALSE)
                  decimal.hbox$packStart(decimal2, expand=FALSE)
                  
                  quote.label <- gtkLabelNew(gettext("Quote"))
                  quote.sep1 <- gtkHSeparatorNew()
                  quote.lab.hbox <- gtkHBoxNew(spacing=4)
                  quote.lab.hbox$packStart(quote.label, expand=FALSE)
                  quote.lab.hbox$packStart(quote.sep1, expand=TRUE, fill=TRUE)                  
                  quote1 <- gtkRadioButtonNewWithLabel(NULL, gettext("Double quote (\")"))
                  quote1$setData("value", "\"")
                  quote2 <- gtkRadioButtonNewWithLabelFromWidget(quote1, gettext("Single quote (')"))
                  quote2$setData("value", "'")
                  quote3 <- gtkRadioButtonNewWithLabelFromWidget(quote1, gettext("None"))
                  quote3$setData("value", "")
                  
                  vbox.options <- gtkVBoxNew(homogeneous=FALSE,spacing=2)
                  vbox.options$setBorderWidth(5)
                  vbox.options$packStart(encode.lab.hbox, expand=FALSE, padding=3)
                  vbox.options$packStart(encode.hbox, expand=FALSE)
                  vbox.options$packStart(header.lab.hbox, expand=FALSE, padding=3)
                  vbox.options$packStart(header.hbox, expand=FALSE)
                  vbox.options$packStart(sepchar.lab.hbox, expand=FALSE, padding=3)
                  vbox.options$packStart(sepchar1, expand=FALSE)
                  vbox.options$packStart(sepchar2, expand=FALSE)
                  vbox.options$packStart(sepchar3, expand=FALSE)
                  vbox.options$packStart(sepchar4, expand=FALSE)
                  vbox.options$packStart(sepchar.hbox, expand=FALSE)
                  vbox.options$packStart(missing.lab.hbox, expand=FALSE, padding=3)
                  vbox.options$packStart(missing1, expand=FALSE)
                  vbox.options$packStart(missing2, expand=FALSE)
                  vbox.options$packStart(missing3, expand=FALSE)
                  vbox.options$packStart(missing4, expand=FALSE)
                  vbox.options$packStart(missing.hbox, expand=FALSE)
                  vbox.options$packStart(decimal.lab.hbox, expand=FALSE, padding=3)
                  vbox.options$packStart(decimal.hbox, expand=FALSE)
                  vbox.options$packStart(quote.lab.hbox, expand=FALSE, padding=3)
                  vbox.options$packStart(quote1, expand=FALSE)
                  vbox.options$packStart(quote2, expand=FALSE)
                  vbox.options$packStart(quote3, expand=FALSE)
                  
                  textview <- gtkTextViewNewWithBuffer()
                  textview$modifyFont(pangoFontDescriptionFromString(rzSettings$getMonospaceFont()))
                  sw1 <- gtkScrolledWindowNew()
                  sw1$setPolicy(GtkPolicyType["automatic"], GtkPolicyType["automatic"])
                  sw1$setShadowType(GtkShadowType["in"])
                  sw1$add(textview)
                  textview.label <- gtkLabelNew(gettext("Input File (first 20 Lines)"))
                  textview.hbox <- gtkHBoxNew()
                  textview.hbox$packStart(textview.label, expand=FALSE)
                  vbox.textview <- gtkVBoxNew(spacing=2)
                  vbox.textview$packStart(textview.hbox, expand=FALSE)
                  vbox.textview$packStart(sw1, expand=TRUE)
                  
                  treeview <- gtkTreeViewNew()
                  sw2 <- gtkScrolledWindowNew()
                  sw2$setPolicy(GtkPolicyType["automatic"], GtkPolicyType["automatic"])
                  sw2$setShadowType(GtkShadowType["in"])
                  sw2$add(treeview)
                  treeview.label <- gtkLabelNew(gettext("Preview (first 20 rows)"))
                  treeview.hbox <- gtkHBoxNew()
                  treeview.hbox$packStart(treeview.label, expand=FALSE)
                  vbox.treeview <- gtkVBoxNew(spacing=2)
                  vbox.treeview$packStart(treeview.hbox, expand=FALSE)
                  vbox.treeview$packStart(sw2, expand=TRUE)
                  
                  paned1 <- gtkVPanedNew()
                  paned1$pack1(vbox.textview)
                  paned1$pack2(vbox.treeview)
                  paned1$setPosition(240)
                  paned1$setBorderWidth(5)
                  
                  paned2 <- gtkHPanedNew()
                  paned2$pack1(vbox.options)
                  paned2$pack2(paned1)
                  
                  widget <<- gtkDialogNewWithButtons(title="Options",
                                                     parent=win,
                                                     flags=5,
                                                     "gtk-ok", GtkResponseType["accept"],
                                                     "gtk-cancel", GtkResponseType["reject"],
                                                     show=FALSE)
                  widget$setSizeRequest(600, -1)
                  widget[["vbox"]]$add(paned2)
                  
                  reload <- function(...) {
                    encode <- localize(encode1$getActiveText())
                    con <- file(filepath, encoding=encode)
                    input <- paste(readLines(con, n=20), collapse="\n")
                    preview <- read.table(con, nrows=20,
                                          header=getActiveData(header1),
                                          sep=getActiveData(sepchar1),
                                          na.strings=getActiveData(missing1),
                                          quote=getActiveData(quote1),
                                          dec=getActiveData(decimal1))
                    textview$getBuffer()$setText(input)
                    lapply(treeview$getColumns(), treeview$removeColumn)
                    treeview$setModel(rGtkDataFrame(preview))
                    renderer <- gtkCellRendererTextNew()
                    for (i in seq_len(ncol(preview))) {
                      column <- gtkTreeViewColumnNewWithAttributes(
                        title=colnames(preview)[i],
                        cell=renderer, text=i-1)
                      treeview$appendColumn(column)
                    }
                  }
                  
                  gSignalConnect(encode1, "changed", reload)
                  gSignalConnect(header1, "toggled", reload)
                  gSignalConnect(header2, "toggled", reload)
                  gSignalConnect(sepchar1, "toggled", reload)
                  gSignalConnect(sepchar2, "toggled", reload)
                  gSignalConnect(sepchar3, "toggled", reload)
                  gSignalConnect(sepchar4, "toggled", reload)
                  gSignalConnect(sepchar5, "toggled", reload)
                  gSignalConnect(missing1, "toggled", reload)
                  gSignalConnect(missing2, "toggled", reload)
                  gSignalConnect(missing3, "toggled", reload)
                  gSignalConnect(missing4, "toggled", reload)
                  gSignalConnect(missing5, "toggled", reload)
                  gSignalConnect(decimal1, "toggled", reload)
                  gSignalConnect(decimal2, "toggled", reload)
                  gSignalConnect(quote1, "toggled", reload)
                  gSignalConnect(quote2, "toggled", reload)
                  gSignalConnect(quote3, "toggled", reload)
                  gSignalConnect(sepchar.entry, "changed", function(...){
                    text <- localize(sepchar.entry$getText())
                    sepchar5$setData("value", text)
                    sepchar5$setActive(TRUE)
                  })
                  gSignalConnect(missing.entry, "changed", function(...){
                    text <- localize(missing.entry$getText())
                    missing5$setData("value", text)
                    missing5$setActive(TRUE)
                  })
                  gSignalConnect(widget, "show", reload)
                  gSignalConnect(widget, "delete-event", function(...){
                    widget$hide()
                    return(TRUE)
                  })
                  gSignalConnect(widget, "response", function(widget, response){
                    widget$hide()
                    if (response == GtkResponseType["accept"]) {
                      spinStart()
                      data <<- read.table(filepath,
                                          header=getActiveData(header1),
                                          sep=getActiveData(sepchar1),
                                          na.strings=getActiveData(missing1),
                                          quote=getActiveData(quote1),
                                          dec=getActiveData(decimal1),
                                          fileEncoding=localize(encode1$getActiveText()))
                    } else {
                      data <<- NULL
                    }
                  })
                  
                },
                
                getActiveData = function(radio) {
                  group <- radio$getGroup()
                  result <- sapply(group,
                                   function(x) if(x$getActive()) x$getData("value"))
                  result <- unlist(result)
                  return(result)
                },
                
                run = function(filepath) {
                  dialog <- gtkFileChooserDialogNew(title=gettext("Select text file"),
                                                    parent=win,
                                                    action=GtkFileChooserAction["open"],
                                                    "gtk-open", GtkResponseType["accept"],
                                                    "gtk-cancel", GtkResponseType["cancel"],
                                                    show=FALSE)
                  response <- dialog$run()
                  dialog$hide()
                  if (response != GtkResponseType["accept"]) {
                    return(NULL)
                  } else {
                    filepath <<- localize(dialog$getFilename())
                    widget$run()
                    return(data)
                  }
                  
                }
                
                )
  )
datasetiotext$accessors("filepath")
