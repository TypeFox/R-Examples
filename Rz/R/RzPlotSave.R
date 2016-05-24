rzplot.save <- 
setRefClass("RzPlotSave",
  fields = c("button.save", "dialog", "file.types", "file.type.list",
             "entry1", "entry2", "combo.unit", "entry.scale", "entry.dpi"),
  methods = list(
    initialize  = function(...) {
      initFields(...)
      # save
      dialog <<- gtkFileChooserDialogNew(title=gettext("Save"), parent=rzTools$getWindow(),
                                         action=GtkFileChooserAction["save"],
                                         "gtk-save", GtkResponseType["accept"],
                                         "gtk-cancel", GtkResponseType["cancel"], 
                                         show=FALSE)
      file.types <<- c(png  = gettext("png (*.png)"),
                       pdf  = gettext("pdf (*.pdf)"),
                       eps  = gettext("eps/ps (*.eps, *.ps)"),
                       wmf  = gettext("wmf(windows only) (*.wmf)"),
                       svg  = gettext("svg (*.svg)"),
                       jpeg = gettext("jpeg (*.jpeg)"),
                       bmp  = gettext("bmp (*.bmp)"),
                       tiff = gettext("tiff (*.tiff)"),
                       tex  = gettext("pictex (*.tex)"))
      file.type.list <<- list(png  = list(name=file.types[["png"]] , pattern="*.png"),
                              pdf  = list(name=file.types[["pdf"]] , pattern="*.pdf"),
                              eps  = list(name=file.types[["eps"]] , pattern=c("*.eps", "*.ps")),
                              wmf  = list(name=file.types[["wmf"]] , pattern="*.wmf"),
                              svg  = list(name=file.types[["svg"]] , pattern="*.svg"),
                              jpeg = list(name=file.types[["jpeg"]], pattern="*.jpeg"),
                              bmp  = list(name=file.types[["bmp"]] , pattern="*.bmp"),
                              tiff = list(name=file.types[["tiff"]], pattern="*.tiff"),
                              tex  = list(name=file.types[["tex"]] , pattern="*.tex")
                              )
      for (i in seq_along(file.type.list)) {
        filter <- gtkFileFilterNew()
        filter$setName(file.type.list[[i]]$name)
        for ( j in file.type.list[[i]]$pattern){ filter$addPattern(j) }
        dialog$addFilter(filter)
      }
      dialogRun <- function(...) dialog$run()
      
      image       <-  gtkImageNewFromStock(GTK_STOCK_SAVE, GtkIconSize["button"])
      button.save <<- gtkButtonNew()
      button.save["sensitive"] <<- TRUE
      button.save$setImage(image)
      gSignalConnect(button.save, "clicked", dialogRun)
      
      label1 <- gtkLabelNew(gettext("width"))
      entry1 <<- gtkEntryNew()
      entry1$setText("(auto)")
      label2 <- gtkLabelNew(gettext("height"))
      entry2 <<- gtkEntryNew()
      entry2$setText("(auto)")
      label.unit <- gtkLabelNew(gettext("unit"))
      combo.unit <<- gtkComboBoxNewText()
      combo.unit$appendText("in")
      combo.unit$appendText("cm")
      combo.unit$appendText("mm")
      combo.unit$setActive(1)
      label.scale <- gtkLabelNew(gettext("scale"))
      entry.scale <<- gtkEntryNew()
      entry.scale$setText("1")
      label.dpi <- gtkLabelNew(gettext("dpi"))
      entry.dpi <<- gtkEntryNew()
      entry.dpi$setText("300")
      
      
      table  <- gtkTableNew(FALSE)
      table["border-width"] <- 5
      table$attach(label1, 0, 1, 0, 1, "shrink", "shrink", 0, 0)
      table$attachDefaults(entry1, 1, 2, 0, 1)
      table$attach(label2, 0, 1, 1, 2, "shrink", "shrink", 0, 0)
      table$attachDefaults(entry2, 1, 2, 1, 2)
      table$attach(label.unit, 0, 1, 2, 3, "shrink", "shrink", 0, 0)
      table$attachDefaults(combo.unit, 1, 2, 2, 3)
      table$attach(label.scale, 2, 3, 0, 1, "shrink", "shrink", 0, 0)
      table$attachDefaults(entry.scale, 3, 4, 0, 1)
      table$attach(label.dpi, 2, 3, 1, 2, "shrink", "shrink", 0, 0)
      table$attachDefaults(entry.dpi, 3, 4, 1, 2)
      table$setColSpacings(5)
      table$setRowSpacings(2)
      
      dialog$setExtraWidget(table)
    },
    
    onSave = function(response.id, p.current, data){
      if (response.id == GtkResponseType["accept"]) {
        dialog$hide()
        
        filename <- localize(dialog$getFilename())
        filetype <- localize(dialog$getFilter()$getName())
        
        match <- NULL
        for (i in seq_along(file.type.list)) {
          match <- c(match, file.type.list[[i]]$name==filetype)
        }
        ind <- which(match)
        ext <- paste(strsplit(file.type.list[[ind]]$pattern[1], "")[[1]][-1], collapse="")
        filename <- ifelse(grepl(paste(file.type.list[[ind]]$pattern, collapse="|"), filename), filename, sprintf(paste("%s", ext, sep=""), filename))
        
        if (file.exists(filename)){
          dialog2 <- gtkMessageDialogNew(dialog, "destroy-with-parent",
                                         GtkMessageType["question"],
                                         GtkButtonsType["ok-cancel"],
                                         gettext("This file already exists. Overwrite it?"))
          response <- dialog2$run()
          dialog2$hide()
          if (response!=GtkResponseType["ok"]){
            dialog$run()
            return()
          }
        }
        
        width  <- localize(entry1$getText())
        width  <- suppressWarnings(as.numeric(width))
        width  <- ifelse(is.na(width), par("din")[1], width)
        height <- localize(entry2$getText())
        height <- suppressWarnings(as.numeric(height))
        height <- ifelse(is.na(height), par("din")[2], height)
        scale  <- localize(entry.scale$getText())
        scale  <- suppressWarnings(as.numeric(scale))
        scale  <- ifelse(is.na(scale), 1, scale)
        dpi    <- localize(entry.dpi$getText())
        dpi    <- suppressWarnings(as.numeric(dpi))
        dpi    <- ifelse(is.na(dpi), 3000, dpi)
        units  <- localize(combo.unit$getActiveText())
        
        if(filetype==file.types[["eps"]]){
          p.current <- p.current + theme(text=element_text(family=""))
          ggsave(filename=filename, plot=p.current, width=width, height=height,
                 family=rzSettings$getPsFont(), units=units, dpi=dpi, scale=scale)
        } else if(filetype==file.types[["pdf"]]){
          p.current <- p.current + theme(text=element_text(family=""))
          ggsave(filename=filename, plot=p.current, width=width, height=height,
                 family=rzSettings$getPdfFont(), units=units, dpi=dpi, scale=scale)
        } else {
          ggsave(filename=filename, plot=p.current, width=width, height=height, units=units, dpi=dpi, scale=scale)
        }
      } else {
        dialog$hide()
      }
      
    }
  )
)
rzplot.save$accessors("button.save", "dialog")
