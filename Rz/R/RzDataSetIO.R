dataSetIO <-
setRefClass("RzDataSetIO",
  methods = list(
    initialize            = function(...) {
      initFields(...)
    },

    save = function(win, data){
      if(is.null(data)) return()

      dialog <- gtkFileChooserDialogNew(title=gettext("Save Data"), parent=win,
                                        action=GtkFileChooserAction["save"],
                                        "gtk-save", GtkResponseType["accept"],
                                        "gtk-cancel", GtkResponseType["cancel"], 
                                        show=FALSE)
      file.types <- c(rzd    = gettext(gettext("Rz Data File (*.rzd)")),
                      spss   = gettext("SPSS Syntax and CSV (*.*)"),
                      stata  = gettext("Stata Do File and CSV (*.*)"),
                      stata2 = gettext("Stata Data File (*.dta)"),
                      csv    = gettext("Comma-Separated Text (*.csv)"),
                      tsv    = gettext("Tab-Separated Text (*.csv)"))
      file.type.list <- list(rzd     = list(name=file.types[["rzd"]]   , pattern="*.rzd"),
                             spss    = list(name=file.types[["spss"]]  , pattern="*.*"),
                             stata   = list(name=file.types[["stata"]] , pattern="*.*"),
                             stata2  = list(name=file.types[["stata2"]], pattern="*.dta"),
                             csv     = list(name=file.types[["csv"]]   , pattern="*.csv"),
                             tsv     = list(name=file.types[["tsv"]]   , pattern="*.csv")
                            )
      for (i in seq_along(file.type.list)) {
        filter <- gtkFileFilterNew()
        filter$setName(file.type.list[[i]]$name)
        for ( j in file.type.list[[i]]$pattern){ filter$addPattern(j) }
        dialog$addFilter(filter)
      }

      # Encoding Selector
      label <- gtkLabel(gettext("Encoding"))
      combo <- gtkComboBoxNewText()
      for(i in iconvlist()) combo$appendText(i)
      index <- which(localeToCharset()[1]==iconvlist()) - 1
      if(length(index)==0) index <- -1
      combo$setActive(index)

      # Header Selector
      header.button <- gtkCheckButtonNewWithLabel(gettext("Hearder(Only for Text File)"))
      header.button$setActive(TRUE)

      # NA Strings
      label.na <- gtkLabel(gettext("NA Strings"))
      entry.na <- gtkEntryNew()
      entry.na$setText("NA")

      hbox <- gtkHBoxNew()
      hbox$packEnd(combo, expand=FALSE, padding=5)
      hbox$packEnd(label, expand=FALSE)
      hbox2 <- gtkHBoxNew()
      hbox2$packEnd(header.button, expand=FALSE, padding=5)
      hbox3 <- gtkHBoxNew()
      hbox3$packEnd(entry.na, expand=FALSE, padding=5)
      hbox3$packEnd(label.na, expand=FALSE)

      vbox <- dialog$getContentArea()
      vbox$packEnd(hbox3, expand=FALSE)
      vbox$packEnd(hbox2, expand=FALSE)
      vbox$packEnd(hbox, expand=FALSE)

      gSignalConnect(dialog, "response", function(dialog, respons.id){
        spinStart()
        if (respons.id == GtkResponseType["accept"]) {
          filename   <- localize(dialog$getFilename())
          filetype   <- localize(dialog$getFilter()$getName())
          encoding   <- localize(combo$getActiveText())
          na.strings <- localize(entry.na$getText())
          header     <- header.button$getActive()

          if (filetype == file.types[["rzd"]]) {
            filename <- ifelse(grepl("*.rzd", filename), filename, sprintf("%s.rzd", filename))
            if(fileCheck(filename, dialog)){
              dialog$hide()
              data$save(file=filename)
            } else { dialog$run() }
          } else if (filetype == file.types[["spss"]]) {
            filename1 <- sprintf("%s.csv", filename)
            filename2 <- sprintf("%s.sps", filename)
            if(fileCheck(filename1, dialog) && fileCheck(filename2, dialog)){
              dialog$hide()
              df <- data$getData.frame()
              varlabels <- data$getVariableLabels()
              write.spss(df, filename1, filename2, varlabels=varlabels)
            } else { dialog$run() }
          } else if (filetype == file.types[["stata"]]) {
            filename1 <- sprintf("%s.csv", filename)
            filename2 <- sprintf("%s.do", filename)
            if(fileCheck(filename1, dialog) && fileCheck(filename2, dialog)){
              dialog$hide()
              df <- data$getData.frame()
              varlabels <- data$getVariableLabels()
              write.stata(df, filename1, filename2, varlabels=varlabels)
            } else { dialog$run() }
          } else if (filetype == file.types[["stata2"]]) {
            filename <- ifelse(grepl("*.dta", filename), filename, sprintf("%s.dta", filename))
            if(fileCheck(filename, dialog)){
              dialog$hide()
              df <- data$getData.frame()
              write.dta(df, filename)
            } else { dialog$run() }
          } else if (filetype == file.types[["csv"]]) {
            filename <- ifelse(grepl("*.csv", filename), filename, sprintf("%s.csv", filename))
            if(fileCheck(filename, dialog)){
              dialog$hide()
              df <- data$getData.frame()
              write.table(df, filename, quote = TRUE, sep = ",",
                          na = na.strings, row.names = FALSE, col.names = header,
                          qmethod = "double", fileEncoding = encoding)
            } else { dialog$run() }
          } else if (filetype == file.types[["tsv"]]) {
            filename <- ifelse(grepl("*.csv", filename), filename, sprintf("%s.csv", filename))
            if(fileCheck(filename, dialog)){
              dialog$hide()
              df <- data$getData.frame()
              write.table(df, filename, quote = TRUE, sep = "\t",
                          na = na.strings, row.names = FALSE, col.names = header,
                          qmethod = "double", fileEncoding = encoding)
            } else { dialog$run() }
          }

        } else {
          dialog$hide()
        }
      })
      dialog$run()
      
    },
    
    open = function(win, type, info.bar){
      # Encoding Selector
      label <- gtkLabel(gettext("Encoding"))
      combo <- gtkComboBoxNewText()
      for(i in iconvlist()) combo$appendText(i)
      index <- which(localeToCharset()[1]==iconvlist()) - 1
      if(length(index)==0) index <- -1
      combo$setActive(index)
      hbox <- gtkHBoxNew()
      hbox$packEnd(combo, expand=FALSE, padding=5)
      hbox$packEnd(label, expand=FALSE)
      
      # File Filter
      file.types <- c(rzdata   = gettext("Rz Data File (*.rzd)"),
                      spss     = gettext("SPSS System File (*.sav)"),
                      spss.por = gettext("SPSS Portable File (*.por)"),
                      stata    = gettext("Stata Data File (*.dta)")
                      )
      file.type.list <- list(rzdata  = list(name=file.types[["rzdata"]]  , pattern="*.rzd"),
                             spss    = list(name=file.types[["spss"]]    , pattern="*.sav"),
                             spss.por= list(name=file.types[["spss.por"]], pattern="*.por"),
                             stata   = list(name=file.types[["stata"]]   , pattern="*.dta")
                            )

      # File Import Dialog
      dialog <- gtkFileChooserDialogFilteredNew(title=gettext("Open"), parent=win, file.type.list=file.type.list[type])

      vbox <- dialog$getContentArea()
      vbox$packEnd(hbox, expand=FALSE)
      
      file       <- dialog$activate()
      encoding   <- localize(combo$getActiveText())
      c.encoding <- localeToCharset()[1]
      data.set   <- NULL
      dialog$hide()
      
      spinStart()
      
      if (is.null(file)){

        return(NULL)
      } else {

        dir  <- dirname(file$filename)
        base <- basename(file$filename)
        
        if (file$filetype == file.types[["rzdata"]]){
        # open rzd
          tmp.env <- new.env()
          load(file=file$filename, envir=tmp.env)
          # for old version
          if(is.null(tmp.env$rzdata)){
            data.set = tmp.env$data@.xData$data.set
          } else {
            rzdata   <- tmp.env$rzdata
            data.set <- rzdata$data.set
          }
          
        } else if (file$filetype == file.types[["spss"]]) {
        # import sav
          importer <- spss.system.file(file$filename)
          data.set   <- as.data.set(importer)
          
        } else if (file$filetype == file.types[["spss.por"]]) {
        # import por
          importer <- spss.portable.file(file$filename)
          data.set   <- as.data.set(importer)
          
        } else if (file$filetype == file.types[["stata"]]) {
        # import dta
          importer <- Stata.file(file$filename)
          data.set   <- as.data.set(importer)
          
        }
        
        # transcoding
        if ( encoding != c.encoding ) {
          info.bar$setText(gettext("Changing encoding may takes several or more minutes. Please wait..."))
          info.bar$setMessageType(GtkMessageType["info"])
          info.bar$show()
          while(gtkEventsPending()) gtkMainIteration()
          names(data.set) <- iconv(names(data.set), from=encoding, to=c.encoding)
          gtkMainIterationDo(FALSE)
          var.labs   <- iconv(description(data.set), from=encoding, to=c.encoding)
          gtkMainIterationDo(FALSE)
          var.labs   <- as.list(var.labs)
          gtkMainIterationDo(FALSE)
          var.labs   <- lapply(var.labs, function(x) gsub("(^')|('$)", "", x))
          gtkMainIterationDo(FALSE)
          labels     <- lapply(data.set, labels)
          gtkMainIterationDo(FALSE)
          labels     <- lapply(labels, function(x){
            x <- iconv(as.character(x), from=encoding, to=c.encoding)
            gtkMainIterationDo(FALSE)
            if(length(x) == 0) return(NULL)
            else return(x)
          })
          for(i in seq_along(var.labs)){
            gtkMainIterationDo(FALSE)
            description(data.set[[i]]) <- var.labs[[i]]
            if( !is.null(labels[[i]]))
              labels(data.set[[i]])@.Data <- labels[[i]]
          }
          info.bar$hide()
        }
        colnames(data.set) <- make.names(colnames(data.set), unique=TRUE)
        data <- new("RzData",
                    file.path=file$filename, original.name=base,
                    data.set=data.set)
        return(data)
      }

    },
    
    importFromText = function(win){
      io <- new("RzDataSetIoText", win=win)
      df <- io$run()
      file.path <- io$getFilepath()
      if (is.null(df)) {
        return(NULL)
      } else {
        data.set   <- data.set(df)
        names(data.set) <- colnames(df)
        encoding <- localeToCharset()[1]
        colnames(data.set) <- make.names(colnames(data.set), unique=TRUE)
        data <- new("RzData",
                    file.path=file.path, original.name=basename(file.path),
                    data.set=data.set)
        return(data)        
      }
    },

    importFromGlobalEnv = function(win){
      combo <- gtkComboBoxNewText()
      dfnames <- unlist(eapply(.GlobalEnv, is.data.frame))
      dfnames <- names(dfnames)[dfnames]
      dsnames <- unlist(eapply(.GlobalEnv, is.data.set))
      dsnames <- names(dsnames)[dsnames]
      importlist <- c(dfnames, dsnames)
      sapply(importlist, combo$appendText)
      dialog <- gtkDialogNewWithButtons(title=gettext("Import from Grobal Environment"), parent=win, flags=c("modal", "destroy-with-parent"),
                                         "gtk-ok", GtkResponseType["ok"], 
                                         "gtk-cancel", GtkResponseType["cancel"],
                                         show=FALSE)
      dialog[["vbox"]]$packStart(combo, expand=FALSE)
      dialog$showAll()
      response <- dialog$run()
      import <- localize(combo$getActiveText())
      dialog$hide()

      spinStart()
      if (response==GtkResponseType["ok"] & length(import)!=0) {
        ds <- get(import, envir=.GlobalEnv)
        if (is.data.frame(ds)) {
          ds2 <- data.set(ds)
          names(ds2) <- colnames(ds)
          ds <- ds2
        }
        data <- new("RzData", file.path=NULL, data.set=ds,
                    original.name=gettextf("%s [from Global Environment]", import))
        return(data)
      } else {
        return(NULL)
      }

    }

))
