rzVVCore <-
setRefClass("RzVVCore",
  fields = c("widget", "vvTreeView", "data", "liststore", "model.selected", "filtered", "vv", "adj",
             "rt.index", "rtg.select", "rt.vars", "rt.var.labs", "rt.val.labs", "rp.msr", "rt.missing",
             "nominalpix", "ordinalpix", "intervalpix", "ratiopix", "rzVarToolBar"),
  methods = list(
    initialize            = function(...) {      
      initFields(...)
      if (class(data)[[1]] == "uninitializedField") return()
      
      rzVarToolBar <<- new("RzVarToolBar")
      toolbar <- rzVarToolBar$getToolBar()
      
      entry <- new("RzSearchEntry")
      
      hbox.vartoolbar <- gtkHBoxNew(spacing=2)
      hbox.vartoolbar$packStart(toolbar, expand=TRUE, fill=TRUE)
      hbox.vartoolbar$packEnd(entry$getEntry.search(), expand=FALSE)
      
      model.selected <<- gtkTreeModelFilterNew(liststore)
      model.selected$setVisibleColumn(column.definition["select"])
      filtered   <<- FALSE
      
      nominalpix    <<- gdkPixbufNewFromFile(file.path(path.package("Rz"), "images/oxygen/cat.png"     ))$retval
      ordinalpix    <<- gdkPixbufNewFromFile(file.path(path.package("Rz"), "images/oxygen/order.png"   ))$retval
      intervalpix   <<- gdkPixbufNewFromFile(file.path(path.package("Rz"), "images/oxygen/interval.png"))$retval
      ratiopix      <<- gdkPixbufNewFromFile(file.path(path.package("Rz"), "images/oxygen/ratio.png"   ))$retval

      vvTreeView <<- gtkTreeViewNewWithModel(liststore)
      vvTreeView$setSearchEntry(entry$getEntry.search())
      vvTreeView$setSearchEqualFunc(entry$searchFunc)
      
      sw   <- gtkScrolledWindowNew()
      sw["shadow-type"] <- GtkShadowType["none"]
      sw$add(vvTreeView)
      sw$setPolicy(GtkPolicyType["automatic"], GtkPolicyType["automatic"])
      if (class(adj)[[1]] == "uninitializedField") adj <<- sw$getVadjustment()
      
      vvTreeView["enable-grid-lines"] <<- GtkTreeViewGridLines["both"]
      vvTreeView["rules-hint"] <<- TRUE
      vvTreeView["has-tooltip"] <<- TRUE
      vvTreeView$modifyFont(pangoFontDescriptionFromString(rzSettings$getVariableViewFont()))
      
      rt.index    <<- gtkCellRendererText()
      rtg.select  <<- gtkCellRendererToggleNew()
      rt.vars     <<- gtkCellRendererText()
      rt.var.labs <<- gtkCellRendererText()
      rp.msr      <<- gtkCellRendererPixbuf()
      rt.val.labs <<- gtkCellRendererText()
      rt.missing  <<- gtkCellRendererText()
      color       <- rt.index["cell-background-gdk"]
      color$red   <- 45000L
      color$green <- 45000L
      color$blue  <- 45000L
      rt.index["cell-background-gdk"] <<- color
      rt.index["xalign"] <<- 0.5
      rt.vars    ["editable"] <<- TRUE
      rt.var.labs["editable"] <<- TRUE
      rt.missing ["editable"] <<- TRUE
      
      columns <- list(
        index   = gtkTreeViewColumnNewWithAttributes(""                     , rt.index   , "text"=column.definition[["index"]]   ),
        select  = gtkTreeViewColumnNewWithAttributes(""                     , rtg.select , "active"=column.definition[["select"]]),
        msr     = gtkTreeViewColumnNewWithAttributes(gettext("Measurement") , rp.msr     , "pixbuf"=column.definition[["msr.image"]]),
        vars    = gtkTreeViewColumnNewWithAttributes(gettext("Names")       , rt.vars    , "text"=column.definition[["vars"]]    ),
        labs    = gtkTreeViewColumnNewWithAttributes(gettext("Labels")      , rt.var.labs, "text"=column.definition[["var.labs"]]),
        val.labs= gtkTreeViewColumnNewWithAttributes(gettext("Value Labels"), rt.val.labs, "text"=column.definition[["val.labs"]]),
        missing = gtkTreeViewColumnNewWithAttributes(gettext("Missing")     , rt.missing , "text"=column.definition[["missing"]] )
      )
      lapply(columns, gtkTreeViewColumnSetSizing   , "fixed")
      lapply(columns, gtkTreeViewColumnSetResizable, TRUE)
      lapply(columns, gtkTreeViewColumnSetSpacing  , 1)
      
      columns$index$setData("attr", c(title="index"))
      columns$index$setMinWidth(30)
      columns$index$setSizing("automatic")
      columns$index$setResizable(FALSE)
      
      columns$select$setData("attr", c(title="select"))
      columns$select$setSizing("automatic")
      columns$select$setResizable(FALSE)
      
      columns$vars$setData("attr", c(title="vars"))
      columns$vars$setFixedWidth(50)
      
      columns$labs$setData("attr", c(title="labs"))
      columns$labs$setFixedWidth(250)
      
      columns$val.labs$setData("attr", c(title="val.labs"))
      columns$val.labs$setFixedWidth(100)
      
      columns$msr$setData("attr", c(title="msr"))
      columns$msr$setFixedWidth(30)
      columns$msr$setMinWidth(30)
      
      columns$missing$setData("attr", c(title="missing"))
      columns$missing$setSizing("automatic")
      columns$missing$setResizable(FALSE)
      
      lapply(columns, function(column) vvTreeView$appendColumn(column))
      
      widget <<- gtkVBoxNew()
      widget$packStart(hbox.vartoolbar, expand=FALSE)
      widget$packStart(sw, expand=TRUE, fill=TRUE)
      
      gSignalConnect(vvTreeView, "row-activated", .self$onRowActivated)
      gSignalConnect(vvTreeView, "query-tooltip", .self$onQueryTooltip)
      gSignalConnect(rtg.select , "toggled", .self$onCelltoggledSelect)
      gSignalConnect(rt.vars    , "edited", .self$onCellEditedVars)
      gSignalConnect(rt.var.labs, "edited", .self$onCellEditedVarLabs)
      gSignalConnect(rt.missing , "edited", .self$onCellEditedMissing)
      gSignalConnect(vvTreeView, "realize", function(...) sw$setVadjustment(adj))
      
      gSignalConnect(rzVarToolBar$getA.select(), "activate", .self$onSelectAll)
      gSignalConnect(rzVarToolBar$getA.unselect(), "activate", .self$onUnselect)
      rzVarToolBar$setSignalId("onlyselec",
                               gSignalConnect(rzVarToolBar$getA.onlyselec(), "toggled",vv$setFiltered))
      gSignalConnect(rzVarToolBar$getA.delete(),    "activate", .self$onDelete)
      gSignalConnect(rzVarToolBar$getA.duplicate(), "activate", .self$onDuplicate)
      
      
      popup <- new("RzPopup")
      gSignalConnect(sw, "button-release-event" , popup$popupmenu)
      gSignalConnect(popup$getA.recode(), "activate", .self$onRecode)
      
    },

    setFiltered = function(active){
      if (active) {
        vvTreeView$setModel(model.selected)
        filtered <<- TRUE
        rzVarToolBar$setFiltered(filtered)
      } else {
        vvTreeView$setModel(liststore)
        filtered <<- FALSE
        rzVarToolBar$setFiltered(filtered)
      }
    },
    
    msrPix = function(msr){
      pix <- NULL
      if(msr=="nominal"){
        pix <- nominalpix
      } else if(msr=="ordinal"){
        pix <- ordinalpix
      } else if(msr=="interval"){
        pix <- intervalpix
      } else if(msr=="ratio"){
        pix <- ratiopix
      }
      return(pix)
    },
    
    getValidIterFromString = function(path) {
      if(!filtered) {
        iter <- liststore$getIterFromString(path)$iter
        return(iter)
      } else {
        iter <- model.selected$getIterFromString(path)$iter
        iter <- model.selected$convertIterToChildIter(iter)$child.iter
        return(iter)
      }
    },
    
    setCell     = function(iter, col, new.value){
      if (liststore$iterIsValid(iter)) {
        liststore$set(iter, col, new.value)
        if(col==column.definition["msr"]){
          liststore$set(iter, column.definition["msr.image"], .self$msrPix(new.value))
        }
      } else {
        iter <- liststore$append()$iter
        liststore$set(iter, col, new.value) 
      }
      row <- unlist(liststore$get(iter, column.definition["index"]))
      vv$setSummary(row)
    },
    
    getCell     = function(path, col){
      if(!filtered) {
        iter <- liststore$getIterFromString(path)$iter
      } else {
        iter <- model.selected$getIterFromString(path)$iter
        iter <- model.selected$convertIterToChildIter(iter)$child.iter
      }
      value <- unlist(liststore$get(iter, col))
      return(value)
    },
    
    
    getSelected = function(){
      iter  <- vvTreeView$getSelection()$getSelected()$iter
      if(filtered) {
        iter  <- model.selected$convertIterToChildIter(iter)$child.iter
      }
      value <- liststore$get(iter, unlist(column.definition))        
      value <- value[-(column.definition["msr.image"]+1)]        
      value <- lapply(value, localize)
      value <- unlist(value)
      names(value) <- names(column.definition[-(column.definition["msr.image"]+1)])
      return(value)
    },
    
    getSelectedRows = function(){
      iter  <- model.selected$getIterFirst()
      inds  <- character(0)
      while(iter$retval){
        value <- model.selected$getValue(iter$iter, column.definition["index"])$value
        inds <- c(inds, value)
        iter$retval <- model.selected$iterNext(iter$iter)
      }
      inds <- sapply(inds, localize)
      inds <- as.numeric(inds)
      return(inds)
    },
    
    onCelltoggledSelect = function(renderer, path){
      active <- renderer$getActive()
      renderer$setActive(!active)
      active <- renderer$getActive()
      .self$setCell(getValidIterFromString(path), column.definition["select"], active)
      
    },
    
    onCellEditedVars    = function(renderer, path, new.text){
      txt     <- localize(new.text)
      txt           <- make.names(txt)
      row           <- as.numeric(getSelected()["index"])
      data.set.name <- data$getData.set.name()
      var.name      <- data$getVariableNames()[row]
      data.set      <- data$getData.set()
      names(data.set)[row] <- txt
      data$setData.set(data.set)
      data$constructVariable(row)
      data$linkDataFrame()
      .self$setCell(getValidIterFromString(path), column.definition["vars"], txt)
    },
    
    onCellEditedVarLabs = function(renderer, path, new.text){
      txt           <- localize(new.text)
      txt           <- sub("^([[:space:]]+)([^[:space:]]+)([[:space:]]+)$", "\\2", txt)
      row           <- as.numeric(getSelected()["index"])
      data.set.name <- data$getData.set.name()
      var.name      <- data$getVariableNames()[row]
      data.set      <- data$getData.set()
      description(data.set[[row]]) <- txt
      data$setData.set(data.set)      
      .self$setCell(getValidIterFromString(path), column.definition["var.labs"], txt)
    },
    
    onEditMsr     = function(renderer, path, new.text){
      data.set.name <- data$getData.set.name()
      selec      <- .self$getSelected()
      var.name   <- selec["vars"]
      var.lab    <- selec["var.labs"]
      data.set   <- data$getData.set()
      var        <- data.set[[var.name]]
      msr        <- measurement(var)
      row <- which(names(data.set)==var.name)
      
      dialog <- gtkWindowNew(show=FALSE)
      dialog$setTransientFor(vvTreeView$getToplevel())
      dialog$setModal(FALSE)
      dialog$setDecorated(FALSE)
      dialog["window-position"] <- GtkWindowPosition["mouse"]
      
      radio1 <- gtkRadioButtonNewWithLabel(label="nominal")
      radio2 <- gtkRadioButtonNewWithLabelFromWidget(group=radio1, label="ordinal")
      radio3 <- gtkRadioButtonNewWithLabelFromWidget(radio1, label="interval")
      radio4 <- gtkRadioButtonNewWithLabelFromWidget(radio1, label="ratio")
      image <- gtkImageNewFromPixbuf(nominalpix)
      radio1$setImage(image)
      image$show()
      image <- gtkImageNewFromPixbuf(ordinalpix)
      radio2$setImage(image)
      image$show()
      image <- gtkImageNewFromPixbuf(intervalpix)
      radio3$setImage(image)
      image$show()
      image <- gtkImageNewFromPixbuf(ratiopix)
      radio4$setImage(image)
      image$show()
      if(msr=="nominal"){
        radio1$setCanDefault(TRUE)
        dialog$setDefault(radio1)
        radio1$setActive(TRUE)
      } else if(msr=="ordinal"){
        radio2$setCanDefault(TRUE)
        dialog$setDefault(radio2)
        radio2$setActive(TRUE)
      } else if(msr=="interval"){
        radio3$setCanDefault(TRUE)
        dialog$setDefault(radio3)
        radio3$setActive(TRUE)
      } else if(msr=="ratio"){
        radio4$setCanDefault(TRUE)
        dialog$setDefault(radio4)
        radio4$setActive(TRUE)
      }
      radio1["draw-indicator"] <- FALSE
      radio2["draw-indicator"] <- FALSE
      radio3["draw-indicator"] <- FALSE
      radio4["draw-indicator"] <- FALSE
      
      onToggled <- function(button){
        if(button$getActive()){
          dialog$hide()
          msr <- localize(button["label"])
          measurement(data.set[[row]]) <- msr
          data$setData.set(data.set)
          data$constructVariable(row)
          data$linkDataFrame()
          cell.row <- as.character(row - 1)
          .self$setCell(liststore$getIterFromString(cell.row)$iter, column.definition["msr"], msr)
        }
      }
      gSignalConnect(radio1, "clicked", onToggled)
      gSignalConnect(radio2, "clicked", onToggled)
      gSignalConnect(radio3, "clicked", onToggled)
      gSignalConnect(radio4, "clicked", onToggled)
      
      vbox <- gtkVBoxNew(spacing=2)
      vbox$packStart(radio1, expand=FALSE)
      vbox$packStart(radio2, expand=FALSE)
      vbox$packStart(radio3, expand=FALSE)
      vbox$packStart(radio4, expand=FALSE)
      frame <- gtkFrameNew()
      frame$setShadowType(GtkShadowType["out"])
      frame$add(vbox)
      dialog$add(frame)
      dialog$showAll()
      
      gSignalConnect(dialog, "leave-notify-event", function(widget, event){
        size <- widget$getSize()
        cond <- c(event$x, event$y ,size$width - event$x, size$height - event$y)
        if (any(cond < 0))  widget$hide()
        return(FALSE)
      })
      
    },
    
    onCellEditedMissing = function(renderer, path, new.text){
      txt           <- localize(new.text)
      txt           <- sub("^([[:space:]]+)([^[:space:]]+)([[:space:]]+)$", "\\2", txt)
      row           <- as.numeric(getSelected()["index"])
      data.set.name <- data$getData.set.name()
      var.name      <- data$getVariableNames()[row]
      data.set      <- data$getData.set()
      result <- try(eval(parse(text=sprintf("c(%s)", txt))), silent=TRUE)
      if (nzchar(txt)&&!is.numeric(result)) return()
      else {
        if (length(result)==2&&grepl("range", txt)) {
          missing.values(data.set[[row]]) <- eval(parse(text=sprintf("list(%s)", txt))) 
        } else if (!nzchar(txt)) {
          missing.values(data.set[[row]]) <- NULL          
        } else {
          missing.values(data.set[[row]]) <- result
        }
      }
      data$setData.set(data.set)
      data$constructVariable(row)
      data$linkDataFrame()
      .self$setCell(getValidIterFromString(path), column.definition["missing"], txt)
    },
    
    onRecode            = function(action){
      onActivatedGTButton <- function(obj, list){
        var.name   <- list$var.name
        liststore2 <- list$liststore2
        textbuffer <- list$textbuffer
        result     <- liststore2$getIterFirst()
        data.tmp   <- NULL
        while(result$retval) {
          data.tmp <- c(data.tmp, lapply(liststore2$get(result$iter, 0, 1), localize))
          result$retval <- liststore2$iterNext(result$iter)
        }
        if(is.null(data.tmp)) return()
        data.tmp <- matrix(data.tmp, byrow=TRUE, ncol=2)
        text <- paste("\"", data.tmp[,2], "\" = ", data.tmp[,1], " <- ", data.tmp[,1], sep="", collapse=",\n")
        text <- sprintf("%s,\n%s", text, "otherwise = \"copy\"")
        textbuffer$setText(text)
      }
      onResponse <- function(dialog, response){
        if (response==GtkResponseType["ok"]) {
          iter <- textbuffer$getBounds()
          text <- localize(textbuffer$getText(iter$start, iter$end))
          if(!nzchar(text)){
            dialog2 <- gtkMessageDialogNew(dialog, "destroy-with-parent",
                                           GtkMessageType["error"],
                                           GtkButtonsType["close"],
                                           gettext("This recode syntax is invalid. Please enter a valid recode syntax."))
            dialog2$run()
            dialog2$hide()
            dialog$run()
            return()
          }
          text <- paste("         ", strsplit(text, "\n")[[1]],sep="", collapse="\n")
          new.var <- localize(entry3$getText())
          new.var.lab <- localize(entry.var.lab$getText())
          new.var <- sub("^([[:space:]]+)([^[:space:]]+)([[:space:]]+)$",
                         "\\2", new.var)
          invalid <- grepl("(^$)|(^[0-9]+)|([]\\[\\^$*?|(){}@!\"#$%&'*+,/:;<=>?~[:space:]-])",
                           new.var)
          if (invalid || length(new.var) == 0 ){
            dialog2 <- gtkMessageDialogNew(dialog, "destroy-with-parent",
                                           GtkMessageType["error"],
                                           GtkButtonsType["close"],
                                           gettext("This variable name is invalid. Please enter a valid variable name."))
            dialog2$run()
            dialog2$hide()
            dialog$run()
            return()
          }
          if( any(new.var==data$getVariableNames()) ){
            dialog2 <- gtkMessageDialogNew(dialog, "destroy-with-parent",
                                           GtkMessageType["question"],
                                           GtkButtonsType["ok-cancel"],
                                           gettext("This Variable already exists. Overwrite it?"))
            response2 <- dialog2$run()
            dialog2$hide()
            if (response2!=GtkResponseType["ok"]){
              dialog$run()
              return()
            }
          }
          text <- sprintf("%s <-\n  recode(%s,\n%s)\n", new.var, var.name, text)
          data.set <- try(within(data.set, eval(parse(text=text))), silent=TRUE)
          if(!is.data.set(data.set)){
            dialog2 <- gtkMessageDialogNew(dialog, "destroy-with-parent",
                                           GtkMessageType["error"],
                                           GtkButtonsType["close"],
                                           gettext("This recode syntax is invalid. Please enter a valid recode syntax."))
            dialog2$run()
            dialog2$hide()
            dialog$run()
            return()
          }
          dialog$hide()
          description(data.set[[new.var]]) <- new.var.lab
          data$setData.set(data.set)
          
          var <- data.set[[new.var]]
          msr <- measurement(var)
          miss.val  <- missing.values(var)
          miss.val  <- ifelse(is.null(miss.val), "", paste(miss.val@filter, collapse=","))
          val.labs  <- labels(var)
          val.labs  <- ifelse(is.null(val.labs), "", paste(val.labs@values, " \"",val.labs@.Data, "\"", sep="", collapse=", "))
          row <- which(names(data.set)==new.var)
          cell.row <- as.character(row - 1)
          data$constructVariable(row)
          data$linkDataFrame()
          .self$setCell(liststore$getIterFromString(cell.row)$iter, column.definition["index"], as.character(row))
          .self$setCell(liststore$getIterFromString(cell.row)$iter, column.definition["select"], filtered)
          .self$setCell(liststore$getIterFromString(cell.row)$iter, column.definition["vars"], new.var)
          .self$setCell(liststore$getIterFromString(cell.row)$iter, column.definition["var.labs"], new.var.lab)
          .self$setCell(liststore$getIterFromString(cell.row)$iter, column.definition["msr"], msr)
          .self$setCell(liststore$getIterFromString(cell.row)$iter, column.definition["msr.image"], .self$msrPix(msr))
          .self$setCell(liststore$getIterFromString(cell.row)$iter, column.definition["val.labs"], val.labs)
          .self$setCell(liststore$getIterFromString(cell.row)$iter, column.definition["missing"] , miss.val)
          
        } else {
          dialog$hide()
        }
      }
      data.set.name <- data$getData.set.name()
      selec      <- .self$getSelected()
      var.name   <- selec["vars"]
      var.lab    <- selec["var.labs"]
      data.set   <- data$getData.set()
      var        <- data.set[[var.name]]
      msr        <- measurement(var)
      labels     <- labels(var)
      liststore2 <- gtkListStoreNew("character", "character")
      if(!is.null(labels)){
        labels <- rbind(labels@values, labels@.Data)
        apply(labels, 2,
              function(x) {
                iter <- liststore2$append()$iter
                liststore2$set(iter, 0, x[1], 1, x[2])
              })
      }
      tw <- gtkTreeViewNewWithModel(liststore2)
      
      rt.val <- gtkCellRendererText()
      rt.lab <- gtkCellRendererText()
      
      col.val <- gtkTreeViewColumnNewWithAttributes(gettext("Values"), rt.val, "text"=0)
      col.lab <- gtkTreeViewColumnNewWithAttributes(gettext("Value Labels"), rt.lab, "text"=1)
      tw$appendColumn(col.val)
      tw$appendColumn(col.lab)
      
      vbox1  <- gtkVBoxNew(spacing=2)
      hbox1  <- gtkHBoxNew(spacing=5)
      vbox2  <- gtkVBoxNew(spacing=2)
      hpaned <- gtkHPanedNew()
      table1  <- gtkTableNew(rows=3, columns=2)
      hpaned$setPosition(160)
      sw1 <- gtkScrolledWindowNew()
      sw1["shadow-type"] <- GtkShadowType["in"]
      
      sw2 <- gtkScrolledWindowNew()
      sw2["shadow-type"] <- GtkShadowType["in"]
      
      sw1$setPolicy("automatic", "automatic")
      sw2$setPolicy("automatic", "automatic")
      label1 <- gtkLabelNew(gettext("Existing Values"))
      entry1 <- gtkEntryNew()
      
      entry1["editable"] <- FALSE
      entry1$setText(paste(sort(unique(na.omit(as.numeric(var)))), collapse=","))
      label2 <- gtkLabelNew(gettext("Source"))
      label3 <- gtkLabelNew(gettext("Destination"))
      label.var.lab <- gtkLabelNew(gettext("Variable Label"))
      label2$setAlignment(0,0.5)
      label3$setAlignment(0,0.5)
      label.var.lab$setAlignment(0,0.5)
      entry2 <- gtkEntryNew()
      entry2$setText(var.name)
      entry2["editable"] <- FALSE
      #      entry2["sensitive"] <- FALSE
      entry3 <- gtkEntryNew()
      entry3$setText(var.name)
      entry.var.lab <- gtkEntryNew()
      entry.var.lab$setText(var.lab)
      entry2$setWidthChars(10)
      entry3$setWidthChars(10)
      entry.var.lab$setWidthChars(10)
      textview <- gtkTextViewNew()
      textview$modifyFont(pangoFontDescriptionFromString(rzSettings$getMonospaceFont()))
      
      textview$setLeftMargin(5)
      textview$setRightMargin(5)
      
      textbuffer <- textview$getBuffer()
      gtbutton   <- gtkButtonNewWithLabel(gettext("Generate Template"))
      
      vbox1$packStart(hbox1, expand=FALSE)
      vbox1$packStart(hpaned)
      
      hbox1$packStart(label1, expand=FALSE)
      hbox1$packStart(entry1)
      
      sw1$add(tw)
      hpaned$add1(sw1)
      hpaned$add2(vbox2)
      
      sw2$add(textview)
      vbox2$packStart(table1, expand=FALSE)
      vbox2$packStart(gtbutton, expand=FALSE, fill=FALSE)
      vbox2$packStart(sw2)
      
      table1$attach(label2, 0, 1, 0, 1, xoptions=GtkAttachOptions["fill"])
      table1$attach(entry2, 1, 2, 0, 1, xpadding = 1, ypadding = 1)
      table1$attach(label3, 0, 1, 1, 2, xoptions=GtkAttachOptions["fill"])
      table1$attach(entry3, 1, 2, 1, 2, xpadding = 1, ypadding = 1)
      table1$attach(label.var.lab, 0, 1, 2, 3, xoptions=GtkAttachOptions["fill"])
      table1$attach(entry.var.lab, 1, 2, 2, 3, xpadding = 1, ypadding = 1)
      table1$setColSpacings(5)
      table1$setRowSpacings(2)
      
      gSignalConnect(gtbutton, "clicked", onActivatedGTButton, list(var.name=var.name, liststore2=liststore2, textbuffer=textbuffer))
      
      dialog <- gtkDialogNewWithButtons(title=gettext("Recode"), parent=vvTreeView$getToplevel(), flags=c("modal", "destroy-with-parent"),
                                        "gtk-ok", GtkResponseType["ok"], 
                                        "gtk-cancel", GtkResponseType["cancel"],
                                        show=FALSE)
      dialog$setDefaultSize(430, 300)
      dialog[["vbox"]]$packStart(vbox1)
      gSignalConnect(dialog, "response", onResponse)
      response <- dialog$run()
    },
    
    onEditValueLabels = function(){
      selec      <- .self$getSelected()
      var.name   <- selec["vars"]
      var.lab    <- selec["var.labs"]
      var.name.label <- gtkLabelNew(gettext("Variable Name"))
      var.lab.label  <- gtkLabelNew(gettext("Variable Label"))
      var.name.entry <- gtkEntryNew()
      var.lab.entry  <- gtkEntryNew()
      var.name.entry$setText(var.name)
      var.lab.entry$setText(var.lab)
      var.name.entry["editable"] <- FALSE
      var.lab.entry["editable"]  <- FALSE
      table <- gtkTableNew(2, 2)
      table$attach(var.name.label, 0, 1, 0, 1, "shrink", "shrink", 5, 2)
      table$attachDefaults(var.name.entry, 1, 2, 0, 1)
      table$attach(var.lab.label , 0, 1, 1, 2, "shrink", "shrink", 5, 2)
      table$attachDefaults(var.lab.entry , 1, 2, 1, 2)
      table$setRowSpacings(2)
      table$setColSpacings(5)
      hbox <- gtkHBoxNew()
      hbox$packStart(table)
      
      data.set   <- data$getData.set()
      var        <- data.set[[var.name]]
      missing.values(var) <- NULL
      labels     <- labels(var)
      deletepix  <- gdkPixbufNewFromFile(file.path(rzSettings$getRzPath(), "images/oxygen/delete.png"))$retval
      addpix     <- gdkPixbufNewFromFile(file.path(rzSettings$getRzPath(), "images/oxygen/add.png"   ))$retval
      liststore2 <- gtkListStoreNew("GdkPixbuf", "character", "character", "logical")
      if(!is.null(labels)){
        labels <- data.frame(values=labels@values, labels=labels@.Data, stringsAsFactors=FALSE)
      } else {
        labels <- data.frame(values=numeric(0), labels=character(0), stringsAsFactors=FALSE)
      }
      
      values <- sort(unique(c(labels[[1]], na.omit(as.numeric(var)))))
      labels <- merge(data.frame(values=values), labels, all=TRUE)
      for(i in seq_len(nrow(labels))) {
        iter <- liststore2$append()$iter
        liststore2$set(iter,
                       0, deletepix,
                       1, labels$values[i],
                       2, ifelse(is.na(labels$labels[i]), "", labels$labels[i]),
                       3, FALSE)
      }
      iter <- liststore2$append()$iter
      liststore2$set(iter, 0, addpix, 1, "", 2, "", 3, TRUE)
      
      rp      <- gtkCellRendererPixbuf()
      rt.val  <- gtkCellRendererText()
      rt.lab  <- gtkCellRendererText()
      rt.val["editable"] <- TRUE          
      rt.lab["editable"] <- TRUE
      
      col.icon <- gtkTreeViewColumnNewWithAttributes("", rp, "pixbuf"=0)
      col.val  <- gtkTreeViewColumnNewWithAttributes(gettext("Values"), rt.val, "text"=1)
      col.lab  <- gtkTreeViewColumnNewWithAttributes(gettext("Value Labels"), rt.lab, "text"=2)
      col.icon$setSizing("fixed")
      col.icon$setFixedWidth(50)
      col.icon$setResizable(FALSE)
      col.val$setSizing("fixed")
      col.val$setFixedWidth(50)
      col.val$setResizable(TRUE)
      col.lab$setSizing("fixed")
      col.lab$setResizable(FALSE)
      
      tw <- gtkTreeViewNewWithModel(liststore2)
      tw["enable-grid-lines"] <- GtkTreeViewGridLines["both"]
      tw["fixed-height-mode"] <- TRUE
      tw["rules-hint"] <- TRUE
      tw$appendColumn(col.icon)
      tw$appendColumn(col.val)
      tw$appendColumn(col.lab)
      
      gSignalConnect(rt.val, "edited", function(object, path, new_text){
        text <- localize(new_text)
        e <- try(as.numeric(text), silent=TRUE)
        if(is.na(e)) return()
        if(nzchar(text)){
          iter <- liststore2$getIterFromString(path)$iter
          add  <- liststore2$getValue(iter, 3)$value
          liststore2$set(iter, 0, deletepix, 1, text, 3, FALSE)
          if(add) {
            iter <- liststore2$append()$iter
            liststore2$set(iter, 0, addpix, 1, "", 2, "", 3, TRUE)
          }
        }
      })
      
      gSignalConnect(rt.lab, "edited", function(object, path, new_text){
        text <- localize(new_text)
        iter <- liststore2$getIterFromString(path)$iter
        liststore2$set(iter, 2, text)
      })
      
      gSignalConnect(tw, "row-activated", function(object, path, column){
        if(column==col.icon){
          iter <- liststore2$getIter(path)$iter
          add  <- liststore2$getValue(iter, 3)$value
          if(!add){
            liststore2$remove(iter)
          }
        }
      })
      
      dialog <- gtkDialogNewWithButtons(title=gettext("Edit Value Labels"), parent=vvTreeView$getToplevel(), flags=c("modal", "destroy-with-parent"),
                                        "gtk-ok", GtkResponseType["ok"], 
                                        "gtk-cancel", GtkResponseType["cancel"],
                                        show=FALSE)
      dialog$setDefaultSize(430, 300)
      sw1 <- gtkScrolledWindowNew()
      sw1["shadow-type"] <- GtkShadowType["in"]
      sw1$setPolicy("automatic", "automatic")
      sw1$add(tw)
      dialog[["vbox"]]$setSpacing(2)
      dialog[["vbox"]]$packStart(hbox, expand=FALSE)
      dialog[["vbox"]]$packStart(sw1)
      gSignalConnect(dialog, "response", function(object, response.id){
        if(response.id==GtkResponseType["ok"]) {
          labels <- numeric(0)
          iter <- liststore2$getIterFirst()
          while(iter$retval){
            val <- liststore2$getValue(iter$iter, 1)$value
            lab <- liststore2$getValue(iter$iter, 2)$value
            val <- localize(val)
            lab <- localize(lab)
            val <- try(as.numeric(val), silent=TRUE)
            if(!is.na(val)&&nzchar(lab)){
              labels[lab] <- val
            }
            iter$retval <- liststore2$iterNext(iter$iter)
          }
          if(any(duplicated(labels))){
            dialog2 <- gtkMessageDialogNew(dialog, "destroy-with-parent",
                                           GtkMessageType["error"],
                                           GtkButtonsType["close"],
                                           gettext("Duplicate label is detected."))
            dialog2$run()
            dialog2$hide()
            dialog$run()
          } else {
            dialog$hide()
            row <- which(names(data.set)==var.name)
            if (length(labels)==0) {
              labels(data.set[[row]]) <- NULL
              labels <- ""
            } else {
              labels <- sort(labels)
              labels(data.set[[row]]) <- labels              
              labels <- paste(labels, " \"", names(labels), "\"", sep="", collapse=", ")
            }
            data$setData.set(data.set)
            data$constructVariable(row)
            data$linkDataFrame()
            cell.row <- as.character(row - 1)
            .self$setCell(liststore$getIterFromString(cell.row)$iter, column.definition["val.labs"], labels)

          }
        } else {
          dialog$hide()
        }
      })
      dialog$run()
    },
    
    onRowActivated      = function(tw, path, column){
      row   <- as.numeric(getCell(path$toString(), column.definition["index"])) 
      col.title <- column$getData("attr")["title"]
      
      if (col.title=="index") {
        data.set <- NULL
        if (data$getSubset.on() & nzchar(data$getSubset.condition())){
          data.set <- data$getData.set.subset()
        } else {
          data.set <- data$getData.set()
        }
        print(codebook(data.set[ row ]))
        
      } else if (col.title=="val.labs") {
        .self$onEditValueLabels()
      } else if (col.title=="msr") {
        .self$onEditMsr()
      }
      
    },
    
    onQueryTooltip      = function(tw, x, y, keyboard_mode, tooltip){
      if(rzSettings$getPopupOff()) return(FALSE)
      path <- tw$getPathAtPos(x, y - 20)$path
      if(is.null(path)) return(FALSE)
      row  <- as.numeric(getCell(path$toString(), column.definition["index"]))
      char <- NULL
      if (data$getSubset.on() & nzchar(data$getSubset.condition())) {
        char  <- vv$getSummaries.subset()[ row ]
      } else {
        char  <- vv$getSummaries()[ row ]
      }
      tooltip$setMarkup(paste("<span font_family=\"", rzSettings$getMonospaceFontFamily(), "\">", char, "</span>", sep="", collapse=""))
      tw$setTooltipRow(tooltip, path)
      return(TRUE)
    },
    
    changeFont = function(){
      vvTreeView$modifyFont(pangoFontDescriptionFromString(rzSettings$getVariableViewFont()))
    },
    
    onSelectAll = function(...){
      liststore$foreach(function(model, path, iter){
        model$setValue(iter, column.definition["select"], TRUE)
        return(FALSE)
      })
    },
    
    onUnselect = function(...){
      ind.env <- new.env()
      model.selected$foreach(function(model, path, iter){
        ind.env$ind <- c(ind.env$ind,
                         model$getValue(iter, column.definition["index"])$value)
        return(FALSE)
      })
      ind <- as.numeric(ind.env$ind)
      for (i in ind) {
        iter <- liststore$getIterFromString(i-1)$iter
        liststore$setValue(iter, column.definition["select"], FALSE)          
      }
    },
    
    onDelete = function(...){
      dialog <- gtkMessageDialogNew(vvTreeView$getToplevel(), "destroy-with-parent",
                                    GtkMessageType["question"], GtkButtonsType["ok-cancel"],
                                    gettext("Are you sure you want to do that?"))
      response <- dialog$run()
      dialog$hide()
      
      if(response==GtkResponseType["ok"]){
        inds <- .self$getSelectedRows()
        if(length(inds)==0) return()
        data$deleteVars(inds)
        vv$reload()
      }
    },
    
    onDuplicate = function(...){
      inds <- .self$getSelectedRows()
      if(length(inds)==0) return()
      data$duplicate(inds)
      vv$reload()
    }
    
    ))
rzVVCore$accessors(c("widget", "adj"))
