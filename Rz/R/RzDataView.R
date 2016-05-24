dataView <-
setRefClass("RzDataView",
  fields = c("RzData", "model"),
  methods = list(
    initialize        = function(...) {
      initFields(...)
      df <- RzData$getData.frame()
      df2 <- data.frame(index=seq_len(nrow(df)))
      
      model <<- rGtkDataFrameNew(df)
      treeView <- gtkTreeViewNewWithModel(model)
      treeView["enable-grid-lines"] <- GtkTreeViewGridLines["both"]
      treeView["fixed-height-mode"] <- TRUE
      treeView["rules-hint"] <- TRUE

      model2 <- rGtkDataFrameNew(df2)
      treeView2 <- gtkTreeViewNewWithModel(model2)
      treeView2["enable-grid-lines"] <- GtkTreeViewGridLines["both"]
      treeView2["can-focus"] <- FALSE
      selec <- treeView2$getSelection()
      gSignalConnect(selec, "changed", function(object) object$unselectAll())
      
      rt.index    <- gtkCellRendererText()
      color       <- rt.index["cell-background-gdk"]
      color$red   <- 45000L
      color$green <- 45000L
      color$blue  <- 45000L
      rt.index["cell-background-gdk"] <- color
      rt.index["xalign"] <- 0.5
      renderer <- gtkCellRendererText()
      
      columns <- lapply(seq_len(dim(model)[2]),function(x){
        gtkTreeViewColumn(title=dimnames(model)[[2]][x], cell=renderer, "text"= x-1)})
      lapply(columns, gtkTreeViewColumnSetSizing, "fixed")
      lapply(columns, gtkTreeViewColumnSetResizable, TRUE)
      lapply(columns, gtkTreeViewColumnSetFixedWidth, 50)
      lapply(columns, treeView$AppendColumn)
      
      column1 <- gtkTreeViewColumn(" ", cell=rt.index, "text"=0)
      column1$setSizing("automatic")
      column1$setResizable(FALSE)
      treeView2$appendColumn(column1)
      
      scrolledWindow <- gtkScrolledWindowNew()
      scrolledWindow["shadow-type"] <- GtkShadowType["in"]
      scrolledWindow$add(treeView)
      
      scrolledWindow2 <- gtkScrolledWindowNew(NULL, scrolledWindow$getVadjustment())
      scrolledWindow2["shadow-type"] <- GtkShadowType["in"]
      scrolledWindow2["hscrollbar-policy"] <- GtkPolicyType["never"]
      scrolledWindow2["vscrollbar-policy"] <- GtkPolicyType["never"]
      scrolledWindow2$add(treeView2)
      
      hbox <- gtkHBoxNew()
      hbox$packStart(scrolledWindow2, expand=FALSE)
      hbox$packStart(scrolledWindow)
      
      
      win <- gtkWindowNew(show=FALSE)
      win$add(hbox)
      win["width-request"] <- 800
      win["height-request"] <- 600
      win["allow-shrink"] <- TRUE
      win["title"] <- gettext("Data")
      win["window-position"] <- GtkWindowPosition["center"]
      
      accel.group <- gtkAccelGroupNew()
      win$addAccelGroup(accel.group)
      gSignalConnect(win, "activate-default", win$destroy)

      accel <- gtkAcceleratorParse("<Ctrl>w")
      win$addAccelerator("activate-default", accel.group, accel$accelerator.key,
                         as.numeric(accel$accelerator.mods), 1)
      accel <- gtkAcceleratorParse("<Ctrl>1")
      win$addAccelerator("activate-default", accel.group, accel$accelerator.key,
                         as.numeric(accel$accelerator.mods), 1)
      
      win$showAll()
      
    }
  )
)
#dataView$accessors(c())
