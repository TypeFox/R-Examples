### R code from vignette source 'ex-RGtk2-UImanager-II.Rnw'

###################################################
### code chunk number 1: helperFUnction
###################################################
require(RGtk2)

##' helper function to bypass lack of cached value in method call
##'
##' @param meth method name
##' @param obj method of object's class
##' @return the method
call_meth <- function(meth, obj) {
  if(exists(meth, obj, inherits=FALSE))
    get(meth, obj)
  else
    methods:::envRefInferField(obj, meth, getClass(class(obj)), obj)
}


###################################################
### code chunk number 2: ex-RGtk2-UImanager-II.Rnw:34-36
###################################################
## Stub UI Manager instance for use with examples
uimanager <- gtkUIManager()


###################################################
### code chunk number 3: ui-xml
###################################################
ui.xml <- readLines(out <- textConnection('
<ui>
  <menubar name="menubar">
    <menu name="FileMenu" action="File">
      <menuitem action="Save"/>
      <menuitem action="SaveAs" />
      <menu name="Export" action="Export">
        <menuitem action="ExportToCSV" />
        <menuitem action="ExportToSaveFile" />
      </menu>
      <separator />
      <menuitem name="FileQuit" action="CloseWindow" />
    </menu>
    <menu action="Edit">
      <menuitem name="EditUndo" action="Undo" />
      <menuitem name="EditRedo" action="Redo" />
      <menuitem action="ChangeColumnName" />
    </menu>
    <menu action="Tools">
      <menuitem action="Filter" />
      <menuitem action="Sort" />
    </menu>
  </menubar>
  <toolbar name="toolbar">
    <toolitem action="Save"/>
    <toolitem action="SaveAs"/>
    <separator />
    <toolitem action="CloseWindow"/>
  </toolbar>
</ui>'), warn=FALSE)
close(out)


###################################################
### code chunk number 4: loadUIFromString
###################################################
id <- uimanager$addUiFromString(ui.xml)


###################################################
### code chunk number 5: ex-RGtk2-UImanager-II.Rnw:93-94
###################################################
fun <- function(...) {}


###################################################
### code chunk number 6: ex-RGtk2-UImanager-II.Rnw:102-118
###################################################
file_list <- 
  list(## name, ID, label, accelerator, tooltip, callback
       list("File",NULL,"_File",NULL,NULL,NULL),
       list("Save", "gtk-save", "Save", "<ctrl>S", 
            "Save data to variable", fun),
       list("SaveAs", "gtk-save", "Save as...", NULL, 
            "Save data to variable", fun),
       list("Export", NULL, "Export", NULL, NULL, NULL),
       list("ExportToCSV", "gtk-export", "Export to CSV", 
            NULL, "Save data to CSV file", fun),
       list("ExportToSaveFile", "gtk-export", 
            "Export to save file", NULL, 
            "Save data to save() file", fun),
       list("CloseWindow", "gtk-close", "Close window", 
            "<ctrl>W", "Close current window", fun)
       )


###################################################
### code chunk number 7: addActionGroup
###################################################
action_group <- gtkActionGroup("FileGroup")
action_group$addActions(file_list)


###################################################
### code chunk number 8: ex-RGtk2-UImanager-II.Rnw:129-130
###################################################
uimanager$insertActionGroup(action_group, 0)


###################################################
### code chunk number 9: GUILayout (eval = FALSE)
###################################################
## window <- gtkWindow(show = FALSE)
## ##
## vbox <- gtkVBox()
## window$add(vbox)
## ##
## menubar <- uimanager$getWidget("/menubar")
## vbox$packStart(menubar, FALSE)
## toolbar <- uimanager$getWidget("/toolbar")
## vbox$packStart(toolbar, FALSE)
## ## ...


###################################################
### code chunk number 10: ex-RGtk2-UImanager-II.Rnw:171-172 (eval = FALSE)
###################################################
## window$addAccelGroup(uimanager$getAccelGroup())


###################################################
### code chunk number 11: ex-RGtk2-UImanager-II.Rnw:193-200
###################################################
Command <- setRefClass("Command",
                       fields = list(
                         receiver="ANY",
                         meth="character",
                         params="list",
                         old_params="list"
                         ))


###################################################
### code chunk number 12: ex-RGtk2-UImanager-II.Rnw:208-218
###################################################
Command$methods(
        initialize = function(receiver, meth, ...) {
          .params <- list(...)
          initFields(receiver = receiver, meth = meth, 
                     params = .params, old_params = .params)
          callSuper()
        },
        execute = function(params) {
          do.call(call_meth(meth, receiver), params)
        })


###################################################
### code chunk number 13: ex-RGtk2-UImanager-II.Rnw:225-233
###################################################
Command$methods(
                do = function() {
                  out <- execute(params)
                  old_params$value <<- out
                },
                undo = function() execute(old_params)
                )



###################################################
### code chunk number 14: illustrateCommand
###################################################
x <- 1
set_x <- function(value) {
  old <- x
  x <<- value
  old
}
cmd <- Command$new(.GlobalEnv, "set_x", value = 2)
cmd$do(); x


###################################################
### code chunk number 15: ex-RGtk2-UImanager-II.Rnw:252-253
###################################################
cmd$undo();


###################################################
### code chunk number 16: ex-RGtk2-UImanager-II.Rnw:255-256
###################################################
x


###################################################
### code chunk number 17: typicalAction (eval = FALSE)
###################################################
## cmd <- Command$new(df_model, "set_col_name", j=j, value=value)
## command_stack$add(cmd)


###################################################
### code chunk number 18: col_name_methods (eval = FALSE)
###################################################
## DfModel$methods(
##                 get_col_name = function(j) varnames[j,1],
##                 get_col_names = function() varnames[ ,1],
##                 set_col_name = function(j, value) {
##                   "Set name, return old"
##                   old_col_name <- get_col_name(j)
##                   varnames[j,1] <<- value
##                   old_col_name
##                 })


###################################################
### code chunk number 19: ensure_type
###################################################
##' S3 generic to ensure we don't change data type when assigning into column
##'
##' @param x column values
##' @param value new value
##' @return coerced new value
ensure_type <- function(x, value) UseMethod("ensure_type")
ensure_type.default <- function(x, value) value
ensure_type.character <- function(x, value) as.character(value)
ensure_type.factor <- function(x, value) {x[length(x) + 1] <- value; tail(x, n=1)}
ensure_type.numeric <- function(x, value) as.numeric(value)
ensure_type.integer <- function(x, value) as.integer(value)
ensure_type.logical <- function(x, value) as.logical(value)


###################################################
### code chunk number 20: DfModel
###################################################
## Define a model to hold the model for an editable data frame
sapply(c("RGtkDataFrame"), setOldClass)
DfModel <- setRefClass("DfModel",
                       fields=list(
                         store="RGtkDataFrame",
                         filtered="ANY",
                         name="character",
                         varnames="RGtkDataFrame",
                         casenames="RGtkDataFrame"
                         ))

## Initialize along with a column for filtering
DfModel$methods(
                initialize=function(DF, nm, ...) {
                  store <<- rGtkDataFrame(cbind(DF, `_visible`=rep(TRUE, nrow(DF))))
                  varnames <<- rGtkDataFrame(data.frame(names(DF), stringsAsFactors=FALSE))
                  casenames <<- rGtkDataFrame(data.frame(rownames(DF), stringsAsFactors=FALSE))
                  if(missing(nm))
                    name <<- deparse(substitute(DF))
                  else
                    name <<- nm
                  filtered <<- store$filter()
                  filtered$setVisibleColumn(length(DF))
                  callSuper()
                })

## Methods to work with the underlying data frame (Get, save, ...)
DfModel$methods(
                get_dataframe=function() {
                  DF <- store[,seq_len(ncol(store)-1)]
                  dimnames(DF) <- list(casenames[,1], varnames[,1])
                  DF
                },
                save=function(nm) {
                  "Save to global workspace"
                  if(!missing(nm))
                    name <<- nm
                  assign(name, get_dataframe(), envir=.GlobalEnv)
                  
                },
                export_to_csv=function(f)  {
                  "Export to csv file"
                  write.csv(get_dataframe(), file=f)
                },
                export_to_save=function(f) {
                  "Export using save()"
                  assign(name, get_dataframe())
                  save(list=name, file=f)
                },
                no_rows=function() dim(store)[1],
                no_cols=function() dim(store)[2] - 1L
                )
## Methods to get and set a cell value. 
DfModel$methods(
                get_cell=function(i,j) {
                  "Return cell value"
                  store[i,j]
                },
                set_cell=function(i, j, value) {
                  "Set cell, return old_value"
                  old <- get_cell(i,j)
                  store[i,j] <<- ensure_type(store[1,j], value)
                  old
                })
## Methods for column names. Similar one for rownames could be implemented, but we
## don't show these in our view. So leave to the reader/
DfModel$methods(
                get_col_name=function(j) varnames[j,1],
                get_col_names=function() varnames[,1],
                set_col_name=function(j, value) {
                  "Set name, return old"
                  old <- get_col_name(j)
                  varnames[j,1] <<- value
                  old
                })
## Code for filtering the display.
DfModel$methods(
                get_filter=function() {
                  "Return logical indicating filter"
                  store[,ncol(store)]
                },
                set_filter=function(value) {
                  "Filter by value. Return old filter value"
                  if(!is.logical(value)) stop("Filter requires a logical variable")
                  ind <- rep(value, length.out=no_rows())
                  old <- get_filter()
                  store[,ncol(store)] <<- value
                  old
                })

## In RGtk2, one can't both sort and filter by proxy. Since R makes sorting easy, 
## we let Gtk handle the filtering and implement sorting below. The "old" value 
## returned by this is what is needed to reverse a sort.
DfModel$methods(
                reorder=function(value) {
                  "Reorder data frame. Return order(value)"
                  perm <- as.integer(value)
                  if(length(perm) != nrow(store)) stop("reorder requires a permutation")
                  if(length(perm) != length(unique(perm))) stop("value has repeated values")
                  if(min(perm) != 1 || max(perm) != nrow(store)) stop("value is not permutation of row indices")

                  store[,] <<- store[perm,]
                  order(perm)  # will revers a[ind][order(ind)] is involution
                })


###################################################
### code chunk number 21: CommandStack
###################################################
## Command Stack
## A list with ptr. delegates call of do or undo to appropriate command
CommandStack <- setRefClass("CommandStack",
                            fields=list(
                              l="list",
                              ptr="integer"
                              ))
## initialize method just sets the list and pointer to a default.
CommandStack$methods(
                     initialize=function() {
                       initFields(l=list(), ptr=0L)
                       callSuper()
                     })
## do method finds the right command then delegates to the commands do method
## undo is similar
## The can_do and can_undo commands are used to check if the command stack allows for
## these operations
CommandStack$methods(
                     do=function() {
                       if(!can_do()) return()
                       cmd <- l[[ptr]]
                       ptr <<- ptr + 1L
                       cmd$do()
                     },
                     undo=function() {
                       if(!can_undo()) return()
                       cmd <- l[[ptr-1]]
                       ptr <<- ptr - 1L
                       cmd$undo()
                     },
                     can_do=function() ptr > 0 && ptr <= length(l),
                     can_undo=function() ptr > 1
                     )
## Methods to add to and clear the command stack
CommandStack$methods(
                     add=function(cmd, call=TRUE) {
                       if(ptr <= 1) {
                         l <<- list(cmd)
                         ptr <<- 1L
                       } else {
                         l <<- l[1:(ptr-1)]
                         l[[length(l) + 1]] <<- cmd
                       }
                       if(call)
                         do()
                     },
                     clear=function(cmd) {
                       l <<- list(); ptr <<- 0L
                     })


###################################################
### code chunk number 22: addCellRenderer
###################################################
## We create our cellrenderers using an S3 generic to dispatch based on the class of the column. This
## works out well, as the view is column based as well. The editable commands have 
## to find a row, a column and a value before make a command to add to the command stack.
## The row comes from the path, but must be "unfiltered" to point to the original data store. 
## The column is passed into the function by the caller.


##' Create an appropriate cell renderer
##'
##' @param x vector to display in column
##' @param nm name of vector for title
##' @param obj a DfModel instance
##' @param view GtkTreeView instance we add the cellrenderer to
##' @param command_stack a CommandStack instance needed for the callback
##' @return NULL
add_cellrenderer_by_class <- function(x, nm, obj, view, j, command_stack) UseMethod("add_cellrenderer_by_class")
add_cellrenderer_by_class.default <- function(x, nm, obj, view, j, command_stack) {
  cr <- gtkCellRendererText()
  cr['editable'] <- TRUE
  gSignalConnect(cr, "edited", f=function(cr, path, newtext) {
    i <- as.numeric(path) + 1
    i <- which(obj$get_filter())[i]     # in regular    
    value <- newtext
    cmd <- Command$new(obj, "set_cell", i=i, j=j, value=value)
    command_stack$add(cmd)
  })
  view$insertColumnWithAttributes(position=-1, 
                                  title=nm,
                                  cell=cr,
                                  text=j-1)
}

add_cellrenderer_by_class.logical <- function(x, nm, obj, view, j, command_stack) {
  cr <- gtkCellRendererToggle()
  cr['activatable'] <- TRUE
  gSignalConnect(cr, "toggled", function(w, path) {
    i <- as.numeric(path) + 1           # in filtered
    i <- which(obj$get_filter())[i]     # in regular
    value <- !obj$get_cell(i,j)
    cmd <- Command$new(obj, "set_cell", i=i, j=j, value=value)
    command_stack$add(cmd)
  })
  view$insertColumnWithAttributes(position=-1, 
                                  title=nm,
                                  cell=cr,
                                  active=j-1)
}

add_cellrenderer_by_class.factor <- function(x, nm, obj, view, j, command_stack) {
  cr <- gtkCellRendererCombo()
  cr_store <- rGtkDataFrame(sort(levels(x)))
  cr['model'] <- cr_store
  cr['text-column'] <- 0
  cr['has-entry'] <- FALSE
  cr['editable'] <- TRUE
  gSignalConnect(cr, "changed", function(w, path, iter, user.data) {
    i <- as.numeric(path) + 1
    i <- which(obj$get_filter())[i]     # in regular    
    value <- cr_store$getValue(iter, 0)$value
    cmd <- Command$new(obj, "set_cell", i=i, j=j, value=value)
    command_stack$add(cmd)
  })
  view$insertColumnWithAttributes(position=-1, 
                                  title=nm,
                                  cell=cr,
                                  text=j-1L)
}


###################################################
### code chunk number 23: EditDataFrame
###################################################
## Main reference class to edit a data frame within a GUI
## The view relies on a DataFrameModel and CommandStack instance, each of which is 
## defined within the initialize method.
EditDataFrame <- setRefClass("EditDataFrame",
                             fields=list(
                               df_model="ANY",
                               command_stack="ANY",
                               actions="list",
                               ## layout
                               mainwindow="ANY",
                               statusbar="ANY",
                               uimanager="ANY", 
                               view="ANY"
                               ))
## The initialize method makes several different calls. Here we initialize the actions into action group.
EditDataFrame$methods(
                      initialize_actions=function(box) {
                        ## our callback. Calls an appropriately named method of this class.
                        fun=function(action) {
                          meth <- action$getName()
                          out <- try(do.call(call_meth(meth, .self), list()), silent=TRUE)
                        }

                        ## Define action groups in a list
                        fileL <- list(## name, ID, label, accelerator, tooltip, callback
                                      list("File",NULL,"_File",NULL,NULL,NULL),
                                      list("Save", "gtk-save", "Save", "<ctrl>S", "Save data to variable", fun),
                                      list("SaveAs", "gtk-save", "Save as...", NULL, "Save data to variable", fun),
                                      list("Export", NULL, "Export", NULL, NULL, NULL),
                                      list("ExportToCSV", "gtk-export", "Export to CSV", NULL, "Save data to CSV file", fun),
                                      list("ExportToSaveFile", "gtk-export", "Export to save() file", NULL, "Save data to save() file", fun),
                                      list("CloseWindow", "gtk-close", "Close window", "<ctrl>W", "Close current window", fun)
                                      )

                        editL <- list(## name, ID, label, accelerator, tooltip, callback
                                      list("Edit", NULL, "_Edit", NULL, NULL, NULL),
                                      list("Undo", "gtk-undo", "Undo", "<ctrl>Z",  "Undo last command", fun),
                                      list("Redo", "gtk-redo", "Redo", "<ctrl>U", "Redo undo command", fun),
                                      list("ChangeColumnName", "gtk-change", "Change column name",
                                                          NULL, "Change a column name", fun)
                                      )
                        
                        toolL <- list(
                                      list("Tools", NULL, "_Tools", NULL, NULL, NULL),
                                      list("Filter", "gtk-filter", "Filter", NULL, "Filter data frame", fun),
                                      list("Sort", "gtk-sort", "Sort", NULL, "Sort data frame by column name", fun)                                      )
                        l <- list(fileL, editL, toolL)

                        ## create UI manager, insert action groups
                        uimanager <<- gtkUIManager()
                        for(i in seq_along(l)) {
                          ag <- gtkActionGroup(sprintf("Group%s",i))
                          ag$addActions(l[[i]])
                          uimanager$insertActionGroup(ag, i-1)
                        }
                      })

## Here we initialize the UI 
EditDataFrame$methods(
                      initialize_ui=function() {
                        ## define xml specifying menu bars and toolbars
                        ui.xml <- readLines(out <- textConnection('
<ui>
  <menubar name="menubar">
    <menu name="FileMenu" action="File">
      <menuitem action="Save"/>
      <menuitem action="SaveAs" />
      <menu name="Export" action="Export">
        <menuitem action="ExportToCSV" />
        <menuitem action="ExportToSaveFile" />
      </menu>
      <separator />
      <menuitem name="FileQuit" action="CloseWindow" />
    </menu>
    <menu action="Edit">
      <menuitem name="EditUndo" action="Undo" />
      <menuitem name="EditRedo" action="Redo" />
      <menuitem action="ChangeColumnName" />
    </menu>
    <menu action="Tools">
      <menuitem action="Filter" />
      <menuitem action="Sort" />
    </menu>
  </menubar>
  <toolbar name="toolbar">
    <toolitem action="Save"/>
    <toolitem action="SaveAs"/>
    <separator />
    <toolitem action="CloseWindow"/>
  </toolbar>
</ui>'), warn=FALSE)
                        close(out)
                        ## specify the UI using XML specification
                        id <- uimanager$addUiFromString(paste(ui.xml, collapse="\n")) 
                      })

## Here we layout the GUI
EditDataFrame$methods(
                      make_gui=function() {
                        DF <- df_model$get_dataframe()
                        nms <- names(DF)
                        view <<- gtkTreeView(df_model$filtered)
                        sapply(seq_len(length(DF)), function(j) {
                          add_cellrenderer_by_class(DF[[j]], nms[j], df_model, view, j, command_stack)
                        })
                        ##
                        ## place into GUI
                        mainwindow <<- w <- gtkWindow(show=FALSE)
                        #
                        vbox <- gtkVBox()
                        w$add(vbox)
                        #
                        menubar <- uimanager$getWidget("/menubar")
                        vbox$packStart(menubar, FALSE)
                        toolbar <- uimanager$getWidget("/toolbar")
                        vbox$packStart(toolbar, FALSE)
                        w$addAccelGroup(uimanager$getAccelGroup())
                        ##
                        sw <- gtkScrolledWindow()
                        sw$add(view)
                        vbox$PackStart(sw, TRUE, TRUE)
                        ##
                        statusbar <<- gtkStatusbar()
                        statusbar$getChildren()[[1]]$setSizeRequest(-1, 25)
                        vbox$PackStart(statusbar, FALSE)
                        w$show()
                      })

## This method call updates the GUI: sets the redo/undo buttons and the status bar.
EditDataFrame$methods(
                      update_UI=function(event="") {
                        ## update actions
                        ## Could save undo/redo actions as we look them up inefficiently below
                        undo <- redo <- NULL
                        for(i in uimanager$getActionGroups()) {
                          tmp <- i$getAction("Redo")
                          if(!is.null(tmp)) redo <- tmp
                          tmp <- i$getAction("Undo")
                          if(!is.null(tmp)) undo <- tmp
                        }
                        undo$setSensitive(command_stack$can_undo())
                        redo$setSensitive(command_stack$can_do())

                        ## update status bar
                        tpl <- "Editing %s. Showing %s lines of %s."
                        statusbar$push(statusbar$getContextId("message"),
                                       sprintf(tpl, df_model$name,
                                               sum(df_model$get_filter()),
                                               df_model$no_rows()))
                        
                      })
## This sets up a callback when parts of the DataFrameModel change. Here is how
## we synchronize column names and why we used a RGtkDataFrame class to hold them in the definition
## of the DataFrameModel class
EditDataFrame$methods(
                      synchronize_view=function() {
                        gSignalConnect(df_model$store, "row-changed", function(model, path, iter) {
                          update_UI()
                        })
                        gSignalConnect(df_model$varnames, "row-changed", function(model, path, iter) {
                          j <- as.numeric(path$toString()) + 1
                          value <- df_model$varnames[j,1]
                          col <- view$getColumn(j-1)
                          col['title'] <- value
                          update_UI()
                        })
                      })
## Finally an initialization method
EditDataFrame$methods(
                      initialize=function(DF) {
                        if(!is.data.frame(DF))
                          stop("Requires a data frame")
                        initFields(df_model=DfModel$new(DF),
                                   command_stack=CommandStack$new())

                        initialize_actions()
                        initialize_ui()
                        make_gui()
                        synchronize_view()
                        update_UI()
                        callSuper()
                      })


###################################################
### code chunk number 24: Actions
###################################################
## Actions are defined here
## Basically we delegate down to data frame model
## We are lazy about some dialogs, so use the gWidgets package
require(gWidgets); options(guiToolkits="RGtk2")
EditDataFrame$methods(
                      Save=function() {
                        df_model$save()
                        command_stack$clear()
                      })

EditDataFrame$methods(
                      Undo=function() {command_stack$undo()},
                      Redo=function() {command_stack$do()}
                      )

EditDataFrame$methods(
                      SaveAs=function() {
                        current_vars <- ls(envir=.GlobalEnv)
                        dlg <- gbasicdialog("Select a variable name...", parent=mainwindow, handler=function(h,...) {
                          var <- svalue(e)
                          if(nchar(var)) {
                            if(exists(var, .GlobalEnv)) {
                              if(!gconfirm(c("Variable exists", "Really overwrite?"), parent=dlg))
                                return()
                            }
                            df_model$save(var)
                            update_UI()
                            command_stack$clear()
                          }
                        })
                        g <- ggroup(cont=dlg, horizontal=FALSE)
                        glabel("Variable name to save as:", cont=g)
                        e <- gedit("", cont=g)
                        e[] <- current_vars
                        visible(dlg, set=TRUE)
                      })

EditDataFrame$methods(
                      ExportToCSV=function() {
                        f <- gfile("Select a filename", type="save")
                        if(!is.na(f))
                          df_model$export_to_csv(f)
                      })

EditDataFrame$methods(
                      ExportToSaveFile=function() {
                        f <- gfile("Select a filename", type="save")
                        if(!is.na(f))
                          df_model$export_to_save(f)
                      })

EditDataFrame$methods(
                      CloseWindow=function() {
                        if(command_stack$can_undo() || command_stack$can_do()) {
                          if(!gconfirm(c("Really quit", "There are pending changes"),
                                       parent=mainwindow))
                            return()
                        }
                        mainwindow$destroy()
                      })

EditDataFrame$methods(
                      ChangeColumnName=function() {
                        j <- NA; value <- character(0)
                        ## get column and new names
                        dlg <- gbasicdialog("Rename a column", parent=mainwindow)
                        g <- ggroup(horizontal=FALSE, cont=dlg)
                        varnames <- df_model$get_col_names()
                        tbl <- gtable(data.frame(Variables=varnames), cont=g, expand=TRUE)
                        size(tbl) <- c(300, 250)
                        l <- glabel("Select a variable", cont=g)
                        e <- gedit("", cont=g); enabled(e) <- FALSE
                        addHandlerClicked(tbl, handler=function(h,...) {
                          val <- svalue(h$obj)
                          assign("j", match(val, varnames), inherits=TRUE)
                          if(!is.na(j)) {
                            svalue(l) <- sprintf("Change %s to:", val)
                            enabled(e) <- TRUE
                          } else {
                            svalue(l) <- "Select a variable"
                            svalue(e) <- ""
                            enabled(e) <- FALSE
                          }
                        })
                        addHandlerKeystroke(e, handler=function(h,...) {
                          assign("value", svalue(h$obj), inherits=TRUE)
                        })
                        ret <- visible(dlg, set=TRUE)
                        if(ret && !is.na(j)) {
                          cmd <- Command$new(df_model, "set_col_name", j=j, value=value)
                          command_stack$add(cmd)
                        }
                      })
                      
EditDataFrame$methods(
                      Filter=function() {
                        ind <- NULL
                        dlg <- gbasicdialog("Enter an expression", parent=mainwindow, handler=function(h,...) {
                          val <- svalue(e)
                          DF <- df_model$get_dataframe()
                          out <- try(eval(parse(text=val), DF), silent=FALSE)
                          if(!inherits(out, "try-error"))
                            assign("ind", out, inherits=TRUE)
                        })
                        g <- ggroup(cont=dlg, horizontal=FALSE)
                        glabel("Enter an expression to filter by:", cont=g)
                        e <- gedit("", cont=g)
                        ret <- visible(dlg, set=TRUE)
                        if(ret && is.logical(ind)) {
                          cmd <- Command$new(df_model, "set_filter", value=ind)
                          command_stack$add(cmd)
                        }
                      })

EditDataFrame$methods(
                      Sort=function() {
                        perm <- integer(0)
                        DF <- df_model$get_dataframe()
                        varnames <- df_model$get_col_names()
                        dlg <- gbasicdialog("Sort by:", parent=mainwindow, handler=function(h,...) {
                          var <- svalue(tbl, index=TRUE)
                          if(length(var) == 0) return()
                          x <- DF[,var]
                          assign("perm", order(x, decreasing=svalue(decreasing)), inherits=TRUE)
                        })
                        g <- ggroup(horizontal=FALSE, cont=dlg)
                        tbl <- gtable(data.frame(Variables=varnames), cont=g)
                        size(tbl) <- c(300, 250)
                        decreasing <- gcheckbox("Decreasing?", checked=FALSE, cont=g)
                        ret <- visible(dlg, set=TRUE)
                        ##
                        if(ret && length(perm)) {
                          cmd <- Command$new(df_model, "reorder",  value=perm)
                          command_stack$add(cmd)
                        }
                      })


###################################################
### code chunk number 25: testItOut
###################################################
## Test it out....
require(MASS)
DF <- Cars93[sample(1:93, 20),c(1, 5, 26)]; DF$American <- DF$Origin == "USA"
a = EditDataFrame$new(DF)


###################################################
### code chunk number 26: oldWay
###################################################
## Old way to add actions, menu bar, For comparison
## not called by the initialize method
EditDataFrame$methods(
                      initialize_actions_old=function() {
                           ## actions. Must have a matching method
                        al <- list()
                        al$save <- gtkAction("Save", "Save", "Save data to variable", "gtk-save")
                        al$saveas <- gtkAction("SaveAs", "Save as...", "Save data to variable", "gtk-save")
                        al$exportAsCSV <- gtkAction("ExportToCSV", "Export to CSV", "Save data to CSV file", "gtk-export")
                        al$exportAsSaveFile <- gtkAction("ExportToSaveFile", "Export to save() file", "Save data to save() file", "gtk-export")
                        al$close <- gtkAction("CloseWindow", "Close window", "Close current window", "gtk-close")
                        ## Edit menu
                        al$undo <- gtkAction("Undo", "Undo", "Undo last command", "gtk-undo")
                        al$redo <- gtkAction("Redo", "Redo", "Redo undo command", "gtk-redo")
                        al$change_column_name <- gtkAction("ChangeColumnName", "Change column name",
                                                          "Change a column name", "gtk-change")
                        ## Tools
                        al$filter <- gtkAction("Filter", "Filter", "Filter data frame", "gtk-filter")
                        al$sort <- gtkAction("Sort", "Sort", "Sort data frame by column name", "gtk-sort")
                        
                        ## stub handler
                        sapply(al, gSignalConnect, "activate", function(action) {
                          meth <- action$getName()
                          out <- try(do.call(get(meth, .self), list()), silent=TRUE)
                          print(out)
                        })
                        
                        actions <<- al
                      },
                      make_menu=function(box) {
                        mb <- gtkMenuBar()

                        fileMenu <- gtkMenu()
                        fileItem <- gtkMenuItem("File")
                        fileItem$setSubmenu(fileMenu)
                        sapply(c("save","saveas", "exportAsCSV","exportAsSaveFile","close"),
                               function(act)
                               fileMenu$append(actions[[act]]$createMenuItem()))

                        editMenu <- gtkMenu()
                        editItem <- gtkMenuItem("Edit")
                        editItem$setSubmenu(editMenu)
                        sapply(c("undo","redo", "change_column_name"),
                               function(act)
                               editMenu$append(actions[[act]]$createMenuItem()))

                        toolsMenu <- gtkMenu()
                        toolsItem <- gtkMenuItem("Tools")
                        toolsItem$setSubmenu(toolsMenu)
                        sapply(c("filter", "sort"),
                               function(act)
                               toolsMenu$append(actions[[act]]$createMenuItem()))

                        sapply(list(fileItem, editItem, toolsItem), mb$append)
                        box$packStart(mb, FALSE)

                      }
)


