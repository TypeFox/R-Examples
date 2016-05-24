################################################################################
# TODO LIST
# TODO: View attributes: Does not work well...


################################################################################
# CHANGE LOG (last 20 changes)
# 06.01.2016: Fixed attributes window bug. Error when close using X.
# 26.10.2015: Fixed attributes window bug.
# 04.10.2015: Added options to limit number of rows, and show attributes.
# 28.08.2015: Added importFrom.
# 01.06.2015: Fixed bug column names not saved. Introduced 02.01.2015 with attributes.
# 11.05.2015: Accepts (the first) column name containing the string 'Sample'
# as alternative to colum name 'Sample.Name'. 'Sample' is case in-sensitive.
# 04.05.2015: Added 'Sample.File.Name' as a defined alternative to Sample.Name.
# 02.01.2015: Copy attribute list to new object upon 'Save As'.
# 11.10.2014: Added 'focus', added 'parent' parameter.
# 28.06.2014: Added help button and moved save gui checkbox.
# 08.05.2014: Implemented 'checkDataset'.
# 02.12.2013: Added parameter 'name' for selection of 'data' in drop menu.
# 30.11.2013: Added info also when 'data' is passed.
# 20.11.2013: Specified package for function 'gtable' -> 'gWidgets::gtable'
# 16.11.2013: New parameter 'data' and 'edit' and implementation.
# 18.07.2013: Check before overwrite object.
# 11.06.2013: Added 'inherits=FALSE' to 'exists'.
# 21.05.2013: Added 'copy to clipboard'
# 17.05.2013: listDataFrames() -> listObjects()
# 09.05.2013: First version.

#' @title Edit or View Data Frames
#'
#' @description
#' GUI to edit and view data frames.
#'
#' @details Select a data frame from the dropdown and view/edit. Optionally
#' save as a new dataframe.
#' @param env environment in wich to search for data frames.
#' @param savegui logical indicating if GUI settings should be saved in the environment.
#' @param data data.frame for instant viewing.
#' @param name character string with the name of the provided dataset.
#' @param edit logical TRUE for enable edit .
#' @param debug logical indicating printing debug information.
#' @param parent widget to get focus when finished.
#' 
#' @export
#' 
#' @importFrom utils help write.table
#' 
#' @return TRUE
#' 
#' @seealso \code{\link{trim_gui}}, \code{\link{cropData_gui}}, \code{\link{combine_gui}}

editData_gui <- function(env=parent.frame(), savegui=NULL, data=NULL,
                         name=NULL, edit=TRUE, debug=FALSE, parent=NULL){
  
  .gData <- data
  .gDataName <- name
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  if(edit){
    guiTitle="Edit or view data frame"
  } else{
    guiTitle="View data frame"
  }
  
  # Create windows.
  w <- gwindow(title=guiTitle, visible=FALSE)
  w_attributes <- gwindow(title="Attributes", visible=FALSE)
  attr_text <- gtext("", container=w_attributes)

  # Runs when window is closed.
  addHandlerDestroy(w, handler = function (h, ...) {
    
    # Save GUI state.
    .saveSettings()
    
    # Focus on parent window.
    if(!is.null(parent)){
      focus(parent)
    }
    
  })
  
  gv <- ggroup(horizontal=FALSE,
               spacing=8,
               use.scrollwindow=FALSE,
               container = w,
               expand=TRUE) 
  
  # Help button group.
  gh <- ggroup(container = gv, expand=FALSE, fill="both")
  
  savegui_chk <- gcheckbox(text="Save GUI settings", checked=FALSE, container=gh)
  
  addSpring(gh)
  
  help_btn <- gbutton(text="Help", container=gh)
  
  addHandlerChanged(help_btn, handler = function(h, ...) {
    
    # Open help page for function.
    print(help("editData_gui", help_type="html"))
    
  })
  
  # FRAME 0 ###################################################################
  
  f0 <- gframe(text = "Datasets",
               horizontal=FALSE,
               spacing = 5,
               container = gv) 
  
  g0 <- glayout(container = f0, spacing = 1)
  
  g0[1,1] <- glabel(text="Select dataset:", container=g0)
  
  g0[1,2] <- dataset_drp <- gdroplist(items=c("<Select data frame>",
                                              listObjects(env=env,
                                                          obj.class="data.frame")), 
                                      selected = 1,
                                      editable = FALSE,
                                      container = g0)
  
  if(!is.null(.gDataName)){
    svalue(dataset_drp) <- .gDataName
  }
  
  g0[1,3] <- g0_samples_lbl <- glabel(text=" 0 samples,", container=g0)
  g0[1,4] <- g0_columns_lbl <- glabel(text=" 0 columns,", container=g0)
  g0[1,5] <- g0_rows_lbl <- glabel(text=" 0 rows", container=g0)
  
  addHandlerChanged(dataset_drp, handler = function (h, ...) {
    
    val_obj <- svalue(dataset_drp)
    
    # Check if suitable.
    ok <- checkDataset(name=val_obj, reqcol=NULL,
                       env=env, parent=w, debug=debug)
    
    if(ok){
      
      # Load or change components.
      .gData <<- get(val_obj, envir=env)
      .gDataName <<- val_obj
      
      if("Sample.Name" %in% names(.gData)){
        samples <- length(unique(.gData$Sample.Name))
      } else  if("Sample.File.Name" %in% names(.gData)){
        samples <- length(unique(.gData$Sample.File.Name))
      } else if(any(grepl("SAMPLE", names(.gData), ignore.case=TRUE))) {
        # Get (first) column name containing "Sample".
        sampleCol <- names(.gData)[grep("SAMPLE", names(.gData), ignore.case=TRUE)[1]]
        # Grab sample names.
        samples <- length(unique(.gData[, sampleCol]))
      } else {
        samples <- "<NA>"
      }
      svalue(g0_samples_lbl) <- paste(" ", samples, "samples,")
      svalue(g0_columns_lbl) <- paste(" ", ncol(.gData), "columns,")
      svalue(g0_rows_lbl) <- paste(" ", nrow(.gData), "rows")
      
      # Refresh info and load table.
      .refreshTbl()
      
    } else {
      svalue(g0_samples_lbl) <- paste(" ", "<NA>", "samples,")
      svalue(g0_columns_lbl) <- paste(" ", "<NA>", "columns,")
      svalue(g0_rows_lbl) <- paste(" ", "<NA>", "rows")
    }
    
  } )  
  
  # FRAME 1 ###################################################################
  
  f1 <- gframe(text = "Options",
               horizontal=FALSE,
               spacing = 5,
               container = gv) 
  
  g1 <- glayout(container = f1, spacing = 1)
  
  g1[1,1] <- f1_show_attr_chk <- gcheckbox(text="Show attributes (separate window)",
                                           checked=FALSE, container=g1)
  
  
  g1[2,1] <- f1_limit_chk <- gcheckbox(text="Limit number of rows to:",
                                       checked=TRUE, container=g1)
  
  
  g1[2,2] <- f1_max_edt <- gedit(text=100, width=8, container=g1)
  
  addHandlerChanged(f1_show_attr_chk, handler = function(h, ...) {

    if(svalue(f1_show_attr_chk)){
      .showAttributes()
    } else {
      if(isExtant(w_attributes)){
        visible(w_attributes) <- FALSE
      }
    }
    
  } )

  addHandlerChanged(f1_limit_chk, handler = function(h, ...) {
    
    .refreshTbl()

  } )
  
  
  # FRAME 2 ###################################################################
  
  f2 <- gframe(text = "Copy | Export | Save",
               horizontal=FALSE,
               spacing = 5,
               container = gv) 
  
  g2 <- glayout(container = f2, spacing = 1)
  
  g2[1,1] <- copy_btn <- gbutton(text="Copy", border=TRUE, container=g2)
  tooltip(copy_btn) <- "Copy to clipboard (NB! large datasets might get truncated)"
  
  g2[1,2] <- export_btn <- gbutton(text="Export", border=TRUE, container=g2)
  tooltip(export_btn) <- "Opens the export dialog"
  
  g2[1,3] <- save_btn <- gbutton(text="Save as", border=TRUE, container=g2)
  tooltip(save_btn) <- "Save as new dataset"
  
  g2[1,4] <- save_txt <- gedit(text=.gDataName,
                               container=g2, width = 50,
                               anchor=c(-1 ,0))
  
  addHandlerChanged(copy_btn, handler = function(h, ...) {
    
    val_tbl <- data_tbl[]
    
    # Change button.
    svalue(copy_btn) <- "Copying..."
    enabled(copy_btn) <- FALSE
    
    # Copy data.
    write.table(val_tbl, "clipboard",
                sep="\t", row.names=FALSE)
    
    # Change button.
    svalue(copy_btn) <- "Copy"
    enabled(copy_btn) <- TRUE
    
  } )
  
  addHandlerChanged(save_btn, handler = function(h, ...) {
    
    val_name <- svalue(save_txt)
    datanew <- data_tbl[]
    
    if(debug){
      print("names(datanew)")
      print(names(datanew))
    }
    
    if (!is.na(val_name) && !is.null(val_name)){
      
      # Copy and add attributes (retain names).
      names(.gData) <- names(datanew)
      attributes(datanew) <- attributes(.gData)
      
      # Change button.
      svalue(save_btn) <- "Saving..."
      enabled(save_btn) <- FALSE
      
      # Save data.
      saveObject(name=val_name, object=datanew, parent=w, env=env, debug=debug)
      
      # Change button.
      svalue(save_btn) <- "Save as"
      enabled(save_btn) <- TRUE
      
    } else {
      
      gmessage(message="A name must be given!",
               title="Error",
               icon = "error")
      
    } 
    
  } )
  
  addHandlerChanged(export_btn, handler = function(h, ...) {
    
    # Open GUI.
    export_gui(env=env, savegui=savegui, debug=debug, parent=parent)
    
  } )
  
  # FRAME 3 ###################################################################
  
  f3 <- gframe(text = "Data frame",
               horizontal=FALSE,
               spacing = 5,
               expand=TRUE,
               container = gv) 
  
  if(is.null(.gData)){
    data_tbl <- gWidgets::gtable(items=data.frame(Data="There is no data"),
                                 container=f3, expand=TRUE)
  } else {
    
    # Get options.
    val_limit <- svalue(f1_limit_chk)
    val_max <- as.numeric(svalue(f1_max_edt))
    val_attr <- svalue(f1_show_attr_chk) 
    
    # Update info.
    if("Sample.Name" %in% names(.gData)){
      samples <- length(unique(.gData$Sample.Name))
    } else if("Sample.File.Name" %in% names(.gData)){
      samples <- length(unique(.gData$Sample.File.Name))
    } else if(any(grepl("SAMPLE", names(.gData), ignore.case=TRUE))) {
      # Get (first) column name containing "Sample".
      sampleCol <- names(.gData)[grep("SAMPLE", names(.gData), ignore.case=TRUE)[1]]
      # Grab sample names.
      samples <- length(unique(.gData[, sampleCol]))
    } else {
      samples <- "<NA>"
    }
    svalue(g0_samples_lbl) <- paste(" ", samples, "samples,")
    svalue(g0_columns_lbl) <- paste(" ", ncol(.gData), "columns,")
    svalue(g0_rows_lbl) <- paste(" ", nrow(.gData), "rows")
    
    # Load data.
    if(edit){
      if(val_limit){
        data_tbl <<- gdf(items=head(.gData, val_max),
                         container=f3, expand=TRUE)
      } else {
        data_tbl <<- gdf(items=.gData, container=f3, expand=TRUE)
      }
    } else {
      if(val_limit){
        data_tbl <<- gWidgets::gtable(items=head(.gData, val_max),
                                      container=f3, expand=TRUE)
      } else {
        data_tbl <<- gWidgets::gtable(items=.gData, container=f3, expand=TRUE)
      }
    }
    
  }
  
  # INTERNAL FUNCTIONS ########################################################
  
  .showAttributes <- function(){
    
    if(debug){
      print(paste("IN:", match.call()[[1]]))
    }
    
    # Get options.
    val_attr <- svalue(f1_show_attr_chk)
    
    if(val_attr & !is.null(.gData)){
      
      if(!isExtant(w_attributes)){

        # Re-create window.        
        w_attributes <<- gwindow(title="Attributes", visible=FALSE)
        attr_text <<- gtext("", container=w_attributes)

      }

      # Get list of attributes.
      attributeList <- attributes(.gData)
      
      # Remove common non-strvalidator attributes (too much to show).
      attributeList$names <- NULL
      attributeList$row.names <- NULL
      attributeList$class <- NULL
      
      # Empty text fiels.
      svalue(attr_text) <- ""

      # Get names of attributes.      
      attrNames <- names(attributeList)

      # Loop over list of attributes and att to text object.
      for(a in seq(along=attrNames)){
        
        # Insert text for current attribute.        
        insert(attr_text, paste(attrNames[a], attributeList[a], sep=": "))
        
      }

      # Show window.
      visible(w_attributes) <- TRUE  
      
    }
    
  }
  
  .refreshTbl <- function(){
    
    if(debug){
      print(paste("IN:", match.call()[[1]]))
    }
    
    # Get options.
    val_limit <- svalue(f1_limit_chk)
    val_max <- as.numeric(svalue(f1_max_edt))
    val_attr <- svalue(f1_show_attr_chk) 
    
    if(!is.null(.gData)){

      if(val_attr){
        .showAttributes()
      }
      
      # Update "save as" with current dataset name.
      svalue(save_txt) <- paste(.gDataName, "_edit", sep="")
      
      # Refresh widget by removing it and...
      delete(f3, data_tbl)
      
      visible(f3) <- FALSE
      # ...creating a new table.
      if(edit){
        if(val_limit){
          data_tbl <<- gdf(items=head(.gData, val_max),
                           container=f3, expand=TRUE)
        } else {
          data_tbl <<- gdf(items=.gData, container=f3, expand=TRUE)
        }
      } else {
        if(val_limit){
          data_tbl <<- gWidgets::gtable(items=head(.gData, val_max),
                                        container=f3, expand=TRUE)
        } else {
          data_tbl <<- gWidgets::gtable(items=.gData, container=f3, expand=TRUE)
        }
      }
      
      visible(f3) <- TRUE
      
    }
    
  }
  
  .loadSavedSettings <- function(){
    
    # First check status of save flag.
    if(!is.null(savegui)){
      svalue(savegui_chk) <- savegui
      enabled(savegui_chk) <- FALSE
      if(debug){
        print("Save GUI status set!")
      }  
    } else {
      # Load save flag.
      if(exists(".strvalidator_editData_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_editData_gui_savegui", envir=env)
      }
      if(debug){
        print("Save GUI status loaded!")
      }
    }
        
    if(debug){
      print(svalue(savegui_chk))
    }

    # Then load settings if true.
    if(svalue(savegui_chk)){
      if(exists(".strvalidator_editData_gui_attr", envir=env, inherits = FALSE)){
        svalue(f1_show_attr_chk) <- get(".strvalidator_editData_gui_attr", envir=env)
      }
      if(exists(".strvalidator_editData_gui_limit", envir=env, inherits = FALSE)){
        svalue(f1_limit_chk) <- get(".strvalidator_editData_gui_limit", envir=env)
      }
      if(exists(".strvalidator_editData_gui_maxrow", envir=env, inherits = FALSE)){
        svalue(f1_max_edt) <- get(".strvalidator_editData_gui_maxrow", envir=env)
      }
      
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
      
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_editData_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_editData_gui_attr", value=svalue(f1_show_attr_chk), envir=env)
      assign(x=".strvalidator_editData_gui_limit", value=svalue(f1_limit_chk), envir=env)
      assign(x=".strvalidator_editData_gui_maxrow", value=svalue(f1_max_edt), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_editData_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_editData_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_editData_gui_attr", envir=env, inherits = FALSE)){
        remove(".strvalidator_editData_gui_attr", envir = env)
      }
      if(exists(".strvalidator_editData_gui_limit", envir=env, inherits = FALSE)){
        remove(".strvalidator_editData_gui_limit", envir = env)
      }
      if(exists(".strvalidator_editData_gui_maxrow", envir=env, inherits = FALSE)){
        remove(".strvalidator_editData_gui_maxrow", envir = env)
      }
      
      if(debug){
        print("Settings cleared!")
      }
    }
    
    if(debug){
      print("Settings saved!")
    }
    
  }
  
  # END GUI ###################################################################
  
  # Load GUI settings.
  .loadSavedSettings()
    
  # Show GUI.
  visible(w) <- TRUE
  focus(w)
  
}
