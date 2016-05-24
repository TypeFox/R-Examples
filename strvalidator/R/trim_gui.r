################################################################################
# TODO LIST
# TODO: ...

# NB! Can't handle Sample.Names as factors?
################################################################################
# CHANGE LOG (last 20 changes)
# 09.01.2016: Added attributes to result.
# 29.08.2015: Added importFrom.
# 23.05.2015: Re-named internal variable 'new' (R function) to 'new_val'.
# 11.05.2015: Accepts (the first) column name containing the string 'Sample'
# as alternative to colum name 'Sample.Name'. All made case in-sensitive.
# 04.05.2015: Implemented 'checkDataset'.
# 07.10.2014: Added 'focus', added 'parent' parameter.
# 28.06.2014: Added help button and moved save gui checkbox.
# 14.01.2014: Removed requirement for column 'Sample.Name'.
# 20.11.2013: Specified package for function 'gtable' -> 'gWidgets::gtable'
# 27.10.2013: Fixed bug when 'samples'=NULL and 'invertS'=TRUE.
# 06.08.2013: Added rows and columns to info.
# 18.07.2013: Check before overwrite object.
# 16.07.2013: Added save GUI settings.
# 11.06.2013: Added 'inherits=FALSE' to 'exists'.
# 06.06.2013: Set initial table height to 200.
# 04.06.2013: Fixed bug in 'missingCol'.
# 24.05.2013: Fixed option 'replace missing values' "NA" -> NA.
# 24.05.2013: Fixed option label for 'replace missing values'.
# 24.05.2013: Improved error message for missing columns.
# 17.05.2013: listDataFrames() -> listObjects()
# 09.05.2013: .result removed, added save as group.
# 27.04.2013: Add selection of dataset in gui. Removed parameter 'data'.

#' @title Trim Data
#'
#' @description
#' GUI wrapper for the \code{\link{trim}} function.
#'
#' @details
#' Simplifies the use of the \code{\link{trim}} function by providing a graphical 
#' user interface to it.
#' 
#' @param env environment in wich to search for data frames and save result.
#' @param savegui logical indicating if GUI settings should be saved in the environment.
#' @param debug logical indicating printing debug information.
#' @param parent widget to get focus when finished.
#' 
#' @return TRUE
#' 
#' @export
#' 
#' @importFrom utils help
#' 
#' @seealso \code{\link{trim}}
#' 


trim_gui <- function(env=parent.frame(), savegui=NULL,
                     debug=FALSE, parent=NULL){

  # Global variables.
  .gData <- data.frame(Sample.Name="NA")
  .gDataName <- NULL
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # Main window.
  w <- gwindow(title="Trim dataset", visible=FALSE)
  
  # Runs when window is closed.
  addHandlerDestroy(w, handler = function (h, ...) {
    
    # Save GUI state.
    .saveSettings()
    
    # Focus on parent window.
    if(!is.null(parent)){
      focus(parent)
    }
    
  })
  
  # Vertical main group.
  gv <- ggroup(horizontal=FALSE,
              spacing=5,
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
    print(help("trim_gui", help_type="html"))
    
  })
  
  # Vertical sub group.
  g0 <- ggroup(horizontal=FALSE,
               spacing=5,
               use.scrollwindow=FALSE,
               container = gv,
               expand=FALSE) 

  # Horizontal sub group.
  g1 <- ggroup(horizontal=TRUE,
              spacing=5,
              use.scrollwindow=FALSE,
              container = gv,
              expand=TRUE) 

  # Vertical sub group.
  g2 <- ggroup(horizontal=FALSE,
              spacing=5,
              use.scrollwindow=FALSE,
              container = gv,
              expand=FALSE) 
  
  
  # DATASET ###################################################################
  
  if(debug){
    print("DATASET")
  }
  
  frame0 <- gframe(text = "Datasets",
                   horizontal=FALSE,
                   spacing = 5,
                   container = g0) 
  
  g0 <- glayout(container = frame0, spacing = 1)
  
  g0[1,1] <- glabel(text="Select dataset:", container=g0)
  
  g0[1,2] <- dataset_drp <- gdroplist(items=c("<Select dataset>",
                                                 listObjects(env=env,
                                                             obj.class="data.frame")),
                                         selected = 1,
                                         editable = FALSE,
                                         container = g0)
  
  g0[1,3] <- g0_samples_lbl <- glabel(text=" 0 samples,", container=g0)
  g0[1,4] <- g0_columns_lbl <- glabel(text=" 0 columns,", container=g0)
  g0[1,5] <- g0_rows_lbl <- glabel(text=" 0 rows", container=g0)
  
  addHandlerChanged(dataset_drp, handler = function (h, ...) {
    
    val_obj <- svalue(dataset_drp)
    
    # Check if suitable.
    requiredCol <- NULL
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       env=env, parent=w, debug=debug)
    
    if(ok){
      
      .gData <<- get(val_obj, envir=env)
      .gDataName <<- val_obj
      
      # Load or change components.
      .refresh_samples_tbl()
      .refresh_columns_tbl()
      # Info.
      if("SAMPLE.NAME" %in% toupper(names(.gData))){
        samples <- length(unique(.gData$Sample.Name))
        svalue(g0_samples_lbl) <- paste(" ", samples, "samples,")
      } else if("SAMPLE.FILE.NAME" %in% toupper(names(.gData))){
        samples <- length(unique(.gData$Sample.File.Name))
        svalue(g0_samples_lbl) <- paste(" ", samples, "samples,")
      } else if(any(grepl("SAMPLE", names(.gData), ignore.case=TRUE))){
        # Get (first) column name containing "Sample".
        sampleCol <- names(.gData)[grep("Sample", names(.gData), ignore.case=TRUE)[1]]
        samples <- length(unique(.gData[sampleCol]))
        svalue(g0_samples_lbl) <- paste(" ", samples, "samples,")
      } else {
        svalue(g0_samples_lbl) <- paste(" ", "<NA>", "samples,")
      }
      svalue(g0_columns_lbl) <- paste(" ", ncol(.gData), "columns,")
      svalue(g0_rows_lbl) <- paste(" ", nrow(.gData), "rows")
      # Result name.
      svalue(f2_save_edt) <- paste(val_obj, "_trim", sep="")
      
    } else {
      
      # Reset components.
      .gData <<- data.frame(Sample.Name="NA")
      .gDataName <<- NULL
      svalue(sample_edt) <- ""
      svalue(column_edt) <- ""
      .refresh_samples_tbl()
      .refresh_columns_tbl()
      svalue(g0_samples_lbl) <- paste(" ", "<NA>", "samples,")
      svalue(g0_columns_lbl) <- paste(" ", "<NA>", "columns,")
      svalue(g0_rows_lbl) <- paste(" ", "<NA>", "rows")
      svalue(f2_save_edt) <- ""
      
    }
  } )
  
  # SAMPLES ###################################################################
  
  if(debug){
    print("SAMPLES")
  }
  
  sample_f <- gframe("Samples", horizontal=FALSE, container=g1, expand=TRUE)
  
  
  sample_opt <- gradio(items=c("Keep","Remove"),
                       selected=1,
                       horizontal=FALSE,
                       container=sample_f)

  sample_lbl <- glabel(text="Selected samples (separate by pipe |):",
                       container=sample_f,
                       anchor=c(-1 ,0))

  sample_edt <- gedit(initial.msg="Doubleklick or drag sample names to list",
                      width = 40,
                      container=sample_f)

  if("SAMPLE.NAME" %in% toupper(names(.gData))){
    df_items <- data.frame(Samples=unique(.gData$Sample.Name),
                           stringsAsFactors=FALSE)
  } else if("SAMPLE.FILE.NAME" %in% toupper(names(.gData))){
    df_items <- data.frame(Samples=unique(.gData$Sample.File.Name),
                           stringsAsFactors=FALSE)
  } else if(any(grepl("SAMPLE", names(.gData), ignore.case=TRUE))){
    # Get (first) column name containing "Sample".
    sampleCol <- names(.gData)[grep("Sample", names(.gData), ignore.case=TRUE)[1]]
    df_items <- data.frame(Samples=unique(.gData[sampleCol]),
                           stringsAsFactors=FALSE)
  }
  sample_tbl <- gWidgets::gtable(items=df_items, container=sample_f, expand=TRUE)
  
  # Set initial size (only height is important here).
  size(sample_tbl) <- c(100,200)

  addDropTarget(sample_edt, handler=function(h,...) {
    # Get values.
    drp_val <- h$dropdata
    sample_val <- svalue(h$obj)
    
    # Add new value to selected.
    new_val <- ifelse(nchar(sample_val) > 0, paste(sample_val, drp_val, sep="|"), drp_val)
    
    # Update text box.
    svalue(h$obj) <- new_val
    
    # Update sample name table.
    tmp_tbl <- sample_tbl[,]  # Get all values.
    tmp_tbl <- tmp_tbl[tmp_tbl!=drp_val]  # Remove value added to selected.
    sample_tbl[,] <- tmp_tbl  # Update table.
    
  })

  
  # COLUMNS ###################################################################
  
  if(debug){
    print("COLUMNS")
  }
  
  column_f <- gframe("Columns", 
                     horizontal=FALSE, 
                     container=g1, 
                     expand=TRUE)
  
  column_opt <- gradio(items=c("Keep","Remove"),
                       selected=1,
                       horizontal=FALSE, 
                       container=column_f)
  
  column_lbl <- glabel(text="Selected columns (separate by pipe |):",
                       container=column_f,
                       anchor=c(-1 ,0))
  
  column_edt <- gedit(initial.msg="Doubleklick or drag column names to list", 
                      width = 40,
                      container=column_f)
  
  
  column_tbl <- gWidgets::gtable(items=names(.gData), 
                       container=column_f,
                       expand=TRUE)

  addDropTarget(column_edt, handler=function(h,...) {
    # Get values.
    drp_val <- h$dropdata
    column_val <- svalue(h$obj)
    
    # Add new value to selected.
    new_val <- ifelse(nchar(column_val) > 0,
                  paste(column_val, drp_val, sep="|"),
                  drp_val)
    
    # Update text box.
    svalue(h$obj) <- new_val
    
    # Update column name table.
    tmp_tbl <- column_tbl[,]  # Get all values.
    tmp_tbl <- tmp_tbl[tmp_tbl!=drp_val]  # Remove value added to selected.
    column_tbl[,] <- tmp_tbl  # Update table.
    
  })
  
  # OPTIONS ###################################################################
  
  if(debug){
    print("OPTIONS")
  }  
  
  option_f <- gframe("Options",
                     horizontal=FALSE, 
                     container=g1, 
                     expand=TRUE)
  
  empty_chk <- gcheckbox(text="Remove empty columns",
                         checked=TRUE,
                         container=option_f)
  
  na_chk <- gcheckbox(text="Remove NA columns",
                      checked=TRUE,
                      container=option_f)
  
  word_chk <- gcheckbox(text="Add word boundaries",
                        checked=FALSE,
                        container=option_f)
  
  case_chk <- gcheckbox(text="Ignore case",
                        checked=TRUE, 
                        container=option_f)
  
  glabel(text="Replace missing values with:",
                  container=option_f)
  na_edt <- gedit(text="NA",
                  container=option_f)

  # FRAME 2 ###################################################################
  
  if(debug){
    print("SAVE")
  }  

  f2 <- gframe(text = "Save as",
               horizontal=TRUE,
               spacing = 5,
               container = g2) 
  
  glabel(text="Name for result:", container=f2)
  
  f2_save_edt <- gedit(text="", container=f2)

  # BUTTON ####################################################################
  
  if(debug){
    print("BUTTON")
  }  
  
  trim_btn <- gbutton(text="Trim dataset",
                      border=TRUE,
                      container=g2)
  
  addHandlerChanged(trim_btn, handler = function(h, ...) {
    
    # Get new dataset name.
    val_name <- svalue(f2_save_edt)

    if(nchar(val_name) > 0) {
      
      # Get values.
      val_data <- .gData
      val_data_name <- .gDataName
      sample_val <- svalue(sample_edt)
      column_val <- svalue(column_edt)
      word_val <- svalue(word_chk)
      case_val <- svalue(case_chk)
      sample_opt_val <- if(svalue(sample_opt, index=TRUE)==1){FALSE}else{TRUE}
      column_opt_val <- if(svalue(column_opt, index=TRUE)==1){FALSE}else{TRUE}
      na_val <- svalue(na_chk)
      empty_val <- svalue(empty_chk)
      na_edt_val <- svalue(na_edt)

      # NA can't be string.
      if(na_edt_val == "NA"){
        na_edt_val <- NA
      }

      # Empty string -> NULL.
      if(sample_val == ""){
        sample_val <- NULL
      }
      
      # Empty string -> NULL.
      if(column_val == ""){
        column_val <- NULL
      }
      
      if(debug){
        print("val_data")
        print(names(val_data))
        print("sample_val")
        print(sample_val)
        print("column_val")
        print(column_val)
        print("word_val")
        print(word_val)
        print("case_val")
        print(case_val)
        print("sample_opt_val")
        print(sample_opt_val)
        print("column_opt_val")
        print(column_opt_val)
        print("na_val")
        print(na_val)
        print("empty_val")
        print(empty_val)
        print("na_edt_val")
        print(na_edt_val)
      }
  
      # Change button.
      svalue(trim_btn) <- "Processing..."
      enabled(trim_btn) <- FALSE
      
      datanew <- trim(data=val_data, samples=sample_val, columns=column_val, 
                   word=word_val, ignore.case=case_val, invert.s=sample_opt_val, invert.c=column_opt_val,
                   rm.na.col=na_val, rm.empty.col=empty_val, missing=na_edt_val, debug=debug)
      
      # Add attributes.
      attr(datanew, which="trim_gui, data") <- val_data_name
      attr(datanew, which="trim_gui, samples") <- sample_val
      attr(datanew, which="trim_gui, columns") <- column_val
      attr(datanew, which="trim_gui, word") <- word_val
      attr(datanew, which="trim_gui, ignore.case") <- case_val
      attr(datanew, which="trim_gui, invert.s") <- sample_opt_val
      attr(datanew, which="trim_gui, invert.c") <- column_opt_val
      attr(datanew, which="trim_gui, rm.na.col") <- na_val
      attr(datanew, which="trim_gui, rm.empty.col") <- empty_val
      attr(datanew, which="trim_gui, missing") <- na_edt_val
      
      # Save data.
      saveObject(name=val_name, object=datanew, parent=w, env=env)
      
      # Close GUI.
      dispose(w)
    
    } else {
      
      gmessage("A file name must be provided!", title="Error",
               icon = "error",
               parent = w) 
    }
    
  } )

  # INTERNAL FUNCTIONS ########################################################
  
  .refresh_samples_tbl <- function(){
    
    if(debug){
      print(paste("IN:", match.call()[[1]]))
    }
    
    if(any(grepl("SAMPLE", names(.gData), ignore.case=TRUE))){

      # Refresh widget by removing it and...
      delete(sample_f, sample_tbl)
      
      if("SAMPLE.NAME" %in% toupper(names(.gData))){
        df_items <- data.frame(Samples=unique(.gData$Sample.Name),
                               stringsAsFactors=FALSE)
      } else if("SAMPLE.FILE.NAME" %in% toupper(names(.gData))){
        df_items <- data.frame(Samples=unique(.gData$Sample.File.Name),
                               stringsAsFactors=FALSE)
      } else if(any(grepl("SAMPLE", names(.gData), ignore.case=TRUE))){
        # Get (first) column name containing "Sample".
        sampleCol <- names(.gData)[grep("Sample", names(.gData), ignore.case=TRUE)[1]]
        df_items <- data.frame(Samples=unique(.gData[sampleCol]),
                               stringsAsFactors=FALSE)
      }
      
      # ...creating a new table.
      sample_tbl <<- gWidgets::gtable(items=df_items, container=sample_f, expand=TRUE)
      
      addDropSource(sample_tbl, handler=function(h,...) svalue(h$obj))
      
      addHandlerDoubleclick(sample_tbl, handler = function(h, ...) {
        
        # Get values.
        tbl_val <- svalue (h$obj)
        sample_val <- svalue(sample_edt)
        
        # Add new value to selected.
        new_val <- ifelse(nchar(sample_val) > 0,
                      paste(sample_val, tbl_val, sep="|"),
                      tbl_val)
        
        # Update text box.
        svalue(sample_edt) <- new_val
        
        
        # Update sample name table.
        tmp_tbl <- sample_tbl[,]  # Get all values.
        tmp_tbl <- tmp_tbl[tmp_tbl!=tbl_val]  # Remove value added to selected.
        sample_tbl[,] <- tmp_tbl  # Update table.
        
      } )
      
    }
    
  }
  
  .refresh_columns_tbl <- function(){
    
    if(debug){
      print(paste("IN:", match.call()[[1]]))
    }
    
    # Refresh widget by removing it and...
    delete(column_f, column_tbl)
    
    # ...creating a new table.
    column_tbl <<- gWidgets::gtable(items=names(.gData), 
                         container=column_f,
                         expand=TRUE)
    
    addDropSource(column_tbl, handler=function(h,...) svalue(h$obj))
    
    addHandlerDoubleclick(column_tbl, handler = function(h, ...) {
    
      # Get values.
      tbl_val <- svalue (h$obj)
      column_val <- svalue(column_edt)
      
      # Add new value to selected.
      new_val <- ifelse(nchar(column_val) > 0,
                    paste(column_val, tbl_val, sep="|"),
                    tbl_val)
      
      # Update text box.
      svalue(column_edt) <- new_val
      
      # Update column name table.
      tmp_tbl <- column_tbl[,]  # Get all values.
      tmp_tbl <- tmp_tbl[tmp_tbl!=tbl_val]  # Remove value added to selected.
      column_tbl[,] <- tmp_tbl  # Update table.
    
    } )
    
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
      if(exists(".strvalidator_trim_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_trim_gui_savegui", envir=env)
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
      if(exists(".strvalidator_trim_gui_sample_option", envir=env, inherits = FALSE)){
        svalue(sample_opt) <- get(".strvalidator_trim_gui_sample_option", envir=env)
      }
      if(exists(".strvalidator_trim_gui_column_option", envir=env, inherits = FALSE)){
        svalue(column_opt) <- get(".strvalidator_trim_gui_column_option", envir=env)
      }
      if(exists(".strvalidator_trim_gui_remove_empty", envir=env, inherits = FALSE)){
        svalue(empty_chk) <- get(".strvalidator_trim_gui_remove_empty", envir=env)
      }
      if(exists(".strvalidator_trim_gui_remove_na", envir=env, inherits = FALSE)){
        svalue(na_chk) <- get(".strvalidator_trim_gui_remove_na", envir=env)
      }
      if(exists(".strvalidator_trim_gui_add_word", envir=env, inherits = FALSE)){
        svalue(word_chk) <- get(".strvalidator_trim_gui_add_word", envir=env)
      }
      if(exists(".strvalidator_trim_gui_ignore_case", envir=env, inherits = FALSE)){
        svalue(case_chk) <- get(".strvalidator_trim_gui_ignore_case", envir=env)
      }
      if(exists(".strvalidator_trim_gui_replace_na", envir=env, inherits = FALSE)){
        svalue(na_edt) <- get(".strvalidator_trim_gui_replace_na", envir=env)
      }
      
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_trim_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_trim_gui_sample_option", value=svalue(sample_opt), envir=env)
      assign(x=".strvalidator_trim_gui_column_option", value=svalue(column_opt), envir=env)
      assign(x=".strvalidator_trim_gui_remove_empty", value=svalue(empty_chk), envir=env)
      assign(x=".strvalidator_trim_gui_remove_na", value=svalue(na_chk), envir=env)
      assign(x=".strvalidator_trim_gui_add_word", value=svalue(word_chk), envir=env)
      assign(x=".strvalidator_trim_gui_ignore_case", value=svalue(case_chk), envir=env)
      assign(x=".strvalidator_trim_gui_replace_na", value=svalue(na_edt), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_trim_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_trim_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_trim_gui_sample_option", envir=env, inherits = FALSE)){
        remove(".strvalidator_trim_gui_sample_option", envir = env)
      }
      if(exists(".strvalidator_trim_gui_column_option", envir=env, inherits = FALSE)){
        remove(".strvalidator_trim_gui_column_option", envir = env)
      }
      if(exists(".strvalidator_trim_gui_remove_empty", envir=env, inherits = FALSE)){
        remove(".strvalidator_trim_gui_remove_empty", envir = env)
      }
      if(exists(".strvalidator_trim_gui_remove_na", envir=env, inherits = FALSE)){
        remove(".strvalidator_trim_gui_remove_na", envir = env)
      }
      if(exists(".strvalidator_trim_gui_add_word", envir=env, inherits = FALSE)){
        remove(".strvalidator_trim_gui_add_word", envir = env)
      }
      if(exists(".strvalidator_trim_gui_ignore_case", envir=env, inherits = FALSE)){
        remove(".strvalidator_trim_gui_ignore_case", envir = env)
      }
      if(exists(".strvalidator_trim_gui_replace_na", envir=env, inherits = FALSE)){
        remove(".strvalidator_trim_gui_replace_na", envir = env)
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
