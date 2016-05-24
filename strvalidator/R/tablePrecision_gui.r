################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 29.08.2015: Added importFrom.
# 05.05.2015: Changed parameter 'ignoreCase' to 'ignore.case' for 'checkSubset' function.
# 11.10.2014: Added 'focus', added 'parent' parameter.
# 28.06.2014: Added help button and moved save gui checkbox.
# 08.05.2014: Implemented 'checkDataset'.
# 06.02.2014: Fixed button locks when error.
# 06.02.2014: Changed name calculatePrecision_gui -> tablePrecision_gui
# 12.01.2014: Replaced 'subset' with native code.
# 08.12.2013: First version.

#' @title Table Precision
#'
#' @description
#' GUI wrapper for the \code{\link{tablePrecision}} function.
#'
#' @details
#' Simplifies the use of the \code{\link{tablePrecision}} function by providing 
#' a graphical user interface.
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
#' @importFrom utils help head str
#' @importFrom graphics title
#' 
#' @seealso \code{\link{tablePrecision}}, \code{\link{checkSubset}}

tablePrecision_gui <- function(env=parent.frame(), savegui=NULL,
                                 debug=FALSE, parent=NULL){
  
  # Global variables.
  .gData <- data.frame(Columns="NA")
  .gRef <- NULL
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # WINDOW ####################################################################
  
  if(debug){
    print("WINDOW")
  }  

  # Main window.
  w <- gwindow(title="Calculate summary statistics for precision", visible=FALSE)

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
    print(help("tablePrecision_gui", help_type="html"))
    
  })
  
  # FRAME 0 ###################################################################
  
  f0 <- gframe(text = "Datasets",
               horizontal=TRUE,
               spacing = 5,
               container = gv) 
  
  g0 <- glayout(container = f0, spacing = 1)

  # Dataset -------------------------------------------------------------------
  
  g0[1,1] <- glabel(text="Select dataset:", container=g0)
  
  dfs <- c("<Select a dataset>", listObjects(env=env, obj.class="data.frame"))
  
  g0[1,2] <- g0_data_drp <- gdroplist(items=dfs, 
                           selected = 1,
                           editable = FALSE,
                           container = g0)
  g0[1,3] <- g0_data_samples_lbl <- glabel(text=" 0 samples", container=g0)
  
  addHandlerChanged(g0_data_drp, handler = function (h, ...) {
    
    val_obj <- svalue(g0_data_drp)
    
    # Check if suitable.
    requiredCol <- c("Sample.Name", "Marker")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       env=env, parent=w, debug=debug)
    
    if(ok){
      
      # Load or change components.
      .gData <<- get(val_obj, envir=env)

      .refresh_key_tbl()
      .refresh_target_tbl()

      svalue(g0_data_samples_lbl) <- paste(length(unique(.gData$Sample.Name)),
                                        "samples.")
      svalue(f4_save_edt) <- paste(val_obj, "_precision_table", sep="")
      
      # Detect kit.
      kitIndex <- detectKit(.gData, index=TRUE)
      # Select in dropdown.
      svalue(f2g2_kit_drp, index=TRUE) <- kitIndex
      
      # Enable buttons.
      enabled(calculate_btn) <- TRUE
        
    } else {

      # Reset components.
      .gData <<- data.frame(Columns="NA")
      svalue(g0_data_drp, index=TRUE) <- 1
      svalue(g0_data_samples_lbl) <- " 0 samples"
      svalue(f4_save_edt) <- ""
      .refresh_key_tbl()
      .refresh_target_tbl()
      
    }
    
  } )  

  # FRAME 2 ###################################################################
  
  f2 <- gframe(text = "Filter",
               horizontal=FALSE,
               spacing = 15,
               container = gv) 
  
  f2_options <- c("Filter by reference dataset",
                  "Filter by kit bins",
                  "Do not filter")
  
  f2_filter_opt <- gradio(items=f2_options,
                             selected=3,
                             horizontal=FALSE,
                             container=f2)
  
  addHandlerChanged(f2_filter_opt, handler = function (h, ...) {
    
    val_opt <- svalue(f2_filter_opt, index=TRUE)
    
    if(val_opt == 1){
      
      enabled(f2g1) <- TRUE
      enabled(f2g2) <- FALSE
      
    } else if(val_opt == 2){
      
      enabled(f2g1) <- FALSE
      enabled(f2g2) <- TRUE
      
    } else {
      
      enabled(f2g1) <- FALSE
      enabled(f2g2) <- FALSE
      
    }
    
  } )  
  
  
  # Reference -----------------------------------------------------------------

  f2g1 <- glayout(container = f2, spacing = 1)
  enabled(f2g1) <- FALSE
  
  f2g1[1,1] <- glabel(text="Select reference dataset:", container=f2g1)
  
  # NB! dfs defined in previous section.
  f2g1[2,1] <- f2g1_ref_drp <- gdroplist(items=dfs, 
                                     selected = 1,
                                     editable = FALSE,
                                     container = f2g1)
  
  f2g1[2,2] <- f2g1_ref_samples_lbl <- glabel(text=" 0 references", container=g0)
  
  addHandlerChanged(f2g1_ref_drp, handler = function (h, ...) {
    
    val_obj <- svalue(f2g1_ref_drp)
    
    if(exists(val_obj, envir=env, inherits = FALSE)){
      
      .gRef <<- get(val_obj, envir=env)
      
      # Check if required columns...
      requiredCol <- c("Sample.Name", "Marker", "Allele")
      slimmed <- sum(grepl("Allele",names(.gRef), fixed=TRUE)) == 1
      
      if(!all(requiredCol %in% colnames(.gRef))){
        
        missingCol <- requiredCol[!requiredCol %in% colnames(.gRef)]
        
        message <- paste("Additional columns required:\n",
                         paste(missingCol, collapse="\n"), sep="")
        
        gmessage(message, title="message",
                 icon = "error",
                 parent = w) 
        
        # Reset components.
        .gRef <<- NULL
        svalue(f2g1_ref_drp, index=TRUE) <- 1
        svalue(f2g1_ref_samples_lbl) <- " 0 references"
        
      } else if (!slimmed) {
        
        message <- paste("The dataset is too fat!\n\n",
                         "There can only be 1 'Allele' column\n",
                         "Slim the dataset in the 'EDIT' tab", sep="")
        
        gmessage(message, title="message",
                 icon = "error",
                 parent = w) 
        
        # Reset components.
        .gRef <<- NULL
        svalue(f2g1_ref_drp, index=TRUE) <- 1
        svalue(f2g1_ref_samples_lbl) <- " 0 references"
        
      }else {
        
        # Load or change components.
        svalue(f2g1_ref_samples_lbl) <- paste(length(unique(.gRef$Sample.Name)),
                                            "samples.")
        
      }
      
    } else {
      
      # Reset components.
      svalue(f2g1_ref_samples_lbl) <- ""
      .gRef <<- NULL
      
    }    
  } )  

  # CHECK ---------------------------------------------------------------------
  
  f2g1[3,1] <- f2g1_check_btn <- gbutton(text="Check subsetting",
                                  border=TRUE,
                                  container=f2g1)
  
  addHandlerChanged(f2g1_check_btn, handler = function(h, ...) {
    
    # Get values.
    val_data <- .gData
    val_ref <- .gRef
    val_ignore <- svalue(f1_ignore_chk)
    val_word <- FALSE
    
    if (!is.null(.gData) || !is.null(.gRef)){
      
      chksubset_w <- gwindow(title = "Check subsetting",
                             visible = FALSE, name=title,
                             width = NULL, height= NULL, parent=w,
                             handler = NULL, action = NULL)
      
      chksubset_txt <- checkSubset(data=val_data,
                                   ref=val_ref,
                                   console=FALSE,
                                   ignore.case=val_ignore,
                                   word=val_word)
      
      gtext (text = chksubset_txt, width = NULL, height = 300, font.attr = NULL, 
             wrap = FALSE, container = chksubset_w)
      
      visible(chksubset_w) <- TRUE
      
    } else {
      
      gmessage(message="Data frame is NULL!\n\n
               Make sure to select a dataset and a reference set",
               title="Error",
               icon = "error")      
      
    } 
    
  } )

  # Kit -------------------------------------------------------------------
  
  f2g2 <- ggroup(horizontal = TRUE, container = f2)
  enabled(f2g2) <- FALSE
  
  glabel(text="Select kit:", container=f2g2)
  
  f2g2_kit_drp <- gdroplist(items=getKit(), 
                       selected = 1,
                       editable = FALSE,
                       container = f2g2) 
  
  f2g2_virtual_chk <- gcheckbox(text="Exclude virtual bins.",
                              checked=TRUE,
                              container=f2g2)
  
  # FRAME 1 ###################################################################
  
  f1 <- gframe(text = "Options",
               horizontal=FALSE,
               spacing = 5,
               expand=TRUE,
               container = gv) 
  
  f1_ignore_chk <- gcheckbox(text="Ignore case",
                             checked=TRUE,
                             container=f1)
  
  # KEY -----------------------------------------------------------------------
  
  f1_key_f <- gframe("Create key from columns", 
                     horizontal=FALSE, 
                     container=f1, 
                     expand=TRUE)

  f1_key_txt <- gedit(initial.msg="Doubleklick or drag column names to list", 
                      width = 40,
                      container=f1_key_f)
  
  f1_key_tbl <- gWidgets::gtable(items=names(.gData), 
                                 container=f1_key_f,
                                 expand=TRUE)
  
  addDropTarget(f1_key_txt, handler=function(h,...) {
    # Get values.
    drp_val <- h$dropdata
    f1_key_val <- svalue(h$obj)
    
    # Add new value to selected.
    new <- ifelse(nchar(f1_key_val) > 0,
                  paste(f1_key_val, drp_val, sep=","),
                  drp_val)
    
    # Update text box.
    svalue(h$obj) <- new
    
    # Update column name table.
    tmp_tbl <- f1_key_tbl[,]  # Get all values.
    tmp_tbl <- tmp_tbl[tmp_tbl!=drp_val]  # Remove value added to selected.
    f1_key_tbl[,] <- tmp_tbl  # Update table.
    
  })
  
  # TARGET --------------------------------------------------------------------

  f1_target_f <- gframe("Calculate precision for target columns", 
                     horizontal=FALSE, 
                     container=f1, 
                     expand=TRUE)
  
  f1_target_txt <- gedit(initial.msg="Doubleklick or drag column names to list", 
                      width = 40,
                      container=f1_target_f)
  
  f1_target_tbl <- gWidgets::gtable(items=names(.gData), 
                                 container=f1_target_f,
                                 expand=TRUE)
  
  addDropTarget(f1_target_txt, handler=function(h,...) {
    # Get values.
    drp_val <- h$dropdata
    f1_target_val <- svalue(h$obj)
    
    # Add new value to selected.
    new <- ifelse(nchar(f1_target_val) > 0,
                  paste(f1_target_val, drp_val, sep=","),
                  drp_val)
    
    # Update text box.
    svalue(h$obj) <- new
    
    # Update column name table.
    tmp_tbl <- f1_target_tbl[,]  # Get all values.
    tmp_tbl <- tmp_tbl[tmp_tbl!=drp_val]  # Remove value added to selected.
    f1_target_tbl[,] <- tmp_tbl  # Update table.
    
  })
  
  # FRAME 4 ###################################################################
  
  f4 <- gframe(text = "Save as",
               horizontal=TRUE,
               spacing = 5,
               container = gv) 
  
  glabel(text="Name for result:", container=f4)
  
  f4_save_edt <- gedit(text="", container=f4)
  
  # BUTTON ####################################################################
  
  calculate_btn <- gbutton(text="Calculate",
                           border=TRUE,
                           container=gv)
  
  addHandlerChanged(calculate_btn, handler = function(h, ...) {
    
    # Get values.
    val_filter <- svalue(f2_filter_opt, index=TRUE)
    val_ignore <- svalue(f1_ignore_chk)
    val_key <- svalue(f1_key_txt)
    val_target <- svalue(f1_target_txt)
    val_data <- .gData
    val_ref <- .gRef
    val_name <- svalue(f4_save_edt)
    val_kit <- svalue(f2g2_kit_drp)
    val_exclude <- svalue(f2g2_virtual_chk)
    
    if(val_filter == 3){
      # Data should not be filtered. Set ref to NA (NULL gives error message.)
      val_ref <- NA
    } else if(val_filter == 2){
      # Filter by kit bins.
      
      # Get markers, bins and flag for virtual bins.
      val_ref <- getKit(kit=val_kit, what="VIRTUAL")
      
      if(val_exclude){
        # Remove virtual bins.
        val_ref <- val_ref[val_ref$Virtual == 0, ]
      }
      
    }
    
    if(debug){
      print("Read Values:")
      print("val_filter")
      print(val_filter)
      print("val_ignore")
      print(val_ignore)
      print("val_target")
      print(val_target)
      print("val_key")
      print(val_key)
      print("val_name")
      print(val_name)
      print("val_data")
      print(head(val_data))
      print("val_ref")
      print(head(val_ref))
    }

    if(!is.null(val_data) & !is.null(val_ref)){

      # Change button.
      svalue(calculate_btn) <- "Processing..."
      enabled(calculate_btn) <- FALSE

      # Filter dataset.
      if(val_filter != 3){
        val_data <- filterProfile(data=val_data, ref=val_ref,
                                add.missing.loci=FALSE, keep.na=FALSE,
                                ignore.case=val_ignore, debug=debug)
      }
      
      # Replace whitespace and split by comma.
      val_key <- gsub("\\s","", val_key)
      val_key <- strsplit(val_key, ",")
      val_key <- unlist(val_key)
      
      # Replace whitespace and split by comma.
      val_target <- gsub("\\s","", val_target)
      val_target <- strsplit(val_target, ",")
      val_target <- unlist(val_target)
      
      if(debug){
        print("Sent Values:")
        print("val_target")
        print(val_target)
        print("val_key")
        print(val_key)
        print("val_data")
        print(head(val_data))
      }
  
      # Calculate precision.      
      datanew <- tablePrecision(data=val_data,
                                    key=val_key,
                                    target=val_target,
                                    debug=debug)
      
      # Save data.
      saveObject(name=val_name, object=datanew, parent=w, env=env)
      
      if(debug){
        print(str(datanew))
        print(paste("EXIT:", match.call()[[1]]))
      }
      
      # Close GUI.
      dispose(w)
      
    } else {

      message <- "A dataset and a reference dataset have to be selected."
      
      gmessage(message, title="Datasets not selected",
               icon = "error",
               parent = w) 
      
    }
    
  } )

  # INTERNAL FUNCTIONS ########################################################

  .refresh_target_tbl <- function(){
    
    if(debug){
      print(paste("IN:", match.call()[[1]]))
    }
    
    # Refresh widget by removing it and...
    delete(f1_target_f, f1_target_tbl)
    
    # ...creating a new table.
    f1_target_tbl <<- gWidgets::gtable(items=names(.gData),
                                    container=f1_target_f,
                                    expand=TRUE)
    
    
    addDropSource(f1_target_tbl, handler=function(h,...) svalue(h$obj))
    
    addHandlerDoubleclick(f1_target_tbl, handler = function(h, ...) {
      
      # Get values.
      tbl_val <- svalue (h$obj)
      target_val <- svalue(f1_target_txt)
      
      # Add new value to selected.
      new <- ifelse(nchar(target_val) > 0,
                    paste(target_val, tbl_val, sep=","),
                    tbl_val)
      
      # Update text box.
      svalue(f1_target_txt) <- new
      
      
      # Update sample name table.
      tmp_tbl <- f1_target_tbl[,]  # Get all values.
      tmp_tbl <- tmp_tbl[tmp_tbl!=tbl_val]  # Remove value added to selected.
      f1_target_tbl[,] <- tmp_tbl  # Update table.
      
    } )
    
  }
  
  .refresh_key_tbl <- function(){
    
    if(debug){
      print(paste("IN:", match.call()[[1]]))
    }
    
    # Refresh widget by removing it and...
    delete(f1_key_f, f1_key_tbl)
    
    # ...creating a new table.
    f1_key_tbl <<- gWidgets::gtable(items=names(.gData), 
                                    container=f1_key_f,
                                    expand=TRUE)
    
    addDropSource(f1_key_tbl, handler=function(h,...) svalue(h$obj))
    
    addHandlerDoubleclick(f1_key_tbl, handler = function(h, ...) {
      
      # Get values.
      tbl_val <- svalue (h$obj)
      key_val <- svalue(f1_key_txt)
      
      # Add new value to selected.
      new <- ifelse(nchar(key_val) > 0,
                    paste(key_val, tbl_val, sep=","),
                    tbl_val)
      
      # Update text box.
      svalue(f1_key_txt) <- new
      
      # Update column name table.
      tmp_tbl <- f1_key_tbl[,]  # Get all values.
      tmp_tbl <- tmp_tbl[tmp_tbl!=tbl_val]  # Remove value added to selected.
      f1_key_tbl[,] <- tmp_tbl  # Update table.
      
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
      if(exists(".strvalidator_tablePrecision_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_tablePrecision_gui_savegui", envir=env)
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
      if(exists(".strvalidator_tablePrecision_gui_filter", envir=env, inherits = FALSE)){
        svalue(f2_filter_opt) <- get(".strvalidator_tablePrecision_gui_filter", envir=env)
      }
      if(exists(".strvalidator_tablePrecision_gui_ignore", envir=env, inherits = FALSE)){
        svalue(f1_ignore_chk) <- get(".strvalidator_tablePrecision_gui_ignore", envir=env)
      }
      if(exists(".strvalidator_tablePrecision_gui_exclude", envir=env, inherits = FALSE)){
        svalue(f2g2_virtual_chk) <- get(".strvalidator_tablePrecision_gui_exclude", envir=env)
      }
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_tablePrecision_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_tablePrecision_gui_filter", value=svalue(f2_filter_opt), envir=env)
      assign(x=".strvalidator_tablePrecision_gui_ignore", value=svalue(f1_ignore_chk), envir=env)
      assign(x=".strvalidator_tablePrecision_gui_exclude", value=svalue(f2g2_virtual_chk), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_tablePrecision_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_tablePrecision_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_tablePrecision_gui_filter", envir=env, inherits = FALSE)){
        remove(".strvalidator_tablePrecision_gui_filter", envir = env)
      }
      if(exists(".strvalidator_tablePrecision_gui_ignore", envir=env, inherits = FALSE)){
        remove(".strvalidator_tablePrecision_gui_ignore", envir = env)
      }
      if(exists(".strvalidator_tablePrecision_gui_exclude", envir=env, inherits = FALSE)){
        remove(".strvalidator_tablePrecision_gui_exclude", envir = env)
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
