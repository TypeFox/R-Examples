################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 26.10.2015: Added attributes.
# 28.08.2015: Added importFrom.
# 05.05.2015: Changed parameter 'ignoreCase' to 'ignore.case' for 'checkSubset' function.
# 05.01.2015: Added kit dropdown and kit attribute to result.
# 07.10.2014: Added 'focus', added 'parent' parameter.
# 03.08.2014: Added detection of kit and add attribute to result.
# 28.06.2014: Added help button and moved save gui checkbox.
# 08.05.2014: Implemented 'checkDataset'.
# 26.07.2013: Changed parameter 'fixed' to 'word' for 'checkSubset' function.
# 18.07.2013: Check before overwrite object.
# 17.07.2013: Added save GUI settings.
# 17.07.2013: 'false' allele checkboxes replaced by gdf table.
# 17.07.2013: Added check subsetting.
# 11.06.2013: Fixed wrong interference passed to 'calculateStutter' (-1).
# 11.06.2013: 'val_replace' and 'val_by' set to NULL if length = 0.
# 04.06.2013: Fixed bug in 'missingCol'.
# 30.05.2013: Added replace 'false' stutters.
# 29.05.2013: Added subset check.
# 28.05.2013: Added warning for additive effects.
# 24.05.2013: Improved error message for missing columns.

#' @title Calculate Stutter
#'
#' @description
#' GUI wrapper for the \code{\link{calculateStutter}} function.
#'
#' @details
#' Simplifies the use of the \code{\link{calculateStutter}} function by providing 
#' a graphical user interface to it.
#' 
#' @param env environment in wich to search for data frames and save result.
#' @param savegui logical indicating if GUI settings should be saved in the environment.
#' @param debug logical indicating printing debug information.
#' @param parent widget to get focus when finished.
#' 
#' @export
#' 
#' @importFrom utils help head
#' @importFrom graphics title
#' 
#' @return TRUE
#' 
#' @seealso \code{\link{calculateStutter}}, \code{\link{checkSubset}}

calculateStutter_gui <- function(env=parent.frame(), savegui=NULL, debug=FALSE, parent=NULL){
  
  # Global variables.
  .gData <- NULL
  .gRef <- NULL
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # Main window.
  w <- gwindow(title="Calculate stutter ratio", visible=FALSE)
  
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
              expand=FALSE) 
  
  # Help button group.
  gh <- ggroup(container = gv, expand=FALSE, fill="both")
  
  savegui_chk <- gcheckbox(text="Save GUI settings", checked=FALSE, container=gh)
  
  addSpring(gh)
  
  help_btn <- gbutton(text="Help", container=gh)
  
  addHandlerChanged(help_btn, handler = function(h, ...) {
    
    # Open help page for function.
    print(help("calculateStutter_gui", help_type="html"))
    
  })
  
  # FRAME 0 ###################################################################
  
  f0 <- gframe(text = "Datasets",
               horizontal=FALSE,
               spacing = 5,
               container = gv) 
  
  f0g0 <- glayout(container = f0, spacing = 1)
  
  # Datasets ------------------------------------------------------------------
  
  f0g0[1,1] <- glabel(text="Select dataset:", container=f0g0)
  
  f0g0[1,2] <- f0_dataset_drp <- gdroplist(items=c("<Select dataset>",
                                              listObjects(env=env,
                                                          obj.class="data.frame")), 
                                      selected = 1,
                                      editable = FALSE,
                                      container = f0g0)
  
  f0g0[1,3] <- f0_samples_lbl <- glabel(text=" 0 samples", container=f0g0)
  
  addHandlerChanged(f0_dataset_drp, handler = function (h, ...) {
    
    val_obj <- svalue(f0_dataset_drp)
    
    # Check if suitable.
    requiredCol <- c("Sample.Name", "Marker", "Allele", "Height")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       env=env, parent=w, debug=debug)
    
    if(ok){
      
      # Load or change components.
      .gData <<- get(val_obj, envir=env)
      samples <- length(unique(.gData$Sample.Name))
      svalue(f0_samples_lbl) <- paste("", samples, "samples")
      svalue(f2_save_edt) <- paste(val_obj, "_stutter", sep="")
      # Detect kit.
      kitIndex <- detectKit(.gData, index=TRUE)
      # Select in dropdown.
      svalue(f2_kit_drp, index=TRUE) <- kitIndex
        
    } else {

      # Reset components.
      .gData <<- NULL
      svalue(f0_dataset_drp, index=TRUE) <- 1
      svalue(f0_samples_lbl) <- " 0 samples"
      svalue(f2_save_edt) <- ""
      
    }
    
  } )  
  
  f0g0[2,1] <- glabel(text="Select reference dataset:", container=f0g0)
  
  f0g0[2,2] <- f0_refset_drp <- gdroplist(items=c("<Select dataset>",
                                             listObjects(env=env,
                                                         obj.class="data.frame")), 
                                     selected = 1,
                                     editable = FALSE,
                                     container = f0g0) 
  
  f0g0[2,3] <- f0_ref_lbl <- glabel(text=" 0 references", container=f0g0)
  
  addHandlerChanged(f0_refset_drp, handler = function (h, ...) {
    
    val_obj <- svalue(f0_refset_drp)
    
    # Check if suitable.
    requiredCol <- c("Sample.Name", "Marker", "Allele")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       env=env, parent=w, debug=debug)
    
    if(ok){
      
      # Load or change components.
      .gRef <<- get(val_obj, envir=env)
      refs <- length(unique(.gRef$Sample.Name))
      svalue(f0_ref_lbl) <- paste("", refs, "references")
        
    } else {
      
      # Reset components.
      .gRef <<- NULL
      svalue(f0_refset_drp, index=TRUE) <- 1
      svalue(f0_ref_lbl) <- " 0 references"
      
    }
    
  } )  

  # CHECK ---------------------------------------------------------------------
  
  if(debug){
    print("CHECK")
  }  
  
  f0g0[3,2] <- f0_check_btn <- gbutton(text="Check subsetting",
                                  border=TRUE,
                                  container=f0g0)
  
  addHandlerChanged(f0_check_btn, handler = function(h, ...) {
    
    # Get values.
    val_data <- .gData
    val_ref <- .gRef
    
    if (!is.null(.gData) || !is.null(.gRef)){
      
      chksubset_w <- gwindow(title = "Check subsetting",
                             visible = FALSE, name=title,
                             width = NULL, height= NULL, parent=w,
                             handler = NULL, action = NULL)
      
      chksubset_txt <- checkSubset(data=val_data,
                                   ref=val_ref,
                                   console=FALSE,
                                   ignore.case=TRUE,
                                   word=FALSE)
      
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
  
  # FRAME 1 ###################################################################
  
  f1 <- gframe("Options", horizontal=FALSE, container=gv)
  
  glabel(text="Calculate stutter ratio within the the following analysis range:",
                      container=f1, anchor=c(-1 ,0))

  f1g1 <- ggroup(horizontal=TRUE,
               spacing=5,
               use.scrollwindow=FALSE,
               container = f1,
               expand=FALSE) 
  
  f1g1_range_b_spb <- gspinbutton (from = 0, to = 3, by = 1,
                            value = 1, digits = 0,
                            container = f1g1) 

  glabel(text="backward stutters to", container=f1g1, anchor=c(-1 ,0))

  f1g1_range_f_spb <- gspinbutton (from = 0, to = 3, by = 1,
                              value = 1, digits = 0,
                              container = f1g1) 
  
  glabel(text="forward stutters.", container=f1g1, anchor=c(-1 ,0))

  glabel(text="NB! Additive effects outside the analysis range cannot be controlled.",
         container=f1, anchor=c(-1 ,0))
  glabel(text="A narrow range like 0 to +1 can be greately affected by neighbouring -1 stutters.",
         container=f1, anchor=c(-1 ,0))
  
  interference_f <- gframe("Level of interference within the given range",
                           horizontal=FALSE, container=gv)
  
  options <- c("no overlap between stutters and alleles",
               "stutter-stutter interference allowed",
               "stutter-allele interference allowed")
  
  interference_opt <- gradio(items=options,
                       selected=1,
                       horizontal=FALSE,
                       container=interference_f)
  
  # FRAME 3 ###################################################################
  
  f3 <- gframe("Replace 'false' stutters", horizontal=FALSE, expand=TRUE, container=gv)
  
  f3g1 <- ggroup(horizontal=FALSE,
                 spacing=5,
                 use.scrollwindow=FALSE,
                 container = f3,
                 expand=TRUE) 

  # Create default data frame.
  f3_replace_val <- c(-1.9, -1.8, -1.7, -0.9, -0.8, -0.7, 0.9, 0.8, 0.7)
  f3_by_val <- c(-1.3, -1.2, -1.1, -0.3, -0.2, -0.1, 0.3, 0.2, 0.1)
  f3_default <- data.frame(False.Stutter=f3_replace_val,
                           True.Stutter=f3_by_val,
                           Replace=TRUE,
                           stringsAsFactors=FALSE)
  
  f3_default_gdf <- gdf(items=f3_default, container=f3g1, expand=TRUE)

  # Set initial minimum size.
  size(f3_default_gdf) <- c(100,100)
  
  # FRAME 2 ###################################################################
  
  f2 <- gframe(text = "Save as",
               horizontal=TRUE,
               spacing = 5,
               container = gv) 
  
  glabel(text="Name for result:", container=f2)
  
  f2_save_edt <- gedit(text="", container=f2)
  
  glabel(text=" Kit attribute:", container=f2)
  
  f2_kit_drp <- gdroplist(items=getKit(), selected = 1,
                          editable = FALSE, container = f2) 
  

  # BUTTON ####################################################################

  if(debug){
    print("BUTTON")
  }  
  
  calculate_btn <- gbutton(text="Calculate",
                      border=TRUE,
                      container=gv)
  
  addHandlerChanged(calculate_btn, handler = function(h, ...) {
    
    # Get values.
    val_back <- svalue(f1g1_range_b_spb)
    val_forward <- svalue(f1g1_range_f_spb)
    val_interference <- svalue(interference_opt, index=TRUE) - 1 # NB! range [0-2]
    val_data <- .gData
    val_ref <- .gRef
    val_name <- svalue(f2_save_edt)
    val_replace_df <- f3_default_gdf[]
    val_chk <- val_replace_df$Replace
    val_replace <- val_replace_df$False.Stutter
    val_by <- val_replace_df$True.Stutter
    val_kit <- svalue(f2_kit_drp)
    
    # Get selected values.
    val_replace <- val_replace[val_chk]
    val_by <- val_by[val_chk]
    if(length(val_replace) == 0){
      val_replace <- NULL      
    }
    if(length(val_by) == 0){
      val_by <- NULL      
    }
    
    if(!is.null(val_data) & !is.null(val_ref)){
        
      if(debug){
        print("val_data")
        print(head(val_data))
        print("val_ref")
        print(head(val_ref))
        print("val_back")
        print(val_back)
        print("val_forward")
        print(val_forward)
        print("val_interference")
        print(val_interference)
        print("val_replace")
        print(val_replace)
        print("val_by")
        print(val_by)
      }
      
      # Change button.
      svalue(calculate_btn) <- "Processing..."
      enabled(calculate_btn) <- FALSE
      
      # Calculate stutter.
      datanew <- calculateStutter(data=val_data, ref=val_ref, 
                                  back=val_back, forward=val_forward,
                                  interference=val_interference,
                                  replace.val=val_replace,
                                  by.val=val_by,
                                  debug=debug)
      
      # Add attributes.
      attr(datanew, which="kit") <- val_kit
      attr(datanew, which="calculateStutter_gui, data") <- svalue(f0_dataset_drp)
      attr(datanew, which="calculateStutter_gui, ref") <- svalue(f0_refset_drp)
      attr(datanew, which="calculateStutter_gui, back") <- val_back
      attr(datanew, which="calculateStutter_gui, forward") <- val_forward
      attr(datanew, which="calculateStutter_gui, interference") <- val_interference
      attr(datanew, which="calculateStutter_gui, replace.val") <- val_replace
      attr(datanew, which="calculateStutter_gui, by.val") <- val_by

      # Save data.
      saveObject(name=val_name, object=datanew, parent=w, env=env)
      
      if(debug){
        print("datanew")
        print(head(datanew))
        print(paste("EXIT:", match.call()[[1]]))
      }
      
      # Close GUI.
      dispose(w)

    } else {
      
      message <- "A dataset and a reference dataset must be selected."
      
      gmessage(message, title="Datasets not selected",
               icon = "error",
               parent = w) 
      
    }
    
  } )

  # INTERNAL FUNCTIONS ########################################################
  
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
      if(exists(".strvalidator_calculateStutter_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_calculateStutter_gui_savegui", envir=env)
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
      if(exists(".strvalidator_calculateStutter_gui_back", envir=env, inherits = FALSE)){
        svalue(f1g1_range_b_spb) <- get(".strvalidator_calculateStutter_gui_back", envir=env)
      }
      if(exists(".strvalidator_calculateStutter_gui_forward", envir=env, inherits = FALSE)){
        svalue(f1g1_range_f_spb) <- get(".strvalidator_calculateStutter_gui_forward", envir=env)
      }
      if(exists(".strvalidator_calculateStutter_gui_interference", envir=env, inherits = FALSE)){
        svalue(interference_opt) <- get(".strvalidator_calculateStutter_gui_interference", envir=env)
      }
      if(exists(".strvalidator_calculateStutter_gui_replace", envir=env, inherits = FALSE)){
        f3_default_gdf[,] <- get(".strvalidator_calculateStutter_gui_replace", envir=env)
      }
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_calculateStutter_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_calculateStutter_gui_back", value=svalue(f1g1_range_b_spb), envir=env)
      assign(x=".strvalidator_calculateStutter_gui_forward", value=svalue(f1g1_range_f_spb), envir=env)
      assign(x=".strvalidator_calculateStutter_gui_interference", value=svalue(interference_opt), envir=env)
      assign(x=".strvalidator_calculateStutter_gui_replace", value=f3_default_gdf[], envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_calculateStutter_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateStutter_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_calculateStutter_gui_back", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateStutter_gui_back", envir = env)
      }
      if(exists(".strvalidator_calculateStutter_gui_forward", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateStutter_gui_forward", envir = env)
      }
      if(exists(".strvalidator_calculateStutter_gui_interference", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateStutter_gui_interference", envir = env)
      }
      if(exists(".strvalidator_calculateStutter_gui_replace", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateStutter_gui_replace", envir = env)
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
