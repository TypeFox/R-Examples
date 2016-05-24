################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 16.12.2015: Implemented option to type a reference sample name.
# 15.12.2015: Implemented 'exact' option.
# 04.12.2015: Removed 'Marker' from required columns.
# 28.08.2015: Added importFrom.
# 05.05.2015: Changed parameter 'ignoreCase' to 'ignore.case' for 'checkSubset' function.
# 11.10.2014: Added 'focus', added 'parent' parameter.
# 28.06.2014: Added help button and moved save gui checkbox.
# 08.05.2014: Implemented 'checkDataset'.
# 25.07.2013: Parameter 'fixed' changed to 'word'.
# 15.07.2013: Added save GUI settings.
# 15.07.2013: Added 'options' group.
# 11.06.2013: Added 'inherits=FALSE' to 'exists'.
# 04.06.2013: Fixed bug in 'missingCol'.
# 24.05.2013: Improved error message for missing columns.
# 17.05.2013: listDataFrames() -> listObjects()
# 27.04.2013: First version.


#' @title Check Subset
#'
#' @description
#' GUI wrapper for the \code{\link{checkSubset}} function.
#'
#' @details
#' Simplifies the use of the \code{\link{checkSubset}} function by providing
#' a graphical user interface to it.
#' 
#' @param env environment in wich to search for data frames.
#' @param savegui logical indicating if GUI settings should be saved in the environment.
#' @param debug logical indicating printing debug information.
#' @param parent widget to get focus when finished.
#' 
#' @export
#' 
#' @importFrom utils help
#' @importFrom graphics title
#' 
#' @return TRUE
#' 
#' @seealso \code{\link{checkSubset}}

checkSubset_gui <- function(env=parent.frame(), savegui=NULL, debug=FALSE, parent=NULL){
  
  # Global variables.
  .gData <- NULL
  .gRef <- NULL

  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # Main window.
  w <- gwindow(title="Check subsetting", visible=FALSE)
  
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
               spacing=15,
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
    print(help("checkSubset_gui", help_type="html"))
    
  })
  
  # DATASET ###################################################################
  
  f0 <- gframe(text = "Datasets",
               horizontal=FALSE,
               spacing = 5,
               container = gv) 

  f0g1 <- glayout(container = f0, spacing = 1)
  
  f0g1[1,1] <- glabel(text="Select dataset:", container=f0g1)
  
  f0g1[1,2] <- dataset_drp <- gdroplist(items=c("<Select dataset>",
                                                 listObjects(env=env,
                                                             obj.class="data.frame")),
                                         selected = 1,
                                         editable = FALSE,
                                         container = f0g1)
  
  f0g1[1,3] <- dataset_samples_lbl <- glabel(text=" 0 samples",
                                              container=f0g1)
  
  addHandlerChanged(dataset_drp, handler = function (h, ...) {
    
    val_obj <- svalue(dataset_drp)
    
    # Check if suitable.
    requiredCol <- c("Sample.Name")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       env=env, parent=w, debug=debug)
    
    if(ok){
      
      # Load or change components.
      .gData <<- get(val_obj, envir=env)
      samples <- length(unique(.gData$Sample.Name))
      svalue(dataset_samples_lbl) <- paste("", samples, "samples")
        
    } else {
      
      # Reset components.
      .gData <<- data.frame(No.Data=NA)
      svalue(dataset_samples_lbl) <- " 0 samples"
      
    }
    
  } )

  f0g1[2,1] <- glabel(text="Select reference set:", container=f0g1)
  
  f0g1[2,2] <- dataset_ref_drp <- gdroplist(items=c("<Select dataset>",
                                                 listObjects(env=env,
                                                             obj.class="data.frame")),
                                         selected = 1,
                                         editable = FALSE,
                                         container = f0g1)
  
  f0g1[2,3] <- dataset_ref_lbl <- glabel(text=" 0 reference samples",
                                              container=f0g1)

  
  f0g1[3,1] <- glabel(text="Or type a reference name:", container=f0g1)
  
  f0g1[3,2] <- dataset_ref_edt <- gedit(container = f0g1)
  
  addHandlerChanged(dataset_ref_drp, handler = function (h, ...) {
    
    val_obj <- svalue(dataset_ref_drp)
    
    # Check if suitable.
    requiredCol <- c("Sample.Name")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       env=env, parent=w, debug=debug)
    
    if(ok){
      
      # Load or change components.
      .gRef <<- get(val_obj, envir=env)
      refs <- length(unique(.gRef$Sample.Name))
      svalue(dataset_ref_lbl) <- paste("", refs, "reference samples")
      
    } else {
      
      # Reset components.
      .gRef <<- data.frame(No.Data=NA)
      svalue(dataset_ref_lbl) <- " 0 reference samples"
      
    }
    
  } )
  
  # FRAME 1 ###################################################################
  
  f1 <- gframe(text = "Options",
               horizontal=FALSE,
               spacing = 5,
               container = gv) 
  
  f1_ignore_case_chk <- gcheckbox(text="Ignore case ('A' will match 'A', 'B-a.2', and 'A2')",
                                  checked = TRUE,
                                  container = f1)

  f1_word_chk <- gcheckbox(text="Add word boundaries ('A' will match 'A', 'B-A.2', and 'A 2' but not 'A2')",
                                  checked = FALSE,
                                  container = f1)
  
  f1_exact_chk <- gcheckbox(text="Exact matching ('A' will match 'A' but not 'B-A.2', 'A 2', or 'A2')",
                           checked = FALSE,
                           container = f1)
  
  # BUTTON ####################################################################

  if(debug){
    print("BUTTON")
  }  
  
  check_btn <- gbutton(text="Subset",
                      border=TRUE,
                      container=gv)
  
  addHandlerChanged(check_btn, handler = function(h, ...) {
    
    # Get values.
    val_data <- .gData
    val_ref <- .gRef
    val_ref_name <- svalue(dataset_ref_edt)
    val_ignore <- svalue(f1_ignore_case_chk)
    val_word <- svalue(f1_word_chk)
    val_exact <- svalue(f1_exact_chk)
    
    if(is.null(.gRef)){
      # If no reference dataset use given name.
      val_ref <- val_ref_name
    }

    # Check that data is available.
    if (!is.null(val_data) && !is.null(val_ref)){
      
      chksubset_w <- gwindow(title = "Check subsetting",
                             visible = FALSE, name=title,
                             width = NULL, height= NULL, parent=w,
                             handler = NULL, action = NULL)
      
      chksubset_txt <- checkSubset(data=val_data,
                                   ref=val_ref,
                                   ignore.case=val_ignore,
                                   word=val_word,
                                   exact=val_exact,
                                   console=FALSE)
      
      gtext (text = chksubset_txt, width = NULL, height = 300, font.attr = NULL, 
             wrap = FALSE, container = chksubset_w)
      
      visible(chksubset_w) <- TRUE
      
    } else {
      
      gmessage(message="Data frame is NULL!\n\n
               Make sure to select a dataset and a reference set or type a reference name",
               title="Error",
               icon = "error")      
      
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
      if(exists(".strvalidator_checkSubset_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_checkSubset_gui_savegui", envir=env)
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
      if(exists(".strvalidator_checkSubset_gui_ignore", envir=env, inherits = FALSE)){
        svalue(f1_ignore_case_chk) <- get(".strvalidator_checkSubset_gui_ignore", envir=env)
      }
      if(exists(".strvalidator_checkSubset_gui_word", envir=env, inherits = FALSE)){
        svalue(f1_word_chk) <- get(".strvalidator_checkSubset_gui_word", envir=env)
      }
      if(exists(".strvalidator_checkSubset_gui_exact", envir=env, inherits = FALSE)){
        svalue(f1_exact_chk) <- get(".strvalidator_checkSubset_gui_exact", envir=env)
      }
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_checkSubset_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_checkSubset_gui_ignore", value=svalue(f1_ignore_case_chk), envir=env)
      assign(x=".strvalidator_checkSubset_gui_word", value=svalue(f1_word_chk), envir=env)
      assign(x=".strvalidator_checkSubset_gui_exact", value=svalue(f1_word_chk), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_checkSubset_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_checkSubset_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_checkSubset_gui_ignore", envir=env, inherits = FALSE)){
        remove(".strvalidator_checkSubset_gui_ignore", envir = env)
      }
      if(exists(".strvalidator_checkSubset_gui_word", envir=env, inherits = FALSE)){
        remove(".strvalidator_checkSubset_gui_word", envir = env)
      }
      if(exists(".strvalidator_checkSubset_gui_exact", envir=env, inherits = FALSE)){
        remove(".strvalidator_checkSubset_gui_exact", envir = env)
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
  
} # End of GUI
