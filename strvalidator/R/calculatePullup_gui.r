################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 28.08.2015: Added importFrom.
# 05.05.2015: Changed parameter 'ignoreCase' to 'ignore.case' for 'checkSubset' function.
# 13.12.2014: Added kit dropdown and kit attribute to result.
# 04.12.2014: First version.

#' @title Calculate Spectral Pull-up
#'
#' @description
#' GUI wrapper for the \code{\link{calculatePullup}} function.
#'
#' @details
#' Simplifies the use of the \code{\link{calculatePullup}} function by
#' providing a graphical user interface.
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
#' @seealso \code{\link{calculatePullup}}, \code{\link{checkSubset}}
#' 


calculatePullup_gui <- function(env=parent.frame(), savegui=NULL,
                                 debug=FALSE, parent=NULL){
  
  # Global variables.
  .gData <- NULL
  .gRef <- NULL
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # WINDOW ####################################################################
  
  if(debug){
    print("WINDOW")
  }  

  # Main window.
  w <- gwindow(title="Calculate spectral pull-up", visible=FALSE)

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
    print(help("calculatePullup_gui", help_type="html"))
    
  })
  
  # FRAME 0 ###################################################################
  
  if(debug){
    print("FRAME 0")
  }  
  
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
    requiredCol <- c("Sample.Name", "Allele", "Marker", "Dye", "Height", "Size", "Data.Point")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       slim=TRUE, slimcol="Height",
                       env=env, parent=w, debug=debug)
    
    if(ok){
      # Load or change components.
      
      # get dataset.
      .gData <<- get(val_obj, envir=env)
      svalue(g0_data_samples_lbl) <- paste(length(unique(.gData$Sample.Name)),
                                           "samples.")
      
      # Suggest a name for result.
      svalue(f4_save_edt) <- paste(val_obj, "_pullup", sep="")
      
      # Detect kit.
      kitIndex <- detectKit(.gData, index=TRUE)
      # Select in dropdown.
      svalue(f4_kit_drp, index=TRUE) <- kitIndex
      
    } else {
      
      # Reset components.
      .gData <<- NULL
      svalue(g0_data_drp, index=TRUE) <- 1
      svalue(g0_data_samples_lbl) <- " 0 samples"
      svalue(f4_save_edt) <- ""
      
    }
    
  } )  
  
  # Reference -----------------------------------------------------------------
  
  g0[2,1] <- glabel(text="Select reference dataset:", container=g0)
  
  # NB! dfs defined in previous section.
  g0[2,2] <- g0_ref_drp <- gdroplist(items=dfs, 
                                     selected = 1,
                                     editable = FALSE,
                                     container = g0)
  
  g0[2,3] <- g0_ref_samples_lbl <- glabel(text=" 0 references", container=g0)
  
  addHandlerChanged(g0_ref_drp, handler = function (h, ...) {
    
    val_obj <- svalue(g0_ref_drp)
    
    # Check if suitable.
    requiredCol <- c("Sample.Name", "Marker", "Allele")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       slim=TRUE, slimcol="Allele",
                       env=env, parent=w, debug=debug)
    
    if(ok){
      # Load or change components.
      
      .gRef <<- get(val_obj, envir=env)
      svalue(g0_ref_samples_lbl) <- paste(length(unique(.gRef$Sample.Name)),
                                          "samples.")
      
    } else {
      
      # Reset components.
      .gRef <<- NULL
      svalue(g0_ref_drp, index=TRUE) <- 1
      svalue(g0_ref_samples_lbl) <- " 0 references"
      
    }    
  } )  
  
  # CHECK ---------------------------------------------------------------------
  
  if(debug){
    print("CHECK")
  }  
  
  g0[3,2] <- g0_check_btn <- gbutton(text="Check subsetting",
                                     border=TRUE,
                                     container=g0)
  
  addHandlerChanged(g0_check_btn, handler = function(h, ...) {
    
    # Get values.
    val_data <- .gData
    val_ref <- .gRef
    val_ignore <- svalue(f1_ignore_chk)
    val_word <- svalue(f1_word_chk)
    
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
  
  # FRAME 1 ###################################################################
  
  if(debug){
    print("FRAME 1")
  }  
  
  f1 <- gframe(text = "Options",
               horizontal=FALSE,
               spacing = 10,
               container = gv) 
  
  f1_ignore_chk <- gcheckbox(text="Ignore case",
                             checked=TRUE,
                             container=f1)
  
  f1_word_chk <- gcheckbox(text="Add word boundaries",
                           checked = FALSE,
                           container = f1)

  f1_ol_chk <- gcheckbox(text="Remove off-ladder peaks",
                           checked = FALSE,
                           container = f1)
  
  # LAYOUT --------------------------------------------------------------------

  f1g1 <- glayout(container = f1, spacing = 1)
  
  f1g1[1,1] <- glabel(text="Pullup analysis range (data points) around known alleles: ", container=f1g1)
  f1g1[1,2] <- f1_pullup_spb <- gspinbutton(from=0, to=1000, by=10, value=6, container=f1g1)

  f1g1[2,1] <- glabel(text="Blocking range (data points) around known alleles: ", container=f1g1)
  f1g1[2,2] <- f1_block_spb <- gspinbutton(from=0, to=1000, by=10, value=70, container=f1g1)
  
  f1_discard_chk <- gcheckbox(text="Discard alleles with no pullup from the result table",
                         checked = FALSE,
                         container = f1)
  
  # FRAME 4 ###################################################################
  
  if(debug){
    print("FRAME 4")
  }  
  
  f4 <- gframe(text = "Save as",
               horizontal=TRUE,
               spacing = 5,
               container = gv) 
  
  glabel(text="Name for result:", container=f4)
  
  f4_save_edt <- gedit(text="", container=f4)

  glabel(text=" Kit attribute:", container=f4)
  
  f4_kit_drp <- gdroplist(items=getKit(), selected = 1,
                          editable = FALSE, container = f4) 
  
  # BUTTON ####################################################################

  if(debug){
    print("BUTTON")
  }  
  
  calculate_btn <- gbutton(text="Calculate",
                      border=TRUE,
                      container=gv)
  
  addHandlerChanged(calculate_btn, handler = function(h, ...) {
    
    # Get values.
    val_data <- .gData
    val_ref <- .gRef
    val_ignore <- svalue(f1_ignore_chk)
    val_word <- svalue(f1_word_chk)
    val_ol <- svalue(f1_ol_chk)
    val_pullup <- svalue(f1_pullup_spb)
    val_block <- svalue(f1_block_spb)
    val_discard <- svalue(f1_discard_chk)
    val_name <- svalue(f4_save_edt)
    val_kit <- svalue(f4_kit_drp)
    
    if(debug){
      print("Read Values:")
      print("val_data")
      print(head(val_data))
      print("val_ref")
      print(head(val_ref))
      print("val_ignore")
      print(val_ignore)
      print("val_word")
      print(val_word)
      print("val_ol")
      print(val_ol)
      print("val_pullup")
      print(val_pullup)
      print("val_block")
      print(val_block)
      print("val_name")
      print(val_name)
    }
    
    # Check if data.
    if(!is.null(.gData) & !is.null(.gRef)){
      
      # Check for NA's in dye column.
      if(!any(is.na(.gData$Dye))){
        
        # Change button.
        svalue(calculate_btn) <- "Processing..."
        enabled(calculate_btn) <- FALSE
        
        datanew <- calculatePullup(data=val_data,
                                   ref=val_ref,
                                   pullup.range=val_pullup,
                                   block.range=val_block,
                                   ol.rm=val_ol,
                                   ignore.case=val_ignore,
                                   word=val_word,
                                   discard=val_discard,
                                   debug=debug)
        
        # Add attribute for detected kit.
        attr(datanew, which="kit") <- val_kit
        
        # Save data.
        saveObject(name=val_name, object=datanew, parent=w, env=env)
        
        if(debug){
          print(str(datanew))
          print(head(datanew))
          print(paste("EXIT:", match.call()[[1]]))
        }
        
        # Close GUI.
        dispose(w)
        
      } else {
        
        message <- "'NA' in 'Dye' column. \nUse add dye function to fix."
        
        gmessage(message, title="NA detected!",
                 icon = "error",
                 parent = w)
        
      }
      
    } else {
      
      message <- "A dataset and a reference dataset have to be selected."
      
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
      if(exists(".strvalidator_calculatePullup_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_calculatePullup_gui_savegui", envir=env)
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
      if(exists(".strvalidator_calculatePullup_gui_window", envir=env, inherits = FALSE)){
        svalue(f1_pullup_spb) <- get(".strvalidator_calculatePullup_gui_window", envir=env)
      }
      if(exists(".strvalidator_calculatePullup_gui_block", envir=env, inherits = FALSE)){
        svalue(f1_block_spb) <- get(".strvalidator_calculatePullup_gui_block", envir=env)
      }
      if(exists(".strvalidator_calculatePullup_gui_ol", envir=env, inherits = FALSE)){
        svalue(f1_ol_chk) <- get(".strvalidator_calculatePullup_gui_ol", envir=env)
      }
      if(exists(".strvalidator_calculatePullup_gui_ignore", envir=env, inherits = FALSE)){
        svalue(f1_ignore_chk) <- get(".strvalidator_calculatePullup_gui_ignore", envir=env)
      }
      if(exists(".strvalidator_calculatePullup_gui_word", envir=env, inherits = FALSE)){
        svalue(f1_word_chk) <- get(".strvalidator_calculatePullup_gui_word", envir=env)
      }
      if(exists(".strvalidator_calculatePullup_gui_discard", envir=env, inherits = FALSE)){
        svalue(f1_discard_chk) <- get(".strvalidator_calculatePullup_gui_discard", envir=env)
      }
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_calculatePullup_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_calculatePullup_gui_window", value=svalue(f1_pullup_spb), envir=env)
      assign(x=".strvalidator_calculatePullup_gui_block", value=svalue(f1_block_spb), envir=env)
      assign(x=".strvalidator_calculatePullup_gui_ol", value=svalue(f1_ol_chk), envir=env)
      assign(x=".strvalidator_calculatePullup_gui_ignore", value=svalue(f1_ignore_chk), envir=env)
      assign(x=".strvalidator_calculatePullup_gui_word", value=svalue(f1_word_chk), envir=env)
      assign(x=".strvalidator_calculatePullup_gui_discard", value=svalue(f1_discard_chk), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_calculatePullup_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculatePullup_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_calculatePullup_gui_window", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculatePullup_gui_window", envir = env)
      }
      if(exists(".strvalidator_calculatePullup_gui_block", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculatePullup_gui_block", envir = env)
      }
      if(exists(".strvalidator_calculatePullup_gui_ol", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculatePullup_gui_ol", envir = env)
      }
      if(exists(".strvalidator_calculatePullup_gui_ignore", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculatePullup_gui_ignore", envir = env)
      }
      if(exists(".strvalidator_calculatePullup_gui_word", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculatePullup_gui_word", envir = env)
      }
      if(exists(".strvalidator_calculatePullup_gui_discard", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculatePullup_gui_discard", envir = env)
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
