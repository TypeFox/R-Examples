################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 12.10.2015: First version.

#' @title Detect Spike
#'
#' @description
#' GUI wrapper for the \code{\link{calculateSpike}} function.
#'
#' @details Simplifies the use of the \code{\link{calculateSpike}} function
#' by providing a graphical user interface.
#' @param env environment in wich to search for data frames.
#' @param savegui logical indicating if GUI settings should be saved in the environment.
#' @param debug logical indicating printing debug information.
#' @param parent widget to get focus when finished.
#' 
#' @return TRUE
#' 
#' @export
#' 
#' @seealso \code{\link{calculateSpike}}

calculateSpike_gui <- function(env = parent.frame(), savegui = NULL, debug = FALSE, parent = NULL){
  
  .gData <- NULL
  .gDataName <- NULL
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  w <- gwindow(title = "Detect spikes", visible = FALSE)
  
  # Runs when window is closed.
  addHandlerDestroy(w, handler = function (h, ...) {
    
    # Save GUI state.
    .saveSettings()
    
    # Focus on parent window.
    if(!is.null(parent)){
      focus(parent)
    }
    
  })
  
  gv <- ggroup(horizontal = FALSE,
               spacing = 8,
               use.scrollwindow = FALSE,
               container = w,
               expand = TRUE) 
  
  # Help button group.
  gh <- ggroup(container = gv, expand = FALSE, fill = "both")
  
  savegui_chk <- gcheckbox(text = "Save GUI settings", checked = FALSE, container = gh)
  
  addSpring(gh)
  
  help_btn <- gbutton(text = "Help", container = gh)
  
  addHandlerChanged(help_btn, handler = function(h, ...) {
    
    # Open help page for function.
    print(help("calculateSpike_gui", help_type = "html"))
    
  })
  
  # FRAME 0 ###################################################################
  
  f0 <- gframe(text = "Dataset",
               horizontal = FALSE,
               spacing = 5,
               container = gv) 
  
  g0 <- glayout(container = f0, spacing = 1)
  
  # Datasets ------------------------------------------------------------------
  
  g0[1,1] <- glabel(text = "Select dataset:", container = g0)
  
  g0[1,2] <- g0_dataset_drp <- gdroplist(items = c("<Select dataset>",
                                                 listObjects(env = env,
                                                             obj.class = "data.frame")), 
                                         selected = 1,
                                         editable = FALSE,
                                         container = g0)
  
  g0[1,3] <- g0_samples_lbl <- glabel(text = " 0 samples", container = g0)
  
  g0[1,4] <- glabel(text = " and the kit used:", container = g0)
  
  g0[1,5] <- kit_drp <- gdroplist(items = getKit(), 
                       selected = 1,
                       editable = FALSE,
                       container = g0) 
  
  addHandlerChanged(g0_dataset_drp, handler = function (h, ...) {
    
    val_obj <- svalue(g0_dataset_drp)
    
    # Check if suitable.
    requiredCol <- c("Sample.Name", "File.Name", "Size")
    ok <- checkDataset(name = val_obj, reqcol = requiredCol,
                       slim = TRUE, slimcol = "Size",
                       env = env, parent = w, debug = debug)
    
    if(ok){
      
      # Load or change components.
      .gData <<- get(val_obj, envir = env)
      .gDataName <<- val_obj
      samples <- length(unique(.gData$Sample.Name))
      svalue(g0_samples_lbl) <- paste("", samples, "samples")
      svalue(f2_save_edt) <- paste(.gDataName, "_spikes", sep = "")
      
      # Detect kit.
      kitIndex <- detectKit(.gData, index=TRUE)
      # Select in dropdown.
      svalue(kit_drp, index = TRUE) <- kitIndex

    } else {
      
      # Reset components.
      .gData <<- NULL
      svalue(g0_dataset_drp, index = TRUE) <- 1
      svalue(g0_samples_lbl) <- " 0 samples"
      svalue(f2_save_edt) <- ""
      
    }
    
  } )
  
  # FRAME 1 ###################################################################
  
  f1 <- gframe(text = "Options",
               horizontal = FALSE,
               spacing = 5,
               container = gv) 
  
  f1g1 <- glayout(container=f1, spacing=1)

  f1g1[1,1] <- glabel(text = "Threshold (number of peaks at similar size):",
                      container = f1g1)  
  f1g1[1,2] <- f1_t_spn <- gspinbutton(from = 1, to = 10, by = 1, value = 3,
                                       container = f1g1)
  
  f1g1[2,1] <- glabel(text = "Round to nearest (bp):", container = f1g1)  
  f1g1[2,2] <- f1_r_spn <- gspinbutton(from = 1, to = 10, by = 1, value = 1,
                                       container = f1g1)
  
  # FRAME 2 ###################################################################
  
  f2 <- gframe(text = "Save as", horizontal = TRUE, spacing = 5, container = gv) 
  
  glabel(text="Name for result:", container = f2)
  
  f2_save_edt <- gedit(text = "", container = f2, expand = TRUE)
  
  # BUTTON ####################################################################
  
  
  calculate_btn <- gbutton(text = "Detect", border = TRUE, container = gv)
  
  addHandlerChanged(calculate_btn, handler = function(h, ...) {
    
    # Get values.
    val_data <- .gData
    val_t <- svalue(f1_t_spn)
    val_r <- svalue(f1_r_spn)
    val_kit <- svalue(kit_drp)
    val_name <- svalue(f2_save_edt)
    
    if(!is.null(val_data)){
      
      # Change button.
      svalue(calculate_btn) <- "Processing..."
      enabled(calculate_btn) <- FALSE
      
      datanew <- calculateSpike(data=val_data,
                                threshold=val_t,
                                round.to=val_r,
                                kit=val_kit,
                                debug=debug)
      
      # Add attributes to result.
      attr(datanew, which="calculateSpike_gui, data") <- svalue(g0_dataset_drp)
      attr(datanew, which="calculateSpike_gui, threshold") <- val_t
      attr(datanew, which="calculateSpike_gui, round.to") <- val_r
      attr(datanew, which="calculateSpike_gui, kit") <- val_kit

      # Save data.
      saveObject(name = val_name, object = datanew, parent = w, env = env)
      
      if(debug){
        print(str(datanew))
        print(paste("EXIT:", match.call()[[1]]))
      }
      
      # Close GUI.
      dispose(w)
      
    } else {
      
      message <- "A dataset must be selected."
      
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
      if(exists(".strvalidator_calculateSpike_gui_savegui", envir = env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_calculateSpike_gui_savegui", envir = env)
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
      if(exists(".strvalidator_calculateSpike_gui_t", envir = env, inherits = FALSE)){
        svalue(f1_t_spn) <- get(".strvalidator_calculateSpike_gui_t", envir = env)
      }
      if(exists(".strvalidator_calculateSpike_gui_r", envir = env, inherits = FALSE)){
        svalue(f1_r_spn) <- get(".strvalidator_calculateSpike_gui_r", envir = env)
      }
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x = ".strvalidator_calculateSpike_gui_savegui",
             value = svalue(savegui_chk), envir = env)
      assign(x = ".strvalidator_calculateSpike_gui_t",
             value = svalue(f1_t_spn), envir = env)
      assign(x = ".strvalidator_calculateSpike_gui_r",
             value = svalue(f1_r_spn), envir = env)

    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_calculateSpike_gui_savegui", envir = env, inherits = FALSE)){
        remove(".strvalidator_calculateSpike_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_calculateSpike_gui_t", envir = env, inherits = FALSE)){
        remove(".strvalidator_calculateSpike_gui_t", envir = env)
      }
      if(exists(".strvalidator_calculateSpike_gui_r", envir = env, inherits = FALSE)){
        remove(".strvalidator_calculateSpike_gui_r", envir = env)
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