################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG
# 11.10.2014: Added 'focus', added 'parent' parameter.
# 28.06.2014: Added help button and moved save gui checkbox.
# 08.05.2014: Implemented 'checkDataset'.
# 23.02.2014: Implemented new option 'ol.rm'.
# 18.07.2013: Check before overwrite object.
# 11.06.2013: Added 'inherits=FALSE' to 'exists'.
# 10.06.2013: New parameter 'savegui'.
# 06.06.2013: Added save GUI settings.
# 04.06.2013: Fixed bug in 'missingCol'.
# 29.05.2013: Disabled button and adding "processing..." after press.
# 24.05.2013: Improved error message for missing columns.
# 21.05.2013: Fixed name on save as.
# 17.05.2013: listDataFrames() -> listObjects()
# 09.05.2013: .result removed, added save as group.
# 04.05.2013: First version.


#' @title Guess Profile
#'
#' @description
#' GUI wrapper for the \code{\link{guessProfile}} function.
#'
#' @details
#' Simplifies the use of the \code{\link{guessProfile}} function by providing
#' a graphical user interface to it.
#' 
#' @param env environment in wich to search for data frames.
#' @param savegui logical indicating if GUI settings should be saved in the environment.
#' @param debug logical indicating printing debug information.
#' @param parent widget to get focus when finished.
#' 
#' @export
#' 
#' @return TRUE
#' 
#' @seealso \code{\link{guessProfile}}, \code{\link{checkSubset}}

guessProfile_gui <- function(env=parent.frame(), savegui=NULL, debug=FALSE, parent=NULL){
  
  # Global variables.
  .gData <- NULL
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  w <- gwindow(title="Guess profile", visible=FALSE)

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
    print(help("guessProfile_gui", help_type="html"))
    
  })
  
  # FRAME 0 ###################################################################
  
  f0 <- gframe(text="Datasets",
                   horizontal=FALSE,
                   spacing = 10,
                   container = gv) 
  
  f0g0 <- glayout(container = f0, spacing = 1)
  
  f0g0[1,1] <- glabel(text="Select dataset:", container=f0g0)
  
  f0g0[1,2] <- f0g0_dataset_drp <- gdroplist(items=c("<Select dataset>",
                                                 listObjects(env=env,
                                                             obj.class="data.frame")),
                                         selected = 1,
                                         editable = FALSE,
                                         container = f0g0)
  
  f0g0[1,3] <- f0g0_samples_lbl <- glabel(text=" 0 samples",
                                              container=f0g0)
  
  addHandlerChanged(f0g0_dataset_drp, handler = function (h, ...) {
    
    val_obj <- svalue(f0g0_dataset_drp)
    
    # Check if suitable.
    requiredCol <- c("Sample.Name", "Marker","Allele","Height")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       env=env, parent=w, debug=debug)
    
    if(ok){
      
      # Load or change components.
      .gData <<- get(val_obj, envir=env)
      samples <- length(unique(.gData$Sample.Name))
      svalue(f0g0_samples_lbl) <- paste(" ", samples, "samples")
      svalue(f2_save_edt) <- paste(val_obj, "_profile", sep="")
        
    } else {

      # Reset components.
      .gData <<- data.frame(No.Data=NA)
      svalue(f0g0_samples_lbl) <- " 0 samples"
      svalue(f2_save_edt) <- ""
      svalue(f0g0_dataset_drp, index=TRUE) <- 1
      
    }
    
  } )
  
  # FRAME 1 ###################################################################
  
  f1 <- gframe(text="Options",
                   horizontal=FALSE,
                   spacing = 10,
                   container = gv) 
  
  f1g1 <- glayout(container = f1, spacing = 5)
  
  f1g1[1,1] <- glabel(text="Accepted ratio >=", container=f1g1)

  f1g1[1,2] <- f1g1_ratio_spb <- gspinbutton(from = 0, to = 1,
                                            by = 0.01, value = 0.6,
                                            container = f1g1)

  f1g1[2,1] <- glabel(text="Accepted peak height >=", container=f1g1)
  
  f1g1[2,2] <- f1g1_height_edt <- gedit(text="100", width=6, container=f1g1)
  
  f1g1[3,1] <- f1g1_na_chk <- gcheckbox(text="Discard NA rows",
                                   checked=FALSE,
                                   container=f1g1)

  f1g1[4,1] <- f1g1_ol_chk <- gcheckbox(text="Ignore off-ladder (OL) alleles",
                                        checked=FALSE,
                                        container=f1g1)
  
  # FRAME 2 ###################################################################
  
  f2 <- gframe(text = "Save as",
               horizontal=TRUE,
               spacing = 5,
               container = gv) 
  
  glabel(text="Name for result:", container=f2)
  
  f2_save_edt <- gedit(text="", width=25, container=f2)

  # BUTTON ####################################################################

  if(debug){
    print("BUTTON")
  }  
  
  check_btn <- gbutton(text="Guess",
                      border=TRUE,
                      container=gv)
  
  addHandlerChanged(check_btn, handler = function(h, ...) {
    
    # Get values.
    val_data <- .gData
    val_ratio <- as.numeric(svalue(f1g1_ratio_spb))
    val_height <- as.numeric(svalue(f1g1_height_edt))
    val_NA <- svalue(f1g1_na_chk)
    val_OL <- svalue(f1g1_ol_chk)
    val_name <- svalue(f2_save_edt)
    
    if(is.na(val_height)){
      val_height <- 0
    }
    
    if (!is.null(.gData)){
      
      # Change button.
      svalue(check_btn) <- "Processing..."
      enabled(check_btn) <- FALSE
      
      datanew <- guessProfile(data=val_data,
                              ratio=val_ratio,
                              height=val_height,
                              na.rm=val_NA,
                              ol.rm=val_OL,
                              debug=debug)
      
      # Save data.
      saveObject(name=val_name, object=datanew, parent=w, env=env)
      
      if(debug){
        print(datanew)
        print(paste("EXIT:", match.call()[[1]]))
      }
      
      # Close GUI.
      dispose(w)
      
    } else {
      
      gmessage(message="Data frame is NULL!\n\n
               Make sure to select a dataset and a reference set",
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
      if(exists(".strvalidator_guessProfile_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_guessProfile_gui_savegui", envir=env)
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
      if(exists(".strvalidator_guessProfile_gui_ratio", envir=env, inherits = FALSE)){
        svalue(f1g1_ratio_spb) <- get(".strvalidator_guessProfile_gui_ratio", envir=env)
      }
      if(exists(".strvalidator_guessProfile_gui_height", envir=env, inherits = FALSE)){
        svalue(f1g1_height_edt) <- get(".strvalidator_guessProfile_gui_height", envir=env)
      }
      if(exists(".strvalidator_guessProfile_gui_na", envir=env, inherits = FALSE)){
        svalue(f1g1_na_chk) <- get(".strvalidator_guessProfile_gui_na", envir=env)
      }
      if(exists(".strvalidator_guessProfile_gui_ol", envir=env, inherits = FALSE)){
        svalue(f1g1_ol_chk) <- get(".strvalidator_guessProfile_gui_ol", envir=env)
      }
      if(debug){
        print("Saved settings loaded!")
      }
    }
      
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_guessProfile_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_guessProfile_gui_ratio", value=svalue(f1g1_ratio_spb), envir=env)
      assign(x=".strvalidator_guessProfile_gui_height", value=svalue(f1g1_height_edt), envir=env)
      assign(x=".strvalidator_guessProfile_gui_na", value=svalue(f1g1_na_chk), envir=env)
      assign(x=".strvalidator_guessProfile_gui_ol", value=svalue(f1g1_ol_chk), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_guessProfile_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_guessProfile_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_guessProfile_gui_ratio", envir=env, inherits = FALSE)){
        remove(".strvalidator_guessProfile_gui_ratio", envir = env)
      }
      if(exists(".strvalidator_guessProfile_gui_height", envir=env, inherits = FALSE)){
        remove(".strvalidator_guessProfile_gui_height", envir = env)
      }
      if(exists(".strvalidator_guessProfile_gui_na", envir=env, inherits = FALSE)){
        remove(".strvalidator_guessProfile_gui_na", envir = env)
      }
      if(exists(".strvalidator_guessProfile_gui_ol", envir=env, inherits = FALSE)){
        remove(".strvalidator_guessProfile_gui_ol", envir = env)
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
