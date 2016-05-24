################################################################################
# TODO LIST
# TODO: Make a general function to add any (selected) kit information?
# TODO: Option to add all uppercase/lower case marker names?

################################################################################
# CHANGE LOG (last 20 changes)
# 09.01.2016: Added more attributes to result.
# 28.08.2015: Added importFrom
# 11.10.2014: Added 'focus', added 'parent' parameter.
# 28.06.2014: Added help button and moved save gui checkbox.
# 06.05.2014: Implemented 'checkDataset'.
# 27.11.2013: First version.


#' @title Add Missing Markers
#'
#' @description
#' GUI wrapper for the \code{\link{addMarker}} function.
#'
#' @details
#' Simplifies the use of the \code{\link{addMarker}} function by providing a graphical 
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
#' @importFrom utils head help
#' 
#' @seealso \code{\link{addMarker}}

addMarker_gui <- function(env=parent.frame(), savegui=NULL, debug=FALSE, parent=NULL){
  
  # Global variables.
  .gData <- data.frame(No.Data=NA)
  .gDataName <- NULL
  .gKit <- 1

  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print(head(.gData))
  }
  
  # Main window.  
  w <- gwindow(title="Add missing markers to dataset", visible=FALSE)

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
    print(help("addMarker_gui", help_type="html"))
    
  })
  
  # DATASET ###################################################################
  
  f0 <- gframe(text = "Dataset and kit",
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
    requiredCol <- c("Sample.Name", "Marker")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       env=env, parent=w, debug=debug)
    
    if(ok){
      # Load or change components.

      # Get data.
      .gData <<- get(val_obj, envir=env)
      
      # Load or change components.
      .gDataName <<- val_obj
      samples <- length(unique(.gData$Sample.Name))
      svalue(dataset_samples_lbl) <- paste(" ", samples, "samples")
      .gKit <<- detectKit(.gData, index=TRUE)
      svalue(kit_drp, index=TRUE) <- .gKit
      svalue(f2_save_edt) <- paste(.gDataName, "_marker", sep="")
      
      if(debug){
        print("Detected kit index")
        print(.gKit)
      }
        
    } else {
      
      # Reset components.
      .gData <<- data.frame(No.Data=NA)
      .gDataName <<- NULL
      svalue(dataset_samples_lbl) <- " 0 samples"
      svalue(f2_save_edt) <- ""
      
    }
    
  } )
  
  # KIT -----------------------------------------------------------------------
  
  f0g1[2,1] <- glabel(text="Kit:", container=f0g1)
  
  kit_drp <- gdroplist(items=getKit(),
                           selected = 1,
                           editable = FALSE,
                           container = f0g1)
  
  f0g1[2,2] <- kit_drp
  
  # FRAME 1 ###################################################################

  # OPTIONS -----------------------------------------------------------------------
  
  f1 <- gframe(text = "Options",
               horizontal=FALSE,
               spacing = 5,
               container = gv) 
  
  f1_ignore_chk <- gcheckbox(text="Ignore case (marker name)",
                             checked=FALSE, container=f1)

  # FRAME 2 ###################################################################
  
  f2 <- gframe(text = "Save as",
               horizontal=TRUE,
               spacing = 5,
               container = gv) 
  
  glabel(text="Name for result:", container=f2)
  
  f2_save_edt <- gedit(text="", container=f2)

  # BUTTON ####################################################################

  if(debug){
    print("BUTTON")
  }  
  
  add_btn <- gbutton(text="Add marker",
                      border=TRUE,
                      container=gv)
  
  addHandlerChanged(add_btn, handler = function(h, ...) {
    
    # Get values.
    val_kit <- svalue(kit_drp)
    val_data <- .gData
    val_data_name <- .gDataName
    val_ignore <- svalue(f1_ignore_chk)
    val_name <- svalue(f2_save_edt)
    
    if(debug){
      print(".gData")
      print(names(.gData))
      print("val_kit")
      print(val_kit)
    }

    # Change button.
    svalue(add_btn) <- "Processing..."
    enabled(add_btn) <- FALSE
    
    datanew <- addMarker(data=val_data,
                         marker=getKit(val_kit, what="Marker"),
                         ignore.case=val_ignore, debug=debug)
    
    # Add attributes.
    attr(datanew, which="addMarker_gui, data") <- val_data_name
    attr(datanew, which="addMarker_gui, kit") <- val_kit
    attr(datanew, which="addMarker_gui, ignore.case") <- val_ignore

    
    # Save data.
    saveObject(name=val_name, object=datanew, parent=w, env=env)
    
    # Close GUI.
    dispose(w)
    
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
      if(exists(".strvalidator_addMarker_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_addMarker_gui_savegui", envir=env)
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
      if(exists(".strvalidator_addMarker_gui_ignore", envir=env, inherits = FALSE)){
        svalue(f1_ignore_chk) <- get(".strvalidator_addMarker_gui_ignore", envir=env)
      }
      
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_addMarker_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_addMarker_gui_ignore", value=svalue(f1_ignore_chk), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_addMarker_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_addMarker_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_addMarker_gui_ignore", envir=env, inherits = FALSE)){
        remove(".strvalidator_addMarker_gui_ignore", envir = env)
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
