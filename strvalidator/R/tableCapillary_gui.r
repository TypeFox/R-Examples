################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 29.08.2015: Added importFrom.
# 11.10.2014: Added 'focus', added 'parent' parameter.
# 10.08.2014: Added scope=RUN.
# 28.06.2014: Added help button and moved save gui checkbox.
# 08.05.2014: Implemented 'checkDataset'.
# 29.10.2013: First version.


#' @title Table Capillary
#'
#' @description
#' GUI wrapper for the \code{\link{tableCapillary}} function.
#'
#' @details
#' Simplifies the use of the \code{\link{tableCapillary}} function by providing a graphical 
#' user interface to it.
#' 
#' @param env environment in wich to search for data frames.
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
#' @seealso \code{\link{tableCapillary}}

tableCapillary_gui <- function(env=parent.frame(), savegui=NULL, debug=FALSE, parent=NULL){
  
  # Global variables.
  .gData <- NULL
  .gDataName <- NULL
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # Main window.
  w <- gwindow(title="Make capillary balance table", visible=FALSE)

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
    print(help("tableCapillary_gui", help_type="html"))
    
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
  
  f0g0[1,3] <- f0g0_samples_lbl <- glabel(text=" 0 rows",
                                              container=f0g0)
  
  addHandlerChanged(f0g0_dataset_drp, handler = function (h, ...) {
    
    val_obj <- svalue(f0g0_dataset_drp)
    
    # Check if suitable.
    requiredCol <- c("Instrument", "Capillary", "Injection", "Well","Mean.Height")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       env=env, parent=w, debug=debug)
    
    if(ok){
      
      # Load or change components.
      .gData <<- get(val_obj, envir=env)
      .gDataName <<- val_obj
      svalue(f0g0_samples_lbl) <- paste(" ", nrow(.gData), "rows")
      svalue(f2_save_edt) <- paste(.gDataName,
                                   "_table_",
                                   svalue(f1g1_scope_opt),
                                   sep="")
        
    } else {

      # Reset components.
      .gData <<- data.frame(No.Data=NA)
      .gDataName <<- NULL
      svalue(f0g0_samples_lbl) <- " 0 rows"
      svalue(f2_save_edt) <- ""
      svalue(f0g0_dataset_drp, index=TRUE) <- 1
      
    }
    
  } )
  
  # FRAME 1 ###################################################################
  
  f1 <- gframe(text="Options",
                   horizontal=FALSE,
                   spacing = 20,
                   container = gv) 
  
  f1g1 <- glayout(container = f1, spacing = 5)
  
  f1g1[1,1] <- glabel(text="Make table by:", container=f1g1)
  
  f1g1[2,1] <- f1g1_scope_opt <- gradio(items=c("capillary","injection",
                                                "plate row", "run",
                                                "instrument"),
                              selected = 1,
                              horizontal = FALSE,
                              container = f1g1)

  addHandlerChanged(f1g1_scope_opt, handler = function (h, ...) {

    svalue(f2_save_edt) <- paste(.gDataName,
                                 "_table_",
                                 svalue(f1g1_scope_opt),
                                 sep="")
    
  })
  
  # FRAME 2 ###################################################################
  
  f2 <- gframe(text = "Save as",
               horizontal=TRUE,
               spacing = 5,
               container = gv) 
  
  glabel(text="Name for result:", container=f2)
  
  f2_save_edt <- gedit(text="", width=45, container=f2)

  # BUTTON ####################################################################

  if(debug){
    print("BUTTON")
  }  
  
  run_btn <- gbutton(text="Make table",
                      border=TRUE,
                      container=gv)
  
  addHandlerChanged(run_btn, handler = function(h, ...) {
    
    # Get values.
    val_data <- .gData
    val_scope <- svalue(f1g1_scope_opt)
    val_name <- svalue(f2_save_edt)

    # Translate into valid scope strings.
    if(val_scope == "capillary"){
      val_scope <- "cap"
    } else if(val_scope == "injection"){
      val_scope <- "inj"
    } else if(val_scope == "plate row"){
      val_scope <- "row"
    } else if(val_scope == "run"){
      val_scope <- "run"
    } else if(val_scope == "instrument"){
      val_scope <- "instr"
    } 
    
    if (!is.null(.gData)){
      
      # Change button.
      svalue(run_btn) <- "Processing..."
      enabled(run_btn) <- FALSE
      
      datanew <- tableCapillary(data=val_data,
                                scope=val_scope,
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

    # Set check state if provided.
    if(!is.null(savegui)){
      svalue(savegui_chk) <- savegui
      enabled(savegui_chk) <- FALSE
      if(debug){
        print("Save GUI status set!")
      }  
    } else {
      # Load save flag.
      if(exists(".strvalidator_tableCapillary_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_tableCapillary_gui_savegui", envir=env)
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
      if(exists(".strvalidator_tableCapillary_gui_scope", envir=env, inherits = FALSE)){
        svalue(f1g1_scope_opt) <- get(".strvalidator_tableCapillary_gui_scope", envir=env)
      }
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){

    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_tableCapillary_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_tableCapillary_gui_scope", value=svalue(f1g1_scope_opt), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_tableCapillary_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_tableCapillary_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_tableCapillary_gui_scope", envir=env, inherits = FALSE)){
        remove(".strvalidator_tableCapillary_gui_scope", envir = env)
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
