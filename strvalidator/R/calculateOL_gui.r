################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 28.08.2015: Added importFrom
# 11.10.2014: Added 'focus', added 'parent' parameter.
# 28.06.2014: Added help button and moved save gui checkbox.
# 21.01.2014: Added parameter 'limit'.
# 06.01.2014: Fixed button name used as 'save as' name.
# 20.11.2013: Fixed result now stored in variable 'datanew' insted of 'val_name'.
# 29.09.2013: First version.


#' @title Analyse Off-ladder Alleles
#'
#' @description
#' GUI wrapper for the \code{\link{calculateOL}} function.
#'
#' @details By analysis of the allelic ladder the risk for getting off-ladder
#' (OL) alleles are calculated. The frequencies from a provided population
#' database is used to calculate the risk per marker and in total for the given
#' kit(s). Virtual alleles can be excluded from the calculation.
#' Small frequencies can be limited to the estimate 5/2N.
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
#' @importFrom utils help str
#' 
#' @seealso \code{\link{calculateOL}}

calculateOL_gui <- function(env=parent.frame(), savegui=NULL, debug=TRUE, parent=NULL){
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # Main window.
  w <- gwindow(title="Analyse off-ladder alleles", visible=FALSE)
  
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
    print(help("calculateOL_gui", help_type="html"))
    
  })
  
  # FRAME 0 ###################################################################
  
  f0 <- gframe(text = "Select kits",
               horizontal=TRUE,
               spacing = 5,
               container = gv) 
  
  kit_checkbox_group <- gcheckboxgroup(items=getKit(),
                                       checked = FALSE,
                                       horizontal = FALSE,
                                       container=f0) 
  
  addHandlerChanged(kit_checkbox_group, handler = function(h, ...) {
    
    val_kits <- svalue(kit_checkbox_group)
    
    if(debug){
      print("val_kits")
      print(val_kits)
    }
    
    # check if any selected kit.
    if(length(val_kits) > 0){
      
      # Enable analyse button.
      enabled(analyse_btn) <- TRUE
      
      # Suggest a save name.
      svalue(f5_save_edt) <- paste(paste(val_kits, collapse = "_"),
                                   "_OL", sep="")
      
    } else {
      
      # Disable analyse button.
      enabled(analyse_btn) <- FALSE
      
      # Empty save name.
      svalue(f5_save_edt) <- ""
      
    }
    
  } )

  # FRAME 1 ###################################################################
  
  f1 <- gframe(text = "Options",
               horizontal=FALSE,
               spacing = 5,
               container = gv) 
  
  glabel(text="Select allele frequency database:",
         anchor=c(-1 ,0), container=f1)
  
  f1_db_names <- getDb()
  
  f1_db_drp <- gdroplist(items=f1_db_names, fill=FALSE,
                         selected = 1, container=f1)

  f1_virtual_chk <- gcheckbox(text="Include virtual bins in analysis",
                              checked=TRUE,
                              container=f1)
  
  f1_msg <- paste("NB! Not all vendors specify which alleles are virtual",
                  "in the bins file.\n",
                  "This can be done manually in the kit.txt file.")
  glabel(text=f1_msg, anchor=c(-1 ,0), container=f1)
  
  f1_limit_chk <- gcheckbox(text="Limit small frequencies to 5/2N",
                              checked=TRUE,
                              container=f1)
  
  # FRAME 5 ###################################################################
  
  f5 <- gframe(text = "Save as",
               horizontal=TRUE,
               spacing = 5,
               container = gv) 
  
  glabel(text="Name for result:", container=f5)
  
  f5_save_edt <- gedit(text="", width = 50, container=f5)
  
  
  # BUTTON ####################################################################
  
  
  analyse_btn <- gbutton(text="Analyse",
                         border=TRUE,
                         container=gv)
  
  addHandlerChanged(analyse_btn, handler = function(h, ...) {
    
    val_name <- svalue(f5_save_edt)
    val_kits <- svalue(kit_checkbox_group)
    val_kitData <- data.frame() # Filled further down.
    val_db_selected <- svalue(f1_db_drp)
    val_db <- NULL  # Filled further down.
    val_virtual <- svalue(f1_virtual_chk)
    val_limit <- svalue(f1_limit_chk)
    
    if(length(val_kits) >0){
    
      # Change button.
      svalue(analyse_btn) <- "Processing..."
      enabled(analyse_btn) <- FALSE
      
      # Get kits.
      for(k in seq(along=val_kits)){
        tmp <- getKit(val_kits[k])
        val_kitData <- rbind(val_kitData, tmp)
      }
      
      # Get allele frequency database.
      val_db <- getDb(val_db_selected)
      
      if(debug){
        print("val_kits")
        print(val_kits)
        print("val_virtual")
        print(val_virtual)
        print("val_limit")
        print(val_limit)
      }
      
      # Analyse bins overlap.
      datanew <- calculateOL(kit=val_kitData,
                             db=val_db,
                             virtual=val_virtual,
                             limit=val_limit,
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
      
      message <- "At least one kit has to be selected."
      
      gmessage(message, title="Not kit selected",
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
      if(exists(".strvalidator_calculateOL_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_calculateOL_gui_savegui", envir=env)
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
      if(exists(".strvalidator_calculateOL_gui_db_name", envir=env, inherits = FALSE)){
        svalue(f1_db_drp) <- get(".strvalidator_calculateOL_gui_db_name", envir=env)
      }
      if(exists(".strvalidator_calculateOL_gui_virtual", envir=env, inherits = FALSE)){
        svalue(f1_virtual_chk) <- get(".strvalidator_calculateOL_gui_virtual", envir=env)
      }
      if(exists(".strvalidator_calculateOL_gui_limit", envir=env, inherits = FALSE)){
        svalue(f1_limit_chk) <- get(".strvalidator_calculateOL_gui_limit", envir=env)
      }
      
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_calculateOL_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_calculateOL_gui_db_name", value=svalue(f1_db_drp), envir=env)
      assign(x=".strvalidator_calculateOL_gui_virtual", value=svalue(f1_virtual_chk), envir=env)
      assign(x=".strvalidator_calculateOL_gui_limit", value=svalue(f1_limit_chk), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_calculateOL_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateOL_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_calculateOL_gui_db_name", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateOL_gui_db_name", envir = env)
      }
      if(exists(".strvalidator_calculateOL_gui_virtual", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateOL_gui_virtual", envir = env)
      }
      if(exists(".strvalidator_calculateOL_gui_limit", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateOL_gui_limit", envir = env)
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
