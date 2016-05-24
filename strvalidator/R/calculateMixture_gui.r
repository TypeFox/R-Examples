################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 28.08.2015: Added importFrom.
# 05.05.2015: Changed parameter 'ignoreCase' to 'ignore.case' for 'checkSubset' function.
# 11.10.2014: Added 'focus', added 'parent' parameter.
# 28.08.2014: Fixed bug in 'Check subsetting' showing extra combinations in many cases.
# 08.07.2014: First version.

#' @title Calculate Mixture
#'
#' @description
#' GUI wrapper for the \code{\link{calculateMixture}} function.
#'
#' @details
#' Simplifies the use of the \code{\link{calculateMixture}} function by
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
#' @seealso \code{\link{calculateMixture}}, \code{\link{checkSubset}}

calculateMixture_gui <- function(env=parent.frame(), savegui=NULL,
                                 debug=FALSE, parent=NULL){
  
  # Global variables.
  .gData <- NULL
  .gRef1 <- NULL
  .gRef2 <- NULL
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # WINDOW ####################################################################
  
  if(debug){
    print("WINDOW")
  }  

  # Main window.
  w <- gwindow(title="Calculate Mixture", visible=FALSE)

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
    print(help("calculateMixture_gui", help_type="html"))
    
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
    requiredCol <- c("Sample.Name", "Marker", "Allele")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       slim=TRUE, slimcol="Allele",
                       env=env, parent=w, debug=debug)
    
    if(ok){
      # Load or change components.
      
      # get dataset.
      .gData <<- get(val_obj, envir=env)
      svalue(g0_data_samples_lbl) <- paste(length(unique(.gData$Sample.Name)),
                                           "samples.")
      svalue(f4_save_edt) <- paste(val_obj, "_mixture", sep="")
      
    } else {
      
      # Reset components.
      .gData <<- NULL
      svalue(g0_data_drp, index=TRUE) <- 1
      svalue(g0_data_samples_lbl) <- " 0 samples"
      svalue(f4_save_edt) <- ""
      
    }
    
  } )  
  
  # Reference 1 ---------------------------------------------------------------
  
  g0[2,1] <- glabel(text="Select reference dataset (major):", container=g0)
  
  # NB! dfs defined in previous section.
  g0[2,2] <- g0_ref1_drp <- gdroplist(items=dfs, 
                                     selected = 1,
                                     editable = FALSE,
                                     container = g0)
  
  g0[2,3] <- g0_ref1_samples_lbl <- glabel(text=" 0 references", container=g0)
  
  addHandlerChanged(g0_ref1_drp, handler = function (h, ...) {
    
    val_obj <- svalue(g0_ref1_drp)
    
    # Check if suitable.
    requiredCol <- c("Sample.Name", "Marker", "Allele")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       slim=TRUE, slimcol="Allele",
                       env=env, parent=w, debug=debug)
    
    if(ok){
      # Load or change components.
      
      .gRef1 <<- get(val_obj, envir=env)
      svalue(g0_ref1_samples_lbl) <- paste(length(unique(.gRef1$Sample.Name)),
                                          "samples.")
      
    } else {
      
      # Reset components.
      .gRef1 <<- NULL
      svalue(g0_ref1_drp, index=TRUE) <- 1
      svalue(g0_ref1_samples_lbl) <- " 0 references"
      
    }    
  } )  

  # Reference 2 ---------------------------------------------------------------
  
  g0[3,1] <- glabel(text="Select reference dataset (minor):", container=g0)
  
  # NB! dfs defined in previous section.
  g0[3,2] <- g0_ref2_drp <- gdroplist(items=dfs, 
                                     selected = 1,
                                     editable = FALSE,
                                     container = g0)
  
  g0[3,3] <- g0_ref2_samples_lbl <- glabel(text=" 0 references", container=g0)
  
  addHandlerChanged(g0_ref2_drp, handler = function (h, ...) {
    
    val_obj <- svalue(g0_ref2_drp)
    
    # Check if suitable.
    requiredCol <- c("Sample.Name", "Marker", "Allele")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       slim=TRUE, slimcol="Allele",
                       env=env, parent=w, debug=debug)
    
    if(ok){
      # Load or change components.
      
      .gRef2 <<- get(val_obj, envir=env)
      svalue(g0_ref2_samples_lbl) <- paste(length(unique(.gRef2$Sample.Name)),
                                          "samples.")
      
    } else {
      
      # Reset components.
      .gRef2 <<- NULL
      svalue(g0_ref2_drp, index=TRUE) <- 1
      svalue(g0_ref2_samples_lbl) <- " 0 references"
      
    }    
  } )  
  
  
  # CHECK ---------------------------------------------------------------------
  
  if(debug){
    print("CHECK")
  }  
  
  g0[4,2] <- g0_check_btn <- gbutton(text="Check subsetting",
                                     border=TRUE,
                                     container=g0)
  
  addHandlerChanged(g0_check_btn, handler = function(h, ...) {
    
    # Get values.
    val_data <- .gData
    val_ref1 <- .gRef1
    val_ref2 <- .gRef2
    val_ignore <- FALSE
    val_word <- FALSE
    
    if (!is.null(.gData) || !is.null(.gRef1) || !is.null(.gRef2)){
      
      chksubset_w <- gwindow(title = "Check subsetting",
                             visible = FALSE, name=title,
                             width = NULL, height= NULL, parent=w,
                             handler = NULL, action = NULL)

      # Create pattern.
      tmp <- paste(".*", unique(.gRef1$Sample.Name), ".*",
                   unique(.gRef2$Sample.Name), ".*", sep="")
      
      # Save as dataframe.
      val_pattern <- data.frame(Sample.Name=tmp, stringsAsFactors=FALSE) 
      
      if(debug){
        print("Pattern")
        print(val_pattern)
      }
      
      chksubset_txt <- checkSubset(data=val_data,
                                   ref=val_pattern,
                                   console=FALSE,
                                   ignore.case=val_ignore,
                                   word=val_word,
                                   debug=debug)
      
      gtext (text = chksubset_txt, width = NULL, height = 300, font.attr = NULL, 
             wrap = FALSE, container = chksubset_w)
      
      visible(chksubset_w) <- TRUE
      
    } else {
      
      gmessage(message="Data frame is NULL!\n\n
               Make sure to select a dataset and two reference sets",
               title="Error",
               icon = "error")      
      
    } 
    
  } )
  
  
  # FRAME 1 ###################################################################
  
  if(debug){
    print("FRAME 1")
  }  
  
  f1 <- gframe(text = "Options", horizontal=FALSE, spacing = 10, container = gv) 

  f1_ol_chk <- gcheckbox(text="Remove off-ladder alleles (affects number of drop-in)",
                         checked = TRUE, container = f1)
  
  f1_drop_chk <- gcheckbox(text="Ignore drop-out (calculate Mx anyway)",
                         checked = TRUE, container = f1)
  
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
  
  # BUTTON ####################################################################

  if(debug){
    print("BUTTON")
  }  
  
  calculate_btn <- gbutton(text="Calculate", border=TRUE, container=gv)
  
  addHandlerChanged(calculate_btn, handler = function(h, ...) {
    
    # Get values.
    val_data <- .gData
    val_ref1 <- .gRef1
    val_ref2 <- .gRef2
    val_name <- svalue(f4_save_edt)
    val_ol <- svalue(f1_ol_chk)
    val_drop <- svalue(f1_drop_chk)
    
    if(debug){
      print("Read Values:")
      print("val_data")
      print(head(val_data))
      print("val_ref1")
      print(head(val_ref1))
      print("val_ref2")
      print(head(val_ref2))
      print("val_ol")
      print(val_ol)
      print("val_drop")
      print(val_drop)
      print("val_name")
      print(val_name)
    }
    
    # Check if data.
    if(!is.null(.gData) & !is.null(.gRef1) & !is.null(.gRef2)){
      
      # Change button.
      svalue(calculate_btn) <- "Processing..."
      enabled(calculate_btn) <- FALSE
      
      datanew <- calculateMixture(data=val_data,
                                  ref1=val_ref1,
                                  ref2=val_ref2,
                                  ol.rm=val_ol,
                                  ignore.dropout=val_drop,
                                  debug=debug)
      
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
      if(exists(".strvalidator_calculateMixture_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_calculateMixture_gui_savegui", envir=env)
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
      if(exists(".strvalidator_calculateMixture_gui_ol", envir=env, inherits = FALSE)){
        svalue(f1_ol_chk) <- get(".strvalidator_calculateMixture_gui_ol", envir=env)
      }
      if(exists(".strvalidator_calculateMixture_gui_dropout", envir=env, inherits = FALSE)){
        svalue(f1_drop_chk) <- get(".strvalidator_calculateMixture_gui_dropout", envir=env)
      }
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_calculateMixture_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_calculateMixture_gui_ol", value=svalue(f1_ol_chk), envir=env)
      assign(x=".strvalidator_calculateMixture_gui_dropout", value=svalue(f1_drop_chk), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_calculateMixture_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateMixture_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_calculateMixture_gui_ol", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateMixture_gui_ol", envir = env)
      }
      if(exists(".strvalidator_calculateMixture_gui_dropout", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateMixture_gui_dropout", envir = env)
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
