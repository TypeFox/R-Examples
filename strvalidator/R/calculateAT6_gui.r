################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 30.11.2015: Added warning.
# 28.08.2015: Added importFrom
# 10.06.2015: Added missing label 'Significance level:'.
# 26.05.2015: First version.

#' @title Calculate Analytical Threshold
#'
#' @description
#' GUI wrapper for the \code{\link{calculateAT6}} function.
#'
#' @details Scores dropouts for a dataset.
#' @param env environment in wich to search for data frames and save result.
#' @param savegui logical indicating if GUI settings should be saved in the environment.
#' @param debug logical indicating printing debug information.
#' @param parent widget to get focus when finished.
#' 
#' @return TRUE
#' 
#' @export
#' 
#' @importFrom utils help head
#' @importFrom graphics title
#' 
#' @seealso \code{\link{calculateAT6}}, \code{\link{calculateAT}},
#'  \code{\link{calculateAT_gui}}, \code{\link{checkSubset}}

calculateAT6_gui <- function(env=parent.frame(), savegui=NULL,
                             debug=FALSE, parent=NULL){
  
  # Global variables.
  .gData <- NULL
  .gRef <- NULL
  .gAm <- NULL
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # Main window.
  w <- gwindow(title="Calculate analytical threshold", visible=FALSE)
  
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
    print(help("calculateAT6_gui", help_type="html"))
    
  })
  
  # FRAME 0 ###################################################################
  
  f0 <- gframe(text = "Datasets",
               horizontal=FALSE,
               spacing = 5,
               container = gv) 
  
  g0 <- glayout(container = f0, spacing = 1)
  
  # Datasets ------------------------------------------------------------------
  
  g0[1,1] <- glabel(text="Select dataset:", container=g0)
  
  g0[1,2] <- dataset_drp <- gdroplist(items=c("<Select dataset>",
                                              listObjects(env=env,
                                                          obj.class="data.frame")), 
                                      selected = 1,
                                      editable = FALSE,
                                      container = g0)
  
  g0[1,3] <- g0_samples_lbl <- glabel(text=" 0 samples", container=g0)
  
  addHandlerChanged(dataset_drp, handler = function (h, ...) {
    
    val_obj <- svalue(dataset_drp)
    
    # Check if suitable.
    requiredCol <- c("Sample.Name", "Marker", "Allele", "Height")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       env=env, parent=w, debug=debug)
    
    if(ok){
      # Load or change components.
      
      .gData <<- get(val_obj, envir=env)
      samples <- length(unique(.gData$Sample.Name))
      svalue(g0_samples_lbl) <- paste("", samples, "samples")
      
      # Suggest a name for result.
      svalue(f2_save_edt) <- paste(val_obj, "_at6", sep="")
      
    } else {
      
      # Reset components.
      .gData <<- NULL
      svalue(dataset_drp, index=TRUE) <- 1
      svalue(g0_samples_lbl) <- " 0 samples"
      svalue(f2_save_edt) <- ""
      
    }
    
  } )  
  
  g0[2,1] <- glabel(text="Select reference dataset:", container=g0)
  
  g0[2,2] <- refset_drp <- gdroplist(items=c("<Select dataset>",
                                             listObjects(env=env,
                                                         obj.class="data.frame")), 
                                     selected = 1,
                                     editable = FALSE,
                                     container = g0) 
  
  g0[2,3] <- g0_ref_lbl <- glabel(text=" 0 references", container=g0)
  
  addHandlerChanged(refset_drp, handler = function (h, ...) {
    
    val_obj <- svalue(refset_drp)
    
    # Check if suitable.
    requiredCol <- c("Sample.Name", "Marker", "Allele")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       env=env, parent=w, debug=debug)
    
    if(ok){
      # Load or change components.
      
      .gRef <<- get(val_obj, envir=env)
      ref <- length(unique(.gRef$Sample.Name))
      svalue(g0_ref_lbl) <- paste("", ref, "references")
      
    } else {
      
      # Reset components.
      .gRef <<- NULL
      svalue(refset_drp, index=TRUE) <- 1
      svalue(g0_ref_lbl) <- " 0 references"
      
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
    val_ignore <- svalue(f1_ignore_case_chk)
    
    if (!is.null(.gData) || !is.null(.gRef)){
      
      chksubset_w <- gwindow(title = "Check subsetting",
                             visible = FALSE, name=title,
                             width = NULL, height= NULL, parent=w,
                             handler = NULL, action = NULL)
      
      chksubset_txt <- checkSubset(data=val_data,
                                   ref=val_ref,
                                   console=FALSE,
                                   ignore.case=val_ignore,
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
  
  # AMOUNT --------------------------------------------------------------------
  
  g0[4,1] <- glabel(text="Select amount dataset:", container=g0)
  
  g0[4,2] <- amset_drp <- gdroplist(items=c("<Select dataset>",
                                            listObjects(env=env,
                                                        obj.class="data.frame")), 
                                    selected = 1,
                                    editable = FALSE,
                                    container = g0) 
  
  g0[4,3] <- g0_am_lbl <- glabel(text=" 0 samples", container=g0)
  
  addHandlerChanged(amset_drp, handler = function (h, ...) {
    
    val_obj <- svalue(amset_drp)
    
    # Check if suitable.
    requiredCol <- c("Sample.Name", "Amount")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       env=env, parent=w, debug=debug)
    
    if(ok){
      # Load or change components.
      
      .gAm <<- get(val_obj, envir=env)
      am <- length(unique(.gAm$Sample.Name))
      svalue(g0_am_lbl) <- paste("", am, "samples")
      
    } else {
      
      # Reset components.
      .gAm <<- NULL
      svalue(amset_drp, index=TRUE) <- 1
      svalue(g0_am_lbl) <- " 0 samples"
      
    }
    
  } )  
  
  # FRAME 1 ###################################################################
  
  f1 <- gframe(text = "Options",
               horizontal=FALSE,
               spacing = 5,
               container = gv) 
  
  glabel(text="NB! This is an indirect method not recommended.",
         anchor=c(-1 ,0), container=f1)
  glabel(text="See 'Help' or reference for limitations.",
         anchor=c(-1 ,0), container=f1)
  
  f1_ignore_case_chk <- gcheckbox(text="Ignore case", checked=TRUE,
                                  container=f1)
  
  f1_items <- c("Linear regression", "Weighted linear regression")
  f1_weighted_opt <- gradio(items=f1_items, selected=2, container=f1)
  
  glabel(text="Significance level:", anchor=c(-1 ,0), container=f1)
  f1_alpha_spn <- gspinbutton(from=0, to=1, by=0.01, value=0.05, container=f1)
  
  # FRAME 2 ###################################################################
  
  f2 <- gframe(text = "Save as",
               horizontal=TRUE,
               spacing = 5,
               container = gv) 
  
  glabel(text="Name for result:", container=f2)
  
  f2_save_edt <- gedit(text="", container=f2)
  
  # BUTTON ####################################################################
  
  calculate_btn <- gbutton(text="Calculate",
                           border=TRUE,
                           container=gv)
  
  addHandlerChanged(calculate_btn, handler = function(h, ...) {
    
    val_ignore_case <- svalue(f1_ignore_case_chk)
    val_weighted <- ifelse(svalue(f1_weighted_opt, index=TRUE)==1, FALSE, TRUE)
    val_alpha <- svalue(f1_alpha_spn)
    val_name <- svalue(f2_save_edt)
    
    if(debug){
      print("GUI options:")
      print("val_ignore_case")
      print(val_ignore_case)
      print("val_weighted")
      print(val_weighted)
      print("val_alpha")
      print(val_alpha)
      print("val_name")
      print(val_name)
    }
    
    if(!is.null(.gData) & !is.null(.gRef)){
      
      # Change button.
      svalue(calculate_btn) <- "Processing..."
      enabled(calculate_btn) <- FALSE
      
      datanew <- calculateAT6(data=.gData,
                              ref=.gRef,
                              amount=.gAm,
                              weighted=val_weighted,
                              alpha=val_alpha,
                              ignore.case=val_ignore_case,
                              debug=debug)
      
      # Save data.
      saveObject(name=val_name, object=datanew, parent=w, env=env)
      
      if(debug){
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
      if(exists(".strvalidator_calculateAT6_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_calculateAT6_gui_savegui", envir=env)
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
      if(exists(".strvalidator_calculateAT6_gui_ignore", envir=env, inherits = FALSE)){
        svalue(f1_ignore_case_chk) <- get(".strvalidator_calculateAT6_gui_ignore", envir=env)
      }
      if(exists(".strvalidator_calculateAT6_gui_weighted", envir=env, inherits = FALSE)){
        svalue(f1_weighted_opt) <- get(".strvalidator_calculateAT6_gui_weighted", envir=env)
      }
      if(exists(".strvalidator_calculateAT6_gui_alpha", envir=env, inherits = FALSE)){
        svalue(f1_alpha_spn) <- get(".strvalidator_calculateAT6_gui_alpha", envir=env)
      }
      
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_calculateAT6_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_calculateAT6_gui_ignore", value=svalue(f1_ignore_case_chk), envir=env)
      assign(x=".strvalidator_calculateAT6_gui_weighted", value=svalue(f1_weighted_opt), envir=env)
      assign(x=".strvalidator_calculateAT6_gui_alpha", value=svalue(f1_alpha_spn), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_calculateAT6_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateAT6_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_calculateAT6_gui_ignore", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateAT6_gui_ignore", envir = env)
      }
      if(exists(".strvalidator_calculateAT6_gui_weighted", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateAT6_gui_weighted", envir = env)
      }
      if(exists(".strvalidator_calculateAT6_gui_alpha", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateAT6_gui_alpha", envir = env)
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