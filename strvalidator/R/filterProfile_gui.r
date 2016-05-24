################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 15.12.2015: Added attributes to result.
# 15.12.2015: Added option 'exact' sample name matching.
# 29.08.2015: Added importFrom.
# 05.05.2015: Changed parameter 'ignoreCase' to 'ignore.case' for 'checkSubset' function.
# 09.04.2015: Added option 'invert' to filter peaks NOT in reference.
# 11.10.2014: Added 'focus', added 'parent' parameter.
# 28.06.2014: Added help button and moved save gui checkbox.
# 08.05.2014: Implemented 'checkDataset'.
# 12.01.2014: Replaced 'subset' with native code.
# 15.12.2013: Fixed filter by kit bins.
# 09.12.2013: Added 'filter by' option.
# 09.12.2013: Added check subset button.
# 18.07.2013: Check before overwrite object.
# 15.07.2013: Added save GUI settings.
# 11.06.2013: Added 'inherits=FALSE' to 'exists'.
# 04.06.2013: Fixed bug in 'missingCol'.
# 24.05.2013: Improved error message for missing columns.
# 17.05.2013: listDataFrames() -> listObjects()
# 09.05.2013: .result removed, added save as group.
# 25.04.2013: First version.

#' @title Filter Profile
#'
#' @description
#' GUI wrapper for the \code{\link{filterProfile}} function.
#'
#' @details Simplifies the use of the \code{\link{filterProfile}} function
#' by providing a graphical user interface to it.
#' All data not matching/matching the reference will be discarded.
#' Useful for filtering stutters and artifacts from raw typing data or
#' to identify drop-ins.
#' 
#' @param env environment in wich to search for data frames.
#' @param savegui logical indicating if GUI settings should be saved in the environment.
#' @param debug logical indicating printing debug information.
#' @param parent widget to get focus when finished.
#' 
#' @export
#' 
#' @importFrom utils help str
#' @importFrom graphics title
#' 
#' @return TRUE
#' 
#' @seealso \code{\link{filterProfile}}, \code{\link{checkSubset}}

filterProfile_gui <- function(env=parent.frame(), savegui=NULL, debug=FALSE, parent=NULL){
  
  .gData <- NULL
  .gDataName <- NULL
  .gRef <- NULL
  .gRefName <- NULL
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  
  w <- gwindow(title="Filter profile", visible=FALSE)
  
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
    print(help("filterProfile_gui", help_type="html"))
    
  })
  
  # FRAME 0 ###################################################################
  
  f0 <- gframe(text = "Datasets",
               horizontal=FALSE,
               spacing = 5,
               container = gv) 
  
  g0 <- glayout(container = f0, spacing = 1)
  
  # Datasets ------------------------------------------------------------------
  
  g0[1,1] <- glabel(text="Select dataset:", container=g0)

  g0[1,2] <- g0_dataset_drp <- gdroplist(items=c("<Select dataset>",
                                   listObjects(env=env,
                                               obj.class="data.frame")), 
                           selected = 1,
                           editable = FALSE,
                           container = g0)
  
  g0[1,3] <- g0_samples_lbl <- glabel(text=" 0 samples", container=g0)
  
  addHandlerChanged(g0_dataset_drp, handler = function (h, ...) {
    
    val_obj <- svalue(g0_dataset_drp)
    
    # Check if suitable.
    requiredCol <- c("Sample.Name", "Marker", "Allele", "Height")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       slim=TRUE, slimcol="Allele",
                       env=env, parent=w, debug=debug)
    
    if(ok){
      
      # Load or change components.
      .gData <<- get(val_obj, envir=env)
      .gDataName <<- val_obj
      samples <- length(unique(.gData$Sample.Name))
      svalue(g0_samples_lbl) <- paste("", samples, "samples")
      svalue(f2_save_edt) <- paste(.gDataName, "_filter", sep="")
      
      # Detect kit.
      kitIndex <- detectKit(.gData, index=TRUE)
      # Select in dropdown.
      svalue(g0_kit_drp, index=TRUE) <- kitIndex
        
    } else {

      # Reset components.
      .gData <<- NULL
      svalue(g0_dataset_drp, index=TRUE) <- 1
      svalue(g0_samples_lbl) <- " 0 samples"
      svalue(f2_save_edt) <- ""
      
    }
    
  } )  
  
  g0[2,1] <- g0_refset_lbl <- glabel(text="Select reference:", container=g0)
  
  g0[2,2] <- g0_refset_drp <- gdroplist(items=c("<Select dataset>",
                                   listObjects(env=env,
                                               obj.class="data.frame")), 
                           selected = 1,
                           editable = FALSE,
                           container = g0) 
  
  g0[2,3] <- g0_ref_lbl <- glabel(text=" 0 references", container=g0)
  
  addHandlerChanged(g0_refset_drp, handler = function (h, ...) {
    
    val_obj <- svalue(g0_refset_drp)
    
    # Check if suitable.
    requiredCol <- c("Sample.Name", "Marker", "Allele")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       slim=TRUE, slimcol="Allele",
                       env=env, parent=w, debug=debug)
    
    if(ok){
      
      # Load or change components.
      .gRef <<- get(val_obj, envir=env)
      .gRefName <<- val_obj
      ref <- length(unique(.gRef$Sample.Name))
      svalue(g0_ref_lbl) <- paste("", ref, "references")
        
    } else {

      # Reset components.
      .gRef <<- NULL
      svalue(g0_refset_drp, index=TRUE) <- 1
      svalue(g0_ref_lbl) <- " 0 references"
      
    }
    
  } )
  
  # CHECK ---------------------------------------------------------------------
  
  g0[3,2] <- g0_check_btn <- gbutton(text="Check subsetting",
                                         border=TRUE,
                                         container=g0)
  
  addHandlerChanged(g0_check_btn, handler = function(h, ...) {
    
    # Get values.
    val_data <- .gData
    val_ref <- .gRef
    val_ignore <- svalue(f1_ignore_case_chk)
    val_exact <- svalue(f1_exact_chk)
    val_word <- FALSE
    
    if (!is.null(.gData) || !is.null(.gRef)){
      
      chksubset_w <- gwindow(title = "Check subsetting",
                             visible = FALSE, name=title,
                             width = NULL, height= NULL, parent=w,
                             handler = NULL, action = NULL)
      
      chksubset_txt <- checkSubset(data=val_data,
                                   ref=val_ref,
                                   console=FALSE,
                                   ignore.case=val_ignore,
                                   exact=val_exact,
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
  
  # Kit -------------------------------------------------------------------
  
  g0[4,1] <- g0_kit_lbl <- glabel(text="Select kit:", container=g0)
  
  g0[4,2] <- g0_kit_drp <- gdroplist(items=getKit(),
                                     selected = 1,
                                     editable = FALSE,
                                     container = g0) 
  
  g0[4,3] <- g0_kit_chk <- gcheckbox(text="Exclude virtual bins.",
                                         checked=TRUE,
                                         container=g0)

  # FRAME 1 ###################################################################
  
  f1 <- gframe(text = "Options",
               horizontal=FALSE,
               spacing = 5,
               container = gv) 

  f1_invert_chk <- gcheckbox(text="Invert (remove peaks matching reference)",
                                       checked = FALSE,
                                       container = f1)
  
  f1_add_missing_loci_chk <- gcheckbox(text="Add missing loci",
                                    checked = TRUE,
                                    container = f1)
  
  f1_keep_na_chk <- gcheckbox(text="Keep loci/sample even if no matching allele",
                           checked = TRUE,
                           container = f1)
  
  f1_ignore_case_chk <- gcheckbox(text="Ignore case ",
                           checked = TRUE,
                           container = f1)
  
  f1_exact_chk <- gcheckbox(text="Exact matching ",
                                  checked = FALSE,
                                  container = f1)
  
  f1_options <- c("Filter by reference dataset",
                  "Filter by kit bins (allelic ladder)")
  
  f1_filter_opt <- gradio(items=f1_options,
                          selected=1,
                          horizontal=FALSE,
                          container=f1)
  
  # HANDLERS ------------------------------------------------------------------
  
  addHandlerChanged(f1_filter_opt, handler = function (h, ...) {
    
    .refreshOptions()
    
  } )  
  

  # FRAME 2 ###################################################################
  
  f2 <- gframe(text = "Save as",
               horizontal=TRUE,
               spacing = 5,
               container = gv) 
  
  glabel(text="Name for result:", container=f2)
  
  f2_save_edt <- gedit(text="", container=f2)

  # BUTTON ####################################################################
  
  
  filter_btn <- gbutton(text="Filter profile",
                        border=TRUE,
                        container=gv)
  
  addHandlerChanged(filter_btn, handler = function(h, ...) {
    
    val_data <- .gData
    val_name_data <- .gDataName
    val_ref <- .gRef
    val_name_ref <- .gRefName
    val_invert <- svalue(f1_invert_chk)
    val_add_missing_loci <- svalue(f1_add_missing_loci_chk)
    val_keep_na <- svalue(f1_keep_na_chk)
    val_ignore_case <- svalue(f1_ignore_case_chk)
    val_exact <- svalue(f1_exact_chk)
    val_name <- svalue(f2_save_edt)
    val_filter <- svalue(f1_filter_opt, index=TRUE)
    val_kit <- svalue(g0_kit_drp)
    val_exclude <- svalue(g0_kit_chk)

    # Check if filter by kit bins.
    if(val_filter == 2){
      
      # Get markers, bins and flag for virtual bins.
      val_ref <- getKit(kit=val_kit, what="VIRTUAL")
      
      if(val_exclude){
        # Remove virtual bins.
        val_ref <- val_ref[val_ref$Virtual == 0, ]
      }
      
    }
      
    if(!is.null(val_data) & !is.null(val_ref)){
      
      # Change button.
      svalue(filter_btn) <- "Processing..."
      enabled(filter_btn) <- FALSE
  
      datanew <- filterProfile(data = val_data,
                               ref = val_ref,
                               add.missing.loci = val_add_missing_loci,
                               keep.na = val_keep_na,
                               ignore.case = val_ignore_case,
                               exact = val_exact,
                               invert = val_invert,
                               debug = debug)

      # Add attributes.
      attr(datanew, which="filterProfile_gui, data") <- val_name_data
      attr(datanew, which="filterProfile_gui, ref") <- val_name_ref
      attr(datanew, which="filterProfile_gui, add.missing.loci") <- val_add_missing_loci
      attr(datanew, which="filterProfile_gui, keep.na") <- val_keep_na
      attr(datanew, which="filterProfile_gui, ignore.case") <- val_ignore_case
      attr(datanew, which="filterProfile_gui, exact") <- val_exact
      attr(datanew, which="filterProfile_gui, invert") <- val_invert

      # Save data.
      saveObject(name=val_name, object=datanew, parent=w, env=env)
      
      if(debug){
        print(str(datanew))
        print(paste("EXIT:", match.call()[[1]]))
      }
      
      # Close GUI.
      dispose(w)
      
    } else {
      
      message <- "A dataset and a reference dataset have to be selected."
      
      gmessage(message, title="Datasets not selected",
               icon = "error",
               parent = w) 
      
    }
    
  } )

  # INTERNAL FUNCTIONS ########################################################
  
  .refreshOptions <- function(){

    val_opt <- svalue(f1_filter_opt, index=TRUE)
    
    if(val_opt == 1){
      
      enabled(g0_refset_lbl) <- TRUE
      enabled(g0_refset_drp) <- TRUE
      enabled(g0_check_btn) <- TRUE
      
      enabled(g0_kit_lbl) <- FALSE
      enabled(g0_kit_drp) <- FALSE
      enabled(g0_kit_chk) <- FALSE
      
    } else if(val_opt == 2){
      
      enabled(g0_refset_lbl) <- FALSE
      enabled(g0_refset_drp) <- FALSE
      enabled(g0_check_btn) <- FALSE
      
      enabled(g0_kit_lbl) <- TRUE
      enabled(g0_kit_drp) <- TRUE
      enabled(g0_kit_chk) <- TRUE
      
    }
    
  }
  
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
      if(exists(".strvalidator_filterProfile_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_filterProfile_gui_savegui", envir=env)
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
      if(exists(".strvalidator_filterProfile_gui_invert", envir=env, inherits = FALSE)){
        svalue(f1_invert_chk) <- get(".strvalidator_filterProfile_gui_invert", envir=env)
      }
      if(exists(".strvalidator_filterProfile_gui_add_loci", envir=env, inherits = FALSE)){
        svalue(f1_add_missing_loci_chk) <- get(".strvalidator_filterProfile_gui_add_loci", envir=env)
      }
      if(exists(".strvalidator_filterProfile_gui_keep_na", envir=env, inherits = FALSE)){
        svalue(f1_keep_na_chk) <- get(".strvalidator_filterProfile_gui_keep_na", envir=env)
      }
      if(exists(".strvalidator_filterProfile_gui_ignore_case", envir=env, inherits = FALSE)){
        svalue(f1_ignore_case_chk) <- get(".strvalidator_filterProfile_gui_ignore_case", envir=env)
      }
      if(exists(".strvalidator_filterProfile_gui_exact", envir=env, inherits = FALSE)){
        svalue(f1_exact_chk) <- get(".strvalidator_filterProfile_gui_exact", envir=env)
      }
      if(exists(".strvalidator_filterProfile_gui_by", envir=env, inherits = FALSE)){
        svalue(f1_filter_opt) <- get(".strvalidator_filterProfile_gui_by", envir=env)
      }
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_filterProfile_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_filterProfile_gui_invert", value=svalue(f1_invert_chk), envir=env)
      assign(x=".strvalidator_filterProfile_gui_add_loci", value=svalue(f1_add_missing_loci_chk), envir=env)
      assign(x=".strvalidator_filterProfile_gui_keep_na", value=svalue(f1_keep_na_chk), envir=env)
      assign(x=".strvalidator_filterProfile_gui_ignore_case", value=svalue(f1_ignore_case_chk), envir=env)
      assign(x=".strvalidator_filterProfile_gui_exact", value=svalue(f1_exact_chk), envir=env)
      assign(x=".strvalidator_filterProfile_gui_by", value=svalue(f1_filter_opt), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_filterProfile_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_filterProfile_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_filterProfile_gui_invert", envir=env, inherits = FALSE)){
        remove(".strvalidator_filterProfile_gui_invert", envir = env)
      }
      if(exists(".strvalidator_filterProfile_gui_add_loci", envir=env, inherits = FALSE)){
        remove(".strvalidator_filterProfile_gui_add_loci", envir = env)
      }
      if(exists(".strvalidator_filterProfile_gui_keep_na", envir=env, inherits = FALSE)){
        remove(".strvalidator_filterProfile_gui_keep_na", envir = env)
      }
      if(exists(".strvalidator_filterProfile_gui_ignore_case", envir=env, inherits = FALSE)){
        remove(".strvalidator_filterProfile_gui_ignore_case", envir = env)
      }
      if(exists(".strvalidator_filterProfile_gui_exact", envir=env, inherits = FALSE)){
        remove(".strvalidator_filterProfile_gui_exact", envir = env)
      }
      if(exists(".strvalidator_filterProfile_gui_by", envir=env, inherits = FALSE)){
        remove(".strvalidator_filterProfile_gui_by", envir = env)
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
  .refreshOptions()
  
  # Show GUI.
  visible(w) <- TRUE
  focus(w)
  
}
