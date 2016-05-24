################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 22.12.2015: First version.

#' @title Calculate Ratio
#'
#' @description
#' GUI wrapper for the \code{\link{calculateRatio}} function.
#'
#' @details
#' Simplifies the use of the \code{\link{calculateRatio}} function
#' by providing a graphical user interface.
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
#' @seealso \code{link{calculateRatio}}, \code{link{checkSubset}}

calculateRatio_gui <- function(env=parent.frame(), savegui=NULL,
                                 debug=FALSE, parent=NULL){
  
  # Global variables.
  .gData <- NULL
  .gRef <- NULL
  .gDataName <- NULL
  .gRefName <- NULL
  .datasetDropDefault <- "<Select dataset>"
  .markerDropDefault <- "<Select marker>"
  .groupDropDefault <- "<Select column>"
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # WINDOW ####################################################################
  
  if(debug){
    print("WINDOW")
  }  

  # Main window.
  w <- gwindow(title="Calculate ratio", visible=FALSE)

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
    print(help("calculateRatio_gui", help_type="html"))
    
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
  
  dfs <- c(.datasetDropDefault, listObjects(env=env, obj.class="data.frame"))
  
  g0[1,2] <- g0_data_drp <- gdroplist(items=dfs, 
                           selected = 1,
                           editable = FALSE,
                           container = g0)
  g0[1,3] <- g0_data_samples_lbl <- glabel(text=" 0 samples", container=g0)
  
  addHandlerChanged(g0_data_drp, handler = function (h, ...) {
    
    val_obj <- svalue(g0_data_drp)
    
    # Check if suitable.
    requiredCol <- c("Sample.Name", "Marker", "Allele", "Height")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       slim=TRUE, slimcol="Height",
                       env=env, parent=w, debug=debug)
    
    if(ok){
      # Load or change components.
      
      # get dataset.
      .gData <<- get(val_obj, envir=env)
      .gDataName <<- val_obj
      svalue(g0_data_samples_lbl) <- paste(length(unique(.gData$Sample.Name)),
                                        "samples.")
      # Update dropdown menues.
      f1_numerator_drp[,] <- unique(c(.markerDropDefault, .gData$Marker))
      f1_denominator_drp[,] <- unique(c(.markerDropDefault, .gData$Marker))
      f1_group_drp[,] <- unique(c(.groupDropDefault, names(.gData)))

      # Suggest a name for the result.
      svalue(f4_save_edt) <- paste(val_obj, "_ratio", sep="")
      
    } else {
      
      # Reset components.
      .gData <<- NULL
      .gDataName <<- NULL
      svalue(g0_data_drp, index=TRUE) <- 1
      svalue(g0_data_samples_lbl) <- " 0 samples"
      svalue(f4_save_edt) <- ""

      # Update dropdown menues.
      f1_numerator_drp[,] <- .markerDropDefault
      f1_denominator_drp[,] <- .markerDropDefault
      f1_group_drp[,] <- .groupDropDefault
      
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
      .gRefName <<- val_obj
      svalue(g0_ref_samples_lbl) <- paste(length(unique(.gRef$Sample.Name)),
                                          "samples.")
        
    } else {
      
      # Reset components.
      .gRef <<- NULL
      .gRefName <<- NULL
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

  f1_ol_chk <- gcheckbox(text="Remove off-ladder alleles", checked=TRUE,
                         container=f1)
  
  f1_ignore_chk <- gcheckbox(text="Ignore case", checked=TRUE,
                         container=f1)
  
  f1_word_chk <- gcheckbox(text="Add word boundaries", checked = FALSE,
                           container = f1)
  
  f1_exact_chk <- gcheckbox(text="Exact matching", checked = FALSE,
                           container = f1)

  f1g1 <- glayout(container = f1)
  
  f1g1[1,1] <- glabel(text = "Select numerator markers:", container = f1g1)
  f1g1[1,2] <- f1_numerator_drp <- gdroplist(items = .markerDropDefault, container=f1g1)
  f1g1[2,1:2] <- f1_numerator_edt <- gedit(text="", container=f1g1)
  
  f1g1[3,1] <- glabel(text = "Select denominator markers:", container = f1g1)
  f1g1[3,2] <- f1_denominator_drp <- gdroplist(items = .markerDropDefault, container=f1g1)
  f1g1[4,1:2] <- f1_denominator_edt <- gedit(text="", container=f1g1)
  
  f1g1[5,1] <- glabel(text = "Group by column:", container = f1g1)
  f1g1[5,2] <- f1_group_drp <- gdroplist(items = .groupDropDefault, container=f1g1)
  
  addHandlerChanged(f1_numerator_drp, handler = function (h, ...) {
    
    val_marker <- svalue(f1_numerator_drp)
    val_value <- svalue(f1_numerator_edt)

    if(!is.null(val_marker)){
      if(val_marker != .markerDropDefault){

        # Add new value to selected.
        if(nchar(val_value) == 0){
          
          svalue(f1_numerator_edt) <- val_marker
          
        } else {
          
          svalue(f1_numerator_edt) <- paste(val_value, val_marker, sep = ",")
          
        }
        
      }
    }

  } )  
  
  addHandlerChanged(f1_denominator_drp, handler = function (h, ...) {
    
    val_marker <- svalue(f1_denominator_drp)
    val_value <- svalue(f1_denominator_edt)

    if(!is.null(val_marker)){
      if(val_marker != .markerDropDefault){

        # Add new value to selected.
        if(nchar(val_value) == 0){
          
          svalue(f1_denominator_edt) <- val_marker
          
        } else {

          svalue(f1_denominator_edt) <- paste(val_value, val_marker, sep = ",")
          
        }
        
      }
    }    
    
  } )  
  

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
  
  calculate_btn <- gbutton(text="Calculate",
                      border=TRUE,
                      container=gv)
  
  addHandlerChanged(calculate_btn, handler = function(h, ...) {
    
    # Get values.
    val_data <- .gData
    val_name_data <- .gDataName
    val_ref <- .gRef
    val_name_ref <- .gRefName
    val_ol <- svalue(f1_ol_chk)
    val_ignore <- svalue(f1_ignore_chk)
    val_word <- svalue(f1_word_chk)
    val_exact <- svalue(f1_exact_chk)
    val_name <- svalue(f4_save_edt)
    val_numerator <- svalue(f1_numerator_edt)
    val_denominator <- svalue(f1_denominator_edt)
    val_group <- svalue(f1_group_drp)
    
    if(val_group == .groupDropDefault){
      val_group <- NULL
    }
    
    if(debug){
      print("Read Values:")
      print("val_data")
      print(head(val_data))
      print("val_ref")
      print(head(val_ref))
      print("val_ol")
      print(val_ol)
      print("val_ignore")
      print(val_ignore)
      print("val_word")
      print(val_word)
      print("val_exact")
      print(val_exact)
      print("val_name")
      print(val_name)
    }
    
    # Check if data.
    if(!is.null(.gData)){

      if(!nchar(val_numerator) > 0){
        
        val_numerator <- NULL

      } else {
      
        val_numerator <- unlist(strsplit(val_numerator, split = ","))
      }

      if(!nchar(val_denominator) > 0){
        
        val_denominator <- NULL
        
      } else {
        
        val_denominator <- unlist(strsplit(val_denominator, split = ","))
        
      }

      if(debug){
        print("Sent Values:")
        print("val_numerator")
        print(val_numerator)
        print("val_denominator")
        print(val_denominator)
        print("val_group")
        print(val_group)
      }
      
      
  
      # Change button.
      svalue(calculate_btn) <- "Processing..."
      enabled(calculate_btn) <- FALSE
      
      datanew <- calculateRatio(data = val_data, ref = val_ref,
                                numerator = val_numerator,
                                denominator = val_denominator,
                                group = val_group, ol.rm = val_ol,
                                ignore.case = val_ignore, word = val_word,
                                exact = val_exact, debug = debug)
      
      # Add attributes.
      attr(datanew, which="calculateRatio_gui, data") <- val_name_data
      attr(datanew, which="calculateRatio_gui, ref") <- val_name_ref
      attr(datanew, which="calculateRatio_gui, numerator") <- val_numerator
      attr(datanew, which="calculateRatio_gui, denominator") <- val_denominator
      attr(datanew, which="calculateRatio_gui, group") <- val_group
      attr(datanew, which="calculateRatio_gui, ol") <- val_ol
      attr(datanew, which="calculateRatio_gui, ignore.case") <- val_ignore
      attr(datanew, which="calculateRatio_gui, word") <- val_word
      attr(datanew, which="calculateRatio_gui, exact") <- val_exact

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
      if(exists(".strvalidator_calculateRatio_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_calculateRatio_gui_savegui", envir=env)
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
      if(exists(".strvalidator_calculateRatio_gui_numerator", envir=env, inherits = FALSE)){
        svalue(f1_numerator_edt) <- get(".strvalidator_calculateRatio_gui_numerator", envir=env)
      }
      if(exists(".strvalidator_calculateRatio_gui_denominator", envir=env, inherits = FALSE)){
        svalue(f1_denominator_edt) <- get(".strvalidator_calculateRatio_gui_denominator", envir=env)
      }
      if(exists(".strvalidator_calculateRatio_gui_ol", envir=env, inherits = FALSE)){
        svalue(f1_ol_chk) <- get(".strvalidator_calculateRatio_gui_ol", envir=env)
      }
      if(exists(".strvalidator_calculateRatio_gui_ignore", envir=env, inherits = FALSE)){
        svalue(f1_ignore_chk) <- get(".strvalidator_calculateRatio_gui_ignore", envir=env)
      }
      if(exists(".strvalidator_calculateRatio_gui_word", envir=env, inherits = FALSE)){
        svalue(f1_word_chk) <- get(".strvalidator_calculateRatio_gui_word", envir=env)
      }
      if(exists(".strvalidator_calculateRatio_gui_exact", envir=env, inherits = FALSE)){
        svalue(f1_exact_chk) <- get(".strvalidator_calculateRatio_gui_exact", envir=env)
      }
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_calculateRatio_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_calculateRatio_gui_numerator", value=svalue(f1_numerator_edt), envir=env)
      assign(x=".strvalidator_calculateRatio_gui_denominator", value=svalue(f1_denominator_edt), envir=env)
      assign(x=".strvalidator_calculateRatio_gui_word", value=svalue(f1_word_chk), envir=env)
      assign(x=".strvalidator_calculateRatio_gui_ignore", value=svalue(f1_ignore_chk), envir=env)
      assign(x=".strvalidator_calculateRatio_gui_exact", value=svalue(f1_exact_chk), envir=env)

    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_calculateRatio_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateRatio_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_calculateRatio_gui_numerator", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateRatio_gui_numerator", envir = env)
      }
      if(exists(".strvalidator_calculateRatio_gui_denominator", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateRatio_gui_denominator", envir = env)
      }
      if(exists(".strvalidator_calculateRatio_gui_ignore", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateRatio_gui_ignore", envir = env)
      }
      if(exists(".strvalidator_calculateRatio_gui_word", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateRatio_gui_word", envir = env)
      }
      if(exists(".strvalidator_calculateRatio_gui_exact", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateRatio_gui_exact", envir = env)
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
