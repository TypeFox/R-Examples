################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 02.12.2016: Fixed options save bug.
# 30.12.2015: First version.

#' @title Calculate Locus Balance
#'
#' @description
#' GUI wrapper for the \code{\link{calculateLb}} function.
#'
#' @details
#' Simplifies the use of the \code{\link{calculateLb}} function
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
#' @seealso \code{link{calculateLb}}, \code{link{checkSubset}}
#' 

calculateLb_gui <- function(env=parent.frame(), savegui=NULL,
                                 debug=FALSE, parent=NULL){
  
  # Global variables.
  .gData <- NULL
  .gRef <- NULL
  .gDataName <- NULL
  .gRefName <- NULL
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # WINDOW ####################################################################
  
  if(debug){
    print("WINDOW")
  }  

  # Main window.
  w <- gwindow(title="Calculate locus balance", visible=FALSE)

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
    print(help("calculateLb_gui", help_type="html"))
    
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
    requiredCol <- c("Sample.Name", "Marker", "Height")
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
      
      # Suggest a name for the result.
      svalue(f4_save_edt) <- paste(val_obj, "_lb", sep="")
      
      # Detect kit.
      kitIndex <- detectKit(data = .gData, index = TRUE)
      # Select in dropdown.
      svalue(kit_drp, index = TRUE) <- kitIndex
      
    } else {
      
      # Reset components.
      .gData <<- NULL
      .gDataName <<- NULL
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
  tooltip(g0_ref_drp) <- "Known alleles will be extracted from data"
  
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
      # Enable checkbox to calculate H.
      enabled(f1_h_chk) <- TRUE
        
    } else {
      
      # Reset components.
      .gRef <<- NULL
      .gRefName <<- NULL
      svalue(g0_ref_drp, index=TRUE) <- 1
      svalue(g0_ref_samples_lbl) <- " 0 references"
      
      # Disable checkbox to calculate H.
      enabled(f1_h_chk) <- FALSE

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
    val_exact <- svalue(f1_exact_chk)
    
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
  
  # Kit -----------------------------------------------------------------------
  
  g0[4,1] <- glabel(text="Select the kit used:", container=g0)
  
  # NB! dfs defined in previous section.
  g0[4,2] <- kit_drp <- gdroplist(items = getKit(), 
                                  selected = 1,
                                  editable = FALSE,
                                  container = g0) 
  

  # FRAME 1 ###################################################################
  
  if(debug){
    print("FRAME 1")
  }  
  
  f1 <- gframe(text = "Options",
               horizontal=FALSE,
               spacing = 10,
               container = gv)
  
  #----------------------------------------------------------------------------
  glabel(text = "Calculate locus balance:", anchor=c(-1 ,0), container = f1)
  
  f1_options_lb <- c("Proportional",
                     "Normalised",
                     "Centred Quantities")
  
  f1_lb_opt <- gradio(items = f1_options_lb,
                      selected = 1,
                      horizontal = TRUE,
                      container = f1)
  
  f1_dye_chk <- gcheckbox(text="Calculate Lb by dye channel", checked = FALSE,
                          container = f1)

  #----------------------------------------------------------------------------
  glabel(text = "Reference sample name matching:", anchor=c(-1 ,0),
         container = f1)

  f1_ignore_chk <- gcheckbox(text="Ignore case", checked=TRUE,
                         container=f1)
  
  f1_word_chk <- gcheckbox(text="Add word boundaries", checked = FALSE,
                           container = f1)
  
  f1_exact_chk <- gcheckbox(text="Exact matching", checked = FALSE,
                           container = f1)
  
  #----------------------------------------------------------------------------
  glabel(text = "Filter data:", anchor=c(-1 ,0), container = f1)
  
  f1_ol_chk <- gcheckbox(text="Remove off-ladder alleles", checked = TRUE,
                         container = f1)
  
  f1_sex_chk <- gcheckbox(text="Remove sex markers", checked = FALSE,
                           container = f1)

  #----------------------------------------------------------------------------
  glabel(text = "Replace missing data with peak height:", anchor=c(-1 ,0),
         container = f1)
  f1_na_edt <- gedit(width = 5, expand = FALSE , container = f1)
  
  #----------------------------------------------------------------------------
  glabel(text = "Post processing:", anchor=c(-1 ,0), container = f1)
  
  f1_h_chk <- gcheckbox(text="Calculate average peak height", checked = TRUE,
                        container = f1)

  # Disable checkbox to calculate H.
  enabled(f1_h_chk) <- FALSE
  
  # FRAME 4 ###################################################################
  
  if(debug){
    print("FRAME 4")
  }  

  f4 <- gframe(text = "Save as",
               horizontal = TRUE,
               spacing = 5,
               container = gv) 
  
  glabel(text="Name for result:", container = f4)
  
  f4_save_edt <- gedit(expand = TRUE, container = f4)
  
  # BUTTON ####################################################################

  if(debug){
    print("BUTTON")
  }  
  
  calculate_btn <- gbutton(text = "Calculate",
                      border = TRUE,
                      container = gv)
  
  addHandlerChanged(calculate_btn, handler = function(h, ...) {
    
    # Get values.
    val_option <- svalue(f1_lb_opt, index=TRUE)
    val_dye <- svalue(f1_dye_chk)
    val_ignore <- svalue(f1_ignore_chk)
    val_word <- svalue(f1_word_chk)
    val_exact <- svalue(f1_exact_chk)
    val_data <- .gData
    val_data_name <- .gDataName
    val_ref <- .gRef
    val_ref_name <- .gRefName
    val_name <- svalue(f4_save_edt)
    val_kit <- svalue(kit_drp)
    val_ol <- svalue(f1_ol_chk)
    val_sex <- svalue(f1_sex_chk)
    val_na <- as.numeric(svalue(f1_na_edt))
    val_h_enabled <- enabled(f1_h_chk)
    val_h <- svalue(f1_h_chk)

    if(debug){
      print("Read Values:")
      print("val_option")
      print(val_option)
      print("val_dye")
      print(val_dye)
      print("val_ol")
      print(val_ol)
      print("val_sex")
      print(val_sex)
      print("val_na")
      print(val_na)
      print("val_ignore")
      print(val_ignore)
      print("val_word")
      print(val_word)
      print("val_exact")
      print(val_exact)
      print("val_name")
      print(val_name)
      print("val_data")
      print(head(val_data))
      print("val_ref")
      print(head(val_ref))
    }
    
    # Check if data.
    if(!is.null(.gData)){

      # Check for NA's in dye column.
      if(!any(is.na(.gData$Dye))){
        
        if(val_option == 1){
          val_option <- "prop"
        } else if(val_option == 2) {
          val_option <- "norm"
        } else if(val_option == 3) {
          val_option <- "cent"
        } else {
          stop("val_option =", val_option, "not implemented!")
        }

        if(is.na(val_na)){
          val_na <- NULL
        }
        
        if(!val_h_enabled){
            val_h <- FALSE
        }
        
        if(debug){
          print("Sent Values:")
          print("val_option")
          print(val_option)
          print("val_na")
          print(val_na)
          print("val_h")
          print(val_h)
          print("val_ignore")
          print(val_ignore)
        }
        
        # Change button.
        svalue(calculate_btn) <- "Processing..."
        enabled(calculate_btn) <- FALSE
        
        datanew <- calculateLb(data = val_data,
                               ref = val_ref,
                               option = val_option,
                               by.dye = val_dye,
                               ol.rm = val_ol,
                               sex.rm = val_sex,
                               na = val_na,
                               kit = val_kit,
                               ignore.case = val_ignore,
                               word = val_word,
                               exact = val_exact,
                               debug=debug)
        
        # Add attributes.
        attr(datanew, which="kit") <- val_kit
        attr(datanew, which="calculateLb_gui, data") <- .gDataName
        attr(datanew, which="calculateLb_gui, ref") <- .gRefName
        attr(datanew, which="calculateLb_gui, option") <- val_option
        attr(datanew, which="calculateLb_gui, by.dye") <- val_dye
        attr(datanew, which="calculateLb_gui, ol.rm") <- val_ol
        attr(datanew, which="calculateLb_gui, sex.rm") <- val_sex
        attr(datanew, which="calculateLb_gui, na") <- val_na
        attr(datanew, which="calculateLb_gui, ignore.case") <- val_ignore
        attr(datanew, which="calculateLb_gui, word") <- val_word
        attr(datanew, which="calculateLb_gui, exact") <- val_exact
        attr(datanew, which="calculateLb_gui, calculate.h") <- val_h
        
        # Calculate and add average peak height.
        if(val_h){
          
          # Heterozygote status is required to calculate 'H'.        
          if(!"Heterozygous" %in% names(val_data)){
            
            if(!"Heterozygous" %in% names(val_ref)){
              
              # Calculate heterozygote indicator for reference set.
              val_ref <- calculateHeterozygous(data=val_ref, debug=debug)
              
              message("Heterozygote indicator calculated for reference set.")
              
            }
            
            # Filter known profile.
            val_data <- filterProfile(data = val_data, ref = val_ref,
                                      add.missing.loci = TRUE, keep.na = TRUE,
                                      ignore.case = val_ignore,
                                      exact = val_exact,
                                      invert = FALSE, debug = debug)
            
            message("Filter known profile from dataset.")
            
            # Add heterozygote indicator to dataset.
            val_data <- addData(data = val_data, new.data = val_ref,
                                by.col = "Sample.Name", then.by.col = "Marker",
                                exact = FALSE, ignore.case = val_ignore,
                                debug = debug)
            
            message("Heterozygote indicator added to dataset.")
            
          }
          
          # Calculate average peak height.
          dfH <- calculateHeight(data = val_data, na = 0, add = FALSE,
                                 exclude = "OL", debug = debug)
          
          message("Average peak height calculated.")
          
          # Add average peak height to dataset.
          datanew <- addData(data = datanew, new.data = dfH,
                             by.col = "Sample.Name", then.by.col=NULL,
                             exact = TRUE, ignore.case=val_ignore,
                             debug=debug)
          
          message("Average peak height added to result.")
          
        }
        
        # Save data.
        saveObject(name = val_name, object = datanew, parent = w, env = env)
        
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

      message <- "A dataset must be selected."
      
      gmessage(message, title = "Datasets not selected",
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
      if(exists(".strvalidator_calculateLb_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_calculateLb_gui_savegui", envir=env)
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
      if(exists(".strvalidator_calculateLb_gui_option", envir=env, inherits = FALSE)){
        svalue(f1_lb_opt) <- get(".strvalidator_calculateLb_gui_option", envir=env)
      }
      if(exists(".strvalidator_calculateLb_gui_dye", envir=env, inherits = FALSE)){
        svalue(f1_dye_chk) <- get(".strvalidator_calculateLb_gui_dye", envir=env)
      }
      if(exists(".strvalidator_calculateLb_gui_ol", envir=env, inherits = FALSE)){
        svalue(f1_ol_chk) <- get(".strvalidator_calculateLb_gui_ol", envir=env)
      }
      if(exists(".strvalidator_calculateLb_gui_sex", envir=env, inherits = FALSE)){
        svalue(f1_sex_chk) <- get(".strvalidator_calculateLb_gui_sex", envir=env)
      }
      if(exists(".strvalidator_calculateLb_gui_na", envir=env, inherits = FALSE)){
        svalue(f1_na_edt) <- get(".strvalidator_calculateLb_gui_na", envir=env)
      }
      if(exists(".strvalidator_calculateLb_gui_ignore", envir=env, inherits = FALSE)){
        svalue(f1_ignore_chk) <- get(".strvalidator_calculateLb_gui_ignore", envir=env)
      }
      if(exists(".strvalidator_calculateLb_gui_word", envir=env, inherits = FALSE)){
        svalue(f1_word_chk) <- get(".strvalidator_calculateLb_gui_word", envir=env)
      }
      if(exists(".strvalidator_calculateLb_gui_exact", envir=env, inherits = FALSE)){
        svalue(f1_exact_chk) <- get(".strvalidator_calculateLb_gui_exact", envir=env)
      }
      if(exists(".strvalidator_calculateLb_gui_h", envir=env, inherits = FALSE)){
        svalue(f1_h_chk) <- get(".strvalidator_calculateLb_gui_h", envir=env)
      }
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_calculateLb_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_calculateLb_gui_option", value=svalue(f1_lb_opt), envir=env)
      assign(x=".strvalidator_calculateLb_gui_dye", value=svalue(f1_dye_chk), envir=env)
      assign(x=".strvalidator_calculateLb_gui_ol", value=svalue(f1_ol_chk), envir=env)
      assign(x=".strvalidator_calculateLb_gui_sex", value=svalue(f1_sex_chk), envir=env)
      assign(x=".strvalidator_calculateLb_gui_na", value=svalue(f1_na_edt), envir=env)
      assign(x=".strvalidator_calculateLb_gui_ignore", value=svalue(f1_ignore_chk), envir=env)
      assign(x=".strvalidator_calculateLb_gui_word", value=svalue(f1_word_chk), envir=env)
      assign(x=".strvalidator_calculateLb_gui_exact", value=svalue(f1_exact_chk), envir=env)
      assign(x=".strvalidator_calculateLb_gui_h", value=svalue(f1_h_chk), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_calculateLb_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateLb_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_calculateLb_gui_option", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateLb_gui_option", envir = env)
      }
      if(exists(".strvalidator_calculateLb_gui_dye", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateLb_gui_dye", envir = env)
      }
      if(exists(".strvalidator_calculateLb_gui_ol", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateLb_gui_ol", envir = env)
      }
      if(exists(".strvalidator_calculateLb_gui_sex", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateLb_gui_sex", envir = env)
      }
      if(exists(".strvalidator_calculateLb_gui_na", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateLb_gui_na", envir = env)
      }
      if(exists(".strvalidator_calculateLb_gui_ignore", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateLb_gui_ignore", envir = env)
      }
      if(exists(".strvalidator_calculateLb_gui_word", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateLb_gui_word", envir = env)
      }
      if(exists(".strvalidator_calculateLb_gui_exact", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateLb_gui_exact", envir = env)
      }
      if(exists(".strvalidator_calculateLb_gui_h", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateLb_gui_h", envir = env)
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
