################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 30.12.2015: Added option for 'exact' matching.
# 30.12.2015: Fixed option 'word' matching not saved.
# 13.11.2015: Added attribute drop.sex.
# 13.11.2015: Added option to calculate Hb as LMW / HMW.
# 08.11.2015: Added automatic calculation of average peak height 'H'.
# 21.10.2015: Added attributes.
# 28.08.2015: Added importFrom
# 17.08.2015: Changed erroneus  to  option 'High peak / low peak' to 'Smaller peak / larger peak'.
# 08.06.2015: Added option to drop sex markers (Fixes issue#9).
# 05.05.2015: Changed parameter 'ignoreCase' to 'ignore.case' for 'checkSubset' function.
# 13.12.2014: Added kit dropdown and kit attribute to result.
# 11.10.2014: Added 'focus', added 'parent' parameter.
# 03.10.2014: Added 'word' parameter (word boundary).
# 07.08.2014: Added check and error message for 'NA' in 'Dye'.
# 28.06.2014: Added help button and moved save gui checkbox.
# 06.05.2014: Implemented 'checkDataset'.
# 23.02.2014: Removed 'perSample' parameter. Added 'OL' check.
# 18.12.2013: Fixed dropdown state not saved (wrong object saved)...
# 27.11.2013: Passed debug to calculateBalance.
# 21.10.2013: Fixed dropdown state not loaded.
# 09.09.2013: Added option 'hb' to specify the definition of Hb.

#' @title Calculate Balance
#'
#' @description
#' GUI wrapper for the \code{\link{calculateBalance}} function.
#'
#' @details
#' Simplifies the use of the \code{\link{calculateBalance}} function
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
#' @seealso \code{link{calculateBalance}}, \code{link{checkSubset}}

calculateBalance_gui <- function(env=parent.frame(), savegui=NULL,
                                 debug=FALSE, parent=NULL){
  
  # Global variables.
  .gData <- NULL
  .gRef <- NULL
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # WINDOW ####################################################################
  
  if(debug){
    print("WINDOW")
  }  

  # Main window.
  w <- gwindow(title="Calculate balance", visible=FALSE)

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
    print(help("calculateBalance_gui", help_type="html"))
    
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
    requiredCol <- c("Sample.Name", "Marker", "Dye", "Height")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       slim=TRUE, slimcol="Height",
                       env=env, parent=w, debug=debug)
    
    if(ok){
      # Load or change components.
      
      # get dataset.
      .gData <<- get(val_obj, envir=env)
      svalue(g0_data_samples_lbl) <- paste(length(unique(.gData$Sample.Name)),
                                        "samples.")
      
      # Suggest a name for the result.
      svalue(f4_save_edt) <- paste(val_obj, "_balance", sep="")
      
      # Detect kit.
      kitIndex <- detectKit(.gData, index=TRUE)
      # Select in dropdown.
      svalue(f4_kit_drp, index=TRUE) <- kitIndex
      
    } else {
      
      # Reset components.
      .gData <<- NULL
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
      svalue(g0_ref_samples_lbl) <- paste(length(unique(.gRef$Sample.Name)),
                                          "samples.")
        
    } else {
      
      # Reset components.
      .gRef <<- NULL
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

  f1_ignore_chk <- gcheckbox(text="Ignore case", checked=TRUE,
                         container=f1)
  
  f1_word_chk <- gcheckbox(text="Add word boundaries", checked = FALSE,
                           container = f1)
  
  f1_drop_chk <- gcheckbox(text="Drop sex markers", checked=TRUE,
                           container=f1)
  
  f1_h_chk <- gcheckbox(text="Calculate average peak height", checked = TRUE,
                        container = f1)

  f1g1 <- ggroup(horizontal = TRUE, spacing = 5, container = f1)
  glabel(text="Calculate balance using:", anchor=c(-1 ,0), container=f1g1)
  f1_methods <- c("High molecular weight / low molecular weight",
                  "Low molecular weight / high molecular weight",
                  "Smaller peak / larger peak")
  f1_method_drp <- gdroplist(items=f1_methods,
                             selected = 1,
                             expand = FALSE,
                             container = f1g1)
  
  f1_options_lb <- c("Calculate proportional locus balance",
                "Calculate normalised locus balance")
  
  f1_lb_opt <- gradio(items=f1_options_lb,
                      selected=1,
                      horizontal=FALSE,
                      container=f1)

  f1_options_perDye <- c("Calculate locus balance within each dye",
                "Calculate locus balance globally across all dyes")
  
  f1_perDye_opt <- gradio(items=f1_options_perDye,
                          selected=2,
                          horizontal=FALSE,
                          container=f1)

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
  
  glabel(text=" Kit attribute:", container=f4)
  
  f4_kit_drp <- gdroplist(items=getKit(), selected = 1,
                       editable = FALSE, container = f4) 
  

  # BUTTON ####################################################################

  if(debug){
    print("BUTTON")
  }  
  
  calculate_btn <- gbutton(text="Calculate",
                      border=TRUE,
                      container=gv)
  
  addHandlerChanged(calculate_btn, handler = function(h, ...) {
    
    # Get values.
    val_method <- svalue(f1_method_drp, index=TRUE)
    val_lb <- svalue(f1_lb_opt, index=TRUE)
    val_perDye <- svalue(f1_perDye_opt, index=TRUE) == 1 # TRUE / FALSE
    val_ignore <- svalue(f1_ignore_chk)
    val_word <- svalue(f1_word_chk)
    val_h <- svalue(f1_h_chk)
    val_data <- .gData
    val_ref <- .gRef
    val_name <- svalue(f4_save_edt)
    val_kit <- svalue(f4_kit_drp)
    val_drop <- svalue(f1_drop_chk)
    
    if(debug){
      print("Read Values:")
      print("val_lb")
      print(val_lb)
      print("val_perDye")
      print(val_perDye)
      print("val_ignore")
      print(val_ignore)
      print("val_word")
      print(val_word)
      print("val_name")
      print(val_name)
      print("val_data")
      print(head(val_data))
      print("val_ref")
      print(head(val_ref))
      print("val_drop")
      print(val_drop)
    }
    
    # Check if data.
    if(!is.null(.gData) & !is.null(.gRef)){

      # Check for NA's in dye column.
      if(!any(is.na(.gData$Dye))){
        
        if(val_lb == 1){
          val_lb <- "prop"
        } else if(val_lb == 2) {
          val_lb <- "norm"
        } else {
          stop("val_lb =", val_lb, "not implemented!")
        }
  
        if(debug){
          print("Sent Values:")
          print("val_lb")
          print(val_lb)
          print("val_perDye")
          print(val_perDye)
          print("val_ignore")
          print(val_ignore)
        }
        
        if("OL" %in% .gData$Allele){
          
          message <- "'OL' alleles might lead to erroneous results!"
          
          gmessage(message, title="'OL' detected in dataset",
                   icon = "warning",
                   parent = w)
          
        }
        
        # Drop sex markers.
        if(val_drop){
          
          # Get sex marker.
          sexMarkers <- getKit(val_kit, what="Sex.Marker")
          
          # Check if sexMarkers was found.
          if(length(sexMarkers) > 0){
            
            # Drop sex markers.
            n0 <- nrow(.gData)
            for(m in seq(along=sexMarkers)){
              .gData <- .gData[.gData$Marker != sexMarkers[m], ]
            }
            n1 <- nrow(.gData)
            message(paste(n1, " rows after removing ", n0-n1, " sex marker rows.", sep=""))
            
          }
          
        }
    
        # Change button.
        svalue(calculate_btn) <- "Processing..."
        enabled(calculate_btn) <- FALSE
        
        datanew <- calculateBalance(data=val_data,
                                    ref=val_ref,
                                    lb=val_lb,
                                    per.dye=val_perDye,
                                    hb=val_method,
                                    ignore.case=val_ignore,
                                    word=val_word,
                                    debug=debug)
        
        # Add attributes.
        attr(datanew, which="kit") <- val_kit
        attr(datanew, which="calculateBalance_gui, data") <- svalue(g0_data_drp)
        attr(datanew, which="calculateBalance_gui, ref") <- svalue(g0_ref_drp)
        attr(datanew, which="calculateBalance_gui, lb") <- val_lb
        attr(datanew, which="calculateBalance_gui, hb") <- val_method
        attr(datanew, which="calculateBalance_gui, ignore.case") <- val_ignore
        attr(datanew, which="calculateBalance_gui, word") <- val_word
        attr(datanew, which="calculateBalance_gui, calculate.h") <- val_h
        attr(datanew, which="calculateBalance_gui, drop.sex") <- val_drop
        
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
            val_data <- filterProfile(data=val_data, ref=val_ref,
                                      add.missing.loci=TRUE, keep.na=TRUE,
                                      ignore.case=val_ignore,
                                      invert=FALSE, debug=debug)
            
            message("Filter known profile from dataset.")
            
            # Add heterozygote indicator to dataset.
            val_data <- addData(data=val_data, new.data=val_ref,
                                by.col="Sample.Name", then.by.col="Marker",
                                exact=FALSE, ignore.case=val_ignore,
                                debug=debug)
            
            message("Heterozygote indicator added to dataset.")
            
          }
          
          # Calculate average peak height.
          dfH <- calculateHeight(data=val_data, na=0, add=FALSE,
                                 exclude="OL", debug=debug)
          
          message("Average peak height calculated.")
          
          # Add average peak height to dataset.
          datanew <- addData(data=datanew, new.data=dfH,
                             by.col="Sample.Name", then.by.col=NULL,
                             exact=TRUE, ignore.case=val_ignore,
                             debug=debug)
          
          message("Average peak height added to result.")
          
        }
        
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
        
        message <- "'NA' in 'Dye' column. \nUse add dye function to fix."
        
        gmessage(message, title="NA detected!",
                 icon = "error",
                 parent = w)
        
      }
      
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
      if(exists(".strvalidator_calculateBalance_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_calculateBalance_gui_savegui", envir=env)
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
      if(exists(".strvalidator_calculateBalance_gui_method", envir=env, inherits = FALSE)){
        svalue(f1_method_drp) <- get(".strvalidator_calculateBalance_gui_method", envir=env)
      }
      if(exists(".strvalidator_calculateBalance_gui_lb", envir=env, inherits = FALSE)){
        svalue(f1_lb_opt) <- get(".strvalidator_calculateBalance_gui_lb", envir=env)
      }
      if(exists(".strvalidator_calculateBalance_gui_perDye", envir=env, inherits = FALSE)){
        svalue(f1_perDye_opt) <- get(".strvalidator_calculateBalance_gui_perDye", envir=env)
      }
      if(exists(".strvalidator_calculateBalance_gui_ignore", envir=env, inherits = FALSE)){
        svalue(f1_ignore_chk) <- get(".strvalidator_calculateBalance_gui_ignore", envir=env)
      }
      if(exists(".strvalidator_calculateBalance_gui_word", envir=env, inherits = FALSE)){
        svalue(f1_word_chk) <- get(".strvalidator_calculateBalance_gui_word", envir=env)
      }
      if(exists(".strvalidator_calculateBalance_gui_sex", envir=env, inherits = FALSE)){
        svalue(f1_drop_chk) <- get(".strvalidator_calculateBalance_gui_sex", envir=env)
      }
      if(exists(".strvalidator_calculateBalance_gui_h", envir=env, inherits = FALSE)){
        svalue(f1_h_chk) <- get(".strvalidator_calculateBalance_gui_h", envir=env)
      }
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_calculateBalance_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_calculateBalance_gui_method", value=svalue(f1_method_drp), envir=env)
      assign(x=".strvalidator_calculateBalance_gui_lb", value=svalue(f1_lb_opt), envir=env)
      assign(x=".strvalidator_calculateBalance_gui_perDye", value=svalue(f1_perDye_opt), envir=env)
      assign(x=".strvalidator_calculateBalance_gui_ignore", value=svalue(f1_ignore_chk), envir=env)
      assign(x=".strvalidator_calculateBalance_gui_word", value=svalue(f1_word_chk), envir=env)
      assign(x=".strvalidator_calculateBalance_gui_sex", value=svalue(f1_drop_chk), envir=env)
      assign(x=".strvalidator_calculateBalance_gui_h", value=svalue(f1_h_chk), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_calculateBalance_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateBalance_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_calculateBalance_gui_method", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateBalance_gui_method", envir = env)
      }
      if(exists(".strvalidator_calculateBalance_gui_lb", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateBalance_gui_lb", envir = env)
      }
      if(exists(".strvalidator_calculateBalance_gui_perDye", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateBalance_gui_perDye", envir = env)
      }
      if(exists(".strvalidator_calculateBalance_gui_ignore", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateBalance_gui_ignore", envir = env)
      }
      if(exists(".strvalidator_calculateBalance_gui_word", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateBalance_gui_word", envir = env)
      }
      if(exists(".strvalidator_calculateBalance_gui_sex", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateBalance_gui_sex", envir = env)
      }
      if(exists(".strvalidator_calculateBalance_gui_h", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateBalance_gui_h", envir = env)
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
