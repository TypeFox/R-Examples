################################################################################
# TODO LIST
# TODO: ...


################################################################################
# CHANGE LOG (last 20 changes)
# 05.01.2016: Fixed 'dev' not find error in ggplot2 2.0.
# 29.08.2015: Added importFrom.
# 11.10.2014: Added 'focus'.
# 06.10.2014: Correct pixel dimensions are now shown.
# 28.06.2014: Added help button and moved save gui checkbox.
# 25.02.2014: Pixel info now update when textbox is changed.
# 09.02.2014: Added info for size in pixel.
# 09.02.2014: Removed unsupported unit 'px'.
# 20.01.2014: First version.

#' @title Save Image
#'
#' @description
#' A simple GUI wrapper for \code{\link{ggsave}}.
#'
#' @details
#' Simple GUI wrapper for ggsave.
#' 
#' @param ggplot plot object.
#' @param name optional string providing a file name.
#' @param env environment where the objects exist.
#' Default is the current environment.
#' @param savegui logical indicating if GUI settings should be saved in the environment.
#' @param debug logical indicating printing debug information.
#' @param parent object specifying the parent widget to center the message box,
#' and to get focus when finished.
#' 
#' @return TRUE
#' 
#' @export
#' 
#' @importFrom ggplot2 ggsave
#' @importFrom utils help
#' @importFrom grDevices dev.cur dev.list dev.size
#' 
#' @seealso \code{\link{ggsave}}

ggsave_gui <- function(ggplot=NULL, name="", env=parent.frame(),
                          savegui=NULL, debug=FALSE, parent=NULL){
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
    print("Current device")
    print(dev.cur())
    print("Device list")
    print(dev.list())
  }
  
  # Constants.
  .separator <- .Platform$file.sep # Platform dependent path separator.
   
  # Main window.
  w <- gwindow(title="Save as image",
               visible=FALSE)
  
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
    print(help("ggsave_gui", help_type="html"))
    
  })
  
  # FRAME 1 ###################################################################
  
  f1 <- gframe(text = "Options",
               horizontal=FALSE,
               spacing = 10,
               container = gv) 
  
  # GRID 1 --------------------------------------------------------------------
  
  f1g1 <- glayout(container = f1, spacing = 2)
  
  f1g1[2,1] <- glabel(text="File name and extension:",
                      container=f1g1,
                      anchor=c(-1 ,0))
  
  f1g1[3,1] <- f1g1_name_edt <- gedit(text=name, width=50, container=f1g1)
  
  f1g1[3,2] <- f1g1_ext_drp <- gdroplist(items=c("eps", "ps", "tex", "pdf",
                                                 "jpeg", "tiff", "png",
                                                 "bmp", "svg", "wmf"),
                                          selected=4,
                                          container=f1g1)
  
  f1g1[4,1] <- f1g1_replace_chk <- gcheckbox(text="Overwrite existing file",
                                             checked = TRUE,
                                             container = f1g1)
  
  f1g1[5,1] <- f1g1_load_chk <- gcheckbox(text="Load size from plot device",
                                            checked=TRUE,
                                            container=f1g1)

  f1g1[6,1] <- f1g1_get_btn <- gbutton(text="Get size", container=f1g1)
  
  addHandlerChanged(f1g1_load_chk, handler = function(h, ...) {
    
    val <- svalue(f1g1_load_chk)
    
    if(val){
      
      # Read size from device.
      .readSize()
      
    } else {
      
      # Could load saved settings...
      
    }
    
  })
  
  addHandlerChanged(f1g1_get_btn, handler = function(h, ...) {
    
    # Read size from device.
    .readSize()
    
  })
  
  # GRID 2 --------------------------------------------------------------------
  
  f1g2 <- glayout(container = f1, spacing = 2)
  
  f1g2[1,1] <- glabel(text="Image settings", container=f1g2, anchor=c(-1 ,0))
  
  f1g2[2,1] <- glabel(text="Unit:",
                      container=f1g2,
                      anchor=c(-1 ,0))
  
  f1g2[2,2] <- f1g2_unit_drp <- gdroplist(items=c("in", "cm"),
                                       selected=2,
                                       container=f1g2)

  addHandlerChanged(f1g2_unit_drp, handler = function(h, ...) {
    
    # Read size from device.
    .readSize()
    
  })
                    
  f1g2[3,1] <- glabel(text="Width:", container=f1g2, anchor=c(-1 ,0))
  
  f1g2[3,2] <- f1g2_width_edt <- gedit(text="",
                                       width=6,
                                       initial.msg="",
                                       container=f1g2)
  
  f1g2[3,3] <- f1g2_width_lbl <- glabel(text=" NA pixels",
                                        container=f1g2, anchor=c(-1 ,0))
  
  f1g2[4,1] <- glabel(text="Height:", container=f1g2, anchor=c(-1 ,0))
  
  f1g2[4,2] <- f1g2_height_edt <- gedit(text="",
                                        width=6,
                                        initial.msg="",
                                        container=f1g2)
  
  f1g2[4,3] <- f1g2_height_lbl <- glabel(text=" NA pixels",
                                         container=f1g2, anchor=c(-1 ,0))

  f1g2[5,1] <- glabel(text="Resolution:", container=f1g2, anchor=c(-1 ,0))
  
  f1g2[5,2] <- f1g2_res_edt <- gedit(text="300",
                                     width=4,
                                     initial.msg="",
                                     container=f1g2)
  
  f1g2[6,1] <- glabel(text="Scaling factor:", container=f1g2, anchor=c(-1 ,0))
  
  f1g2[6,2] <- f1g2_scale_edt <- gedit(text="1",
                                     width=4,
                                     initial.msg="",
                                     container=f1g2)
  
  addHandlerKeystroke(f1g2_width_edt, handler = function(h, ...) {
      
    # Get values.    
    val_w <- as.numeric(svalue(f1g2_width_edt))
    val_u <- svalue(f1g2_unit_drp)
    val_d <- as.numeric(svalue(f1g2_res_edt))
    val_s <- as.numeric(svalue(f1g2_scale_edt))
    
    # Convert to pixel.
    pixels <- .toPixel(unit=val_u, val=val_w, dpi=val_d, scale=val_s)

    # Update label.
    svalue(f1g2_width_lbl) <- paste(" ", pixels,"pixels")
    
  })

  addHandlerKeystroke(f1g2_height_edt, handler = function(h, ...) {
    
    # Get values.    
    val_h <- as.numeric(svalue(f1g2_height_edt))
    val_u <- svalue(f1g2_unit_drp)
    val_d <- as.numeric(svalue(f1g2_res_edt))
    val_s <- as.numeric(svalue(f1g2_scale_edt))
    
    # Convert to pixel.
    pixels <- .toPixel(unit=val_u, val=val_h, dpi=val_d, scale=val_s)
    
    svalue(f1g2_height_lbl) <- paste(" ", pixels,"pixels")
    
  })

  addHandlerChanged(f1g2_width_edt, handler = function(h, ...) {
    
    # Get values.    
    val_w <- as.numeric(svalue(f1g2_width_edt))
    val_u <- svalue(f1g2_unit_drp)
    val_d <- as.numeric(svalue(f1g2_res_edt))
    val_s <- as.numeric(svalue(f1g2_scale_edt))
    
    # Convert to pixel.
    pixels <- .toPixel(unit=val_u, val=val_w, dpi=val_d, scale=val_s)
    
    # Update label.
    svalue(f1g2_width_lbl) <- paste(" ", pixels,"pixels")
    
  })
  
  addHandlerChanged(f1g2_height_edt, handler = function(h, ...) {
    
    val_h <- as.numeric(svalue(f1g2_height_edt))
    val_u <- svalue(f1g2_unit_drp)
    val_d <- as.numeric(svalue(f1g2_res_edt))
    val_s <- as.numeric(svalue(f1g2_scale_edt))
    
    # Convert to pixel.
    pixels <- .toPixel(unit=val_u, val=val_h, dpi=val_d, scale=val_s)
    
    svalue(f1g2_height_lbl) <- paste(" ", pixels,"pixels")
    
  })
  
  addHandlerKeystroke(f1g2_res_edt, handler = function(h, ...) {
    
    # Get values.    
    val_w <- as.numeric(svalue(f1g2_width_edt))
    val_h <- as.numeric(svalue(f1g2_height_edt))
    val_u <- svalue(f1g2_unit_drp)
    val_d <- as.numeric(svalue(f1g2_res_edt))
    val_s <- as.numeric(svalue(f1g2_scale_edt))
    
    # Convert to pixel.
    pixels_w <- .toPixel(unit=val_u, val=val_w, dpi=val_d, scale=val_s)
    pixels_h <- .toPixel(unit=val_u, val=val_h, dpi=val_d, scale=val_s)
    
    # Update label.
    svalue(f1g2_width_lbl) <- paste(" ", pixels_w,"pixels")
    svalue(f1g2_height_lbl) <- paste(" ", pixels_h,"pixels")
    
  })
  
  addHandlerChanged(f1g2_res_edt, handler = function(h, ...) {
    
    # Get values.    
    val_w <- as.numeric(svalue(f1g2_width_edt))
    val_h <- as.numeric(svalue(f1g2_height_edt))
    val_u <- svalue(f1g2_unit_drp)
    val_d <- as.numeric(svalue(f1g2_res_edt))
    val_s <- as.numeric(svalue(f1g2_scale_edt))
    
    # Convert to pixel.
    pixels_w <- .toPixel(unit=val_u, val=val_w, dpi=val_d, scale=val_s)
    pixels_h <- .toPixel(unit=val_u, val=val_h, dpi=val_d, scale=val_s)
    
    # Update label.
    svalue(f1g2_width_lbl) <- paste(" ", pixels_w,"pixels")
    svalue(f1g2_height_lbl) <- paste(" ", pixels_h,"pixels")
    
  })

  addHandlerKeystroke(f1g2_scale_edt, handler = function(h, ...) {
    
    # Get values.    
    val_w <- as.numeric(svalue(f1g2_width_edt))
    val_h <- as.numeric(svalue(f1g2_height_edt))
    val_u <- svalue(f1g2_unit_drp)
    val_d <- as.numeric(svalue(f1g2_res_edt))
    val_s <- as.numeric(svalue(f1g2_scale_edt))
    
    # Convert to pixel.
    pixels_w <- .toPixel(unit=val_u, val=val_w, dpi=val_d, scale=val_s)
    pixels_h <- .toPixel(unit=val_u, val=val_h, dpi=val_d, scale=val_s)
    
    # Update label.
    svalue(f1g2_width_lbl) <- paste(" ", pixels_w,"pixels")
    svalue(f1g2_height_lbl) <- paste(" ", pixels_h,"pixels")
    
  })
  
  addHandlerChanged(f1g2_scale_edt, handler = function(h, ...) {
    
    # Get values.    
    val_w <- as.numeric(svalue(f1g2_width_edt))
    val_h <- as.numeric(svalue(f1g2_height_edt))
    val_u <- svalue(f1g2_unit_drp)
    val_d <- as.numeric(svalue(f1g2_res_edt))
    val_s <- as.numeric(svalue(f1g2_scale_edt))
    
    # Convert to pixel.
    pixels_w <- .toPixel(unit=val_u, val=val_w, dpi=val_d, scale=val_s)
    pixels_h <- .toPixel(unit=val_u, val=val_h, dpi=val_d, scale=val_s)
    
    # Update label.
    svalue(f1g2_width_lbl) <- paste(" ", pixels_w,"pixels")
    svalue(f1g2_height_lbl) <- paste(" ", pixels_h,"pixels")
    
  })
  
  # GRID 3 --------------------------------------------------------------------
  
  f1g3 <- glayout(container = f1, spacing = 5)
  
  f1g3[1,1] <- glabel(text="File path:",
                      container=f1g3,
                      anchor=c(-1 ,0))
  
  f1g3[2,1:2] <- f1g3_save_brw <- gfilebrowse(text=getwd(),
                                              quote=FALSE,
                                              type="selectdir",
                                              container=f1g3)
  
  # BUTTON ####################################################################
  
  g_save_btn <- gbutton(text="Save",
                          border=TRUE,
                          container=gv) 
  
  # HANDLERS ##################################################################
  

  
  addHandlerChanged(g_save_btn, handler = function(h, ...) {
    
    # Get values.
    val_name <- svalue(f1g1_name_edt)
    val_ggplot <- ggplot
    val_ext <- paste(".", svalue(f1g1_ext_drp), sep="")
    val_scale <- as.numeric(svalue(f1g2_scale_edt))
    val_unit <- svalue(f1g2_unit_drp)
    val_replace <- svalue(f1g1_replace_chk)
    val_w <- as.numeric(svalue(f1g2_width_edt))
    val_h <- as.numeric(svalue(f1g2_height_edt))
    val_r <- as.numeric(svalue(f1g2_res_edt))
    val_path <- svalue(f1g3_save_brw)
    
    # Check file name.
    if(nchar(val_name) == 0){
      val_name <- NA
    }
    
    # Check path.
    if(nchar(val_path) == 0){
      val_path <- NA
    }
    
    if(debug){
      print("val_name")
      print(val_name)
      print("val_ext")
      print(val_ext)
      print("val_replace")
      print(val_replace)
      print("val_w")
      print(val_w)
      print("val_h")
      print(val_h)
      print("val_r")
      print(val_r)
      print("val_path")
      print(val_path)
    }
    
    # Check for file name and path.
    ok <- !is.na(val_name) && !is.na(val_path) && !is.null(val_ggplot)
    
    if(ok){
      
      svalue(g_save_btn) <- "Processing..."

      # Add trailing path separator if not present.
      if(substr(val_path, nchar(val_path), nchar(val_path)+1) != .separator){
        val_path <- paste(val_path, .separator, sep="")
      }
      
      # Repeat until saved or cancel.
      okToSave <- FALSE
      cancel <- FALSE
      repeat{
        
        # Construct complete file name.
        fullFileName <- paste(val_path, val_name, val_ext, sep="")
        
        if(val_replace){
          # Ok to overwrite.
          okToSave <- TRUE
          
          if(debug){
            print("Replace=TRUE. Ok to save!")
          }
          
        } else {
          # Not ok to overwrite.
          
          if(debug){
            print("Replace=FALSE. Check if file exist!")
          }
          
          # Check if file exist.
          if(file.exists(fullFileName)){
            
            if(debug){
              
              print(paste("file '", name, "' already exist!", sep=""))
              
            }
            
            # Create dialog.
            dialog <- gbasicdialog(title="Save error", parent=w,
                                   do.buttons=FALSE, width=200,
                                   height=200, horizontal=FALSE)
            
            glabel(text="The file already exist!",
                   anchor=c(-1 ,0), container=dialog)
            
            glabel(text="Chose to cancel, overwrite or give a new name.",
                   anchor=c(-1 ,0), container=dialog)
            
            # Edit box for new name.
            newName <- gedit(container=dialog)
            
            # Container for buttons.
            gg <- ggroup(container=dialog) 
            
            btn_cancel <- gbutton("Cancel", container = gg, handler = function(h, ...) {
              cancel <<- TRUE
              dispose(dialog)
            })
            
            btn_replace <- gbutton("Overwrite", container = gg, handler = function(h, ...) {
              val_replace <<- TRUE
              dispose(dialog)
            })
            
            btn_retry <- gbutton("Retry", container = gg, handler = function(h, ...) {
              val_name <<- svalue(newName)
              if(debug){
                print("val_name")
                print(val_name)
              }
              dispose(dialog)
            })

            # Show dialog.
            visible(dialog, set=TRUE)
            
          } else {
            okToSave <- TRUE
          }
          
        }
        
        if(cancel){
          # Chose to cancel.
          
          if(debug){
            print("Chose to cancel!")
          }
          
          break ## EXIT REPEAT.
          
        }
        
        if(okToSave){
          
          # Save plot device as image.
          ggsave(filename = paste(val_name, val_ext, sep=""),
                 plot = val_ggplot,
                 path = val_path,
                 scale = val_scale,
                 width = val_w, height = val_h,
                 units = val_unit, dpi = val_r)
          
          
          if(debug){
            print("Image saved!")
          }
          
          break ## EXIT REPEAT.
          
        }
        
      }  ## END REPEAT.
      
      # Close GUI.
      dispose(w)
      
    } else {
      
      gmessage(message="Plot object, file name and path must be provided.",
               title="Error",
               parent=w,
               icon = "error")      
    }    
  } )
  
  # INTERNAL FUNCTIONS ########################################################
  
  .toPixel <- function(unit, val, dpi=72, scale=1){
    
    # Convert to pixel.
    if(unit == "cm"){
      pixels <- (val / (2.54 / dpi)) * scale
    } else if(unit == "in"){
      pixels <- val * dpi * scale
    } else {
      pixels <- NA
    }
    
    return(round(pixels,0))
    
  }
  
  .readSize <- function() {
    
    # Get values.
    val_unit <- svalue(f1g2_unit_drp)
    val_size <- round(dev.size(val_unit),2)
    val_px <- dev.size("px")
    
    # Update.
    svalue(f1g2_width_edt) <- val_size[1]
    svalue(f1g2_height_edt) <- val_size[2]
    
#     svalue(f1g2_width_lbl) <- paste(" ", val_px[1],"pixels")
#     svalue(f1g2_height_lbl) <- paste(" ", val_px[2],"pixels")
    
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
      if(exists(".strvalidator_ggsave_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_ggsave_gui_savegui", envir=env)
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
      if(exists(".strvalidator_ggsave_gui_ext", envir=env, inherits = FALSE)){
        svalue(f1g1_ext_drp) <- get(".strvalidator_ggsave_gui_ext", envir=env)
      }
      if(exists(".strvalidator_ggsave_gui_replace", envir=env, inherits = FALSE)){
        svalue(f1g1_replace_chk) <- get(".strvalidator_ggsave_gui_replace", envir=env)
      }
      if(exists(".strvalidator_ggsave_gui_load", envir=env, inherits = FALSE)){
        svalue(f1g1_load_chk) <- get(".strvalidator_ggsave_gui_load", envir=env)
      }
      if(exists(".strvalidator_ggsave_gui_unit", envir=env, inherits = FALSE)){
        svalue(f1g2_unit_drp) <- get(".strvalidator_ggsave_gui_unit", envir=env)
      }
      if(exists(".strvalidator_ggsave_gui_width", envir=env, inherits = FALSE)){
        svalue(f1g2_width_edt) <- get(".strvalidator_ggsave_gui_width", envir=env)
      }
      if(exists(".strvalidator_ggsave_gui_height", envir=env, inherits = FALSE)){
        svalue(f1g2_height_edt) <- get(".strvalidator_ggsave_gui_height", envir=env)
      }
      if(exists(".strvalidator_ggsave_gui_res", envir=env, inherits = FALSE)){
        svalue(f1g2_res_edt) <- get(".strvalidator_ggsave_gui_res", envir=env)
      }
      if(exists(".strvalidator_ggsave_gui_scale", envir=env, inherits = FALSE)){
        svalue(f1g2_scale_edt) <- get(".strvalidator_ggsave_gui_scale", envir=env)
      }
#       if(exists(".strvalidator_ggsave_gui_path", envir=env, inherits = FALSE)){
#         svalue(f1g3_save_brw) <- get(".strvalidator_ggsave_gui_path", envir=env)
#       }
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_ggsave_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_ggsave_gui_ext", value=svalue(f1g1_ext_drp), envir=env)
      assign(x=".strvalidator_ggsave_gui_replace", value=svalue(f1g1_replace_chk), envir=env)
      assign(x=".strvalidator_ggsave_gui_load", value=svalue(f1g1_load_chk), envir=env)
      assign(x=".strvalidator_ggsave_gui_unit", value=svalue(f1g2_unit_drp), envir=env)
      assign(x=".strvalidator_ggsave_gui_width", value=svalue(f1g2_width_edt), envir=env)
      assign(x=".strvalidator_ggsave_gui_height", value=svalue(f1g2_height_edt), envir=env)
      assign(x=".strvalidator_ggsave_gui_res", value=svalue(f1g2_res_edt), envir=env)
      assign(x=".strvalidator_ggsave_gui_scale", value=svalue(f1g2_scale_edt), envir=env)
#       assign(x=".strvalidator_ggsave_gui_path", value=svalue(f1g3_save_brw), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_ggsave_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_ggsave_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_ggsave_gui_ext", envir=env, inherits = FALSE)){
        remove(".strvalidator_ggsave_gui_ext", envir = env)
      }
      if(exists(".strvalidator_ggsave_gui_replace", envir=env, inherits = FALSE)){
        remove(".strvalidator_ggsave_gui_replace", envir = env)
      }
      if(exists(".strvalidator_ggsave_gui_load", envir=env, inherits = FALSE)){
        remove(".strvalidator_ggsave_gui_load", envir = env)
      }
      if(exists(".strvalidator_ggsave_gui_unit", envir=env, inherits = FALSE)){
        remove(".strvalidator_ggsave_gui_unit", envir = env)
      }
      if(exists(".strvalidator_ggsave_gui_width", envir=env, inherits = FALSE)){
        remove(".strvalidator_ggsave_gui_width", envir = env)
      }
      if(exists(".strvalidator_ggsave_gui_height", envir=env, inherits = FALSE)){
        remove(".strvalidator_ggsave_gui_height", envir = env)
      }
      if(exists(".strvalidator_ggsave_gui_res", envir=env, inherits = FALSE)){
        remove(".strvalidator_ggsave_gui_res", envir = env)
      }
      if(exists(".strvalidator_ggsave_gui_scale", envir=env, inherits = FALSE)){
        remove(".strvalidator_ggsave_gui_scale", envir = env)
      }
#       if(exists(".strvalidator_ggsave_gui_path", envir=env, inherits = FALSE)){
#         remove(".strvalidator_ggsave_gui_path", envir = env)
#       }
      
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
  
  # Read size.
  if(svalue(f1g1_load_chk)){
    .readSize()
  }
  
  # Show GUI.
  visible(w) <- TRUE
  focus(w)
  
}
