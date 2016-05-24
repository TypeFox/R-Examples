################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 29.08.2015: Added importFrom.
# 11.10.2014: Added 'focus', added 'parent' parameter.
# 28.06.2014: Added help button and moved save gui checkbox.
# 20.11.2013: Specified package for function 'gtable' -> 'gWidgets::gtable'
# 27.10.2013: Added warning when no object selected.
# 15.07.2013: Added save GUI settings.
# 10.07.2013: First version.

#' @title Export
#'
#' @description
#' GUI wrapper for the \code{\link{export}} function.
#'
#' @details
#' Simplifies the use of the \code{\link{export}} function by providing a graphical 
#' user interface to it.
#' 
#' @param env environment where the objects exist.
#' Default is the current environment.
#' @param savegui logical indicating if GUI settings should be saved in the environment.
#' @param debug logical indicating printing debug information.
#' @param parent widget to get focus when finished.
#' 
#' @return TRUE
#' 
#' @importFrom utils help
#' 
#' @seealso \code{\link{export}}


export_gui <- function(env=parent.frame(), savegui=NULL, debug=FALSE, parent=NULL){
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # Main window.
  w <- gwindow(title="Export objects as files or images",
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
    print(help("export_gui", help_type="html"))
    
  })
  
  # FRAME 0 ###################################################################
  
  f0 <- gframe(text = "Objects",
               horizontal=TRUE,
               spacing = 5,
               expand=TRUE,
               container = gv) 
  
  # Create list of objects.
  itemList <- listObjects(env=env, obj.class="data.frame")
  itemList <- c(itemList, listObjects(env=env, obj.class="ggplot"))
  
  f0_object_tbl <- gWidgets::gtable(items=itemList, multiple = TRUE,
                          expand=TRUE, container = f0)

  # Set initial minimal size.
  size(f0_object_tbl) <- c(100,150)

  
  # FRAME 1 ###################################################################
  
  f1 <- gframe(text = "File name",
               horizontal=TRUE,
               spacing = 5,
               container = gv) 
  
  # GRID 1 --------------------------------------------------------------------
  
  f1g1 <- glayout(container = f1)
  
  f1g1[1,1] <- f1g1_name_chk <- gcheckbox(text="Use object names",
                                          checked = TRUE,
                                          container = f1g1)
  
  f1g1[2,1] <- glabel(text="File name (separated by | ):",
                      container=f1g1,
                      anchor=c(-1 ,0))
  
  f1g1[3,1] <- f1g1_name_edt <- gedit(text="", width=50, container=f1g1)

  # Defult is disabled.
  enabled(f1g1_name_edt) <- FALSE
  
  f1g1[4,1] <- f1g1_replace_chk <- gcheckbox(text="Overwrite existing files",
                                          checked = TRUE,
                                          container = f1g1)
  
  # FRAME 2 ###################################################################
  
  f2 <- gframe(text = "Options",
               horizontal=FALSE,
               spacing = 15,
               container = gv) 
  
  # GRID 1 --------------------------------------------------------------------

  f2g1 <- glayout(container = f2)
  
  f2g1[1,1] <- glabel(text="File extension:",
                      container=f2g1,
                      anchor=c(-1 ,0))
  
  f2g1[1,2] <- f2g1_ext_drp <- gdroplist(items=c("auto", ".RData"),
                                         selected = 1,
                                         editable = FALSE,
                                         container = f2g1)
  
  f2g1[2,1] <- glabel(text="Delimiter:",
                      container=f2g1,
                      anchor=c(-1 ,0))
  
  f2g1[2,2] <- f2g1_del_drp <- gdroplist(items=c("TAB","SPACE","COMMA"),
                                         selected = 1,
                                         editable = FALSE,
                                         container = f2g1)

  # FRAME 3 ###################################################################

  f3 <- gframe(text="Image settings",
               horizontal=FALSE, spacing = 10, container = gv)
  
  # GRID 1 --------------------------------------------------------------------
  
  f3g1 <- glayout(container = f3, spacing = 5)
  
  f3g1[1,1] <- glabel(text="Width:", container=f3g1, anchor=c(-1 ,0))
  
  f3g1[1,2] <- f3g1_width_edt <- gedit(text="3000",
                                       width=4,
                                       initial.msg="",
                                       container=f3g1)
  
  f3g1[1,3] <- glabel(text="Height:", container=f3g1, anchor=c(-1 ,0))
  
  f3g1[1,4] <- f3g1_height_edt <- gedit(text="2000",
                                        width=4,
                                        initial.msg="",
                                        container=f3g1)
  
  f3g1[1,5] <- glabel(text="Resolution:", container=f3g1, anchor=c(-1 ,0))
  
  f3g1[1,6] <- f3g1_res_edt <- gedit(text="250",
                                     width=4,
                                     initial.msg="",
                                     container=f3g1)
  
  
  # FRAME 4 ###################################################################
  
  f4 <- gframe(text="Location",
               horizontal=FALSE, spacing = 10, container = gv)
  
  # GRID 1 --------------------------------------------------------------------
  
  f4g1 <- glayout(container = f4, spacing = 5)
  
  f4g1[1,1] <- glabel(text="File path:",
                      container=f4g1,
                      anchor=c(-1 ,0))

  expDefText <- "Select folder..."
  f4g1[2,1:2] <- f4g1_save_brw <- gfilebrowse(text=expDefText,
                                              quote=FALSE,
                                              type="selectdir",
                                              container=f4g1)
  
  # BUTTON ####################################################################
  
  g_export_btn <- gbutton(text="Export",
                           border=TRUE,
                           container=gv) 
  
  # HANDLERS ##################################################################
  
  addHandlerChanged(f1g1_name_chk, handler = function(h, ...) {
    
    # Get values.
    val <- svalue(f1g1_name_chk)
    
    if(val){
      enabled(f1g1_name_edt) <- FALSE
    } else {
      enabled(f1g1_name_edt) <- TRUE
    }
    
  } )

  addHandlerChanged(f2g1_ext_drp, handler = function(h, ...) {
    
    # Get values.
    val <- svalue(f2g1_ext_drp)
    
    if(val == ".RData"){
      enabled(f2g1_del_drp) <- FALSE
      enabled(f3g1) <- FALSE
    } else {
      enabled(f2g1_del_drp) <- TRUE
      enabled(f3g1) <- TRUE
    }
    
  } )
  
  addHandlerChanged(g_export_btn, handler = function(h, ...) {
    
    # Get values.
    val_object <- svalue(f0_object_tbl)
    val_use_obj <- svalue(f1g1_name_chk)
    val_name <- svalue(f1g1_name_edt)
    val_replace <- svalue(f1g1_replace_chk)
    val_ext <- svalue(f2g1_ext_drp)
    val_del <- svalue(f2g1_del_drp, index=TRUE)
    val_w <- as.numeric(svalue(f3g1_width_edt))
    val_h <- as.numeric(svalue(f3g1_height_edt))
    val_r <- as.numeric(svalue(f3g1_res_edt))
    val_path <- svalue(f4g1_save_brw)
    
    # Assign a delimiter character.
    if(val_del == 1){
      val_delimiter <- "\t"   
    } else if(val_del == 2){
      val_delimiter <- " "
    } else if(val_del == 3){
      val_delimiter <- ","
    } 

    # Check file name.
    if(nchar(val_name) == 0){
      val_name <- NA
    }

    # Check path.
    if(nchar(val_path) == 0 || val_path == expDefText){
      val_path <- NA
    }

    if(debug){
      print("val_object")
      print(val_object)
      print("val_use_obj")
      print(val_use_obj)
      print("val_name")
      print(val_name)
      print("val_ext")
      print(val_ext)
      print("val_del")
      print(val_del)
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
    ok <- val_use_obj || !is.na(val_name)
    ok <- ok && !is.na(val_path)
    
    # Check if any objects have been selected.
    if(length(val_object) == 0){
      ok <- FALSE
    }

    if(ok){
      
      svalue(g_export_btn) <- "Processing..."
      
      repeat{
      
        fail <- export(object=val_object, name=val_name, use.object.name=val_use_obj,
                        env=env, path=val_path, 
                        ext=val_ext, delim=val_delimiter, 
                        width=val_w, height=val_h, res=val_r,
                        overwrite=val_replace, debug=debug)
        
        if(is.data.frame(fail)){
          
          if(debug){
            print("The following objects failed:")
            print(fail)
          }
          
          dialog <- gbasicdialog(title="Export failed!", parent=w, do.buttons=FALSE, 
                                 width=200, height=200, horizontal=FALSE)
          
          msgtxt <- paste("\nThe following objects were not saved because the file names existed.\n",
                          "Make sure to exit the last edited cell before continuing.")
          
          msg <- glabel(text=msgtxt,
                        anchor=c(-1 ,0),
                        container=dialog)
          
          tbl <- gdf(items = fail, row.names=FALSE, container = dialog)
          size(tbl) <- c(100,200)
          
          gg <- ggroup(container=dialog) 
          btn_cancel <- gbutton("Cancel", container = gg, handler = function(h, ...) {
            fail <<- NULL
            dispose(dialog)
          })
          
          btn_replace <- gbutton("Overwrite", container = gg, handler = function(h, ...) {
            val_replace <<- TRUE
            dispose(dialog)
          })
          
          btn_retry <- gbutton("Retry", container = gg, handler = function(h, ...) {
            val_object <<- tbl[ , "Object"]
            val_name <<- tbl[ , "New.Name"]
            val_use_obj <<- FALSE

            if(debug){
              print("val_object")
              print(val_object)
              print("val_name")
              print(val_name)
            }
            
            dispose(dialog)
          })
          
          visible(dialog, set=TRUE)
          
        } else {
          
           break
        }
        
        if(is.null(fail)){
          break
        }
        
      }

      svalue(g_export_btn) <- "Export"
      
    } else {
      
      gmessage(message=paste("At least one object must be selected.",
                             "\nFile name and path must be provided."),
               title="Error",
               parent=w,
               icon = "error")      
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
      if(exists(".strvalidator_export_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_export_gui_savegui", envir=env)
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
      if(exists(".strvalidator_export_gui_objName", envir=env, inherits = FALSE)){
        svalue(f1g1_name_chk) <- get(".strvalidator_export_gui_objName", envir=env)
      }
      if(exists(".strvalidator_export_gui_replace", envir=env, inherits = FALSE)){
        svalue(f1g1_replace_chk) <- get(".strvalidator_export_gui_replace", envir=env)
      }
      if(exists(".strvalidator_export_gui_ext", envir=env, inherits = FALSE)){
        svalue(f2g1_ext_drp) <- get(".strvalidator_export_gui_ext", envir=env)
      }
      if(exists(".strvalidator_export_gui_del", envir=env, inherits = FALSE)){
        svalue(f2g1_del_drp) <- get(".strvalidator_export_gui_del", envir=env)
      }
      if(exists(".strvalidator_export_gui_width", envir=env, inherits = FALSE)){
        svalue(f3g1_width_edt) <- get(".strvalidator_export_gui_width", envir=env)
      }
      if(exists(".strvalidator_export_gui_height", envir=env, inherits = FALSE)){
        svalue(f3g1_height_edt) <- get(".strvalidator_export_gui_height", envir=env)
      }
      if(exists(".strvalidator_export_gui_res", envir=env, inherits = FALSE)){
        svalue(f3g1_res_edt) <- get(".strvalidator_export_gui_res", envir=env)
      }
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_export_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_export_gui_objName", value=svalue(f1g1_name_chk), envir=env)
      assign(x=".strvalidator_export_gui_replace", value=svalue(f1g1_replace_chk), envir=env)
      assign(x=".strvalidator_export_gui_ext", value=svalue(f2g1_ext_drp), envir=env)
      assign(x=".strvalidator_export_gui_del", value=svalue(f2g1_del_drp), envir=env)
      assign(x=".strvalidator_export_gui_width", value=svalue(f3g1_width_edt), envir=env)
      assign(x=".strvalidator_export_gui_height", value=svalue(f3g1_height_edt), envir=env)
      assign(x=".strvalidator_export_gui_res", value=svalue(f3g1_res_edt), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_export_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_export_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_export_gui_objName", envir=env, inherits = FALSE)){
        remove(".strvalidator_export_gui_objName", envir = env)
      }
      if(exists(".strvalidator_export_gui_replace", envir=env, inherits = FALSE)){
        remove(".strvalidator_export_gui_replace", envir = env)
      }
      if(exists(".strvalidator_export_gui_ext", envir=env, inherits = FALSE)){
        remove(".strvalidator_export_gui_ext", envir = env)
      }
      if(exists(".strvalidator_export_gui_del", envir=env, inherits = FALSE)){
        remove(".strvalidator_export_gui_del", envir = env)
      }
      if(exists(".strvalidator_export_gui_width", envir=env, inherits = FALSE)){
        remove(".strvalidator_export_gui_width", envir = env)
      }
      if(exists(".strvalidator_export_gui_height", envir=env, inherits = FALSE)){
        remove(".strvalidator_export_gui_height", envir = env)
      }
      if(exists(".strvalidator_export_gui_res", envir=env, inherits = FALSE)){
        remove(".strvalidator_export_gui_res", envir = env)
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
