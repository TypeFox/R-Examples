################################################################################
# TODO LIST
# TODO: Check folder DOES NOT WORK BECAUSE \ IS ESCAPE CHARACTER.
# TODO: Update when 'import' is changed to use 'fread'.

################################################################################
# CHANGE LOG (last 20 changes)
# 15.12.2015: Removed "0" from the default 'na.strings'.
# 04.12.2015: Implemented new parameter 'na.strings'.
# 05.10.2015: Added attributes.
# 29.08.2015: Added importFrom.
# 23.05.2015: Added new options available in 'import'.
# 11.10.2014: Added 'focus', added 'parent' parameter.
# 28.06.2014: Added help button and moved save gui checkbox.
# 20.01.2014: Remove redundant "overwrite?" message dialog.
# 13.01.2014: Handle empty dataframe by stay in gui and show message.
# 10.12.2013: Updated with new parameter names in function 'import'.
# 12.11.2013: Pass debug to function.
# 18.07.2013: Check before overwrite object.
# 15.07.2013: Added save GUI settings.
# 11.06.2013: Fixed 'exists' added 'inherits=FALSE'. Added parameter 'debug'.
# 16.04.2013: Added object name check.

#' @title Import Data
#'
#' @description
#' GUI wrapper for the \code{\link{import}} function.
#'
#' @details
#' Simplifies the use of the \code{\link{import}} function by providing a graphical 
#' user interface to it.
#' 
#' @param env environment into which the object will be saved.
#' Default is the current environment.
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
#' @seealso \code{\link{import}}


import_gui <- function(env=parent.frame(), savegui=NULL, debug=FALSE, parent=NULL){
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # Define variables.
  defaultDir <- "Select a directory..."
  defaultFile <- "Select a file..."
  
  
  # Add new parameter , settings=FALSE  
  #  # Load settings.
  #   if(settings){
  #     if(exists(".strvalidator_import_gui_file")){
  #       defaultFile <- .strvalidator_import_gui_file
  #     }
  #     if(exists(".strvalidator_import_gui_dir")){
  #       defaultDir <- .strvalidator_import_gui_dir
  #     }
  #     
  #   }
  
  # Main window.  
  w <- gwindow(title="Import from files", 
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
    print(help("import_gui", help_type="html"))
    
  })
  
  # GUI #######################################################################
  
  options <- c("Import multiple files from a directory into one dataset", 
               "Import a single file")
  
  import_opt <- gradio(items=options, selected=2,
                       horizontal=FALSE, container=gv)
  
  addHandlerChanged(import_opt, handler = function (h, ...) {
    
    .refresh()
    
  })
  
  import_file <- gfilebrowse(text=defaultFile,
                             initial.filename = defaultFile, # Not implemented in current version?
                             type="open",
                             quote = FALSE,
                             container=gv)
  
  
  import_folder <- gfilebrowse(text=defaultDir, 
                               initial.dir = defaultDir, # Not implemented in current version?
                               type="selectdir",
                               quote = FALSE,
                               container=gv)
  
  enabled(import_folder) <- FALSE
  
  # OPTIONS -------------------------------------------------------------------
  
  opt_frm <- gframe(text="Options", pos=0, horizontal=FALSE, container=gv)
  
  opt_file_chk <- gcheckbox(text="Save file name", checked = TRUE,
                            container=opt_frm)
  
  opt_time_chk <- gcheckbox(text="Save file time stamp", checked = TRUE,
                            container=opt_frm)
  
  glabel(text="Delimiter:", container=opt_frm, anchor=c(-1 ,0))
  opt_sep_drp <- gdroplist(items=c("TAB","SPACE","COMMA","SEMICOLON"),
                           selected=1, editable=FALSE, container=opt_frm)
  
  glabel(text="NA strings (separated by comma):",
         container=opt_frm, anchor=c(-1 ,0))
  opt_na_edt <- gedit(text="NA,,", container=opt_frm)

  opt_trim_chk <- gcheckbox(text="Auto trim samples", checked = FALSE,
                            container=opt_frm)
  
  opt_slim_chk <- gcheckbox(text="Auto slim repeated columns",
                            checked = FALSE, container=opt_frm)
  
  addHandlerChanged(opt_trim_chk, handler = function (h, ...) {
    
    .refresh()
    
  })
  
  addHandlerChanged(opt_slim_chk, handler = function (h, ...) {
    
    .refresh()
    
  })
  
  
  # MULTIPLE FILES OPTIONS ----------------------------------------------------
  
  multi_frm <- gexpandgroup(text="Multiple files options",
                            horizontal=FALSE, container=gv)
  
  enabled(multi_frm) <- FALSE
  
  multi_case_chk <- gcheckbox(text="Ignore case", checked = TRUE,
                              container=multi_frm)
  
  glabel(text="Prefix:", container=multi_frm, anchor=c(-1 ,0))
  multi_pre_edt <- gedit(initial.msg="", width = 25,
                         container=multi_frm, expand=TRUE)
  
  glabel(text="Suffix:", container=multi_frm, anchor=c(-1 ,0))
  multi_suf_edt <- gedit(initial.msg="", width = 25,
                         container=multi_frm, expand=TRUE)
  
  glabel(text="Extension:", container=multi_frm, anchor=c(-1 ,0))
  multi_ext_edt <- gedit(text="txt", width = 25,
                         container=multi_frm, expand=TRUE)
  
  # TRIM ----------------------------------------------------------------------
  
  trim_frm <- gexpandgroup(text="Trim options",
                           horizontal=FALSE, container=gv)
  
  glabel(text="Trim samples containing the word (separate by pipe |):",
         container=trim_frm, anchor=c(-1 ,0))
  
  trim_samples_edt <- gedit(text="pos|neg|ladder",
                            width = 25, container=trim_frm, expand=TRUE)
  
  trim_invert_chk <- gcheckbox(text="Invert (remove matching samples)", checked=TRUE,
                               container=trim_frm)
  
  # SLIM ----------------------------------------------------------------------
  
  slim_frm <- gexpandgroup(text="Slim options",
                           horizontal=FALSE, container=gv)
  
  slim_fix_chk <- gcheckbox(text="Keep all fixed (keep a row even if no data)",
                            checked=TRUE, container=slim_frm)
  
  # SAVE --------------------------------------------------------------------
  
  save_frm <- gframe(text="Save options", pos=0,
                     horizontal=FALSE, container=gv)
  
  glabel(text="Name:", container=save_frm, anchor=c(-1 ,0))
  import_edt <- gedit(initial.msg="Name for new dataset",
                      width = 25, container=save_frm, expand=TRUE)
  
  # IMPORT --------------------------------------------------------------------
  
  import_btn <- gbutton(text="Import", border=TRUE, container=gv)
  
  
  addHandlerChanged(import_btn, handler = function(h, ...) {
    
    # Get values.
    file_val <- svalue(import_file)
    folder_val <- svalue(import_folder)
    ignore_val <- svalue(multi_case_chk)
    prefix_val <- svalue(multi_pre_edt)
    suffix_val <- svalue(multi_suf_edt)
    extension_val <- svalue(multi_ext_edt)
    folder_opt_val <- if(svalue(import_opt, index=TRUE)==1){TRUE}else{FALSE}
    val_name <- svalue(import_edt)
    get_file_val <- svalue(opt_file_chk)
    get_time_val <- svalue(opt_time_chk)
    del_val <- svalue(opt_sep_drp, index=TRUE)
    na_val <- svalue(opt_na_edt)
    trim_val <- svalue(opt_trim_chk)
    trim_what_val <- svalue(trim_samples_edt)
    trim_invert_val <- svalue(trim_invert_chk)
    slim_val <- svalue(opt_slim_chk)
    slim_fix_val <- svalue(slim_fix_chk)
    
    # Assign a delimiter character.
    if(del_val == 1){
      val_delimiter <- "\t"   
    } else if(del_val == 2){
      val_delimiter <- " "
    } else if(del_val == 3){
      val_delimiter <- ","
    } else if(del_val == 4){
      val_delimiter <- ";"
    } 
    
    # Convert to character vector.
    val_na <- unlist(strsplit(na_val,","))

    # Initiate variable.  
    ok <- TRUE

    # Check that a name has been provided for the new data object.
    if(nchar(val_name) == 0){
      
      gmessage("A name for the dataset must be provided.",
               title="Error", icon="error", parent=w)
      
      ok <- FALSE
      
    }
    
    #     # TODO: DOES NOT WORK BECAUSE \ IS ESCAPE CHARACTER.
    #     # Check that folder exist.
    #     if(!file.exists(folder_val)){
    # 
    #       ok <- FALSE
    #       
    #       gmessage("The provided folder does not exist or is not accessible.",
    #                title="Error", icon="error", parent=w)
    #       
    #     }
    
    
    # Check if ok to import data to 'env'.
    if(ok){
      
      # Set arguments.
      if(folder_opt_val){
        file_val <- NA
      } else {
        folder_val <- NA
      }
      
      if(!nchar(prefix_val)>0){
        prefix_val <- NA
      }
      
      if(!nchar(suffix_val)>0){
        suffix_val <- NA
      }
      
      if(debug){
        print("val_name")
        print(val_name)
        print("ignore_val")
        print(ignore_val)
        print("prefix_val")
        print(prefix_val)
        print("suffix_val")
        print(suffix_val)
        print("folder_opt_val")
        print(folder_opt_val)
        print("file_val")
        print(file_val)
        print("folder_val")
        print(folder_val)
        print("get_file_val")
        print(get_file_val)
        print("get_time_val")
        print(get_time_val)
        print("del_val")
        print(del_val)
        print("na_val")
        print(na_val)
        print("val_na")
        print(val_na)
        print("trim_val")
        print(trim_val)
        print("trim_what_val")
        print(trim_what_val)
        print("trim_invert_val")
        print(trim_invert_val)
        print("slim_val")
        print(slim_val)
        print("slim_fix_val")
        print(slim_fix_val)
      }
      
      # Change button.
      svalue(import_btn) <- "Processing..."
      enabled(import_btn) <- FALSE
      
      # Call function.
      datanew <- import(folder=folder_opt_val,
                        extension=extension_val,
                        suffix=suffix_val,
                        prefix=prefix_val,
                        import.file=file_val,
                        folder.name=folder_val,
                        file.name=get_file_val,
                        time.stamp=get_time_val,
                        separator=val_delimiter,
                        ignore.case=ignore_val,
                        auto.trim=trim_val,
                        trim.samples=trim_what_val,
                        trim.invert=trim_invert_val,
                        auto.slim=slim_val,
                        slim.na=slim_fix_val,
                        na.strings=val_na,
                        debug=debug)
      
      if(length(datanew) == 0){
        
        # Show warning.
        gmessage(message="Dataset empty!\nCheck your file filter.",
                 title="Error",
                 icon = "error",
                 parent = w)
        
        # Change button.
        svalue(import_btn) <- "Import"
        enabled(import_btn) <- TRUE
        
      } else {
        
        # Add attributes.
        attr(datanew, which="import_gui, strvalidator") <- as.character(utils::packageVersion("strvalidator"))
        attr(datanew, which="import_gui, folder") <- folder_opt_val
        attr(datanew, which="import_gui, extension") <- extension_val
        attr(datanew, which="import_gui, suffix") <- suffix_val
        attr(datanew, which="import_gui, prefix") <- prefix_val
        attr(datanew, which="import_gui, import.file") <- file_val
        attr(datanew, which="import_gui, folder.name") <- folder_val
        attr(datanew, which="import_gui, file.name") <- get_file_val
        attr(datanew, which="import_gui, time.stamp") <- get_time_val
        attr(datanew, which="import_gui, ignore.case") <- ignore_val
        attr(datanew, which="import_gui, auto.trim") <- trim_val
        attr(datanew, which="import_gui, trim.samples") <- trim_what_val
        attr(datanew, which="import_gui, trim.invert") <- trim_invert_val
        attr(datanew, which="import_gui, auto.slim") <- slim_val
        attr(datanew, which="import_gui, slim.na") <- slim_fix_val
        attr(datanew, which="import_gui, separator") <- val_delimiter
        attr(datanew, which="import_gui, na.strings") <- val_na
        
        # Save data.
        saveObject(name=val_name, object=datanew, parent=w, env=env)
        
        # Close GUI.
        dispose(w)
        
      }
      
    }
    
  } )
  
  # INTERNAL FUNCTIONS ########################################################
  
  .refresh <- function(){
    
    # Get values.
    val_trim <- svalue(opt_trim_chk)
    val_slim <- svalue(opt_slim_chk)
    folder_opt_val <- if(svalue(import_opt, index=TRUE)==1){TRUE}else{FALSE}
    
    
    if(val_trim){
      enabled(trim_samples_edt) <- TRUE
      enabled(trim_invert_chk) <- TRUE
    } else {
      enabled(trim_samples_edt) <- FALSE
      enabled(trim_invert_chk) <- FALSE
    }
    
    if(val_slim){
      enabled(slim_fix_chk) <- TRUE
    } else {
      enabled(slim_fix_chk) <- FALSE
    }
    
    if(folder_opt_val){
      enabled(multi_frm) <- TRUE
      enabled(import_folder) <- TRUE
      enabled(import_file) <- FALSE
    } else {
      enabled(multi_frm) <- FALSE
      enabled(import_folder) <- FALSE
      enabled(import_file) <- TRUE
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
      if(exists(".strvalidator_import_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_import_gui_savegui", envir=env)
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
      if(exists(".strvalidator_import_gui_import_opt", envir=env, inherits = FALSE)){
        svalue(import_opt) <- get(".strvalidator_import_gui_import_opt", envir=env)
      }
      if(exists(".strvalidator_import_gui_file", envir=env, inherits = FALSE)){
        svalue(opt_file_chk) <- get(".strvalidator_import_gui_file", envir=env)
      }
      if(exists(".strvalidator_import_gui_time", envir=env, inherits = FALSE)){
        svalue(opt_time_chk) <- get(".strvalidator_import_gui_time", envir=env)
      }
      if(exists(".strvalidator_import_gui_sep", envir=env, inherits = FALSE)){
        svalue(opt_sep_drp) <- get(".strvalidator_import_gui_sep", envir=env)
      }
      if(exists(".strvalidator_import_gui_na", envir=env, inherits = FALSE)){
        svalue(opt_na_edt) <- get(".strvalidator_import_gui_na", envir=env)
      }
      if(exists(".strvalidator_import_gui_ignore", envir=env, inherits = FALSE)){
        svalue(multi_case_chk) <- get(".strvalidator_import_gui_ignore", envir=env)
      }
      if(exists(".strvalidator_import_gui_prefix", envir=env, inherits = FALSE)){
        svalue(multi_pre_edt) <- get(".strvalidator_import_gui_prefix", envir=env)
      }
      if(exists(".strvalidator_import_gui_suffix", envir=env, inherits = FALSE)){
        svalue(multi_suf_edt) <- get(".strvalidator_import_gui_suffix", envir=env)
      }
      if(exists(".strvalidator_import_gui_extension", envir=env, inherits = FALSE)){
        svalue(multi_ext_edt) <- get(".strvalidator_import_gui_extension", envir=env)
      }
      if(exists(".strvalidator_import_gui_trim", envir=env, inherits = FALSE)){
        svalue(opt_trim_chk) <- get(".strvalidator_import_gui_trim", envir=env)
      }
      if(exists(".strvalidator_import_gui_trim_samples", envir=env, inherits = FALSE)){
        svalue(trim_samples_edt) <- get(".strvalidator_import_gui_trim_samples", envir=env)
      }
      if(exists(".strvalidator_import_gui_trim_invert", envir=env, inherits = FALSE)){
        svalue(trim_invert_chk) <- get(".strvalidator_import_gui_trim_invert", envir=env)
      }
      if(exists(".strvalidator_import_gui_slim", envir=env, inherits = FALSE)){
        svalue(opt_slim_chk) <- get(".strvalidator_import_gui_slim", envir=env)
      }
      if(exists(".strvalidator_import_gui_slim_fix", envir=env, inherits = FALSE)){
        svalue(slim_fix_chk) <- get(".strvalidator_import_gui_slim_fix", envir=env)
      }
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_import_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_import_gui_import_opt", value=svalue(import_opt), envir=env)
      assign(x=".strvalidator_import_gui_file", value=svalue(opt_file_chk), envir=env)
      assign(x=".strvalidator_import_gui_time", value=svalue(opt_time_chk), envir=env)
      assign(x=".strvalidator_import_gui_sep", value=svalue(opt_sep_drp), envir=env)
      assign(x=".strvalidator_import_gui_na", value=svalue(opt_na_edt), envir=env)
      assign(x=".strvalidator_import_gui_ignore", value=svalue(multi_case_chk), envir=env)
      assign(x=".strvalidator_import_gui_prefix", value=svalue(multi_pre_edt), envir=env)
      assign(x=".strvalidator_import_gui_suffix", value=svalue(multi_suf_edt), envir=env)
      assign(x=".strvalidator_import_gui_extension", value=svalue(multi_ext_edt), envir=env)
      assign(x=".strvalidator_import_gui_trim", value=svalue(opt_trim_chk), envir=env)
      assign(x=".strvalidator_import_gui_trim_samples", value=svalue(trim_samples_edt), envir=env)
      assign(x=".strvalidator_import_gui_trim_invert", value=svalue(trim_invert_chk), envir=env)
      assign(x=".strvalidator_import_gui_slim", value=svalue(opt_slim_chk), envir=env)
      assign(x=".strvalidator_import_gui_slim_fix", value=svalue(slim_fix_chk), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_import_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_import_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_import_gui_import_opt", envir=env, inherits = FALSE)){
        remove(".strvalidator_import_gui_import_opt", envir = env)
      }
      if(exists(".strvalidator_import_gui_file", envir=env, inherits = FALSE)){
        remove(".strvalidator_import_gui_file", envir = env)
      }
      if(exists(".strvalidator_import_gui_time", envir=env, inherits = FALSE)){
        remove(".strvalidator_import_gui_time", envir = env)
      }
      if(exists(".strvalidator_import_gui_sep", envir=env, inherits = FALSE)){
        remove(".strvalidator_import_gui_sep", envir = env)
      }
      if(exists(".strvalidator_import_gui_na", envir=env, inherits = FALSE)){
        remove(".strvalidator_import_gui_na", envir = env)
      }
      if(exists(".strvalidator_import_gui_ignore", envir=env, inherits = FALSE)){
        remove(".strvalidator_import_gui_ignore", envir = env)
      }
      if(exists(".strvalidator_import_gui_prefix", envir=env, inherits = FALSE)){
        remove(".strvalidator_import_gui_prefix", envir = env)
      }
      if(exists(".strvalidator_import_gui_suffix", envir=env, inherits = FALSE)){
        remove(".strvalidator_import_gui_suffix", envir = env)
      }
      if(exists(".strvalidator_import_gui_extension", envir=env, inherits = FALSE)){
        remove(".strvalidator_import_gui_extension", envir = env)
      }
      if(exists(".strvalidator_import_gui_trim", envir=env, inherits = FALSE)){
        remove(".strvalidator_import_gui_trim", envir = env)
      }
      if(exists(".strvalidator_import_gui_trim_samples", envir=env, inherits = FALSE)){
        remove(".strvalidator_import_gui_trim_samples", envir = env)
      }
      if(exists(".strvalidator_import_gui_trim_invert", envir=env, inherits = FALSE)){
        remove(".strvalidator_import_gui_trim_invert", envir = env)
      }
      if(exists(".strvalidator_import_gui_slim", envir=env, inherits = FALSE)){
        remove(".strvalidator_import_gui_slim", envir = env)
      }
      if(exists(".strvalidator_import_gui_slim_fix", envir=env, inherits = FALSE)){
        remove(".strvalidator_import_gui_slim_fix", envir = env)
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
  .refresh()
  
  # Show GUI.
  visible(w) <- TRUE
  focus(w)
  
}