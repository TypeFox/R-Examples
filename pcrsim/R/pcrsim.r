################################################################################
# TODO LIST
# TODO: Create parameter file from strvalidator kit definition file (Short name must be identical).
# TODO: Rewrite compact/strvalidator::compact to use data.table?

# IMPORTANT: To manually run R CMD check in RStudio all packages must be installed in
# both the 32 and 64 bit version. Make sure it is possible to start manually
# (GTK+ must be installed by clicking 'OK' on the message box).

# See http://r-pkgs.had.co.nz/release.html for advice on release.
# IMPORTANT: Use build_win() to test on current R and R-dev 'library(devtools)'.
# IMPORTANT: Use devtools::release() to submitt to CRAN.

# Versioning convention (x.yy.z[.9###]):
# Increment x on major change.
# Increment yy on new features.
# Increment z on minor changes and bug fixes.
# [optional]Increment ### on development versions.
# NB! Write changes in NEWS for x.yy.z.9### versions, but move changes to NEWS under x.yy.z upon release official version.

# NOTE:
# \u00B5 is the unicode for my 

################################################################################
# CHANGE LOG (10 last changes)
# 14.04.2016: Version 1.0.0 released.

#' @title GUI for PCRsim
#'
#' @description
#' A graphical user interface for simulation of the entire DNA process.
#'
#' @details
#' This graphical user interface give access to parameters and functions
#' for simulation of the forensic DNA process.
#' Detailes are entered for each sub process organised into tabs.
#' Simulation is performed and the result can be viewed within the GUI,
#' plotted as an electropherogram (EPG), or saved to a text file.
#' The EPG can be saved as an image.
#' 
#' @param debug logical, indicating if debug information should be printed.
#' 
#' @export
#' 
#' @import gWidgets
#' @importFrom strvalidator listObjects editData_gui getDb trim slim generateEPG import_gui getKit
#' @importFrom utils head tail str packageVersion
#' 
#' @examples
#' \dontrun{
#' # Open the graphical user interface.
#' pcrsim()
#' }

pcrsim <- function(debug=FALSE){
  
  # Global variables.
  .pcrsim_env <- new.env()
  .separator <- .Platform$file.sep # Platform dependent path separator.
  .start_tab_name <- "Welcome"
  .project_tab_name <- "Projects"
  .file_tab_name <- "Workspace"
  .profile_tab_name <- "Profile"
  .sample_tab_name <- "Sample"
  .degradation_tab_name <- "Degradation"
  .extraction_tab_name <- "Extraction"
  .dilution_tab_name <- "Normalization"
  .amplification_tab_name <- "PCR"
  .analysis_tab_name <- "CE"
  .epg_tab_name <- "EPG"
  .simulation_tab_name <- "Simulation"
  .mixture_tab_name <- "Mixtures"
  
  .project_description_variable <- ".pcrsim_project_description"
  .project_tmp_env <- new.env()
  .project_name_list <- NULL
  .project_path_list <- NULL

  .simDataRes <- NULL # Current simulation result for the complete DNA process.
  .simDataPh <- NULL # Current simulation result (after CE) with peak heights.
  .simEPG <- NULL # ggplot2 object.
  
  .ws_name_variable <- ".pcrsim_project_name"
  .ws_path_variable <- ".pcrsim_project_path"
  
  # Main window.
  w <- gwindow(title=paste("PCRsim",packageVersion("pcrsim")), visible = FALSE)
  
  g <- ggroup(horizontal = TRUE, spacing = 5, use.scrollwindow = FALSE,
              container = w, expand = TRUE)
  
  nb <- gnotebook(closebuttons = FALSE, dontCloseThese = NULL,
                  container = g, expand = TRUE)
  
  
  # Define groups.
  start_gf <- ggroup(horizontal = FALSE, spacing=10, use.scrollwindow=FALSE,
                     container = nb, label=.start_tab_name, expand=TRUE)
  
  project_gf <- ggroup(horizontal = FALSE, spacing=10, use.scrollwindow=FALSE,
                        container = nb, label=.project_tab_name, expand=TRUE)

  file_gf <- ggroup(horizontal = FALSE, spacing=10, use.scrollwindow=FALSE,
                    container = nb, label=.file_tab_name, expand=TRUE)
  
  profile_gf <- ggroup(horizontal = FALSE, spacing=10, use.scrollwindow=FALSE,
                       container = nb, label=.profile_tab_name, expand=TRUE)
  
  sample_gf <- ggroup(horizontal = FALSE, spacing=5, use.scrollwindow=FALSE,
                      container = nb, label=.sample_tab_name, expand=TRUE)
  
  deg_gf <- ggroup(horizontal = FALSE, spacing=5, use.scrollwindow=FALSE,
                   container = nb, label=.degradation_tab_name, expand=TRUE)
  
  ex_gf <- ggroup(horizontal = FALSE, spacing=5, use.scrollwindow=FALSE,
                  container = nb, label=.extraction_tab_name, expand=TRUE)
  
  dil_gf <- ggroup(horizontal = FALSE, spacing=5, use.scrollwindow=FALSE,
                   container = nb, label=.dilution_tab_name, expand=TRUE)
  
  amp_gf <- ggroup(horizontal = FALSE, spacing=5, use.scrollwindow=FALSE,
                   container = nb, label=.amplification_tab_name, expand=TRUE)
  
  ce_gf <- ggroup(horizontal = FALSE, spacing=5, use.scrollwindow=FALSE,
                  container = nb, label=.analysis_tab_name, expand=TRUE)

  epg_gf <- ggroup(horizontal = FALSE, spacing=5, use.scrollwindow=FALSE,
                  container = nb, label=.epg_tab_name, expand=TRUE)
  
  sim_gf <- ggroup(horizontal = FALSE, spacing=5, use.scrollwindow=FALSE,
                   container = nb, label=.simulation_tab_name, expand=TRUE)

  # START #####################################################################
  
  glabel("", container=start_gf) # Adds some space.
  
  # STR TYPING KIT ------------------------------------------------------------
  
  start_f1 <- gframe(text = "PCR sim", markup = FALSE, pos = 0, horizontal = TRUE,
                     container = start_gf, expand = TRUE) 
  
  about_txt <- paste("PCR sim is a package for simulation of the forensic ",
                     "DNA process. This graphical user interface give access ",
                     "to function parameters. Parameters are organised ",
                     "into tabs for the respective subprocess. ",
                     "Simulation is performed and the result for each step can ",
                     "be viewed as a table within the GUI or plotted as an ",
                     "electropherogram (EPG).\n\n",
                     "Effort has been made to mimic each step of the real process as closely ",
                     "as possible - consequently, due to performance, it is not ",
                     "well suitable for very large simulations.\n\n",
                     "The simulator must be calibrated to the quantification ",
                     "method used and for each capillary electrophoresis ",
                     "instrument for realistic simulations.\n\n",
                     "General information and tutorials:\n",
                     "https://sites.google.com/site/forensicapps/pcrsim\n\n",
                     "Please report bugs to:\n",
                     "https://github.com/OskarHansson/pcrsim/issues\n\n",
                     "The source is hosted at GitHub:\n",
                     "https://github.com/OskarHansson/pcrsim", sep="")
  
  gtext(text=about_txt, width = NULL, height = 300, font.attr = NULL, 
        wrap = TRUE, expand=TRUE, container = start_f1) 
  
  start_f2 <- gframe(text = "License", markup = FALSE, pos = 0, horizontal=TRUE,
                     container = start_gf, expand=TRUE) 
  
  license_txt <- paste("Copyright (C) 2013 Oskar Hansson\n\n",
                       "This program is free software; you can redistribute it and/or ",
                       "modify it under the terms of the GNU General Public License ",
                       "as published by the Free Software Foundation; either version 2 ",
                       "of the License, or (at your option) any later version.\n\n",
                       "This program is distributed in the hope that it will be useful, ",
                       "but WITHOUT ANY WARRANTY; without even the implied warranty of ",
                       "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the ",
                       "GNU General Public License for more details.\n\n",
                       "You should have received a copy of the GNU General Public License ",
                       "along with this program; if not, write to the Free Software ",
                       "Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, ",
                       "MA  02110-1301, USA.", sep="")
  
  gtext(text=license_txt, width = NULL, height = 300, font.attr = NULL, 
        wrap = TRUE, expand=TRUE, container = start_f2) 
  
  # PROJECT MANAGER ###########################################################
  
  # Vertical main group.
  project_f1 <- ggroup(horizontal=FALSE, container = project_gf,
                       spacing=10, expand=TRUE)
  
  # FOLDER --------------------------------------------------------------------
  
  glabel(text="Folder:", anchor=c(-1 ,0), container=project_f1)
  
  project_fb <- gfilebrowse(text=getwd(), type="selectdir", quote = FALSE,
                            container=project_f1)
  
  
  addHandlerChanged(project_fb, handler = function (h, ...) {
    
    .updateProjectList()
    
  } )  
  
  addHandlerBlur(project_fb, handler = function (h, ...) {
    
    .updateProjectList()
    
  } )  
  
  # PROJECTS ------------------------------------------------------------------
  
  # Horizontal main group.
  project_f2 <- gframe(text="Projects", horizontal=TRUE, spacing=10,
                       container = project_f1, expand=TRUE)
  
  # Button group.
  project_g1 <- ggroup(horizontal=FALSE, spacing=10, container = project_f2,
                       expand=FALSE)
  
  glabel("", container=project_g1) # Adds some space.
  
  project_open_btn <- gbutton(text="Open", border=TRUE, container = project_g1)
  tooltip(project_open_btn) <- "Open selected project"
  
  project_add_btn <- gbutton(text="Add", border=TRUE, container = project_g1)
  tooltip(project_add_btn) <- "Merge with current project"
  
  project_delete_btn <- gbutton(text="Delete", border=TRUE,
                                container = project_g1)
  tooltip(project_delete_btn) <- "Delete selected project from the file system"
  
  addSpring(project_g1)

  # HANDLERS ..................................................................
  
  addHandlerChanged(project_open_btn, handler = function (h, ...) {
    
    # Get selected projects file name.
    val_name <- svalue(project_tbl)
    val_id <- as.numeric(project_tbl[svalue(project_tbl, index=TRUE), "Id"])
    val_prj <- .project_path_list[val_id]
    val_env <- .pcrsim_env
    
    if(debug){
      print(paste("Selected path", val_prj))
      print(paste("Selected project", val_name))
      print(paste("Selected index", val_id))
    }
    
    # Check if file exist.
    if(length(val_prj) > 0){
      
      if(file.exists(val_prj)){
        
        # Clear environment.
        remove(list=ls(envir=val_env, all.names=TRUE), 
               envir=val_env, inherits=FALSE)
        
        # Load project to workspace.
        load(file=val_prj, envir=val_env, verbose=FALSE)
        
        # Select to workspace tab.
        svalue(nb) <- match(.file_tab_name, names(nb))
        
      }
      
    }
    
  } )
  
  addHandlerChanged(project_add_btn, handler = function (h, ...) {
    
    # Get selected projects file name.
    val_name <- svalue(project_tbl)
    val_id <- as.numeric(project_tbl[svalue(project_tbl, index=TRUE), "Id"])
    val_prj <- .project_path_list[val_id]
    val_env <- .pcrsim_env
    
    if(debug){
      print(paste("Selected path", val_prj))
      print(paste("Selected project", val_name))
      print(paste("Selected index", val_id))
    }
    
    # Check if file exist.
    if(length(val_prj) > 0){
      
      if(file.exists(val_prj)){
        
        # Load project to workspace.
        load(file=val_prj, envir=val_env, verbose=FALSE)
        message(paste("Loaded", val_prj))
        
      }
      
    }
    
  } )
  
  addHandlerChanged(project_delete_btn, handler = function (h, ...) {
    
    # Get selected projects file name.
    val_name <- svalue(project_tbl)
    val_id <- as.numeric(project_tbl[svalue(project_tbl, index=TRUE), "Id"])
    val_prj <- .project_path_list[val_id]
    
    if(debug){
      print(paste("Selected path", val_prj))
      print(paste("Selected project", val_name))
      print(paste("Selected index", val_id))
    }
    
    # Check if file exist.
    if(length(val_prj) > 0){
      
      if(file.exists(val_prj)){
        
        # Delete project file and update list.
        file.remove(val_prj)
        message(paste("Deleted", val_prj))
        .updateProjectList()
        
        # Clear description box.
        svalue(proj_info_lbl) <- paste("Project:")
        svalue(proj_info_txt) <- ""
        
      }
      
    }
    
  } )
  
  # Projects group.
  project_g2 <- ggroup(horizontal=FALSE,
                       use.scrollwindow=FALSE,
                       container = project_f2,
                       expand=TRUE) 
  
  # Projects list.
  project_tbl <- gWidgets::gtable(items=data.frame(Name="", Date="",
                                                   Size="", Id="",
                                                   stringsAsFactors=FALSE), 
                                  multiple = TRUE,
                                  chosencol = 1,
                                  expand = TRUE,
                                  container = project_g2) 
  
  addHandlerClicked(project_tbl, handler = function (h, ...) {
    
    # Get selected projects file name.
    val_name <- svalue(project_tbl)
    val_id <- as.numeric(project_tbl[svalue(project_tbl, index=TRUE), "Id"])
    val_prj <- .project_path_list[val_id]
    val_obj <- .project_description_variable
    val_env <- .project_tmp_env
    
    if(debug){
      print(paste("In addHandlerClicked(project_tbl"))
      print(paste("Selected path", val_prj))
      print(paste("Selected project", val_name))
      print(paste("Selected index", val_id))
    }
    
    # Enable possibly disabled save button upon changed selectioin.
    enabled(project_save_btn) <- TRUE
    
    # Clear environment.
    remove(list=ls(envir=val_env, all.names=TRUE), envir=val_env, inherits=FALSE)
    
    # Check if file exist.
    if(length(val_prj) > 0){
      
      if(file.exists(val_prj)){
        
        # Load project in temporary environment.
        load(file=val_prj, envir=val_env, verbose=FALSE)
        if(exists(x=val_obj, envir=val_env, inherits=FALSE)){
          description <- get(x=val_obj, envir=val_env, inherits=FALSE)
        } else {
          description <- "Write a project description here!"
        }
        
        # Load description.
        svalue(proj_info_lbl) <- paste("Project:", val_name)
        svalue(proj_info_txt) <- description
        
      }
      
    }
    
  } )
  
  # DESCRIPTION ---------------------------------------------------------------
  
  # Horizontal main group.
  project_f3 <- gframe(text="Description", horizontal=TRUE, spacing=10,
                       container = project_f1, expand=TRUE)
  
  # Button group.
  project_g3 <- ggroup(horizontal=FALSE, spacing=10, container = project_f3,
                       expand=FALSE)
  
  project_save_btn <- gbutton(text="Save", border=TRUE, container = project_g3)
  tooltip(project_save_btn) <- "Save project description"
  
  glabel("", container=project_g1) # Adds some space.
  
  glabel("", container=project_g1) # Adds some space.
  
  addHandlerChanged(project_save_btn, handler = function (h, ...) {
    
    enabled(project_save_btn) <- FALSE
    
    # Get selected projects file name.
    val_name <- svalue(project_tbl)
    val_id <- project_tbl[svalue(project_tbl, index=TRUE),"Id"]
    val_id <- as.numeric(val_id)
    val_prj <- .project_path_list[val_id]
    val_obj <- .project_description_variable
    val_env <- .project_tmp_env
    val_description <- svalue(proj_info_txt)
    
    message("Assign: ", val_obj)
    message("Save: ", val_prj)
    
    # Save project description and write to disc.
    assign(x=val_obj, value=val_description, envir=val_env, inherits=FALSE)
    save(file=val_prj, list=ls(envir=val_env, all.names=TRUE), envir=val_env)
    
    enabled(project_save_btn) <- TRUE
    
  } )  
  
  # Button group.
  project_g4 <- ggroup(horizontal=FALSE, spacing=10, container=project_f3,
                       expand=TRUE)
  
  # Project description window.
  proj_info_lbl <- glabel(text="Project:", anchor=c(-1 ,0),
                          container=project_g4)
  proj_info_txt <- gtext(text=NULL, height=300, expand=TRUE,
                         wrap=TRUE, container=project_g4) 
  
    
  # FILE ######################################################################
  
  # LOADED DATASETS -----------------------------------------------------------
  
  if(debug){
    print("LOADED DATASETS")
  }
  
  # Variables and constants ...................................................
  
  # Get dataframes in workspace.
  file_list <- strvalidator::listObjects(env=.pcrsim_env,
                                          obj.class=c("data.frame","ggplot","environment"))
  
  
  # GUI .......................................................................
  
  file_f1 <- gframe(text = "Project workspace", markup = FALSE, pos = 0,
                    horizontal=TRUE, container = file_gf, expand=TRUE)
  
  file_f1g1 <- ggroup(horizontal=FALSE, container = file_f1, expand=FALSE)
  
  file_f1g1_load_btn <- gbutton(text = "Load", border = TRUE,
                                container = file_f1g1)
  tooltip(file_f1g1_load_btn) <- "Load project"
  
  file_f1g1_ws_save_btn <- gbutton(text = "Save", border = TRUE,
                                   container = file_f1g1)
  tooltip(file_f1g1_ws_save_btn) <- "Save project"

  file_f1g1_ws_saveas_btn <- gbutton(text = "Save As", border = TRUE,
                                     container = file_f1g1)
  tooltip(file_f1g1_ws_saveas_btn) <- "Save project in specified location"
  
  file_f1g1_import_btn <- gbutton(text = "Import", border = TRUE,
                                  container = file_f1g1)
  tooltip(file_f1g1_import_btn) <- "Import from file"
  
  file_f1g1_refresh_btn <- gbutton(text = "Refresh", border = TRUE,
                                   container = file_f1g1) 
  tooltip(file_f1g1_refresh_btn) <- "Refresh project workspace"
  
  file_f1g1_view_btn <- gbutton(text = "View", border = TRUE,
                                container = file_f1g1)
  tooltip(file_f1g1_view_btn) <- "View dataset or profile"
  
  file_f1g1_rm_btn <- gbutton(text = "Remove", border = TRUE,
                              container = file_f1g1) 
  tooltip(file_f1g1_rm_btn) <- "Remove selected object"
  
  file_f1g1_rename_btn <- gbutton(text = "Rename", border = TRUE,
                                  container = file_f1g1)
  tooltip(file_f1g1_rename_btn) <- "Rename selected object"
  
  file_f1g1_save_chk <- gcheckbox(text="Save GUI settings",
                                  checked = TRUE,
                                  container = file_f1g1)
  
  # To show debug messages in simulation.
  file_f1g1_debug_chk <- gcheckbox(text="Debug",
                                   checked = FALSE,
                                   container = file_f1g1)
  
  # To show debug messages in functions called from simulation.
  file_f1g1_ext_debug_chk <- gcheckbox(text="Extended debug",
                                       checked = FALSE,
                                       container = file_f1g1)

  glabel(text="\nSimulation template", container = file_f1g1)  
  file_f1g1_save_temp_btn <- gbutton(text = "Save", border = TRUE,
                                     container = file_f1g1)
  tooltip(file_f1g1_save_temp_btn) <- "Save simulation settings"
  
  file_f1g1_load_temp_btn <- gbutton(text = "Load", border = TRUE,
                                     container = file_f1g1)
  tooltip(file_f1g1_load_temp_btn) <- "Load simulation settings"
  
  file_f1g1_tbl <- gtable(items = file_list, multiple = TRUE,
                          expand = TRUE, container = file_f1) 
  
  # HANDLERS ..................................................................
  
  addHandlerChanged(file_f1g1_rename_btn, handler = function (h, ...) {
    
    
    object <- svalue(file_f1g1_tbl)
    
    if(length(object) == 1){
      
      newName_inp <- ginput(message="Enter new name",
                            title="Input",
                            icon = "info",
                            parent=w)
      
      newName_inp <- make.names(newName_inp)
      
      if(debug){
        print("newName_inp")
        print(newName_inp)
      }
      
      assign(newName_inp,
             get(object, envir = .pcrsim_env),
             envir = .pcrsim_env)
      
      remove(list=object, envir=.pcrsim_env)
      
      .refreshLoaded()
      
      
    } else {
      gmessage(message="You can only rename one object at a time!",
               title="Error",
               icon = "error",
               parent = w) 
    }
    
  } )
  
  addHandlerChanged(file_f1g1_load_btn, handler = function (h, ...) {
    
    
    ws_path <- gfile(text="Select a saved project", type="open",
                     filter = list("R files" = list(patterns = c("*.R","*.Rdata"))),
                     multi=FALSE)
    
    if(!is.na(ws_path)){
      
      if(file.exists(ws_path)){
        
        # Load project and refresh list.
        load(file=ws_path, envir = .pcrsim_env)
        .refreshLoaded()
        
        # Load saved settings.
        if(svalue(file_f1g1_save_chk)){
          
          .loadSavedSettings(.pcrsim_env)
        }
        
      } else {
        
        gmessage(message="The project file was not found",
                 title="File not found",
                 icon = "error",
                 parent = w) 
      }
    }    
    
  } )
  
  addHandlerChanged(file_f1g1_import_btn, handler = function (h, ...) {
    
    import_gui(env=.pcrsim_env)
    .refreshLoaded()
    
  } )  
  
  addHandlerChanged(file_f1g1_refresh_btn, handler = function (h, ...) {
    
    .refreshLoaded()
    
  })
  
  addHandlerChanged(file_f1g1_view_btn, handler = function(h, ...) {
    
    # Get selected dataset name(s).
    val_obj <- svalue(file_f1g1_tbl)
    val_data <- get(val_obj, envir=.pcrsim_env)
    val_class <- class(val_data)
    
    if (!is.null(val_obj) && !is.na(val_obj) && length(val_obj) > 0){
      
      if("data.frame" %in% val_class){
        
        # Open GUI.
        strvalidator::editData_gui(env=.pcrsim_env, data=val_data,
                                   name=val_obj, edit=FALSE, debug=debug)
        
      } else if("ggplot" %in% val_class){
        
        # Plot object.
        print(val_data)
        
      } else {
        message(paste("Object of type", val_class, "not supported!"))
      }
      
    }
    
  } )
  
  addHandlerChanged(file_f1g1_rm_btn, handler = function(h, ...) {
    
    # Get selected dataset name(s).
    val_obj <- svalue(file_f1g1_tbl)
    
    if(debug){
      print(paste("IN:", match.call()[[1]]))
      print("Changed, file_f1g1_rm_btn")
      print(val_obj)
    }
    
    if (!is.null(val_obj) && !is.na(val_obj)){
      
      # Get active reference data frame.
      remove(list=val_obj, envir=.pcrsim_env)
      
      .refreshLoaded()
      
      
    } 
  } )
  
  addHandlerChanged(file_f1g1_ws_save_btn, handler = function (h, ...) {

    # Initiate flag.
    ok <- TRUE
    
    # Get project name if available.
    if(exists(.ws_name_variable, envir=.pcrsim_env)){
      ws_name <- get(.ws_name_variable, envir=.pcrsim_env,
                     inherits=FALSE)
    } else {
      ok <- FALSE
    }
    
    # Get project path if available.
    if(exists(.ws_path_variable, envir=.pcrsim_env)){
      ws_save_path <- get(.ws_path_variable, envir=.pcrsim_env,
                          inherits=FALSE)
    } else {
      ok <- FALSE
    }
    
    if(ok){
      if(!is.na(ws_name) && !ws_name==""){
        
        ws_full_name <- paste(ws_save_path, .separator, ws_name, ".RData", sep="")
        
        if(debug){
          print(ws_full_name)
        }
        
        if(file.exists(ws_save_path)){
          
          .saveSettings()
          
          save(file=ws_full_name, 
               list=ls(envir = .pcrsim_env, all.names = TRUE),
               envir = .pcrsim_env)
          
          gmessage(message=paste("Project saved!\n\n", ws_full_name),
                   title="PCRsim",
                   icon ="info",
                   parent=w)
          
        } else {
          
          gmessage(message="The project directory was not found",
                   title="Directory not found",
                   icon = "error",
                   parent = w) 
        }
        
      } else {
        gmessage(message="A file name must be given",
                 title="File name required",
                 icon = "error",
                 parent = w) 
      }
      
    } else {
      
      gmessage(message="No project name or path!\nUse 'Save As' instead.",
               title="Property not found",
               icon = "error",
               parent = w) 
      
    }  
    
  } )
  
  addHandlerChanged(file_f1g1_ws_saveas_btn, handler = function (h, ...) {
    
    # Initiate flag.
    ok <- TRUE
    
    # Pick save location.
    ws_save_path <- gfile(text="Select a directory to save project in",
                          type="selectdir",
                          filter = list("R files" = list(patterns = c("*.R","*.Rdata"))),
                          multi=FALSE)
    
    # Ask for project name.
    ws_name <- ginput(message="Input project name",
                      text="",
                      title="Save as",
                      icon ="info",
                      parent=w)
    
    # Check if valid name.
    if(!is.na(ws_name) && !ws_name==""){
      
      # Create complete path.
      ws_full_name <- paste(ws_save_path, .separator, ws_name, ".RData", sep="")
      
      if(debug){
        print(ws_full_name)
      }
      
      # Check if file exist.
      if(file.exists(ws_full_name)){
        
        # Ask if overwrite.
        ok <- gconfirm(message=paste(ws_full_name,
                                     "\nalready exist!\n\n Overwrite?"),
                       title="Confirm", icon="question", parent=w)
        
      }
      
      # Check if ok to overwrite.
      if(ok){
        
        # Save project.
        if(file.exists(ws_save_path)){
          
          # Save project variables in workspace.
          assign(x=.ws_name_variable, value=ws_name, envir=.pcrsim_env)
          assign(x=.ws_path_variable, value=ws_save_path, envir=.pcrsim_env)
          
          # Save settings.        
          .saveSettings()
          
          # Save project. 
          save(file=ws_full_name, 
               list=ls(envir = .pcrsim_env, all.names = TRUE),
               envir = .pcrsim_env)
          
          gmessage(message=paste("Project saved!\n\n", ws_full_name),
                   title="STR-validator",
                   icon ="info",
                   parent=w)
          
        } else {
          
          gmessage(message="The project directory was not found",
                   title="Directory not found",
                   icon = "error",
                   parent = w)
          
        }
        
      } else {
        gmessage(message="Project was not saved!",
                 title="Info",
                 icon = "info",
                 parent = w) 
      }
      
      
    } else {
      gmessage(message="A file name must be given",
               title="File name required",
               icon = "error",
               parent = w) 
    }
    
  } )
  
  addHandlerChanged(file_f1g1_save_temp_btn, handler = function(h, ...) {
    
    # Ask for name.
    valName <- ginput(message="Name for simulation settings template", text="",
                      title="Input", icon = "info", parent=w)
    
    if(!is.na(valName)){
      
      # Create new environment.
      template <- new.env()
      
      # Save settings in environment.
      .saveSettings(env=template)
      
      # Save 'template' (environment) in workspace.
      assign(x=valName, value=template, envir=.pcrsim_env, inherits=FALSE)
      
      # Confirm with message.
      message(paste("Simulation settings template saved as", valName))
      
    }
    
    # Update workspace.
    .refreshLoaded()
    
  } )
  
  addHandlerChanged(file_f1g1_load_temp_btn, handler = function(h, ...) {
    
    # Get selected template name.
    valObj <- svalue(file_f1g1_tbl)
    
    # Check if any selcted.
    if (!is.na(valObj) && !is.null(valObj)){
      
      # Load template into workspace.
      .loadSavedSettings(valObj)
      
      # Confirm with message.
      message(paste("Simulation settings template", valObj , "loaded"))
      
    } 
    
  } )
  
  
  # DATASETS ------------------------------------------------------------------  
  
  if(debug){
    print("DATASETS")
  }
  
  
  file_f2 <- gframe(text = "Load dataframe from R workspace",
                    markup = FALSE,
                    pos = 0,
                    horizontal=TRUE,
                    container = file_gf,
                    expand=FALSE)
  
  file_f2g1 <- ggroup(horizontal=FALSE,
                      container = file_f2,
                      expand=FALSE)
  
  glabel("", container=file_f2g1) # Adds some space.
  
  file_f2g1_refresh_btn <- gbutton(text="Refresh",
                                   border=TRUE,
                                   container = file_f2g1)
  tooltip(file_f2g1_refresh_btn) <- "Refresh dropdown"
  
  file_f2g1_load_btn <- gbutton(text="Load",
                                border=TRUE,
                                container = file_f2g1) 
  tooltip(file_f2g1_load_btn) <- "Load dataset"
  
  glabel("", container=file_f2g1) # Adds some space.
  
  file_f2g1_drp <- gdroplist(items=c("<Select dataframe>", 
                                     strvalidator::listObjects(env=.pcrsim_env,
                                                               obj.class="data.frame")), 
                             selected = 1,
                             editable = FALSE,
                             container = file_f2g1) 
  
  # HANDLERS ..................................................................

  addHandlerChanged(file_f2g1_refresh_btn, handler = function (h, ...) {
    
    .refreshWs()
  } )
  
  addHandlerChanged(file_f2g1_load_btn, handler = function(h, ...) {
    
    # Get selected dataset name.
    val_obj <- svalue(file_f2g1_drp)
    
    if (!is.null(val_obj) && !is.na(val_obj)){
      
      # Load dataset.
      assign(val_obj, get(val_obj), envir=.pcrsim_env)
      
      # Update list.
      .refreshLoaded()
      
    } 
  } )
  
  # PROFILE PARAMETERS ########################################################
  
  glabel("", container=profile_gf) # Adds some space.
  
  # STR TYPING KIT ------------------------------------------------------------
  
  profile_f1 <- gframe(text = "STR typing kit",
                       markup = FALSE,
                       pos = 0,
                       horizontal = TRUE,
                       container = profile_gf,
                       expand = FALSE) 
  
  profile_f1g1 <- glayout(container = profile_f1)
  
  profile_f1g1[1,1] <- glabel("", container=profile_f1g1) # Adds some space.
  
  profile_f1g1[2,1] <- glabel(text="The following kit will be used:",
                              container=profile_f1g1, anchor=c(-1 ,0))
  
  profile_f1g1[2,2] <- profile_kit_drp <- gdroplist(items=getParameter(), 
                                                     selected = 1,
                                                     editable = FALSE,
                                                     container = profile_f1g1) 

  profile_f1g1[2,3] <- glabel(text="with settins from method:",
                              container=profile_f1g1, anchor=c(-1 ,0))
  
  profile_f1g1[2,4] <- profile_method_drp <- gdroplist(items=getParameter(what="Methods"), 
                                                     selected = 1,
                                                     editable = FALSE,
                                                     container = profile_f1g1) 
  
  
  profile_f1g1[3,1] <- glabel("", container=profile_f1g1) # Adds some space.
  
  # HANDLERS ..................................................................

  addHandlerChanged(profile_kit_drp, handler = function(h, ...) {
    val <- svalue (h$obj)
    
    if(is.null(val)){
      tmpkit <- NA
    } else {
      tmpkit <- val
    }
    markerList <- getKit(kit=tmpkit, what="Marker")
    
    dnaProfile <- data.frame(Marker=markerList,
                             Allele.1="NA",
                             Allele.2="NA",
                             stringsAsFactors=FALSE)
    
    profile_tbl[,] <- dnaProfile
    
    # Update 
    .refreshProfileTab()
    .refreshCETab()
    
  } )

  addHandlerChanged(profile_method_drp, handler = function(h, ...) {
    val <- svalue (h$obj)

    if(!is.null(val)){
      # Update 
      .refreshCETab()
    }
    
  } )

  
  # DNA PROFILE ---------------------------------------------------------------
  
  # Variables and constants ...................................................
  
  profile_options <- c("Manually enter a profile",
                       "Select a data frame from workspace",
                       "Import a profile from file", 
                       "Generate a random profile", 
                       "Generate a random profile for each simulation")
  
  profileText <- "Select for import..."
  
  # Get names of allele frequency databases.
  profileDb <- strvalidator::getDb()
  
  wsObj <- strvalidator::listObjects(env=.GlobalEnv, obj.class="data.frame")
  
  # GUI .......................................................................
  
  profile_f2 <- gframe(text = "Create profile", markup = FALSE, pos = 0,
                       horizontal = TRUE, container = profile_gf,
                       expand = FALSE) 
  
  profile_f2g1 <- glayout(container = profile_f2)
  
  profile_f2g1[2,1:2] <- profile_opt <- gradio(items=profile_options,
                                               selected = 1,
                                               horizontal = FALSE,
                                               container = profile_f2g1)
  
  profile_f2g1[3,1] <- glabel(text="Load from R workspace:",
                              container=profile_f2g1, anchor=c(-1 ,0))
  
  profile_f2g1[3,2] <- profile_ws_drp <- gdroplist(items=wsObj,
                                                   selected = 1,
                                                   editable = FALSE,
                                                   container = profile_f2g1) 
  enabled(profile_ws_drp) <- FALSE        
  
  profile_f2g1[3,3] <- profile_ws_btn <- gbutton(text = "Select",
                                                 border=TRUE, container = profile_f2g1)
  enabled(profile_ws_btn) <- FALSE      
  
  profile_f2g1[4,1] <- glabel(text="Import from file:",
                              container=profile_f2g1, anchor=c(-1 ,0))
  
  profile_f2g1[4,2:3] <- profile_file_browser <- gfilebrowse(text=profileText,
                                                           quote=FALSE,
                                                           type="open",
                                                           container=profile_f2g1)
  enabled(profile_file_browser) <- FALSE      

  profile_f2g1[5,1] <- glabel(text="Population database:",
                              container=profile_f2g1, anchor=c(-1 ,0))
  
  profile_f2g1[5,2] <- profile_pop_drp <- gdroplist(items=profileDb,
                                                    selected = 1,
                                                    editable = FALSE,
                                                    container = profile_f2g1) 
  enabled(profile_pop_drp) <- FALSE        
  
  profile_f2g1[5,3] <- profile_generate_btn <- gbutton(text = "Generate",
                                                       border=TRUE,
                                                       container = profile_f2g1)
  
  enabled(profile_generate_btn) <- FALSE      
  
  
  # HANDLERS ..................................................................
  
  addHandlerChanged(profile_opt, handler = function(h, ...) {
    val <- svalue (h$obj, index=TRUE)
    
    # Disable all.
    enabled(profile_ws_btn) <- FALSE      
    enabled(profile_ws_drp) <- FALSE      
    enabled(profile_file_browser) <- FALSE      
    enabled(profile_pop_drp) <- FALSE      
    enabled(profile_generate_btn) <- FALSE      
    enabled(profile_tbl) <- FALSE
    
    # Enable allowed buttons.    
    if (val == 1 ){
      enabled(profile_tbl) <- TRUE
      
    } else if (val == 2 ){
      enabled(profile_ws_drp) <- TRUE      
      enabled(profile_ws_btn) <- TRUE
      wsObj <- strvalidator::listObjects(env=.GlobalEnv, obj.class="data.frame")
      profile_ws_drp[] <- wsObj
      
    } else if (val == 3 ){
      enabled(profile_file_browser) <- TRUE

    } else if ( val == 4) {
      enabled(profile_pop_drp) <- TRUE      
      enabled(profile_generate_btn) <- TRUE
      
    } else if ( val == 5) {
      enabled(profile_pop_drp) <- TRUE      
      
    }
    
  } )
  
  addHandlerChanged(profile_ws_btn, handler = function(h, ...) {
    val_obj <- svalue(profile_ws_drp)
    valExtDebug <- svalue(file_f1g1_ext_debug_chk)

    if (!is.null(val_obj) && !is.na(val_obj)){

      dnaProfile <- get(val_obj)
      
      if(all(c("Marker", "Allele.1", "Allele.2") %in% names(dnaProfile))){

        if("Sample.Name" %in% names(dnaProfile)){
          # Find number of profiles and update number of contributors.
          svalue(sim_sim_edt) <- length(unique(dnaProfile$Sample.Name))
        } else {
          message("Imported profiles contain no 'Sample.Name' column")
          message("Number of contributors not changed (check manually in tab 'Simulation')")
        }
        
        # Get Marker and Allele columns.
        dnaProfile <- strvalidator::trim(data=dnaProfile, samples=NULL,
                                         columns="Marker|Allele", ignore.case=TRUE,
                                         invert.s=FALSE, invert.c=FALSE, rm.na.col=TRUE,
                                         rm.empty.col=TRUE, missing=NA, debug=valExtDebug)

        if(!is.null(dnaProfile)){

          # Update profile table.
          profile_tbl[,] <- dnaProfile
          message("Profile table updated!")
          
        } else {
          
          message("Dataset must contain columns 'Marker', 'Allele.1' and 'Allele.2'.")
          
        }
        
      } else {
        
        message("Dataset must contain columns 'Marker', 'Allele.1' and 'Allele.2'.")
        
      }
      
    }

  } )
  
  addHandlerChanged(profile_file_browser, handler = function(h, ...) {
    val <- svalue(profile_file_browser)
    valExtDebug <- svalue(file_f1g1_ext_debug_chk)
    
    if(enabled(profile_file_browser)){

      if (file.exists(val)){
        
        dnaProfile <- strvalidator::import(import.file=val, time.stamp=FALSE,
                                           file.name=FALSE, debug=valExtDebug)
        
        if(all(c("Marker", "Allele.1", "Allele.2") %in% names(dnaProfile))){
          
          if("Sample.Name" %in% names(dnaProfile)){
            # Find number of profiles and update number of contributors.
            svalue(sim_sim_edt) <- length(unique(dnaProfile$Sample.Name))
          } else {
            message("Imported profiles contain no 'Sample.Name' column")
            message("Number of contributors not changed (check manually in tab 'Simulation')")
          }
          
          # Get Marker and Allele columns.
          dnaProfile <- strvalidator::trim(data=dnaProfile, samples=NULL,
                                           columns="Marker|Allele", ignore.case=TRUE,
                                           invert.s=FALSE, invert.c=FALSE, rm.na.col=TRUE,
                                           rm.empty.col=TRUE, missing=NA, debug=valExtDebug)
          
          # Update profile table.        
          profile_tbl[,] <- dnaProfile
          message("Profile table updated!")
          
        } else {
          
          message("Dataset must contain columns 'Marker', 'Allele.1' and 'Allele.2'.")
        }
        
      } else {
        
        gmessage(message="File not found!", title="Error", icon = "error")      
        
      } 
      
    } else {
      
      message(paste("Please select the option '", profile_opt[3],
                    "' before choosing a file.", sep=""))
      
    }
    
    
  } )

  addHandlerChanged(profile_generate_btn, handler = function(h, ...) {
    
    # Get values.
    val_kit <- svalue(profile_kit_drp)
    val_db <- svalue(profile_pop_drp)
    valExtDebug <- svalue(file_f1g1_ext_debug_chk)
    
    if (!is.na(val_db) && !is.null(val_kit)){
      
      # Load database file
      profile_db <- strvalidator::getDb(db.name.or.index=val_db, debug=valExtDebug)
      
      if(debug){
        print("Database used to generate profile.")
        str(profile_db)
        head(profile_db)
        tail(profile_db)
      }
      
      # Simulate sample.
      generatedData <- simProfile(data=NULL, sim=1, kit=val_kit,
                                  name=NULL, db=profile_db,
                                  debug=valExtDebug)
      
      if(debug){
        print("generatedData:")
        print(head(generatedData))
      }
      
      # Check result.
      if(nrow(generatedData) > 0){
        
        # Extract relevant columns.
        tmp <- generatedData[, c("Marker","Allele")]
        
        # Get rows for allele 1 and allele 2.
        a1 <- seq(from=1, to=nrow(tmp), by=2)
        a2 <- seq(from=2, to=nrow(tmp), by=2)
        
        # Create dataframe.
        tmpdf <- data.frame(Marker=tmp$Marker[a1],
                            Allele.1=tmp$Allele[a1],
                            Allele.2=tmp$Allele[a2],
                            stringsAsFactors=FALSE)
        
        # Update profile table.        
        profile_tbl[,] <- tmpdf
        message("Profile table updated!")
        
      } else {
        
        warning("Function 'simProfile' returned no data!")
        
      }
      
    } 
    
  } )
  
  # EDIT PROFILE ----------------------------------------------------------------
  
  profile_f3 <- gframe(text = "View and edit profile", markup = FALSE, pos = 0,
                       horizontal = TRUE, container = profile_gf, expand = FALSE) 
  
  profile_f3g1 <- glayout(container = profile_f3)
  
  if(is.null(svalue(profile_kit_drp))){
    tmpkit <- NA
  } else {
    tmpkit <- svalue(profile_kit_drp)
  }
  markerList <- getKit(kit=tmpkit, what="Marker")
  
  profile_profile_item <- data.frame(Marker=markerList,
                                     Allele.1="NA",
                                     Allele.2="NA",
                                     stringsAsFactors=FALSE)
  
  profile_f3g1[1,1] <- profile_tbl <- gdf(items=profile_profile_item,
                                          container = profile_f3g1)
  
  
  size(profile_tbl) <- c(600, 300)
  
  # SAMPLE PARAMETERS ###########################################################
  
  glabel("", container=sample_gf) # Adds some space.
  
  # SAMPLE NAME ---------------------------------------------------------------
  
  sample_f1_txt_width <- 20
  
  sample_f1 <- gframe(text = "General",
                      markup = FALSE,
                      pos = 0,
                      horizontal = TRUE,
                      container = sample_gf,
                      expand = FALSE) 
  
  
  sample_f1g1 <- glayout(container = sample_f1)
  
  sample_f1g1[1,1] <- glabel("", container=sample_f1g1) # Adds some space.
  
  
  sample_f1g1[2,1] <- glabel(text="Sample name:",
                             container=sample_f1g1, anchor=c(-1 ,0))
  
  sample_f1g1[2,2] <- sample_name_edt <- gedit(text="",
                                               width=sample_f1_txt_width,
                                               container=sample_f1g1)
  
  sample_f1g1[3,1:2] <- sample_mix_chk <- gcheckbox(text="", checked = FALSE,
                              container=sample_f1g1, anchor=c(-1 ,0))

  msgtmp <- paste("To simulate a mixture different values (separated by comma)\n",
                  "can be given per sample e.g. number of cells '1000,500' and haplotype flag 'T,F'.")

  sample_f1g1[4,1:2] <- glabel(text=msgtmp, container=sample_f1g1)

  glabel("", container=sample_f1) # Adds some space.

  # AMOUNT OF DNA ---------------------------------------------------------------
  
  sample_f2_txt_width <- 10
  
  sample_f2 <- gframe(text = "Amount of DNA",
                      markup = FALSE,
                      pos = 0,
                      horizontal=FALSE,
                      container = sample_gf,
                      expand=FALSE) 
  
  glabel("", container=sample_f2) # Adds some space.
  
  glabel("Estimate the amount of DNA by:", container=sample_f2, anchor=c(-1 ,0))
  
  sample_by_items <- c("Number of cells",
                       "DNA concentration",
                       "Regression (don't simulate degradation if this is used)")
  
  sample_by_opt <- gradio(items=sample_by_items, selected = 1,
                          horizontal=FALSE, container=sample_f2, anchor=c(-1 ,0))
  
  
  sample_f2g1 <- glayout(container = sample_f2)
  
  
  sample_f2g1[2,1] <- glabel(text="Number of cells:",
                             container=sample_f2g1, anchor=c(-1 ,0))
  
  sample_f2g1[2,2] <- sample_ncells_edt <- gedit(text="",
                                                 width=sample_f2_txt_width,
                                                 container=sample_f2g1)
  
  sample_f2g1[2,3] <- glabel(text="Standard deviation:",
                             container=sample_f2g1, anchor=c(-1 ,0))
  
  sample_f2g1[2,4] <- sample_ncells_sd_edt <- gedit(text="0",
                                                    width=sample_f2_txt_width,
                                                    container=sample_f2g1)
  
  sample_f2g1[3,1] <- glabel(text="DNA concentration (ng/\u00B5l):",
                             container=sample_f2g1, anchor=c(-1 ,0))
  
  sample_f2g1[3,2] <- sample_conc_edt <- gedit(text="",
                                               width=sample_f2_txt_width,
                                               container=sample_f2g1)
  
  sample_f2g1[3,3] <- glabel(text="Standard deviation:",
                             container=sample_f2g1, anchor=c(-1 ,0))
  
  sample_f2g1[3,4] <- sample_conc_sd_edt <- gedit(text="0",
                                                  width=sample_f2_txt_width,
                                                  container=sample_f2g1)
  
  
  sample_f2g1[4,1] <- glabel(text="DNA per diploid cell (pg):",
                             container=sample_f2g1, anchor=c(-1 ,0))
  
  sample_f2g1[4,2] <- sample_celldna_edt <- gedit(text="6",
                                                  width=sample_f2_txt_width,
                                                  container=sample_f2g1)
  
  
  sample_f2g1[5,1] <- glabel(text="Sample volume (\u00B5l):",
                             container=sample_f2g1, anchor=c(-1 ,0))
  
  sample_f2g1[5,2] <- sample_vol_edt <- gedit(text="",
                                              width=sample_f2_txt_width,
                                              container=sample_f2g1)
  
  sample_f2g1[5,3] <- glabel(text="Standard deviation:",
                             container=sample_f2g1, anchor=c(-1 ,0))
  
  sample_f2g1[5,4] <- sample_vol_sd_edt <- gedit(text="0",
                                                 width=sample_f2_txt_width,
                                                 container=sample_f2g1)
  
  sample_f2g1[6,1] <- glabel(text="Slope:",
                             container=sample_f2g1, anchor=c(-1 ,0))
  
  sample_f2g1[6,2] <- sample_slope_edt <- gedit(text="-0.008",
                                                width=sample_f2_txt_width,
                                                container=sample_f2g1)
  
  sample_f2g1[6,3] <- glabel(text="Intercept:",
                             container=sample_f2g1, anchor=c(-1 ,0))
  
  sample_f2g1[6,4] <- sample_intercept_edt <- gedit(text="2.5",
                                                    width=sample_f2_txt_width,
                                                    container=sample_f2g1)
  
  sample_f2g1[7,1] <- glabel(text="Concentration @ 100bp (ng/\u00B5l):",
                             container=sample_f2g1, anchor=c(-1 ,0))
  
  sample_f2g1[7,2] <- sample_ex1_edt <- gedit(text="", width=sample_f2_txt_width,
                                              container=sample_f2g1)
  enabled(sample_ex1_edt) <- FALSE
  
  
  sample_f2g1[8,1] <- glabel(text="Concentration @ 400bp (ng/\u00B5l):",
                             container=sample_f2g1, anchor=c(-1 ,0))
  
  sample_f2g1[8,2] <- sample_ex2_edt <- gedit(text="", width=sample_f2_txt_width,
                                              container=sample_f2g1)
  enabled(sample_ex2_edt) <- FALSE


  sample_f2g1[9,1:2] <- glabel("Haplotype flag (FALSE=diploid, TRUE=haploid):",
                               container=sample_f2g1, anchor=c(-1 ,0))
  
  sampleTypeItems <- c("FALSE", "TRUE")
  
  sample_f2g1[9,3:4] <- sample_type_cmb <- gcombobox(items=sampleTypeItems,
                                                   selected = 1, editable = TRUE,
                                                   container=sample_f2g1)
  
  
  glabel("", container=sample_f2) # Adds some space.

  # HANDLERS ..................................................................
  
  addHandlerChanged(sample_by_opt, handler = function(h, ...) {
    
    .refreshSampleTab()
    
  } )
  
  addHandlerKeystroke(sample_slope_edt, handler = function(h, ...) {
    
    .refreshSampleTab()
    
  } )
  
  addHandlerKeystroke(sample_intercept_edt, handler = function(h, ...) {
    
    .refreshSampleTab()
    
  } )
  

  # DEGRADATION PARAMETERS ####################################################

  glabel("", container=deg_gf) # Adds some space.
  
  # DEGRADATION -----------------------------------------------------------------
  
  deg_f1_txt_w <- 10
  
  deg_f1 <- gframe(text = "Degradation", markup = FALSE, pos = 0,
                   horizontal=TRUE, container = deg_gf, expand=FALSE) 
  
  deg_f1g1 <- glayout(container = deg_f1)
  
  deg_f1g1[1,1] <- glabel(text = "Current P(deg)=", container = deg_f1g1,
                          anchor = c(1 ,0))
  deg_f1g1[1,2] <- deg_pam_lbl <- glabel(text = "NA", container = deg_f1g1)

  deg_f1g1[2,1] <- glabel(text = "Quantification amplicon size (bp):",
                          container=deg_f1g1, anchor = c(1 ,0))
  deg_f1g1[2,2] <- deg_target_edt <- gedit(text = "80", container = deg_f1g1,
                                           width = deg_f1_txt_w)
  
  
  deg_tmp_items <- c("Specify manually",
                     "Calculate from concentrations")
  
  deg_f1g1[3:3,1] <- deg_opt <- gradio(items = deg_tmp_items, selected = 1,
                                       horizontal = FALSE, container = deg_f1g1)
  
  deg_f1g1[4,1] <- glabel(text = "Degradation parameter (probability per bp):",
                          container = deg_f1g1, anchor = c(-1 ,0))
  
  deg_f1g1[4,2] <- deg_spn <- gspinbutton(from = 0, to = 1, by = 0.001,
                                          value = 0.01, digits = 6,
                                          container = deg_f1g1)
  tooltip(deg_spn) <- "NB! Use the spin buttons, or type manually AND hit enter to confirm!"
  
  deg_f1g1[5,1] <- glabel(text = "Concentration in ng/\u00B5l (separated by comma):",
                          container = deg_f1g1, anchor = c(-1 ,0))
  
  deg_f1g1[5,2] <- deg_conc_edt <- gedit(text = "", width = deg_f1_txt_w,
                                         container = deg_f1g1)
  
  deg_f1g1[6,1] <- glabel(text = "Size in base pair (separated by comma):",
                          container = deg_f1g1, anchor = c(-1 ,0))
  
  deg_f1g1[6,2] <- deg_size_edt <- gedit(text = "", width = deg_f1_txt_w,
                                         container = deg_f1g1)

  deg_f1g1[7,1] <- glabel(text = "", container = deg_f1g1) # Adds some space.
  
  # HANDLERS ..................................................................
  
  addHandlerChanged(deg_opt, handler = function(h, ...) {
    
    .refreshDegradationTab()
    
  } )

  addHandlerChanged(deg_spn, handler = function(h, ...) {
    
    .refreshDegradationTab()
    
  } )

  addHandlerKeystroke(deg_spn, handler = function(h, ...) {
    
    .refreshDegradationTab()
    
  } )
  
  addHandlerKeystroke(deg_conc_edt, handler = function(h, ...) {
    
    .refreshDegradationTab()
    
  } )
  
  addHandlerKeystroke(deg_size_edt, handler = function(h, ...) {
    
    .refreshDegradationTab()
    
  } )
  
  # EXTRACTION PARAMETERS #######################################################
  
  glabel("", container=ex_gf) # Adds some space.
  
  # EFFICIENCY ------------------------------------------------------------------
  
  ex_frame_1_txt_width <- 6 
  
  ex_frame_1 <- gframe(text = "Extraction efficiency", markup = FALSE, pos = 0,
                       horizontal=TRUE, container = ex_gf, expand=FALSE) 
  
  ex_grid_1 <- glayout(container = ex_frame_1)
  
  ex_grid_1[1,1] <- glabel("", container=ex_grid_1) # Adds some space.

  ex_grid_1[2,1] <- glabel(text="Extraction efficiency:",
                           container=ex_grid_1, anchor=c(-1 ,0))
  
  ex_grid_1[2,2] <- ex_eff_spn <- gspinbutton(from = 0, to = 1, by = 0.05,
                                              value = 0.3, digits = 2,
                                              container = ex_grid_1)
  tooltip(ex_eff_spn) <- "NB! Use the spin buttons, or type manually AND hit enter to confirm!"
  
  
  ex_grid_1[2,3] <- glabel(text="Standard deviation:",
                           container=ex_grid_1, anchor=c(-1 ,0))
  
  ex_grid_1[2,4] <- ex_eff_sd_edt <- gedit(text="0",
                                           width=ex_frame_1_txt_width,
                                           container=ex_grid_1)
  
  ex_grid_1[3,1] <- glabel("", container=ex_grid_1) # Adds some space.
  
  # VOLUME ----------------------------------------------------------------------
  
  ex_frame_2_txt_width <- 6 
  
  ex_frame_2 <- gframe(text = "Extraction volume",
                       markup = FALSE,
                       pos = 0,
                       horizontal=TRUE,
                       container = ex_gf,
                       expand=FALSE) 
  
  
  ex_grid_2 <- glayout(container = ex_frame_2)
  
  ex_grid_2[1,1] <- glabel("", container=ex_grid_2) # Adds some space.
  
  
  ex_grid_2[2,1] <- glabel(text="Extraction volume (\u00B5l):",
                           container=ex_grid_2, anchor=c(-1 ,0))
  
  ex_grid_2[2,2] <- ex_vol_edt <- gedit(text="100",
                                        width=ex_frame_2_txt_width,
                                        container=ex_grid_2)
  
  ex_grid_2[2,3] <- glabel(text="Standard deviation:",
                           container=ex_grid_2, anchor=c(-1 ,0))
  
  ex_grid_2[2,4] <- ex_vol_sd_edt <- gedit(text="0",
                                           width=ex_frame_2_txt_width,
                                           container=ex_grid_2)
  
  ex_grid_2[3,1] <- glabel("", container=ex_grid_2) # Adds some space.
  
  
  # NORMALIZATION PARAMETERS #######################################################
  
  glabel("", container=dil_gf) # Adds some space.
  
  # AMOUNT --------------------------------------------------------------------
  
  dil_f1_txt_width <- 6 
  
  dil_f1 <- gframe(text = "Current settings",
                   markup = FALSE,
                   pos = 0,
                   horizontal=TRUE,
                   container = dil_gf,
                   expand=FALSE) 
  
  dil_f1g1 <- glayout(container = dil_f1)
  
  dil_f1g1[1,1] <- glabel("Optimal amount for PCR is set to:",
                          container=dil_f1g1, anchor=c(-1 ,0))
  dil_f1g1[1,2] <- dil_amount_lbl <- glabel(text="", container=dil_f1g1, anchor=c(-1 ,0))
  
  dil_f1g1[1,3] <- glabel(text="ng", container=dil_f1g1, anchor=c(-1 ,0))
  
  dil_f1g1[2,1] <- glabel("Aliquot for PCR is set to:",
                          container=dil_f1g1, anchor=c(-1 ,0))
  dil_f1g1[2,2] <- dil_aliquot_lbl <- glabel(text="", container=dil_f1g1, anchor=c(-1 ,0))
  
  dil_f1g1[2,3] <- glabel(text="\u00B5l", container=dil_f1g1, anchor=c(-1 ,0))
  
  dil_f1g1[3,1] <- glabel("", container=dil_f1g1) # Make some space.
  
  # CONCENTRATION -------------------------------------------------------------
  
  dil_f2_txt_width <- 6 
  
  dil_f2 <- gframe(text = "Target concentration",
                   markup = FALSE,
                   pos = 0,
                   horizontal=TRUE,
                   container = dil_gf,
                   expand=FALSE) 
  
  dil_f2g1 <- glayout(container = dil_f2)
  
  dil_tmp_items <- c("Calculated from settings:",
                     "Specify manually:")
  
  dil_f2g1[1:2,1] <- dil_target_opt <- gradio(items=dil_tmp_items, selected = 1,
                                              horizontal = FALSE, container = dil_f2g1)
  
  dil_f2g1[1,2] <- dil_target_lbl <- glabel(text="", container=dil_f2g1)
  
  dil_f2g1[1,3] <- glabel(text="ng/\u00B5l", container=dil_f2g1)
  
  dil_f2g1[2,2] <- dil_target_edt <- gedit(text="",
                                           width=dil_f2_txt_width,
                                           container=dil_f2g1)
  
  
  dil_f2g1[2,3] <- glabel(text="ng/\u00B5l", container=dil_f2g1)
  
  dil_f2g1[3,1] <- glabel(text="Tolerance for target:",
                          container=dil_f2g1, anchor=c(-1 ,0))
  
  dil_f2g1[3,2] <- dil_tolerance_edt <- gedit(text="0",
                                              width=dil_f2_txt_width,
                                              container=dil_f2g1)
  
  dil_f2g1[3,3] <- dil_max_lbl <- glabel(text="", container=dil_f2g1)
  
  dil_f2g1[4,1] <- glabel("", container=dil_f2g1) # Adds some space.
  
  addHandlerKeystroke(dil_tolerance_edt, handler = function(h, ...) {
    
    .refreshDilutionTab()
    
  } )
  
  
  addHandlerChanged(dil_target_opt, handler = function(h, ...) {
    
    .refreshDilutionTab()
    
  } )
  
  # ACCURACY ------------------------------------------------------------------
  
  dil_f3_txt_width <- 4 
  
  dil_f3 <- gframe(text = "Pipetting accuracy",
                   markup = FALSE,
                   pos = 0,
                   horizontal=TRUE,
                   container = dil_gf,
                   expand=FALSE) 
  
  dil_f3g1 <- glayout(container = dil_f3)
  
  dil_f3g1[1,1] <- glabel(text="Minimum pipetting volume:",
                          container=dil_f3g1)
  
  
  dil_f3g1[1,2] <- dil_accuracy_edt <- gedit(text="1",
                                             width=dil_f3_txt_width,
                                             container=dil_f3g1)
  
  dil_f3g1[1,3] <- glabel(text="ng/\u00B5l", container=dil_f3g1)
  
  dil_f3g1[2,1] <- glabel("", container=dil_f3g1) # Adds some space.
  
  # VOLUME --------------------------------------------------------------------
  
  dil_f4_txt_width <- 6 
  
  dil_f4 <- gframe(text = "Final volume",
                   markup = FALSE,
                   pos = 0,
                   horizontal=TRUE,
                   container = dil_gf,
                   expand=FALSE) 
  
  dil_f4g1 <- glayout(container = dil_f4)
  
  dil_tmp_items <- c("Use extraction volume:",
                     "Specify manually:")
  
  dil_f4g1[1:2,1] <- dil_volume_opt <- gradio(items=dil_tmp_items, selected = 1,
                                              horizontal = FALSE, container = dil_f4g1)
  
  
  dil_f4g1[1,2] <- dil_volume_lbl <- glabel(text="", container=dil_f4g1)
  
  dil_f4g1[1,3] <- glabel(text="\u00B5l", container=dil_f4g1)
  
  dil_f4g1[2,2] <- dil_volume_edt <- gedit(text="",
                                           width=dil_f4_txt_width,
                                           container=dil_f4g1)
  
  dil_f4g1[2,3] <- glabel(text="\u00B5l", container=dil_f4g1)
  
  dil_f4g1[3,1] <- glabel("", container=dil_f4g1) # Adds some space.
  
  addHandlerChanged(dil_volume_opt, handler = function(h, ...) {
    
    .refreshDilutionTab()
    
  } )
  
  # AMPLIFICATION PARAMETERS ####################################################
  
  glabel("", container=amp_gf) # Adds some space.
  
  # ALIQUOT --------------------------------------------------------------------
  
  amp_f1_txt_width <- 6 
  
  amp_f1 <- gframe(text = "Aliquot",
                   markup = FALSE,
                   pos = 0,
                   horizontal=TRUE,
                   container = amp_gf,
                   expand=FALSE) 
  
  
  amp_f1g1 <- glayout(container = amp_f1)
  
  amp_f1g1[1,1] <- glabel("", container=amp_f1g1) # Adds some space.
  
  
  
  amp_f1g1[2,1] <- glabel(text="Aliquot for PCR (\u00B5l):",
                          container=amp_f1g1, anchor=c(-1 ,0))
  
  amp_f1g1[2,2] <- amp_aliquot_edt <- gedit(text="10",
                                            width=amp_f1_txt_width,
                                            container=amp_f1g1)
  
  amp_f1g1[2,3] <- glabel(text="Standard deviation:",
                          container=amp_f1g1, anchor=c(-1 ,0))
  
  amp_f1g1[2,4] <- amp_aliquot_sd_edt <- gedit(text="0",
                                               width=amp_f1_txt_width,
                                               container=amp_f1g1)
  
  amp_f1g1[3,1] <- glabel("", container=amp_f1g1) # Adds some space.
  
  # PCR -------------------------------------------------------------------------
  
  amp_f2_txt_width <- 6 
  
  amp_f2 <- gframe(text = "PCR",
                   markup = FALSE,
                   pos = 0,
                   horizontal=TRUE,
                   container = amp_gf,
                   expand=FALSE) 
  
  
  amp_f2g1 <- glayout(container = amp_f2)
  
  amp_f2g1[1,1] <- glabel("", container=amp_f2g1) # Adds some space.
  
  
  amp_f2g1[2,1] <- glabel(text="Optimal amount of DNA:",
                          container=amp_f2g1, anchor=c(-1 ,0))
  
  amp_f2g1[2,2] <- amp_amount_edt <- gedit(text="0.5",
                                           width=amp_f2_txt_width,
                                           container=amp_f2g1)
  
  amp_f2g1[3,1] <- glabel(text="Total amplification volume:",
                          container=amp_f2g1, anchor=c(-1 ,0))
  
  amp_f2g1[3,2] <- amp_total_vol_edt <- gedit(text="25",
                                              width=amp_f2_txt_width,
                                              container=amp_f2g1)
  
  amp_f2g1[3,3] <- glabel(text="Standard deviation:",
                          container=amp_f2g1, anchor=c(-1 ,0))
  
  amp_f2g1[3,4] <- amp_total_vol_sd_edt <- gedit(text="0",
                                                 width=amp_f2_txt_width,
                                                 container=amp_f2g1)
  
  amp_f2g1[4,1] <- glabel("Number of PCR cycles:",
                          container=amp_f2g1, anchor=c(-1 ,0))

  amp_f2g1[4,2] <- cyc_sb <- gspinbutton(from=1, to = 50, by =1,
                                         value=30, container = amp_f2g1)
  tooltip(cyc_sb) <- "NB! Use the spin buttons, or type manually AND hit enter to confirm!"

  amp_f2g1[5,1] <- glabel("", container=amp_f2g1) # Adds some space.
  
  # EFFICIENCY ----------------------------------------------------------------
  
  amp_f3_txt_width <- 6 
  
  amp_f3 <- gframe(text = "PCR efficiency",
                   markup = FALSE,
                   pos = 0,
                   horizontal = FALSE,
                   container = amp_gf,
                   expand = TRUE) 
  
  amp_current_lbl <- glabel("", container=amp_f3, anchor=c(-1 ,0))

  amp_f3_items <- c("Get marker specific PCR efficiency from kit information.",
                    "Use custom PCR efficiency for all markers:")
  
  amp_f3_eff_opt <- gradio(items=amp_f3_items, selected = 1,
                           horizontal = FALSE, container = amp_f3)
  
  
  amp_f3g1 <- glayout(container = amp_f3)
  
  amp_f3g1[1,1] <- glabel("", container=amp_f3g1) # Adds some space.
  
  amp_f3g1[2,1] <- glabel(text="Amplification probability:",
                          container=amp_f3g1, anchor=c(-1 ,0))
  
  amp_f3_eff_spn <- gspinbutton(from = 0, to = 1, by = 0.001,
                                value = 0.90, digits = 3, container=amp_f3g1)
  amp_f3g1[2,2] <- amp_f3_eff_spn
  
  amp_f3g1[2,3] <- glabel(text="Standard deviation:",
                          container=amp_f3g1, anchor=c(-1 ,0))
  
  amp_f3g1[2,4] <- amp_f3_eff_sd_edt <- gedit(text="0",
                                              width=amp_f3_txt_width,
                                              container=amp_f3g1)
  
  amp_f3g1[3,1] <- glabel(text="Stutter probability:",
                          container=amp_f3g1, anchor=c(-1 ,0))
  
  amp_eff_stutt_spn <- gspinbutton(from = 0, to = 1, by = 0.001,
                                   value = 0.005, digits = 3, container=amp_f3g1)
  amp_f3g1[3,2] <- amp_eff_stutt_spn
  
  amp_f3g1[3,3] <- glabel(text="Standard deviation:", container=amp_f3g1,
                          anchor=c(-1 ,0))
  
  amp_eff_stutt_sd_edt <- gedit(text="0", width=amp_f3_txt_width,
                                container=amp_f3g1)
  amp_f3g1[3,4] <- amp_eff_stutt_sd_edt
  
    
  glabel("", container=amp_f3) # Adds some space.
  
  amp_stutter_chk <- gcheckbox(text="Simulate stutters",
                               checked=TRUE, container = amp_f3)
  
  addHandlerChanged(amp_f3_eff_opt, handler = function(h, ...) {
    
    .refreshAmplificationTab()
    
  } )
  
  addHandlerChanged(amp_stutter_chk, handler = function(h, ...) {
    
    .refreshAmplificationTab()
    
  } )
  
  
  # ANALYSIS PARAMETERS #######################################################
  
  glabel("", container=ce_gf) # Adds some space.
  
  # CE ------------------------------------------------------------------------
  
  ce_f2_txt_width <- 8 
  
  ce_f2 <- gframe(text = "Capillary electrophoresis",
                  markup = FALSE,
                  pos = 0,
                  horizontal=TRUE,
                  container = ce_gf,
                  expand=FALSE) 
  
  
  ce_f2g1 <- glayout(container = ce_f2)
  
  ce_f2g1[1,1] <- glabel("", container=ce_f2g1) # Adds some space.
  
  
  ce_f2g1[2,1] <- glabel(text="Volume PCR product analysed (\u00B5l):",
                         container=ce_f2g1, anchor=c(-1 ,0))
  
  ce_f2g1[2,2] <- ce_aliquot_edt <- gedit(text="1",
                                          width=ce_f2_txt_width,
                                          initial.msg="",
                                          container=ce_f2g1)
  
  ce_f2g1[2,3] <- glabel(text="Standard deviation:",
                         container=ce_f2g1, anchor=c(-1 ,0))
  
  ce_f2g1[2,4] <- ce_aliquot_sd_edt <- gedit(text="0",
                                             width=ce_f2_txt_width,
                                             initial.msg="",
                                             container=ce_f2g1)
  
  ce_f2g1[3,1] <- glabel("", container=ce_f2g1) # Adds some space.
  
  # LOAD FROM PARAMETERS ------------------------------------------------------

  ce_f4_txt_width <- 8 
  
  ce_f4 <- gframe(text = "Load parameters",
                  markup = FALSE,
                  pos = 0,
                  horizontal = FALSE,
                  container = ce_gf,
                  expand = FALSE) 
  
  ce_current_lbl <- glabel("", container=ce_f4, anchor=c(-1 ,0))
  
  ce_settings_opt <- gradio(items = c("Load threshold and scaling from parameter file",
                                      "Use custom settings for threshold and scaling"),
                            selected = 1, container = ce_f4)
  
  
  addHandlerChanged(ce_settings_opt, handler = function(h, ...) {
    
    .refreshCETab()

  } )
  
  # DETECTION THRESHOLD -------------------------------------------------------
  
  ce_f3_txt_width <- 8 
  
  ce_f3 <- gframe(text = "Detection threshold",
                  markup = FALSE,
                  pos = 0,
                  horizontal=TRUE,
                  container = ce_gf,
                  expand=FALSE) 
  
  
  ce_f3g1 <- glayout(container = ce_f3)
  
  ce_f3g1[1,1:2] <- glabel("lm(log(M) ~ log(H))", container=ce_f3g1)
  
  ce_f3g1[2,1] <- glabel(text="Threshold intercept:",
                         container=deg_f1g1, anchor=c(-1 ,0))
  
  ce_f3g1[2,2] <- ce_t_intercept_edt <- gedit(text="13.9880",
                                              width=ce_f3_txt_width,
                                              container=ce_f3g1)
  
  ce_f3g1[3,1] <- glabel(text="Threshold slope:",
                         container=ce_f3g1, anchor=c(-1 ,0))
  
  ce_f3g1[3,2] <- ce_t_slope_edt <- gedit(text="0.9047",
                                          width=ce_f3_txt_width,
                                          container=ce_f3g1)
  
  ce_f3g1[4,1] <- glabel(text="Residual standard error:",
                         container=ce_f3g1, anchor=c(-1 ,0))
  
  ce_f3g1[4,2] <- ce_t_residual_edt <- gedit(text="0.5951",
                                             width=ce_f3_txt_width,
                                             container=ce_f3g1)
  
  
  ce_f3g1[5,1] <- glabel("", container=ce_f3g1) # Adds some space.
  
  # SCALING -------------------------------------------------------------------
  
  ce_f1_txt_width <- 8 
  
  ce_f1 <- gframe(text = "Peak height scaling", markup = FALSE, pos = 0,
                  horizontal=TRUE, container = ce_gf, expand=FALSE) 
  
  ce_f1g1 <- glayout(container = ce_f1)
  
  ce_f1g1[1,1:2] <- glabel("lm(log(H) ~ log(M))", container=ce_f1g1)
  
  
  ce_f1g1[2,1] <- glabel(text="Scaling intercept:",
                         container=deg_f1g1, anchor=c(-1 ,0))
  
  ce_f1g1[2,2] <- ce_intercept_edt <- gedit(text="-10.4812",
                                            width=ce_f1_txt_width,
                                            container=ce_f1g1)
  
  ce_f1g1[3,1] <- glabel(text="Scaling slope:",
                         container=ce_f1g1, anchor=c(-1 ,0))
  
  ce_f1g1[3,2] <- ce_slope_edt <- gedit(text="0.8590",
                                        width=ce_f1_txt_width,
                                        container=ce_f1g1)
  
  ce_f1g1[4,1] <- glabel(text="Residual standard error:",
                         container=ce_f1g1, anchor=c(-1 ,0))
  
  ce_f1g1[4,2] <- ce_residual_edt <- gedit(text="0.5798",
                                           width=ce_f1_txt_width,
                                           container=ce_f1g1)
  
  
  ce_f1g1[5,1] <- glabel("", container=ce_f1g1) # Adds some space.
  
  # ELECTROPHEROGRAM #######################################################
  
  glabel("", container=epg_gf) # Adds some space.
  gv <- epg_gf # As used in the copied code below.
  
  # NB! Below code directly copied from strvalidator::generateEPG_gui()
  # NB! Comment out the plot title field.
  #// START OF COPIED CODE.
  # FRAME 1 ###################################################################
  
  f1 <- gframe(text = "Options", horizontal=FALSE, spacing = 10, container = gv)
  
  #glabel(text="Plot title:", anchor=c(-1 ,0), container=f1)
  #f1_title_edt <- gedit(text="", width=25, container=f1)
  
  # Layout --------------------------------------------------------------------
  f1g1 <- glayout(container = f1, spacing = 1)
  
  f1g1[1,1] <- glabel(text = "Axis scales:   ", anchor=c(-1 ,0), container = f1g1)  
  f1g1[2:3,1] <- f1_scale_opt <- gradio(items = c("free", "free_y", "free_x"),
                                        selected = 2, horizontal = FALSE, container = f1g1)
  
  f1g1[1,2] <- glabel(text="Allele label text size:", container=f1g1)
  f1g1[1,3] <- f1_size_spb <- gspinbutton(from=0, to=10, by=1, value=2,
                                          container=f1g1)
  
  f1g1[1,4] <- glabel(text="Vertical justification:", container=f1g1)
  f1g1[1,5] <- f1_vjust_spb <- gspinbutton(from=0, to=1, by=0.5, value=1,
                                           container=f1g1)
  
  f1g1[2,2] <- glabel(text="Allele label angle:", container=f1g1)
  f1g1[2,3] <- f1_angle_spb <- gspinbutton(from=0, to=360, by=15, value=0,
                                           container=f1g1)
  
  f1g1[2,4] <- glabel(text="Horizontal justification:", container=f1g1)
  f1g1[2,5] <- f1_hjust_spb <- gspinbutton(from=0, to=1, by=0.5, value=0.5,
                                           container=f1g1)
  
  f1g1[3,2] <- glabel(text="Plot area expansion:", container=f1g1)
  f1g1[3,3] <- f1_expand_spb <- gspinbutton(from=0, to=1, by=0.05, value=0.10,
                                            container=f1g1)
  
  f1g1[3,4] <- glabel(text="Analytical threshold:", container=f1g1)
  f1g1[3,5] <- f1_at_spb <- gspinbutton(from=0, to=1000, by=10, value=0,
                                        container=f1g1)
  
  f1_ignore_chk <- gcheckbox(text="Ignore case in marker names",
                             checked=TRUE, container=f1)
  
  f1_wrap_chk <- gcheckbox(text="Wrap by dye and add marker ranges and allele names",
                           checked=TRUE, container=f1)
  
  f1_fix_chk <- gcheckbox(text="Fix x-axis to size range",
                          checked=TRUE, container=f1)
  
  f1_collapse_chk <- gcheckbox(text="Collapse (add peak heights of identical alleles. Discards OL)",
                               checked=TRUE, container=f1)
  
  f1_box_chk <- gcheckbox(text="Plot peak height distribution (boxplot)",
                          checked=FALSE, container=f1)
  
  f1_peaks_chk <- gcheckbox(text="Plot mean peak height for distributions",
                            checked=TRUE, container=f1)
  
  
  addHandlerChanged(f1_collapse_chk, handler = function(h, ...) {
    
    val_collapse <- svalue(f1_collapse_chk)
    
    if(val_collapse){
      
      enabled(f1_box_chk) <- TRUE
      enabled(f1_peaks_chk) <- TRUE
      
    } else {
      
      enabled(f1_box_chk) <- FALSE
      enabled(f1_peaks_chk) <- FALSE
      
    }
    
  } )
  
  #// END OF COPIED CODE.
  
  # SIMULATION PARAMETERS #####################################################
  
  glabel("", container=sim_gf) # Adds some space.
  
  # SIMULATION ----------------------------------------------------------------
  
  # sim_f1_txt_width <- 8 
  
  sim_f1 <- gframe(text = "Options", markup = FALSE, pos = 0, horizontal=FALSE,
                   container = sim_gf, expand=FALSE) 
  
  sim_f1g1 <- glayout(container = sim_f1)
  
  sim_f1g1[1,1] <- glabel("", container=sim_f1g1) # Adds some space.
  
  
  sim_f1g1[2,1] <- glabel(text="Number of simulations or contributors:",
                          container=sim_f1g1, anchor=c(-1 ,0))
  
  sim_f1g1[2,2] <- sim_sim_edt <- gedit(text="1", width=10,
                                        container=sim_f1g1)
  
  sim_f1g1[3,1] <- glabel(text="Name for result:",
                          container=sim_f1g1, anchor=c(-1 ,0))
  
  sim_f1g1[3,2] <- sim_name_edt <- gedit(text="sim", width=25,
                                         container=sim_f1g1)
  
  sim_f1g1[4,1] <- glabel(text="Plot title:", 
                          container=sim_f1g1, anchor=c(-1 ,0))
  
  sim_f1g1[4,2] <- sim_epg_title_edt <- gedit(text="Simulated data",
                                              width=25, container=sim_f1g1)
  
  
  sim_f1g1[5,1] <- sim_rnd_seed_chk <- gcheckbox(text = "Set random seed:",
                                                 checked = FALSE,
                                                 container = sim_f1g1) 
  sim_f1g1[5,2] <- sim_rnd_seed_edt <- gedit(text="", width=6,
                                             container=sim_f1g1)

  sim_update_chk <- gcheckbox(text="Auto update EPG",
                              checked=FALSE, container = sim_f1) 
  
  sim_save_chk <- gcheckbox(text="Save in workspace, create name from:",
                            checked=FALSE, container = sim_f1) 
  
  sim_save_s_chk <- gcheckbox(text="Sample Name",
                              checked=FALSE, container = sim_f1) 
  sim_save_t_chk <- gcheckbox(text="Time stamp",
                              checked=FALSE, container = sim_f1) 
  sim_save_n_chk <- gcheckbox(text="Result name",
                              checked=FALSE, container = sim_f1) 
  
  addHandlerChanged(sim_rnd_seed_chk, handler = function(h, ...) {
    
    .refreshSimulationTab()
    
  } )
  addHandlerChanged(sim_save_chk, handler = function(h, ...) {
    
    .refreshSimulationTab()
    
  } )

  sim_f1g1[5,1] <- glabel("", container=sim_f1g1) # Adds some space.
  
  # Simulation options --------------------------------------------------------
  
  sim_f2 <- gframe(text = "Simulate", markup = FALSE, pos = 0, horizontal=FALSE,
                   container = sim_gf, expand=FALSE) 
  
  sim_profile_chk <- gcheckbox(text="Profile",
                                  checked = FALSE, container = sim_f2)
  
  sim_sample_chk <- gcheckbox(text="Sample",
                                 checked = FALSE, container = sim_f2)
  
  sim_deg_chk <- gcheckbox(text="Degradation",
                                      checked = FALSE, container = sim_f2)
  
  sim_extr_chk <- gcheckbox(text="Extraction",
                                     checked = FALSE, container = sim_f2)
  
  sim_norm_chk <- gcheckbox(text="Normalization",
                                   checked = FALSE, container = sim_f2)
  
  sim_pcr_chk <- gcheckbox(text="PCR Amplification",
                              checked = FALSE, container = sim_f2)
  
  sim_ce_chk <- gcheckbox(text="Capillary Electrophoresis",
                             checked = FALSE, container = sim_f2)
  
  addHandlerChanged(sim_profile_chk, handler = function(h, ...) {
    
    .refreshSimulationTab()
    
  } )

  addHandlerChanged(sim_sample_chk, handler = function(h, ...) {
    
    .refreshSimulationTab()
    
  } )

  addHandlerChanged(sim_deg_chk, handler = function(h, ...) {
    
    .refreshSimulationTab()
    
  } )

  addHandlerChanged(sim_extr_chk, handler = function(h, ...) {
    
    .refreshSimulationTab()
    
  } )
  
  addHandlerChanged(sim_norm_chk, handler = function(h, ...) {
    
    .refreshSimulationTab()
    
  } )
  
  addHandlerChanged(sim_pcr_chk, handler = function(h, ...) {
    
    .refreshSimulationTab()
    
  } )
  
  addHandlerChanged(sim_ce_chk, handler = function(h, ...) {
    
    .refreshSimulationTab()
    
  } )
  
  # Button ------------------------------------------------------------------  
  
  sim_sim_btn <- gbutton(text = "Simulate", border=TRUE, container=sim_gf)
  
  sim_view_btn <- gbutton(text="View result of the DNA process",
                          border=TRUE, container=sim_gf)
  
  sim_view_ph_btn <- gbutton(text="View result after capillary electrophoresis",
                             border=TRUE, container=sim_gf)
  
  sim_epg_btn <- gbutton(text = "Generate EPG",
                         border=TRUE, container=sim_gf) 
  
  
  addHandlerChanged(sim_sim_btn, handler = function(h, ...) {
    
    val_ok <- TRUE
    
    # Disable button until simulation is finished.
    enabled(sim_sim_btn) <- FALSE
    
    # LOAD ARGUMENTS ----------------------------------------------------------
    
    message("****** LOADING PARAMETERS FROM GUI")
    
    # File tab.
    val_debug <- svalue(file_f1g1_debug_chk)
    valExtDebug <- svalue(file_f1g1_ext_debug_chk)
    
    if(val_debug){
      print("** TAB FILE")
      print("val_debug:")
      print(val_debug)
      print("valExtDebug:")
      print(valExtDebug)
    }
    
    # Profile tab.
    val_kit <- svalue(profile_kit_drp)
    val_method <- svalue(profile_method_drp)
    val_profile_opt <- svalue(profile_opt, index=TRUE)
    val_db <- svalue(profile_pop_drp)
    val_profile <- profile_tbl[,]
    
    if(val_debug){
      print("** TAB PROFILE")
      print("val_kit:")
      print(val_kit)
      print("val_method:")
      print(val_method)
      print("val_profile_opt:")
      print(val_profile_opt)
      print("val_db:")
      print(val_db)
      print("val_profile:")
      print(head(val_profile))
    }
    
    # Sample tab.
    val_smpl_name <- svalue(sample_name_edt)
    val_smpl_mix <- svalue(sample_mix_chk)
    val_smpl_opt <- svalue(sample_by_opt, index=TRUE)
    val_smpl_vol <- svalue(sample_vol_edt)
    val_smpl_vol_sd <- svalue(sample_vol_sd_edt)
    val_smpl_conc <- svalue(sample_conc_edt)
    val_smpl_conc_sd <- svalue(sample_conc_sd_edt)
    val_smpl_ncells <- svalue(sample_ncells_edt)
    val_smpl_ncells_sd <- svalue(sample_ncells_sd_edt)
    val_smpl_celldna <- as.numeric(svalue(sample_celldna_edt))
    val_smpl_haploid <- svalue(sample_type_cmb)
    val_smpl_slope <- as.numeric(svalue(sample_slope_edt))
    val_smpl_intercept <- as.numeric(svalue(sample_intercept_edt))

    # Check and create vectors.
    val_smpl_vol <- .stringToVector(data=val_smpl_vol)
    val_smpl_vol_sd <- .stringToVector(data=val_smpl_vol_sd)
    val_smpl_conc <- .stringToVector(data=val_smpl_conc)
    val_smpl_conc_sd <- .stringToVector(data=val_smpl_conc_sd)
    val_smpl_ncells <- .stringToVector(data=val_smpl_ncells)
    val_smpl_ncells_sd <- .stringToVector(data=val_smpl_ncells_sd)
    val_smpl_vol <- .stringToVector(data=val_smpl_vol)
    val_smpl_vol_sd <- .stringToVector(data=val_smpl_vol_sd)
    val_smpl_haploid <- .stringToVector(data=val_smpl_haploid, type="logical")
    
    
    if(val_debug){
      print("** TAB SAMPLE")
      print("val_smpl_name:")
      print(val_smpl_name)
      print("val_smpl_opt:")
      print(val_smpl_opt)
      print("val_smpl_vol:")
      print(val_smpl_vol)
      print("val_smpl_vol_sd:")
      print(val_smpl_vol_sd)
      print("val_smpl_conc:")
      print(val_smpl_conc)
      print("val_smpl_conc_sd:")
      print(val_smpl_conc_sd)
      print("val_smpl_ncells:")
      print(val_smpl_ncells)
      print("val_smpl_ncells_sd:")
      print(val_smpl_ncells_sd)
      print("val_smpl_celldna:")
      print(val_smpl_celldna)
      print("val_smpl_haploid:")
      print(val_smpl_haploid)
      print("val_smpl_slope:")
      print(val_smpl_slope)
      print("val_smpl_intercept:")
      print(val_smpl_intercept)
    }
    
    # Extraction tab.
    val_exprob <- as.numeric(svalue(ex_eff_spn))
    val_exprob_sd <- as.numeric(svalue(ex_eff_sd_edt))
    val_exvol <- as.numeric(svalue(ex_vol_edt))
    val_exvol_sd <- as.numeric(svalue(ex_vol_sd_edt))
    
    val_title <- svalue(sim_epg_title_edt)
    
    if(val_debug){
      print("** TAB EXTRACTION")
      print("val_exprob:")
      print(val_exprob)
      print("val_exprob_sd:")
      print(val_exprob_sd)
      print("val_exvol:")
      print(val_exvol)
      print("val_exvol_sd:")
      print(val_exvol_sd)
    }
    
    # Dilution tab.
    val_dil_amount_t <- as.numeric(svalue(dil_tolerance_edt))
    val_dil_target_opt <- svalue(dil_target_opt, index=TRUE)
    val_dil_target_edt <- as.numeric(svalue(dil_target_edt))
    val_dil_target_lbl <- as.numeric(svalue(dil_target_lbl))
    val_dil_accuracy <- as.numeric(svalue(dil_accuracy_edt))
    val_dil_volume_opt <- svalue(dil_volume_opt, index=TRUE)
    val_dil_volume_edt <- as.numeric(svalue(dil_volume_edt))
    val_dil_volume_lbl <- as.numeric(svalue(dil_volume_lbl))
    
    if(val_debug){
      print("** TAB DILUTION")
      print("val_dil_amount_t:")
      print(val_dil_amount_t)
      print("val_dil_target_opt:")
      print(val_dil_target_opt)
      print("val_dil_target_edt:")
      print(val_dil_target_edt)
      print("val_dil_target_lbl:")
      print(val_dil_target_lbl)
      print("val_dil_accuracy:")
      print(val_dil_accuracy)
      print("val_dil_volume_opt:")
      print(val_dil_volume_opt)
      print("val_dil_volume_edt:")
      print(val_dil_volume_edt)
      print("val_dil_volume_lbl:")
      print(val_dil_volume_lbl)
    }
    
    # Degradation tab.
    val_deg_deg <- as.numeric(svalue(deg_pam_lbl))
    val_deg_size <- as.numeric(svalue(deg_target_edt))
    
    if(val_debug){
      print("** TAB DEGRADATION")
      print("val_deg_deg:")
      print(val_deg_deg)
      print("val_deg_size:")
      print(val_deg_size)
    }
    
    # Amplification tab.
    val_amp_alq <- as.numeric(svalue(amp_aliquot_edt))
    val_amp_alq_sd <- as.numeric(svalue(amp_aliquot_sd_edt))
    val_amp_amount <- as.numeric(svalue(amp_amount_edt))
    val_amp_tvol <- as.numeric(svalue(amp_total_vol_edt))
    val_amp_tvol_sd <- as.numeric(svalue(amp_total_vol_sd_edt))
    val_amp_cyc <- as.numeric(svalue(cyc_sb))
    val_amp_opt <- svalue(amp_f3_eff_opt, index=TRUE)
    val_amp_eff <- svalue(amp_f3_eff_spn)
    val_amp_eff_sd <- as.numeric(svalue(amp_f3_eff_sd_edt))
    val_amp_eff_stutt <- svalue(amp_eff_stutt_spn)
    val_amp_eff_stutt_sd <- as.numeric(svalue(amp_eff_stutt_sd_edt))
    val_amp_stutter <- svalue(amp_stutter_chk)
    
    if(val_debug){
      print("** TAB PCR")
      print("val_amp_alq:")
      print(val_amp_alq)
      print("val_amp_alq_sd:")
      print(val_amp_alq_sd)
      print("val_amp_amount:")
      print(val_amp_amount)
      print("val_amp_tvol:")
      print(val_amp_tvol)
      print("val_amp_tvol_sd:")
      print(val_amp_tvol_sd)
      print("val_amp_cyc:")
      print(val_amp_cyc)
      print("val_amp_opt:")
      print(val_amp_opt)
      print("val_amp_eff:")
      print(val_amp_eff)
      print("val_amp_eff_sd:")
      print(val_amp_eff_sd)
      print("val_amp_eff_stutt:")
      print(val_amp_eff_stutt)
      print("val_amp_eff_stutt_sd:")
      print(val_amp_eff_stutt_sd)
      print("val_amp_stutter:")
      print(val_amp_stutter)
    }
    
    # Capillary electrophoresis tab.
    val_ce_aliq <- as.numeric(svalue(ce_aliquot_edt))
    val_ce_aliq_sd <- as.numeric(svalue(ce_aliquot_sd_edt))
    val_ce_opt <- svalue(ce_settings_opt, index = TRUE)
    val_ce_intercept <- as.numeric(svalue(ce_intercept_edt))
    val_ce_slope <- as.numeric(svalue(ce_slope_edt))
    val_ce_residual <- as.numeric(svalue(ce_residual_edt))
    val_ce_t_intercept <- as.numeric(svalue(ce_t_intercept_edt))
    val_ce_t_slope <- as.numeric(svalue(ce_t_slope_edt))
    val_ce_t_residual <- as.numeric(svalue(ce_t_residual_edt))
    
    if(val_debug){
      print("** TAB CE")
      print("val_ce_aliq:")
      print(val_ce_aliq)
      print("val_ce_aliq_sd:")
      print(val_ce_aliq_sd)
      print("val_ce_opt:")
      print(val_ce_opt)
      print("val_ce_intercept:")
      print(val_ce_intercept)
      print("val_ce_slope:")
      print(val_ce_slope)
      print("val_ce_residual:")
      print(val_ce_residual)
      print("val_ce_t_intercept:")
      print(val_ce_t_intercept)
      print("val_ce_t_slope:")
      print(val_ce_t_slope)
      print("val_ce_t_residual:")
      print(val_ce_t_residual)
    }

    # Electropherogram tab.
    val_epg_wrap <- svalue(f1_wrap_chk)
    val_epg_box <- svalue(f1_box_chk)
    val_epg_peaks <- svalue(f1_peaks_chk)
    val_epg_collapse <- svalue(f1_collapse_chk)
    val_epg_ignore <- svalue(f1_ignore_chk)
    val_epg_at <- as.numeric(svalue(f1_at_spb))
    val_epg_scale <- svalue(f1_scale_opt)
    val_epg_fix <- svalue(f1_fix_chk)
    val_epg_size <- as.numeric(svalue(f1_size_spb))
    val_epg_angle <- as.numeric(svalue(f1_angle_spb))
    val_epg_vjust <- as.numeric(svalue(f1_vjust_spb))
    val_epg_hjust <- as.numeric(svalue(f1_hjust_spb))
    val_epg_expand <- as.numeric(svalue(f1_expand_spb))
    
    if(val_debug){
      print("** TAB EPG")
      print("val_epg_wrap:")
      print(val_epg_wrap)
      print("val_epg_box:")
      print(val_epg_box)
      print("val_epg_peaks:")
      print(val_epg_peaks)
      print("val_epg_collapse:")
      print(val_epg_collapse)
      print("val_epg_ignore:")
      print(val_epg_ignore)
      print("val_epg_at:")
      print(val_epg_at)
      print("val_epg_scale:")
      print(val_epg_scale)
      print("val_epg_fix:")
      print(val_epg_fix)
      print("val_epg_size:")
      print(val_epg_size)
      print("val_epg_angle:")
      print(val_epg_angle)
      print("val_epg_vjust:")
      print(val_epg_vjust)
      print("val_epg_hjust:")
      print(val_epg_hjust)
      print("val_epg_expand:")
      print(val_epg_expand)
    }
    
    # Simulation tab.
    val_simulations <- as.numeric(svalue(sim_sim_edt))
    val_sim_name <- svalue(sim_name_edt)
    val_sim_rnd_chk <- svalue(sim_rnd_seed_chk)
    val_sim_rnd_seed <- as.numeric(svalue(sim_rnd_seed_edt))
    val_sim_save <- svalue(sim_save_chk)
    val_sim_save_s <- svalue(sim_save_s_chk)
    val_sim_save_t <- svalue(sim_save_t_chk)
    val_sim_save_n <- svalue(sim_save_n_chk)
    val_sim_profile <- svalue(sim_profile_chk)
    val_sim_sample <- svalue(sim_sample_chk)
    val_sim_sample_enabled <- enabled(sim_sample_chk)
    val_sim_extraction <- svalue(sim_extr_chk)
    val_sim_extraction_enabled <- enabled(sim_extr_chk)
    val_sim_degradation <- svalue(sim_deg_chk)
    val_sim_degradation_enabled <- enabled(sim_deg_chk)
    val_sim_dilution <- svalue(sim_norm_chk)
    val_sim_dilution_enabled <- enabled(sim_norm_chk)
    val_sim_pcr <- svalue(sim_pcr_chk)
    val_sim_pcr_enabled <- enabled(sim_pcr_chk)
    val_sim_ce <- svalue(sim_ce_chk)
    val_sim_ce_enabled <- enabled(sim_ce_chk)
    
    if(val_debug){
      print("** TAB SIMULATION")
      print("val_simulations:")
      print(val_simulations)
      print("val_sim_name:")
      print(val_sim_name)
      print("val_sim_save:")
      print(val_sim_save)
      print("val_sim_save_s:")
      print(val_sim_save_s)
      print("val_sim_save_t:")
      print(val_sim_save_t)
      print("val_sim_save_n:")
      print(val_sim_save_n)
      print("val_sim_profile:")
      print(val_sim_profile)
      print("val_sim_sample:")
      print(val_sim_sample)
      print("val_sim_extraction:")
      print(val_sim_extraction)
      print("val_sim_degradation:")
      print(val_sim_degradation)
      print("val_sim_dilution:")
      print(val_sim_dilution)
      print("val_sim_pcr:")
      print(val_sim_pcr)
      print("val_sim_ce:")
      print(val_sim_ce)
    }

    # CHECK SIM SETTINGS ------------------------------------------------------
    
    # Sample is required for downstream simulations.
    if(!val_sim_sample_enabled){
      val_sim_degradation <- FALSE
      val_sim_extraction <- FALSE
      val_sim_dilution <- FALSE
      val_sim_pcr <- FALSE
      val_sim_ce <- FALSE
    }

    # Extraction is required for downstream simulations.
    if(!val_sim_extraction_enabled){
      val_sim_dilution <- FALSE
      val_sim_pcr <- FALSE
      val_sim_ce <- FALSE
    }

    # PCR is required for downstream simulations.
    if(!val_sim_pcr_enabled){
      val_sim_ce <- FALSE
    }
    
    # PREPARE -----------------------------------------------------------------

    # Shared variables must always be prepared! ...............................

    if(val_debug){
      print("Prepare essential parameters.")
    }

    # Set random seed.
    if(val_sim_rnd_chk){
      set.seed(val_sim_rnd_seed)
    }
    
    # Reset simulation data.
    simData <- NULL
    
    # Convert from pico grams to nano grams.
    valCellDNA <- val_smpl_celldna / 1000
    
    # Get sample name.
    valSampleName <- val_smpl_name
    
    # Get mixture simulation flag.
    valMixture <- val_smpl_mix
    

    # Sample ..................................................................

    if(val_sim_sample){
      # Simulate sample.

      if(val_debug){
        print("Prepare sample parameters.")
      }
      
      # Estimate the amount of DNA by..
      if(val_smpl_opt == 1){
        # Use number of cells.
        val_smpl_conc <- NULL # Set concentration to null.
        val_smpl_slope <- NULL # Set slope to null.
        val_smpl_intercept <- NULL # Set intercept to null.
        if(val_debug){
          print("Parameter updated.")
          print("val_smpl_conc:")
          print(val_smpl_conc)
          print("val_smpl_slope:")
          print(val_smpl_slope)
          print("val_smpl_intercept:")
          print(val_smpl_intercept)
        }
      } else if(val_smpl_opt == 2){
        # Use concentration.
        val_smpl_ncells <- NULL # Set number of cells to null.
        val_smpl_slope <- NULL # Set slope to null.
        val_smpl_intercept <- NULL # Set intercept to null.
        if(val_debug){
          print("Parameter updated.")
          print("val_smpl_ncells:")
          print(val_smpl_ncells)
          print("val_smpl_slope:")
          print(val_smpl_slope)
          print("val_smpl_intercept:")
          print(val_smpl_intercept)
        }
      } else if(val_smpl_opt == 3){
        # Use slope and intercept.
        val_smpl_conc <- NULL # Set concentration to null.
        val_smpl_ncells <- NULL # Set number of cells to null.
        if(val_debug){
          print("Parameter updated.")
          print("val_smpl_ncells:")
          print(val_smpl_ncells)
          print("val_smpl_conc:")
          print(val_smpl_conc)
        }
      } else {
        warning(paste("val_smpl_opt =", val_smpl_opt, "not supported!"))
      }
      
    }

    # Save ....................................................................
    
    # Check if save result and create name for result.
    if(val_sim_save){
      
      if(val_debug){
        print("Prepare save parameters.")
      }
      
      # Default name if nothing.
      val_res_name <- "sim"
      
      # Add simulation name.
      if(val_sim_save_n){
        if(val_sim_name != ""){
          val_res_name <- val_sim_name
        }
      }
      
      # Add sample name.
      if(val_sim_save_s){
        val_res_name <- paste(val_res_name, valSampleName, sep="_")
      }
      
      # Add time-stamp.
      if(val_sim_save_t){
        val_res_name <- paste(val_res_name, format(Sys.time(), "%d.%m.%Y_%H.%M.%S"), sep="_")
      }
      
      if(val_debug){
        print("Base name for result created.")
        print("val_res_name:")
        print(val_res_name)
      }
      
    }
    
    # Profile .................................................................

    if(val_sim_profile){
      # Simulate profile.

      if(val_debug){
        print("Prepare profile parameters.")
      }
      
      # Check if fat data.
      if(length(grep("Allele", names(val_profile), fixed=TRUE)) > 1){
        # Slim data.
        val_profile <- strvalidator::slim(data=val_profile,
                                          fix = "Marker", stack = "Allele",
                                          keep.na = TRUE, debug = valExtDebug)
        if(val_debug){
          print("Fat data slimmed!")
        }
      }
      
      # Replace character "NA" with NA.
      val_profile$Allele[val_profile$Allele == "NA"] <- NA
      
      # Store profile in simData.
      simData <- val_profile
      
      # Get profile name.
      profileName <- valSampleName
      # Check arguments.
      if(profileName == ""){
        # If no sample name 'name' must be NULL.
        profileName <- NULL
      }
      
      if(val_profile_opt == 5){
        # If options 5 then 'data' must be NULL.
        simData <- NULL
        # Load database file.
        profile_db <- strvalidator::getDb(db.name.or.index=val_db, debug=valExtDebug)
      } else {
        # If options 1-4 then 'db' must be NULL.
        profile_db <- NULL
      }
      
    
    }
    
    # Dilution ................................................................

    if(val_sim_dilution){
      
      if(val_debug){
        print("Prepare dilution parameters.")
      }
      
      # Initiate additional variables.
      val_dil_conc <- NULL # Target concentration.
      val_dil_vol <- NULL # Final dilution volume.
      
      # Get target concentration.
      if(val_dil_target_opt == 1){
        val_dil_conc <- val_dil_target_lbl
      } else if (val_dil_target_opt == 2){
        val_dil_conc <- val_dil_target_edt 
      } else {
        stop(paste("val_dil_target_opt=", val_dil_target_opt, "not supported!"))
      }
      if(val_debug){
        print("val_dil_conc:")
        print(val_dil_conc)
      }
      
      # Get final dilution volume.
      if(val_dil_volume_opt == 1){
        val_dil_vol <- val_dil_volume_lbl
      } else if (val_dil_volume_opt == 2){
        val_dil_vol <- val_dil_volume_edt 
      } else {
        stop(paste("val_dil_volume_opt=", val_dil_volume_opt, "not supported!"))
      }
      if(val_debug){
        print("val_dil_vol:")
        print(val_dil_vol)
      }
      
    }
    
    # SIMULATE PROFILE --------------------------------------------------------
    
    if(val_sim_profile){
      # Simulate profile.
      
      if(val_debug){
        print("****** SIMULATE PROFILE")
      }
      
      # Simulate profile.
      simData <- simProfile(data=simData, kit=NULL, sim=val_simulations,
                            name=profileName, db=profile_db,
                            debug=valExtDebug)
      
    }
    
    # SIMULATE SAMPLE ---------------------------------------------------------
    
    if(val_sim_sample){
      
      if(val_debug){
        print("****** SIMULATE SAMPLE")
        print("simData")
        print(head(simData))
      }
      
      # Simulate sample.
      simData <- simSample(data=simData,
                           cells=val_smpl_ncells, sd.cells=val_smpl_ncells_sd,
                           conc=val_smpl_conc, sd.conc=val_smpl_conc_sd,
                           vol=val_smpl_vol, sd.vol=val_smpl_vol_sd,
                           cell.dna=valCellDNA, haploid=val_smpl_haploid,
                           kit=val_kit, slope=val_smpl_slope, intercept=val_smpl_intercept,
                           debug=valExtDebug)
      
    }

    # SIMULATE DEGRADATION ----------------------------------------------------
    
    if(val_sim_degradation){
      
      if(val_debug){
        print("****** SIMULATE DEGRADATION")
      }
      
      # Simulate degradation.
      simData <- simDegradation(data=simData, kit=val_kit, deg=val_deg_deg,
                                quant.target=val_deg_size, debug=valExtDebug)
      
    }
    
    # SIMULATE EXTRACTION -----------------------------------------------------
    
    if(val_sim_extraction){
      
      if(val_debug){
        print("****** SIMULATE EXTRACTION")
      }
      
      # Simulate extraction.
      simData <- simExtraction(data=simData, 
                               vol.ex=val_exvol, sd.vol=val_exvol_sd,
                               prob.ex=val_exprob, sd.prob=val_exprob_sd,
                               cell.dna=valCellDNA, debug=valExtDebug)
      
    }
    
    # SIMULATE DILUTION -------------------------------------------------------
    
    if(val_sim_dilution){
      
      if(val_debug){
        print("****** SIMULATE NORMALIZATION")
      }
      
      # Simulate dilution.
      simData <- simNormalize(data=simData, volume=val_dil_vol,
                             accuracy=val_dil_accuracy, target=val_dil_conc,
                             tolerance=val_dil_amount_t, multiple=FALSE,
                             debug=valExtDebug)
      
    }
    
    # SIMULATE AMPLIFICATION --------------------------------------------------
    
    if(val_sim_pcr){
      
      if(val_debug){
        print("****** SIMULATE PCR")
      }
      
      # Arguments.
      pcr_kit <- val_kit
      pcr_method <- val_method
      
      # Check arguments:
      if(val_amp_opt == 1){
        
        # If automatic, method and kit must be provided (i.e. no change).
        
      } else if(val_amp_opt == 2){
        
        # If custom pcr efficiency 'kit' must be NULL and method is not used.
        pcr_kit <- NULL
        pcr_method <- "DEFAULT"
        
      } else {
        
        stop("PCR efficiency option = ", val_amp_opt, " not implemented!")
        
      }
      
      # Simulate PCR amplification.
      simData <- simPCR(data=simData, kit=pcr_kit, method=pcr_method, pcr.cyc=val_amp_cyc,
                        pcr.prob=val_amp_eff, sd.pcr.prob=val_amp_eff_sd,
                        stutter.prob=val_amp_eff_stutt, sd.stutter.prob=val_amp_eff_stutt_sd,
                        vol.aliq=val_amp_alq, sd.vol.aliq=val_amp_alq_sd,
                        vol.pcr=val_amp_tvol, sd.vol.pcr=val_amp_tvol_sd,
                        stutter=val_amp_stutter, debug=valExtDebug)
      
    }
    
    # ATTACH ATTRIBUTES -------------------------------------------------------

    simData <- .addAttributes(simData)
    
    # SAVE FIRST RESULT -------------------------------------------------------
    
    # Save result in project.
    if(val_sim_save){
      assign(x=val_res_name, value=simData, envir=.pcrsim_env)
      message(paste("Simulation result saved as", val_res_name))
    }
    
    # Save current simulation in varaiable.
    .simDataRes <<- simData

    # COMPACT DATA ------------------------------------------------------------

    # Compact stutters if simulated.
    if(any(grepl("PCR.Stutter.", names(simData), fixed = TRUE))){
      
      simData <- compactStutter(data=simData, stutterin="PCR.Stutter.",
                                   targetcol="PCR.Amplicon")
      
    }

    # Always compact alleles.
    if(any(grepl("PCR.Amplicon", names(simData)))){
      
      simCompact <- compact(data=simData, per.sample=!valMixture,
                            col="PCR.Amplicon", sim=TRUE)
      
    }

    # SIMULATE CAPILLARY ELECTROPHORESIS --------------------------------------
    
    if(val_sim_ce){
      
      if(val_debug){
        print("****** SIMULATE CE")
      }

      # Simulate capillary electrophoresis.
      simPh <- simCE(data=simCompact, vol=val_ce_aliq, sd.vol=val_ce_aliq_sd,
                       intercept=val_ce_intercept, slope=val_ce_slope, sigma=val_ce_residual,
                       t.intercept=val_ce_t_intercept, t.slope=val_ce_t_slope, t.sigma=val_ce_t_residual,
                       debug=valExtDebug)

      # ATTACH ATTRIBUTES -------------------------------------------------------
      
      simPh <- .addAttributes(simPh)

      # SAVE SECOND RESULT ------------------------------------------------------
      
      # Save result in project.
      if(val_sim_save){
        tmp_name <- paste(val_res_name,"_ce", sep="")
        assign(x=tmp_name, value=simPh, envir=.pcrsim_env)
        message(paste("Capillary electrophoresis result saved as", tmp_name))
      }

      # Save current simulation in varaiable.
      .simDataPh <<- simPh

    }
    
    # PLOT EPG ----------------------------------------------------------------
    
    if(svalue(sim_update_chk) && !is.null(.simDataPh)){
      # Generate EPG.      
      .simEPG <<- strvalidator::generateEPG(data=.simDataPh, 
                                            kit=val_kit, 
                                            title=val_title,
                                            wrap=val_epg_wrap,
                                            boxplot=val_epg_box,
                                            peaks=val_epg_peaks,
                                            collapse=val_epg_collapse,
                                            silent=FALSE,
                                            ignore.case=val_epg_ignore,
                                            at=val_epg_at,
                                            scale=val_epg_scale,
                                            limit.x=val_epg_fix,
                                            label.size=val_epg_size,
                                            label.angle=val_epg_angle,
                                            label.vjust=val_epg_vjust,
                                            label.hjust=val_epg_hjust,
                                            expand=val_epg_expand,
                                            debug=valExtDebug)
      
      # Save result in project.
      if(val_sim_save){
        tmp_name <- paste(val_res_name,"_epg", sep="")
        assign(x=tmp_name, value=.simEPG, envir=.pcrsim_env)
        message(paste("Electropherogram saved as", tmp_name))
      }
      
    }
    
    # Enable buttons.
    enabled(sim_sim_btn) <- TRUE
    
  } )
  
  addHandlerChanged(sim_view_btn, handler = function(h, ...) {
    
    if(!is.null(.simDataRes)){
      
      strvalidator::editData_gui(data = .simDataRes,
                                 edit = FALSE, env = .pcrsim_env)

    } else {
      
      gmessage(message="Currently no data!", title="Error",
               icon = "error", parent = w) 
      
    }

  } )
  
  addHandlerChanged(sim_view_ph_btn, handler = function(h, ...) {

    if(!is.null(.simDataPh)){

      strvalidator::editData_gui(data = .simDataPh,
                                 edit = FALSE, env = .pcrsim_env)
      
  } else {
    
    gmessage(message="Currently no data!", title="Error",
             icon = "error", parent = w) 
    
  }
  
    
  } )
  
  addHandlerChanged(sim_epg_btn, handler = function(h, ...) {

    if(!is.null(.simDataPh)){
      
      # Load all settings.
      val_title <- svalue(sim_epg_title_edt)
      val_kit <- svalue(profile_kit_drp)
      # Electropherogram tab.
      val_epg_wrap <- svalue(f1_wrap_chk)
      val_epg_box <- svalue(f1_box_chk)
      val_epg_peaks <- svalue(f1_peaks_chk)
      val_epg_collapse <- svalue(f1_collapse_chk)
      val_epg_ignore <- svalue(f1_ignore_chk)
      val_epg_at <- as.numeric(svalue(f1_at_spb))
      val_epg_scale <- svalue(f1_scale_opt)
      val_epg_fix <- svalue(f1_fix_chk)
      val_epg_size <- as.numeric(svalue(f1_size_spb))
      val_epg_angle <- as.numeric(svalue(f1_angle_spb))
      val_epg_vjust <- as.numeric(svalue(f1_vjust_spb))
      val_epg_hjust <- as.numeric(svalue(f1_hjust_spb))
      val_epg_expand <- as.numeric(svalue(f1_expand_spb))
      valExtDebug <- svalue(file_f1g1_ext_debug_chk)
      
      # Generate EPG.      
      .simEPG <<- strvalidator::generateEPG(data=.simDataPh, 
                                            kit=val_kit, 
                                            title=val_title,
                                            wrap=val_epg_wrap,
                                            boxplot=val_epg_box,
                                            peaks=val_epg_peaks,
                                            collapse=val_epg_collapse,
                                            silent=FALSE,
                                            ignore.case=val_epg_ignore,
                                            at=val_epg_at,
                                            scale=val_epg_scale,
                                            limit.x=val_epg_fix,
                                            label.size=val_epg_size,
                                            label.angle=val_epg_angle,
                                            label.vjust=val_epg_vjust,
                                            label.hjust=val_epg_hjust,
                                            expand=val_epg_expand,
                                            debug=valExtDebug)
      
    
    } else {
      
      gmessage(message="Currently no data!", title="Error", icon = "error", parent = w) 
      
    }
    
  } )
  
  # MAIN EVENT HANDLERS #########################################################
  
  addHandlerChanged(nb, handler = function (h, ...) {
    
    if(debug){
      print("NOTEBOOK CHANGED")
      print(if(is.null(h$pageno)) svalue(h$obj) else h$pageno)
    }
    
    # Refresh depending on active tab.
    #tab <- svalue(nb)
    tab <- if(is.null(h$pageno)) svalue(h$obj) else h$pageno
    tabName <- names(nb)[tab]
    
    if(tabName == .file_tab_name){
      
      .refreshLoaded()
      .refreshWs()
      
    }
    
    if(tabName == .sample_tab_name){
      
      .refreshSampleTab()
      
    }
    
    if(tabName == .degradation_tab_name){
      
      .refreshDegradationTab()
      
    }
    
    if(tabName == .dilution_tab_name){
      
      .refreshDilutionTab()
      
    }
    
    if(tabName == .amplification_tab_name){
      
      .refreshAmplificationTab()
      
    }
    
    if(tabName == .simulation_tab_name){
      
      .refreshSimulationTab()
      
    }
    
  })
  
  addHandlerFocus(w, handler = function (h, ...) {
    
    if(debug){
      print(paste("IN:", match.call()[[1]]))
      print("FOCUS")
    }
    
    # Refresh depending on active tab.
    tab <- svalue(nb)
    tabName <- names(nb)[tab]
    
    if(tabName == .file_tab_name){
      
      .refreshLoaded()
      .refreshWs()
      
    }
    
  })

  # MAIN EVENT HANDLERS #########################################################
  addHandlerChanged(nb, handler = function (h, ...) {
    
    if(debug){
      print("NOTEBOOK CHANGED")
      print(if(is.null(h$pageno)) svalue(h$obj) else h$pageno)
    }
    
    # Refresh depending on active tab.
    tab <- if(is.null(h$pageno)) svalue(h$obj) else h$pageno
    tabName <- names(nb)[tab]
    
    # Check if a tab name exist and then perform tasks.
    if(length(tabName) != 0){
      
      if(tabName == .profile_tab_name){
        
        .refreshProfileTab()
        
      } else if(tabName == .analysis_tab_name){
        
        .refreshCETab()
        
      } else if(tabName == .amplification_tab_name){
        
        .refreshAmplificationTab()
        
      } else if(tabName == .simulation_tab_name){
        
        .refreshSimulationTab()
        
      }
      
    } # End check.
    
  })
  
  addHandlerFocus(w, handler = function (h, ...) {
    
    if(debug){
      print(paste("IN:", match.call()[[1]]))
      print("FOCUS")
    }
    
    # Refresh depending on active tab.
    tab <- svalue(nb)
    tabName <- names(nb)[tab]
    
    # Check if a tab name exist and then perform tasks.
    if(length(tabName) != 0){
      
      if(tabName == .profile_tab_name){
        
        .refreshProfileTab()
        
      } else if(tabName == .analysis_tab_name){
        
        .refreshCETab()

      } else if(tabName == .amplification_tab_name){
        
        .refreshAmplificationTab()
        
      }
      
    } # End check.
    
  })
  
  
# INTERNAL FUNCTIONS ##########################################################

  .addAttributes <- function(data){
    
    if(!is.null(data)){

      val_rnd_chk <- svalue(sim_rnd_seed_chk)
      val_rnd_seed <- as.numeric(svalue(sim_rnd_seed_edt))
      
      # Add attributes to result.
      attr(data, which="pcrsim_version") <- as.character(utils::packageVersion("pcrsim"))
      attr(data, which="pcrsim_date") <- date()
      if(val_rnd_chk){
        attr(data, which="pcrsim_random_seed") <- val_rnd_seed
      }
      
    }
    
    return(data)
    
  }
  
  .updateProjectList <- function(){
    
    # Get project folder.
    projectdir <- svalue(project_fb)
    
    # Create filter for only 'RData' files.
    fileFilter <- paste(".*","\\.", "RData", sep="") 
    
    # Get list of result files.
    .project_path_list <<- list.files(path = projectdir, pattern = fileFilter,
                                      full.names = TRUE, recursive = FALSE,
                                      ignore.case = TRUE, include.dirs = FALSE)
    
    .project_name_list <<- list.files(path = projectdir, pattern = fileFilter,
                                      full.names = FALSE, recursive = FALSE,
                                      ignore.case = TRUE, include.dirs = FALSE)
    
    df <- file.info(.project_path_list)
    
    # Check if any project in list.
    if(length(.project_name_list) > 0){
      
      # Update projects list.    
      project_tbl[,] <<- data.frame(Project=.project_name_list, 
                                    Date=paste(df$mtime),
                                    Size=df$size,
                                    Id=seq(length(.project_name_list)),
                                    stringsAsFactors=FALSE)
      
    }
    
    if(debug){
      print("Project list updated!")
    }
    
  }

  .stringToVector <- function(data, split="[, |_]", type="numeric", empty=NULL, nothing=NULL){
    # Separate a string to a vector by character in 'split'.
    # Convert the vector to 'type'.
    # Remove NA's
    
    # data - vector to convert.
    # split - regex expression.
    # type - convert to this data type.
    # empty - if data is empty string return 'empty'
    # nothing - if data is of length zero return 'nothing'
    
    # Original data.
    valuesBefore <- data

    if(is.character(data)){
      
      if(length(data) == 0) {
        
        data <- nothing
        
      } else if(nchar(data) == 0){
        
        data <- empty
        
      } else {     
        
        # Split.
        data <- unlist(strsplit(x=data, split=split, fixed=FALSE))
        
        if(type == "numeric"){
          data <- as.numeric(data)
        } else if(type == "logical"){
          data <- as.logical(data)
        } else if(type == "character"){
          data <- as.character(data)
        } else {
          warning(paste("Type",type,"not handled, no conversion performed!"))
        }
        
        if(any(is.na(data))){
          data <- data[!is.na(data)]
          warning(paste("Removed invalid values from", valuesBefore, "\nReturn", data))
        }
        
      }
      
    } else {
      
      # message(paste("Not a string (data is", class(data),")! Return unchanged:", data))
      
    }
    
    return(data)
    
  }

  .refreshProfileTab <- function(){
    
    # Get values.    
    val_kit <- svalue(profile_kit_drp)
    val_method <- svalue(profile_method_drp)

    # Update method dropdown.
    profile_method_drp[,] <- getParameter(kit=svalue(profile_kit_drp), what="method")
    
    # Restore selection.
    if(val_method %in% profile_method_drp[]){
      svalue(profile_method_drp) <- val_method
    } else {
      svalue(profile_method_drp, index=TRUE) <- 1
    }

  }

  .refreshSampleTab <- function(){
    
    # Get values.
    valOpt <- svalue (sample_by_opt, index=TRUE)
    valSlope <- as.numeric(svalue(sample_slope_edt))
    valIntercept <- as.numeric(svalue(sample_intercept_edt))
    valMix <- svalue(sample_mix_chk)
    valSim <- svalue(sim_sim_edt)
    
    # Enable/disable fields.
    if (valOpt == 1){
      # Enable.
      enabled(sample_ncells_edt) <- TRUE
      enabled(sample_ncells_sd_edt) <- TRUE
      
      # Disable.
      enabled(sample_conc_edt) <- FALSE
      enabled(sample_conc_sd_edt) <- FALSE
      
      enabled(sample_vol_edt) <- FALSE
      enabled(sample_vol_sd_edt) <- FALSE
      
      enabled(sample_celldna_edt) <- FALSE
      
      enabled(sample_slope_edt) <- FALSE
      enabled(sample_intercept_edt) <- FALSE
      
    } else if (valOpt ==2) {
      
      # Enable.
      enabled(sample_conc_edt) <- TRUE
      enabled(sample_conc_sd_edt) <- TRUE
      
      enabled(sample_vol_edt) <- TRUE
      enabled(sample_vol_sd_edt) <- TRUE
      
      enabled(sample_celldna_edt) <- TRUE
      
      # Disable.
      enabled(sample_ncells_edt) <- FALSE
      enabled(sample_ncells_sd_edt) <- FALSE
      
      enabled(sample_slope_edt) <- FALSE
      enabled(sample_intercept_edt) <- FALSE
      
    } else if (valOpt ==3) {
      
      # Enable.
      enabled(sample_slope_edt) <- TRUE
      enabled(sample_intercept_edt) <- TRUE
      
      # Disable.
      enabled(sample_conc_edt) <- FALSE
      enabled(sample_conc_sd_edt) <- FALSE
      
      enabled(sample_vol_edt) <- FALSE
      enabled(sample_vol_sd_edt) <- FALSE
      
      enabled(sample_celldna_edt) <- FALSE
      
      enabled(sample_ncells_edt) <- FALSE
      enabled(sample_ncells_sd_edt) <- FALSE
      
    } else {
      
      warning(paste("Option=", valOpt, "in tab", .sample_tab_name, "not handled!"))
      
    }
    
    # Calculate example concentrations.
    if(!is.na(valSlope) && !is.na(valIntercept)){
      if (valSlope != "" && valIntercept!= ""){
        # Concentration = 10^(slope * size + intercept)
        conc1 <- 10^(valSlope * 100 + valIntercept)
        conc2 <- 10^(valSlope * 400 + valIntercept)
        # Update values.
        svalue(sample_ex1_edt) <- round(conc1, 4)
        svalue(sample_ex2_edt) <- round(conc2, 4)
      }
    }
    
    # Update mixture checkbox label.
    labTmp <- paste("Simulate a mixture with", valSim,
                    "contributors (set in", .simulation_tab_name, "tab)")
    sample_mix_chk[] <- labTmp 
    
  }
  
  .refreshDegradationTab <- function(){
    
    # Get values.
    valOpt <- svalue(deg_opt, index=TRUE)
    valDeg <- as.numeric(svalue(deg_spn))
    valConc <- svalue(deg_conc_edt)
    valSize <- svalue(deg_size_edt)
    valExtDebug <- svalue(file_f1g1_ext_debug_chk)
    
    # Check and convert to vector.
    valConc <- .stringToVector(data=valConc, empty=NA, nothing=NA)
    valSize <- .stringToVector(data=valSize, empty=NA, nothing=NA)
    
    # Enable/disable widgets depending on selected option.
    if (valOpt == 1){
      # Manual.
      
      # Enable widgets.
      enabled(deg_spn) <- TRUE
      
      # Disable widgets.
      enabled(deg_conc_edt) <- FALSE
      enabled(deg_size_edt) <- FALSE

      # Update degradation label.
      svalue(deg_pam_lbl) <- valDeg

      if(debug){
        print(paste("Degradation parameter =", svalue(deg_pam_lbl)))
      }
      
    } else if (valOpt == 2) {
      # Concentrations.
      
      # Enable widgets.
      enabled(deg_conc_edt) <- TRUE
      enabled(deg_size_edt) <- TRUE
      
      # Disable widgets.
      enabled(deg_spn) <- FALSE

      # Calculate degradation parameter.
      if (!is.na(valConc) && !is.na(valSize)){
        if (length(valConc) == length(valSize)){

          # Calculate degradation parameter.
          val_deg <- calculateDegradation(conc=valConc, size=valSize)
                                         
          # Update degradation label.
          svalue(deg_pam_lbl) <- val_deg
          
          if(debug){
            print(paste("Degradation parameter =", svalue(deg_pam_lbl)))
          }
          
        }
      }
      
    } else {
      warning(paste("Option=", valOpt, "in tab", .degradation_tab_name, "not handled!"))
    }
    
  }
  
  
  .refreshDilutionTab <- function(){
    
    # Declare temporary variable.
    val_tmp <- NULL
    
    # Get values.
    val_conc <- NULL # Target conc calculated below.
    val_amount <- as.numeric(svalue(amp_amount_edt)) # Optimal amount.
    val_volume <- as.numeric(svalue(ex_vol_edt)) # Final dilution volume.
    val_aliquot <- as.numeric(svalue(amp_aliquot_edt)) # Final dilution volume.
    val_tolerance <- as.numeric(svalue(dil_tolerance_edt)) # Tolerance around target conc.
    val_target_opt <- svalue(dil_target_opt, index=TRUE)
    val_volume_opt <- svalue(dil_volume_opt, index=TRUE)
    
    # Enable/disable widget group depending on selected option.
    if (val_volume_opt == 1){
      enabled(dil_volume_lbl) <- TRUE
      enabled(dil_volume_edt) <- FALSE
    } else if (val_volume_opt == 2) {
      enabled(dil_volume_lbl) <- FALSE
      enabled(dil_volume_edt) <- TRUE
    }
    
    # Enable/disable widget group depending on selected option.
    if (val_target_opt == 1){
      enabled(dil_target_lbl) <- TRUE
      enabled(dil_target_edt) <- FALSE
    } else if (val_target_opt == 2) {
      enabled(dil_target_lbl) <- FALSE
      enabled(dil_target_edt) <- TRUE
    }
    
    # Update amount, volume, and aliquot labels:
    svalue(dil_amount_lbl) <- val_amount
    svalue(dil_volume_lbl) <- val_volume
    svalue(dil_aliquot_lbl) <- val_aliquot
    
    # Update target concentration label.
    if(!is.na(val_aliquot) && !is.na(val_amount)){
      val_conc <- val_amount / val_aliquot
      svalue(dil_target_lbl) <- round(val_conc, 4)
    }
    
    # Update maximum concentration label.
    if (!is.na(val_tolerance)){
      # Calculate min and max amount.
      conc_max <- round((val_conc + (val_tolerance * val_conc)), 4)
      # Update label.
      svalue(dil_max_lbl) <- paste("(dilute when concentration >",
                                   conc_max, " ng/\u00B5l)", sep="")
    } else {
      # Update label.
      svalue(dil_max_lbl) <- ""
    }
    
    
  }
  
  .refreshAmplificationTab <- function(){
    # Get values.    
    val <- svalue (amp_f3_eff_opt, index=TRUE)
    val_kit <- svalue(profile_kit_drp)
    val_method <- svalue(profile_method_drp)

    # Update label.
    svalue(amp_current_lbl) <- paste("Current kit: ", val_kit,
                                     " and method: ", val_method,
                                     " (set in '", .profile_tab_name,
                                     "' tab)", sep="")
    
    # Enable/disable widget group depending on selected option.
    if (val == 1){
      enabled(amp_f3g1) <- FALSE
    } else if (val == 2) {
      enabled(amp_f3g1) <- TRUE
    }
    
    # Check state of widget group.
    if (enabled(amp_f3g1)) {
      
      # Check state of stutter checkbox.
      if(svalue(amp_stutter_chk)) {
        enabled(amp_eff_stutt_spn) <- TRUE
        enabled(amp_eff_stutt_sd_edt) <- TRUE
      } else {
        enabled(amp_eff_stutt_spn) <- FALSE
        enabled(amp_eff_stutt_sd_edt) <- FALSE
      }
      
    }
    
  }

  .refreshCETab <- function(){
    
    # Get values.    
    val_kit <- svalue(profile_kit_drp)
    val_method <- svalue(profile_method_drp)
    #val_method <- svalue(ce_method_drp)
    val_settings <- svalue(ce_settings_opt, index = TRUE)
    val_debug <- svalue(file_f1g1_debug_chk)
    valExtDebug <- svalue(file_f1g1_ext_debug_chk)
    
    svalue(ce_current_lbl) <- paste("Current kit: ", val_kit,
                                    " and method: ", val_method,
                                    " (set in '", .profile_tab_name,
                                    "' tab)", sep="")

    if(val_settings == 1){
      # Settings from method.
      
      # Inactivate gui.
      enabled(ce_f1) <- FALSE
      enabled(ce_f3) <- FALSE
      
      if(is.null(val_method) || nchar(val_method) == 0){
        val_method <- "DEFAULT"
      }
        
      # Get parameters.    
      val_threshold <- getParameter(kit=val_kit, what="THRESHOLD", method=val_method)
      val_scaling <- getParameter(kit=val_kit, what="SCALING", method=val_method)
      
      if(val_debug){
        print(".refreshCETab")
        print("val_method:")
        print(val_method)
        print("val_threshold:")
        print(val_threshold)
        print("val_scaling:")
        print(val_scaling)
      }
      
      # Update gui.
      svalue(ce_t_intercept_edt) <- val_threshold$Intercept
      svalue(ce_t_slope_edt) <- val_threshold$Slope
      svalue(ce_t_residual_edt) <- val_threshold$Sigma
      svalue(ce_intercept_edt) <- val_scaling$Intercept
      svalue(ce_slope_edt) <- val_scaling$Slope
      svalue(ce_residual_edt) <- val_scaling$Sigma

    } else if(val_settings == 2) {
      # Custom settings.
      
      # Activate gui.
      enabled(ce_f1) <- TRUE
      enabled(ce_f3) <- TRUE
      
    } else {
      
      stop("CE seetings option = ", val_settings, " not implemented!")
      
    }
    
  }

  .refreshSimulationTab <- function(){

    # Get values.    
    val_rnd_chk <- svalue(sim_rnd_seed_chk)
    val_save <- svalue(sim_save_chk)
    val_profile <- svalue(sim_profile_chk)
    val_sample <- svalue(sim_sample_chk)
    val_deg <- svalue(sim_deg_chk)
    val_extr <- svalue(sim_extr_chk)
    val_norm <- svalue(sim_norm_chk)
    val_pcr <- svalue(sim_pcr_chk)
    val_ce <- svalue(sim_ce_chk)
    #val_debug <- svalue(file_f1g1_debug_chk)
    #valExtDebug <- svalue(file_f1g1_ext_debug_chk)
    
    # Enable buttons.
    enabled(sim_sim_btn) <- TRUE
    enabled(sim_view_btn) <- TRUE
    enabled(sim_view_ph_btn) <- TRUE
    enabled(sim_epg_btn) <- TRUE
    
    # Change gui depending on simulation settings.
    if(val_rnd_chk){
      enabled(sim_rnd_seed_edt) <- TRUE
    } else {
      enabled(sim_rnd_seed_edt) <- FALSE
    }
    
    # Change gui depending on simulation settings.
    if(val_profile){
      
      # Activate in gui.
      enabled(sim_sample_chk) <- TRUE

    } else {
      
      # Inactivate gui.
      enabled(sim_sample_chk) <- FALSE

    }
    
    # Change gui depending on simulation settings.
    if(val_sample){
      
      # Activate in gui.
      enabled(sim_deg_chk) <- TRUE
      enabled(sim_extr_chk) <- TRUE

    } else {
      
      # Inactivate gui.
      enabled(sim_deg_chk) <- FALSE
      enabled(sim_extr_chk) <- FALSE
      
    }
    
    # Change gui depending on simulation settings.
    if(val_extr){
      
      # Activate in gui.
      enabled(sim_norm_chk) <- TRUE
      enabled(sim_pcr_chk) <- TRUE
      
    } else {
      
      # Inactivate gui.
      enabled(sim_norm_chk) <- FALSE
      enabled(sim_pcr_chk) <- FALSE
      
    }
    
    # Change gui depending on simulation settings.
    if(val_pcr){
      
      # Activate in gui.
      enabled(sim_ce_chk) <- TRUE

    } else {
      
      # Inactivate gui.
      enabled(sim_ce_chk) <- FALSE

    }
    
    # Change gui depending on save value.
    if(val_save){

      # Activate in gui.
      enabled(sim_save_s_chk) <- TRUE
      enabled(sim_save_t_chk) <- TRUE
      enabled(sim_save_n_chk) <- TRUE
      
    } else {
      
      # Inactivate in gui.
      enabled(sim_save_s_chk) <- FALSE
      enabled(sim_save_t_chk) <- FALSE
      enabled(sim_save_n_chk) <- FALSE
      
    }
    
  }
  
  .refreshWs <- function(){
    
    # Get data frames in global workspace.
    dfs <- strvalidator::listObjects(env=.GlobalEnv, obj.class="data.frame")
    
    if(!is.null(dfs)){
      
      blockHandler(file_f2g1_drp)
      
      # Populate drop list.
      file_f2g1_drp[] <- c("<Select dataframe>", dfs)
      
      # Select first item.
      svalue(file_f2g1_drp, index=TRUE) <- 1 
      
      unblockHandler(file_f2g1_drp)
      
    }
  }
  
  .refreshLoaded <- function(){
    
    if(debug){
      print(paste("IN:", match.call()[[1]]))
    }
    
    # Get data frames.
    dfs <- strvalidator::listObjects(env=.pcrsim_env, obj.class=c("data.frame","ggplot",
                                                                        "environment"))
    
    if(!is.null(dfs)){
      
      #blockHandler(file_f1g1_tbl)
      
      #delete(file_f1g1_f, file_f1g1_tbl)
      #file_f1g1_tbl <<- gtable(items=dfs, 
      #                        multiple = TRUE,
      #                        container = file_f1g1_f) 
      
      # Populate table.
      file_f1g1_tbl[] <- dfs
      
      #unblockHandler(file_f1g1_tbl)
      
    }
    
    if(debug){
      print(paste("EXIT:", match.call()[[1]]))
    }
  }
  
  .loadSavedSettings <- function(env){
    # Argument 'env' is either the name of the environment or the environment itself.
    
    # Check argument.
    if(class(env) == "environment"){
      # Ok, env is environemnt.
    } else if(class(env) == "character"){
      # Get environment that varible name refer to.
      env <- get(x=env, envir=.pcrsim_env, inherits=FALSE)
    } else {
      message(paste("Argument of class", class(env), "not handled!"))
    }
    
    if(debug){
      print(paste("IN:", match.call()[[1]]))
      print(paste("Environment:", env))
    }
    
    # TAB FILE:
    if(debug){
      print("TAB FILE")
    }
    if(exists(".pcrsim_file_debug_chk", envir=env, inherits=FALSE)){
      svalue(file_f1g1_debug_chk) <- get(".pcrsim_file_debug_chk", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_file_ext_debug_chk", envir=env, inherits=FALSE)){
      svalue(file_f1g1_ext_debug_chk) <- get(".pcrsim_file_ext_debug_chk", envir=env, inherits=FALSE)
    }
    
    # TAB PROFILE:
    if(debug){
      print("TAB PROFILE")
    }
    if(exists(".pcrsim_profile_kit", envir=env, inherits=FALSE)){
      svalue(profile_kit_drp) <- get(".pcrsim_profile_kit", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_profile_method", envir=env, inherits=FALSE)){
      svalue(profile_method_drp) <- get(".pcrsim_profile_method", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_profile_opt", envir=env, inherits=FALSE)){
      svalue(profile_opt) <- get(".pcrsim_profile_opt", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_profile_pop", envir=env, inherits=FALSE)){
      svalue(profile_pop_drp) <- get(".pcrsim_profile_pop", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_profile_profile", envir=env, inherits=FALSE)){
      profile_tbl[,] <- get(".pcrsim_profile_profile", envir=env, inherits=FALSE)
    }
    
    # TAB SAMPLE:
    if(debug){
      print("TAB SAMPLE")
    }
    if(exists(".pcrsim_sample_name", envir=env, inherits=FALSE)){
      svalue(sample_name_edt) <- get(".pcrsim_sample_name", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sample_mix", envir=env, inherits=FALSE)){
      svalue(sample_mix_chk) <- get(".pcrsim_sample_mix", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sample_by", envir=env, inherits=FALSE)){
      svalue(sample_by_opt) <- get(".pcrsim_sample_by", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sample_ncells", envir=env, inherits=FALSE)){
      svalue(sample_ncells_edt) <- get(".pcrsim_sample_ncells", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sample_ncells_sd", envir=env, inherits=FALSE)){
      svalue(sample_ncells_sd_edt) <- get(".pcrsim_sample_ncells_sd", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sample_conc", envir=env, inherits=FALSE)){
      svalue(sample_conc_edt) <- get(".pcrsim_sample_conc", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sample_conc_sd", envir=env, inherits=FALSE)){
      svalue(sample_conc_sd_edt) <- get(".pcrsim_sample_conc_sd", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sample_celldna", envir=env, inherits=FALSE)){
      svalue(sample_celldna_edt) <- get(".pcrsim_sample_celldna", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sample_vol", envir=env, inherits=FALSE)){
      svalue(sample_vol_edt) <- get(".pcrsim_sample_vol", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sample_vol_sd", envir=env, inherits=FALSE)){
      svalue(sample_vol_sd_edt) <- get(".pcrsim_sample_vol_sd", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sample_slope", envir=env, inherits=FALSE)){
      svalue(sample_slope_edt) <- get(".pcrsim_sample_slope", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sample_intercept", envir=env, inherits=FALSE)){
      svalue(sample_intercept_edt) <- get(".pcrsim_sample_intercept", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sample_type", envir=env, inherits=FALSE)){
      svalue(sample_type_cmb) <- get(".pcrsim_sample_type", envir=env, inherits=FALSE)
    }
    
    # TAB DEGRADATION:
    if(debug){
      print("TAB DEGRADATION")
    }
    if(exists(".pcrsim_deg_lbl", envir=env, inherits=FALSE)){
      svalue(deg_pam_lbl) <- get(".pcrsim_deg_lbl", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_deg_target", envir=env, inherits=FALSE)){
      svalue(deg_target_edt) <- get(".pcrsim_deg_target", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_deg_opt", envir=env, inherits=FALSE)){
      svalue(deg_opt) <- get(".pcrsim_deg_opt", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_deg_spn", envir=env, inherits=FALSE)){
      svalue(deg_spn) <- get(".pcrsim_deg_spn", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_deg_conc", envir=env, inherits=FALSE)){
      svalue(deg_conc_edt) <- get(".pcrsim_deg_conc", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_deg_size", envir=env, inherits=FALSE)){
      svalue(deg_size_edt) <- get(".pcrsim_deg_size", envir=env, inherits=FALSE)
    }

    # TAB EXTRACTION:
    if(debug){
      print("TAB EXTRACTION")
    }
    if(exists(".pcrsim_ex_eff", envir=env, inherits=FALSE)){
      svalue(ex_eff_spn) <- get(".pcrsim_ex_eff", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_ex_eff_sd", envir=env, inherits=FALSE)){
      svalue(ex_eff_sd_edt) <- get(".pcrsim_ex_eff_sd", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_ex_vol", envir=env, inherits=FALSE)){
      svalue(ex_vol_edt) <- get(".pcrsim_ex_vol", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_ex_vol_sd", envir=env, inherits=FALSE)){
      svalue(ex_vol_sd_edt) <- get(".pcrsim_ex_vol_sd", envir=env, inherits=FALSE)
    }
    
    # TAB DILUTION:
    if(debug){
      print("TAB DILUTION")
    }
    if(exists(".pcrsim_dil_tolerance", envir=env, inherits=FALSE)){
      svalue(dil_tolerance_edt) <- get(".pcrsim_dil_tolerance", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_dil_target", envir=env, inherits=FALSE)){
      svalue(dil_target_edt) <- get(".pcrsim_dil_target", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_dil_target_opt", envir=env, inherits=FALSE)){
      svalue(dil_target_opt) <- get(".pcrsim_dil_target_opt", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_dil_accuracy", envir=env, inherits=FALSE)){
      svalue(dil_accuracy_edt) <- get(".pcrsim_dil_accuracy", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_dil_volume", envir=env, inherits=FALSE)){
      svalue(dil_volume_edt) <- get(".pcrsim_dil_volume", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_dil_volume_opt", envir=env, inherits=FALSE)){
      svalue(dil_volume_opt) <- get(".pcrsim_dil_volume_opt", envir=env, inherits=FALSE)
    }
    
    # TAB PCR AMPLIFICATION:
    if(debug){
      print("TAB PCR")
    }
    if(exists(".pcrsim_amp_aliquot", envir=env, inherits=FALSE)){
      svalue(amp_aliquot_edt) <- get(".pcrsim_amp_aliquot", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_amp_aliquot_sd", envir=env, inherits=FALSE)){
      svalue(amp_aliquot_sd_edt) <- get(".pcrsim_amp_aliquot_sd", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_amp_amount", envir=env, inherits=FALSE)){
      svalue(amp_amount_edt) <- get(".pcrsim_amp_amount", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_amp_total_vol", envir=env, inherits=FALSE)){
      svalue(amp_total_vol_edt) <- get(".pcrsim_amp_total_vol", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_amp_total_vol_sd", envir=env, inherits=FALSE)){
      svalue(amp_total_vol_sd_edt) <- get(".pcrsim_amp_total_vol_sd", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_amp_cyc_sb", envir=env, inherits=FALSE)){
      svalue(cyc_sb) <- get(".pcrsim_amp_cyc_sb", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_amp_eff_opt", envir=env, inherits=FALSE)){
      svalue(amp_f3_eff_opt) <- get(".pcrsim_amp_eff_opt", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_amp_eff", envir=env, inherits=FALSE)){
      svalue(amp_f3_eff_spn) <- get(".pcrsim_amp_eff", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_amp_eff_sd", envir=env, inherits=FALSE)){
      svalue(amp_f3_eff_sd_edt) <- get(".pcrsim_amp_eff_sd", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_amp_stutt_eff", envir=env, inherits=FALSE)){
      svalue(amp_eff_stutt_spn) <- get(".pcrsim_amp_stutt_eff", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_amp_stutt_eff_sd", envir=env, inherits=FALSE)){
      svalue(amp_eff_stutt_sd_edt) <- get(".pcrsim_amp_stutt_eff_sd", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_amp_sim_stutt", envir=env, inherits=FALSE)){
      svalue(amp_stutter_chk) <- get(".pcrsim_amp_sim_stutt", envir=env, inherits=FALSE)
    }
    
    # TAB CE ANALYSIS:
    if(debug){
      print("TAB CE")
    }
    if(exists(".pcrsim_ce_aliquot", envir=env, inherits=FALSE)){
      svalue(ce_aliquot_edt) <- get(".pcrsim_ce_aliquot", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_ce_aliquot_sd", envir=env, inherits=FALSE)){
      svalue(ce_aliquot_sd_edt) <- get(".pcrsim_ce_aliquot_sd", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_ce_t_intercept", envir=env, inherits=FALSE)){
      svalue(ce_t_intercept_edt) <- get(".pcrsim_ce_t_intercept", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_ce_t_slope", envir=env, inherits=FALSE)){
      svalue(ce_t_slope_edt) <- get(".pcrsim_ce_t_slope", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_ce_t_residual", envir=env, inherits=FALSE)){
      svalue(ce_t_residual_edt) <- get(".pcrsim_ce_t_residual", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_ce_intercept", envir=env, inherits=FALSE)){
      svalue(ce_intercept_edt) <- get(".pcrsim_ce_intercept", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_ce_slope", envir=env, inherits=FALSE)){
      svalue(ce_slope_edt) <- get(".pcrsim_ce_slope", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_ce_residual", envir=env, inherits=FALSE)){
      svalue(ce_residual_edt) <- get(".pcrsim_ce_residual", envir=env, inherits=FALSE)
    }

    # TAB EPG:
    if(debug){
      print("TAB EPG")
    }
    # NB! Naming follow code copied from strvalidator::generateEPG_gui()
    if(exists(".pcrsim_epg_scale", envir=env, inherits=FALSE)){
      svalue(f1_scale_opt) <- get(".pcrsim_epg_scale", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_epg_size", envir=env, inherits=FALSE)){
      svalue(f1_size_spb) <- get(".pcrsim_epg_size", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_epg_vjust", envir=env, inherits=FALSE)){
      svalue(f1_vjust_spb) <- get(".pcrsim_epg_vjust", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_epg_angle", envir=env, inherits=FALSE)){
      svalue(f1_angle_spb) <- get(".pcrsim_epg_angle", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_epg_hjust", envir=env, inherits=FALSE)){
      svalue(f1_hjust_spb) <- get(".pcrsim_epg_hjust", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_epg_expand", envir=env, inherits=FALSE)){
      svalue(f1_expand_spb) <- get(".pcrsim_epg_expand", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_epg_at", envir=env, inherits=FALSE)){
      svalue(f1_at_spb) <- get(".pcrsim_epg_at", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_epg_ignore", envir=env, inherits=FALSE)){
      svalue(f1_ignore_chk) <- get(".pcrsim_epg_ignore", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_epg_wrap", envir=env, inherits=FALSE)){
      svalue(f1_wrap_chk) <- get(".pcrsim_epg_wrap", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_epg_fix", envir=env, inherits=FALSE)){
      svalue(f1_fix_chk) <- get(".pcrsim_epg_fix", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_epg_collapse", envir=env, inherits=FALSE)){
      svalue(f1_collapse_chk) <- get(".pcrsim_epg_collapse", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_epg_box", envir=env, inherits=FALSE)){
      svalue(f1_box_chk) <- get(".pcrsim_epg_box", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_epg_peaks", envir=env, inherits=FALSE)){
      svalue(f1_peaks_chk) <- get(".pcrsim_epg_peaks", envir=env, inherits=FALSE)
    }
    
    # TAB SIMULATION:
    if(debug){
      print("TAB SIMULATION")
    }
    if(exists(".pcrsim_sim_sim", envir=env, inherits=FALSE)){
      svalue(sim_sim_edt) <- get(".pcrsim_sim_sim", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sim_name", envir=env, inherits=FALSE)){
      svalue(sim_name_edt) <- get(".pcrsim_sim_name", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sim_title", envir=env, inherits=FALSE)){
      svalue(sim_epg_title_edt) <- get(".pcrsim_sim_title", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sim_rnd_chk", envir=env, inherits=FALSE)){
      svalue(sim_rnd_seed_chk) <- get(".pcrsim_sim_rnd_chk", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sim_rnd_edt", envir=env, inherits=FALSE)){
      svalue(sim_rnd_seed_edt) <- get(".pcrsim_sim_rnd_edt", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sim_update_chk", envir=env, inherits=FALSE)){
      svalue(sim_update_chk) <- get(".pcrsim_sim_update_chk", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sim_ce_save_chk", envir=env, inherits=FALSE)){
      svalue(sim_save_chk) <- get(".pcrsim_sim_ce_save_chk", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sim_ce_save_s_chk", envir=env, inherits=FALSE)){
      svalue(sim_save_s_chk) <- get(".pcrsim_sim_ce_save_s_chk", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sim_ce_save_t_chk", envir=env, inherits=FALSE)){
      svalue(sim_save_t_chk) <- get(".pcrsim_sim_ce_save_t_chk", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sim_ce_save_n_chk", envir=env, inherits=FALSE)){
      svalue(sim_save_n_chk) <- get(".pcrsim_sim_ce_save_n_chk", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sim_profile_chk", envir=env, inherits=FALSE)){
      svalue(sim_profile_chk) <- get(".pcrsim_sim_profile_chk", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sim_sample_chk", envir=env, inherits=FALSE)){
      svalue(sim_sample_chk) <- get(".pcrsim_sim_sample_chk", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sim_extraction_chk", envir=env, inherits=FALSE)){
      svalue(sim_extr_chk) <- get(".pcrsim_sim_extraction_chk", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sim_degradation_chk", envir=env, inherits=FALSE)){
      svalue(sim_deg_chk) <- get(".pcrsim_sim_degradation_chk", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sim_dilution_chk", envir=env, inherits=FALSE)){
      svalue(sim_norm_chk) <- get(".pcrsim_sim_dilution_chk", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sim_pcr_chk", envir=env, inherits=FALSE)){
      svalue(sim_pcr_chk) <- get(".pcrsim_sim_pcr_chk", envir=env, inherits=FALSE)
    }
    if(exists(".pcrsim_sim_ce_chk", envir=env, inherits=FALSE)){
      svalue(sim_ce_chk) <- get(".pcrsim_sim_ce_chk", envir=env, inherits=FALSE)
    }
    
    if(debug){
      print("Saved settings loaded!")
      print(paste("EXIT:", match.call()[[1]]))
    }
    
  }
  
  .saveSettings <- function(env = .pcrsim_env){
    # Argument 'env' is either the name of the environment or the environment itself.
    
    # Check argument.
    if(class(env) == "environment"){
      # Ok, env is environemnt.
    } else if(class(env) == "character"){
      # Get environment that varible name refer to.
      env <- get(x=env, envir=.pcrsim_env, inherits=FALSE)
    } else {
      message(paste("Argument of class", class(env), "not handled!"))
    }
    
    if(debug){
      print(paste("IN:", match.call()[[1]]))
      print(paste("Environment:", env))
    }
    
    # Then save settings if true.
    if(svalue(file_f1g1_save_chk)){
      
      # TAB FILE:
      assign(x=".pcrsim_file_savegui", value=svalue(file_f1g1_save_chk), envir=env, inherits=FALSE)
      assign(x=".pcrsim_file_debug_chk", value=svalue(file_f1g1_debug_chk), envir=env, inherits=FALSE)
      assign(x=".pcrsim_file_ext_debug_chk", value=svalue(file_f1g1_ext_debug_chk), envir=env, inherits=FALSE)
      
      # TAB PROFILE:
      assign(x=".pcrsim_profile_kit", value=svalue(profile_kit_drp), envir=env, inherits=FALSE)
      assign(x=".pcrsim_profile_method", value=svalue(profile_method_drp), envir=env, inherits=FALSE)
      assign(x=".pcrsim_profile_opt", value=svalue(profile_opt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_profile_pop", value=svalue(profile_pop_drp), envir=env, inherits=FALSE)
      assign(x=".pcrsim_profile_profile", value=profile_tbl[,], envir=env, inherits=FALSE)
      
      # TAB SAMPLE:
      assign(x=".pcrsim_sample_name", value=svalue(sample_name_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sample_mix", value=svalue(sample_mix_chk), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sample_by", value=svalue(sample_by_opt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sample_ncells", value=svalue(sample_ncells_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sample_ncells_sd", value=svalue(sample_ncells_sd_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sample_conc", value=svalue(sample_conc_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sample_conc_sd", value=svalue(sample_conc_sd_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sample_celldna", value=svalue(sample_celldna_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sample_vol", value=svalue(sample_vol_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sample_vol_sd", value=svalue(sample_vol_sd_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sample_slope", value=svalue(sample_slope_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sample_intercept", value=svalue(sample_intercept_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sample_type", value=svalue(sample_type_cmb), envir=env, inherits=FALSE)
      
      # TAB DEGRADATION:
      assign(x=".pcrsim_deg_lbl", value=svalue(deg_pam_lbl), envir=env, inherits=FALSE)
      assign(x=".pcrsim_deg_target", value=svalue(deg_target_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_deg_opt", value=svalue(deg_opt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_deg_spn", value=svalue(deg_spn), envir=env, inherits=FALSE)
      assign(x=".pcrsim_deg_conc", value=svalue(deg_conc_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_deg_size", value=svalue(deg_size_edt), envir=env, inherits=FALSE)

      # TAB EXTRACTION:
      assign(x=".pcrsim_ex_eff", value=svalue(ex_eff_spn), envir=env, inherits=FALSE)
      assign(x=".pcrsim_ex_eff_sd", value=svalue(ex_eff_sd_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_ex_vol", value=svalue(ex_vol_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_ex_vol_sd", value=svalue(ex_vol_sd_edt), envir=env, inherits=FALSE)
      
      # TAB DILUTION:
      assign(x=".pcrsim_dil_tolerance", value=svalue(dil_tolerance_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_dil_target", value=svalue(dil_target_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_dil_target_opt", value=svalue(dil_target_opt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_dil_accuracy", value=svalue(dil_accuracy_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_dil_volume", value=svalue(dil_volume_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_dil_volume_opt", value=svalue(dil_volume_opt), envir=env, inherits=FALSE)
      
      # TAB AMPLIFICATION:
      assign(x=".pcrsim_amp_aliquot", value=svalue(amp_aliquot_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_amp_aliquot_sd", value=svalue(amp_aliquot_sd_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_amp_amount", value=svalue(amp_amount_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_amp_total_vol", value=svalue(amp_total_vol_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_amp_total_vol_sd", value=svalue(amp_total_vol_sd_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_amp_cyc_sb", value=svalue(cyc_sb), envir=env, inherits=FALSE)
      assign(x=".pcrsim_amp_eff_opt", value=svalue(amp_f3_eff_opt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_amp_eff", value=svalue(amp_f3_eff_spn), envir=env, inherits=FALSE)
      assign(x=".pcrsim_amp_eff_sd", value=svalue(amp_f3_eff_sd_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_amp_stutt_eff", value=svalue(amp_eff_stutt_spn), envir=env, inherits=FALSE)
      assign(x=".pcrsim_amp_stutt_eff_sd", value=svalue(amp_eff_stutt_sd_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_amp_sim_stutt", value=svalue(amp_stutter_chk), envir=env, inherits=FALSE)
      
      # TAB ANALYSIS:
      assign(x=".pcrsim_ce_aliquot", value=svalue(ce_aliquot_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_ce_aliquot_sd", value=svalue(ce_aliquot_sd_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_ce_t_intercept", value=svalue(ce_t_intercept_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_ce_t_slope", value=svalue(ce_t_slope_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_ce_t_residual", value=svalue(ce_t_residual_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_ce_intercept", value=svalue(ce_intercept_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_ce_slope", value=svalue(ce_slope_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_ce_residual", value=svalue(ce_residual_edt), envir=env, inherits=FALSE)

      # TAB EPG:
      # NB! Naming follow code copied from strvalidator::generateEPG_gui()
      assign(x=".pcrsim_epg_scale", value=svalue(f1_scale_opt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_epg_size", value=svalue(f1_size_spb), envir=env, inherits=FALSE)
      assign(x=".pcrsim_epg_vjust", value=svalue(f1_vjust_spb), envir=env, inherits=FALSE)
      assign(x=".pcrsim_epg_angle", value=svalue(f1_angle_spb), envir=env, inherits=FALSE)
      assign(x=".pcrsim_epg_hjust", value=svalue(f1_hjust_spb), envir=env, inherits=FALSE)
      assign(x=".pcrsim_epg_expand", value=svalue(f1_expand_spb), envir=env, inherits=FALSE)
      assign(x=".pcrsim_epg_at", value=svalue(f1_at_spb), envir=env, inherits=FALSE)
      assign(x=".pcrsim_epg_ignore", value=svalue(f1_ignore_chk), envir=env, inherits=FALSE)
      assign(x=".pcrsim_epg_wrap", value=svalue(f1_wrap_chk), envir=env, inherits=FALSE)
      assign(x=".pcrsim_epg_fix", value=svalue(f1_fix_chk), envir=env, inherits=FALSE)
      assign(x=".pcrsim_epg_collapse", value=svalue(f1_collapse_chk), envir=env, inherits=FALSE)
      assign(x=".pcrsim_epg_box", value=svalue(f1_box_chk), envir=env, inherits=FALSE)
      assign(x=".pcrsim_epg_peaks", value=svalue(f1_peaks_chk), envir=env, inherits=FALSE)

      # TAB SIMULATION:
      assign(x=".pcrsim_sim_sim", value=svalue(sim_sim_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sim_name", value=svalue(sim_name_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sim_title", value=svalue(sim_epg_title_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sim_rnd_chk", value=svalue(sim_rnd_seed_chk), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sim_rnd_edt", value=svalue(sim_rnd_seed_edt), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sim_update_chk", value=svalue(sim_update_chk), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sim_ce_save_chk", value=svalue(sim_save_chk), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sim_ce_save_s_chk", value=svalue(sim_save_s_chk), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sim_ce_save_t_chk", value=svalue(sim_save_t_chk), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sim_ce_save_n_chk", value=svalue(sim_save_n_chk), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sim_profile_chk", value=svalue(sim_profile_chk), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sim_sample_chk", value=svalue(sim_sample_chk), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sim_extraction_chk", value=svalue(sim_extr_chk), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sim_degradation_chk", value=svalue(sim_deg_chk), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sim_dilution_chk", value=svalue(sim_norm_chk), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sim_pcr_chk", value=svalue(sim_pcr_chk), envir=env, inherits=FALSE)
      assign(x=".pcrsim_sim_ce_chk", value=svalue(sim_ce_chk), envir=env, inherits=FALSE)
      
    } else { # or remove all saved values if false.
      
      # TAB FILE:
      if(debug){
        print("TAB FILE")
      }
      if(exists(".pcrsim_file_savegui", envir=env, inherits=FALSE)){
        remove(".pcrsim_file_savegui", envir = env)
      }
      if(exists(".pcrsim_file_debug_chk", envir=env, inherits=FALSE)){
        remove(".pcrsim_file_debug_chk", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_file_ext_debug_chk", envir=env, inherits=FALSE)){
        remove(".pcrsim_file_ext_debug_chk", envir=env, inherits=FALSE)
      }
      
      # TAB PROFILE:
      if(debug){
        print("TAB PROFILE")
      }
      if(exists(".pcrsim_profile_kit", envir=env, inherits=FALSE)){
        remove(".pcrsim_profile_kit", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_profile_method", envir=env, inherits=FALSE)){
        remove(".pcrsim_profile_method", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_profile_opt", envir=env, inherits=FALSE)){
        remove(".pcrsim_profile_opt", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_profile_pop", envir=env, inherits=FALSE)){
        remove(".pcrsim_profile_pop", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_profile_profile", envir=env, inherits=FALSE)){
        remove(".pcrsim_profile_profile", envir=env, inherits=FALSE)
      }
      
      # TAB SAMPLE:
      if(debug){
        print("TAB SAMPLE")
      }
      if(exists(".pcrsim_sample_name", envir=env, inherits=FALSE)){
        remove(".pcrsim_sample_name", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sample_mix", envir=env, inherits=FALSE)){
        remove(".pcrsim_sample_mix", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sample_by", envir=env, inherits=FALSE)){
        remove(".pcrsim_sample_by", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sample_ncells", envir=env, inherits=FALSE)){
        remove(".pcrsim_sample_ncells", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sample_ncells_sd", envir=env, inherits=FALSE)){
        remove(".pcrsim_sample_ncells_sd", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sample_conc", envir=env, inherits=FALSE)){
        remove(".pcrsim_sample_conc", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sample_conc_sd", envir=env, inherits=FALSE)){
        remove(".pcrsim_sample_conc_sd", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sample_celldna", envir=env, inherits=FALSE)){
        remove(".pcrsim_sample_celldna", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sample_vol", envir=env, inherits=FALSE)){
        remove(".pcrsim_sample_vol", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sample_vol_sd", envir=env, inherits=FALSE)){
        remove(".pcrsim_sample_vol_sd", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sample_slope", envir=env, inherits=FALSE)){
        remove(".pcrsim_sample_slope", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sample_intercept", envir=env, inherits=FALSE)){
        remove(".pcrsim_sample_intercept", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sample_type", envir=env, inherits=FALSE)){
        remove(".pcrsim_sample_type", envir=env, inherits=FALSE)
      }
      
      # TAB DEGRADATION:
      if(debug){
        print("TAB DEGRADATION")
      }
      if(exists(".pcrsim_deg_lbl", envir=env, inherits=FALSE)){
        remove(".pcrsim_deg_lbl", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_deg_target", envir=env, inherits=FALSE)){
        remove(".pcrsim_deg_target", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_deg_opt", envir=env, inherits=FALSE)){
        remove(".pcrsim_deg_opt", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_deg_spn", envir=env, inherits=FALSE)){
        remove(".pcrsim_deg_spn", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_deg_conc", envir=env, inherits=FALSE)){
        remove(".pcrsim_deg_conc", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_deg_size", envir=env, inherits=FALSE)){
        remove(".pcrsim_deg_size", envir=env, inherits=FALSE)
      }

      # TAB EXTRACTION:
      if(debug){
        print("TAB EXTRACTION")
      }
      if(exists(".pcrsim_ex_eff", envir=env, inherits=FALSE)){
        remove(".pcrsim_ex_eff", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_ex_eff_sd", envir=env, inherits=FALSE)){
        remove(".pcrsim_ex_eff_sd", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_ex_vol", envir=env, inherits=FALSE)){
        remove(".pcrsim_ex_vol", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_ex_vol_sd", envir=env, inherits=FALSE)){
        remove(".pcrsim_ex_vol_sd", envir=env, inherits=FALSE)
      }
      
      # TAB DILUTION:
      if(debug){
        print("TAB DILUTION")
      }
      if(exists(".pcrsim_dil_tolerance", envir=env, inherits=FALSE)){
        remove(".pcrsim_dil_tolerance", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_dil_target", envir=env, inherits=FALSE)){
        remove(".pcrsim_dil_target", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_dil_target_opt", envir=env, inherits=FALSE)){
        remove(".pcrsim_dil_target_opt", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_dil_accuracy", envir=env, inherits=FALSE)){
        remove(".pcrsim_dil_accuracy", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_dil_volume", envir=env, inherits=FALSE)){
        remove(".pcrsim_dil_volume", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_dil_volume_opt", envir=env, inherits=FALSE)){
        remove(".pcrsim_dil_volume_opt", envir=env, inherits=FALSE)
      }
      
      # TAB PCR AMPLIFICATION:
      if(debug){
        print("TAB PCR")
      }
      if(exists(".pcrsim_amp_aliquot", envir=env, inherits=FALSE)){
        remove(".pcrsim_amp_aliquot", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_amp_aliquot_sd", envir=env, inherits=FALSE)){
        remove(".pcrsim_amp_aliquot_sd", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_amp_amount", envir=env, inherits=FALSE)){
        remove(".pcrsim_amp_amount", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_amp_total_vol", envir=env, inherits=FALSE)){
        remove(".pcrsim_amp_total_vol", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_amp_total_vol_sd", envir=env, inherits=FALSE)){
        remove(".pcrsim_amp_total_vol_sd", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_amp_cyc_sb", envir=env, inherits=FALSE)){
        remove(".pcrsim_amp_cyc_sb", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_amp_eff_opt", envir=env, inherits=FALSE)){
        remove(".pcrsim_amp_eff_opt", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_amp_eff", envir=env, inherits=FALSE)){
        remove(".pcrsim_amp_eff", envir=env, inherits=FALSE)
      }
      
      # TAB CE ANALYSIS:
      if(debug){
        print("TAB CE")
      }
      if(exists(".pcrsim_ce_aliquot", envir=env, inherits=FALSE)){
        remove(".pcrsim_ce_aliquot", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_ce_aliquot_sd", envir=env, inherits=FALSE)){
        remove(".pcrsim_ce_aliquot_sd", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_ce_t_intercept", envir=env, inherits=FALSE)){
        remove(".pcrsim_ce_t_intercept", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_ce_t_slope", envir=env, inherits=FALSE)){
        remove(".pcrsim_ce_t_slope", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_ce_t_residual", envir=env, inherits=FALSE)){
        remove(".pcrsim_ce_t_residual", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_ce_intercept", envir=env, inherits=FALSE)){
        remove(".pcrsim_ce_intercept", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_ce_slope", envir=env, inherits=FALSE)){
        remove(".pcrsim_ce_slope", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_ce_residual", envir=env, inherits=FALSE)){
        remove(".pcrsim_ce_residual", envir=env, inherits=FALSE)
      }
      
      # TAB EPG:
      if(debug){
        print("TAB EPG")
      }
      if(exists(".pcrsim_epg_size", envir=env, inherits=FALSE)){
        remove(".pcrsim_epg_size", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_epg_vjust", envir=env, inherits=FALSE)){
        remove(".pcrsim_epg_vjust", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_epg_angle", envir=env, inherits=FALSE)){
        remove(".pcrsim_epg_angle", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_epg_hjust", envir=env, inherits=FALSE)){
        remove(".pcrsim_epg_hjust", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_epg_expand", envir=env, inherits=FALSE)){
        remove(".pcrsim_epg_expand", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_epg_scale", envir=env, inherits=FALSE)){
        remove(".pcrsim_epg_scale", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_epg_type", envir=env, inherits=FALSE)){
        remove(".pcrsim_epg_type", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_epg_fix", envir=env, inherits=FALSE)){
        remove(".pcrsim_epg_fix", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_epg_peaks", envir=env, inherits=FALSE)){
        remove(".pcrsim_epg_peaks", envir=env, inherits=FALSE)
      }
      
      # TAB SIMULATION:
      if(debug){
        print("TAB SIMULATION")
      }
      if(exists(".pcrsim_sim_sim", envir=env, inherits=FALSE)){
        remove(".pcrsim_sim_sim", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sim_name", envir=env, inherits=FALSE)){
        remove(".pcrsim_sim_name", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sim_title", envir=env, inherits=FALSE)){
        remove(".pcrsim_sim_title", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sim_rnd_chk", envir=env, inherits=FALSE)){
        remove(".pcrsim_sim_rnd_chk", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sim_rnd_edt", envir=env, inherits=FALSE)){
        remove(".pcrsim_sim_rnd_edt", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sim_update_chk", envir=env, inherits=FALSE)){
        remove(".pcrsim_sim_update_chk", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sim_ce_save_chk", envir=env, inherits=FALSE)){
        remove(".pcrsim_sim_ce_save_chk", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sim_ce_save_s_chk", envir=env, inherits=FALSE)){
        remove(".pcrsim_sim_ce_save_s_chk", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sim_ce_save_t_chk", envir=env, inherits=FALSE)){
        remove(".pcrsim_sim_ce_save_t_chk", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sim_ce_save_n_chk", envir=env, inherits=FALSE)){
        remove(".pcrsim_sim_ce_save_n_chk", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sim_profile_chk", envir=env, inherits=FALSE)){
        remove(".pcrsim_sim_profile_chk", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sim_sample_chk", envir=env, inherits=FALSE)){
        remove(".pcrsim_sim_sample_chk", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sim_extraction_chk", envir=env, inherits=FALSE)){
        remove(".pcrsim_sim_extraction_chk", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sim_degradation_chk", envir=env, inherits=FALSE)){
        remove(".pcrsim_sim_degradation_chk", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sim_dilution_chk", envir=env, inherits=FALSE)){
        remove(".pcrsim_sim_dilution_chk", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sim_pcr_chk", envir=env, inherits=FALSE)){
        remove(".pcrsim_sim_pcr_chk", envir=env, inherits=FALSE)
      }
      if(exists(".pcrsim_sim_ce_chk", envir=env, inherits=FALSE)){
        remove(".pcrsim_sim_ce_chk", envir=env, inherits=FALSE)
      }
      
      if(debug){
        print("Settings cleared!")
      }
    }
    
    if(debug){
      print("Settings saved!")
      print(paste("EXIT:", match.call()[[1]]))
    }
    
  }
  
  # END GUI ###################################################################
  
  
  # Show GUI, with tab one.
  svalue(nb) <- 1
  visible(w) <- TRUE
  message("PCRsim graphical user interface loaded!")
  
  return(NULL)
  
}
