################################################################################
# TODO LIST
# TODO: Object size not sorted correct (seem to sort as character)
# TODO: Migrate to gWidgets2.
# TODO: Save .importPath in ws for last used path (only in coming gWidgets2 ??)
# TODO: Multiple selection not working.
# TODO: USe viwweports instead of grid.arrange in complex plots?
# http://www.imachordata.com/extra-extra-get-your-gridextra/#comment-146

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

# Attributes:
# Put these in 'base functions':
# attr(dataDrop, which="[function], strvalidator") <- as.character(utils::packageVersion("strvalidator"))
# attr(dataDrop, which="[function], call") <- match.call()
# attr(at.rank, which="[function], date") <- date()
# Add additional attributes as:
# attr(datanew, which="[function], [attribute]") <- parameter
# Except for 'global' parameters used by other functions:
# attr(datanew, which="kit") <- kit

# NOTE:
# NOTE: Can't import data frame named 'drop'
# NOTE: Buttons named 'Plot' will show up 'plot'.
# NOTE: Some button names will change due to locale.

################################################################################
# CHANGE LOG (last 20 changes)
# 30.12.2015: Added button to new function 'calculateLb' in the 'Balance' tab.
# 22.12.2015: Added new group 'Marker ratio' in the 'Balance' tab.
# 18.12.2015: Added tooltips in the 'Workspace' tab.
# 30.11.2015: Moved updated description. Previous information moved to 'strvalidator-package'.
# 12.10.2015: Added 'Calculate' and 'Filter' button in 'Result' tab.
# 12.10.2015: Added new group 'Drop-in tools' in 'Result' tab.
# 29.08.2015: Added importFrom.
# 01.06.2015: Added 'Calculate' and 'Plot' (AT6) button in 'AT' tab.
# 24.05.2015: Added 'Columns' button in 'Tools' tab.
# 04.05.2015: Added 'AT' tab.
# 01.01.2015: Fixed error in 'Workspace' tab when no selection and 'Delete' is pressed.
# 19.12.2014: Added 'EPG' button in 'Tools' tab.
# 12.12.2014: Re-named 'Edit' tab to 'Tools' and change name on some buttons.
# 04.12.2014: Added 'Pull-up' tab.
# 28.10.2014: Fixed "Error in if (tabName == .file_tab_name) { : argument is of length zero"
# 06.10.2014: Added ggplot support for 'View' button in 'Workspace' tab.
# 28.08.2014: Fixed error message when no projects in projects folder.
# 08.07.2014: Added 'Mixture' tab.
# 04.07.2014: Ask to overwrite project file if exist.
# 04.07.2014: Added new button 'Add' and 'Save As' to 'Workspace' tab.

#' @title Graphical User Interface For The STR-validator Package
#'
#' @description
#' GUI simplifying the use of the STR-validator package.
#'
#' @details The graphical user interface give easy access to all graphical
#' versions of the functions available in the strvalidator package. It connects
#' functions 'under the hood' to allow a degree of automation not available
#' using the command based functions. In addition it provides a project based
#' workflow.\cr\cr
#' Click \code{Index} at the bottom of the help page to see a complete list
#' of functions.
#' 
#' @param debug logical indicating printing debug information.
#' 
#' @return TRUE
#' 
# @import ggplot2
#' @import gWidgets
#' @import gWidgetsRGtk2
#' @import RGtk2
# @import data.table
# @import gridExtra
#' @importFrom utils packageVersion help object.size
#' @importFrom graphics title
#' 
#' @export
#' 
#' @examples
#' # To start the graphical user interface.
#' \dontrun{
#' strvalidator()
#' }



strvalidator <- function(debug=FALSE){
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # Specify toolkit.
  options("guiToolkit"="RGtk2")
  
  # Global variables.
  .strvalidator_env <- new.env()
  .separator <- .Platform$file.sep # Platform dependent path separator.
  .save_gui <- TRUE
  .start_tab_name <- "Welcome"
  .file_tab_name <- "Workspace"
  .project_tab_name <- "Projects"
  .drylab_tab_name <- "DryLab"
  .edit_tab_name <- "Tools"
  .at_tab_name <- "AT"
  .stutter_tab_name <- "Stutter"
  .balance_tab_name <- "Balance"
  .concordance_tab_name <- "Concordance"
  .drop_tab_name <- "Dropout"
  .mixture_tab_name <- "Mixture"
  .result_tab_name <- "Result"
  .precision_tab_name <- "Precision"
  .pullup_tab_name <- "Pull-up"
  .object_classes_view <- c("data.frame", "ggplot")
  .object_classes_import <- c("data.frame", "ggplot")
  .project_description_variable <- ".strvalidator_project_description"
  .project_tmp_env <- new.env()
  .project_name_list <- NULL
  .project_path_list <- NULL
  .ws_name_variable <- ".strvalidator_project_name"
  .ws_path_variable <- ".strvalidator_project_path"
  
  
  # MAIN WINDOW  ##############################################################
  
  # Main window.
  w <- gwindow(title=paste("STR-validator",packageVersion("strvalidator"),
                           " - a forensic validation toolbox"),
               visible = FALSE,
               name=title)
  
  # Vertical main group.
  gv <- ggroup(horizontal=FALSE,
               use.scrollwindow=FALSE,
               container = w,
               expand=TRUE) 
  
  # Help button group.
  gh <- ggroup(container = gv, expand=FALSE, fill="both")
  
  savegui_chk <- gcheckbox(text="Save GUI settings", checked=TRUE, container=gh)
  
  addHandlerChanged(savegui_chk, handler = function(h, ...) {
    
    # Update variable.
    .save_gui <<- svalue(savegui_chk)
    
  })
  
  addSpring(gh)
  
  help_btn <- gbutton(text="Help", container=gh)
  
  addHandlerChanged(help_btn, handler = function(h, ...) {
    
    # Open help page for function.
    print(help("strvalidator", help_type="html"))
    
  })
  
  # Main client area.
  nb <- gnotebook(closebuttons = FALSE,
                  dontCloseThese = NULL,
                  container = gv)
  
  
  # NOTEBOOK ##################################################################
  
  # Define groups.
  start_tab <- ggroup(horizontal = FALSE,
                      spacing=10,
                      use.scrollwindow=FALSE,
                      container = nb,
                      label=.start_tab_name,
                      expand=TRUE)
  
  project_tab <- ggroup(horizontal = FALSE,
                        spacing=10,
                        use.scrollwindow=FALSE,
                        container = nb,
                        label=.project_tab_name,
                        expand=TRUE)
  
  file_tab <- ggroup(horizontal = FALSE,
                     spacing=10,
                     use.scrollwindow=FALSE,
                     container = nb,
                     label=.file_tab_name,
                     expand=TRUE)
  
  drylab_tab <- ggroup(horizontal = FALSE,
                       spacing=10,
                       use.scrollwindow=FALSE,
                       container = nb,
                       label=.drylab_tab_name,
                       expand=TRUE)
  
  edit_tab <- ggroup(horizontal = FALSE,
                     spacing=5,
                     use.scrollwindow=FALSE,
                     container = nb,
                     label=.edit_tab_name,
                     expand=TRUE)
  
  at_tab <- ggroup(horizontal = FALSE,
                   spacing=10,
                   use.scrollwindow=FALSE,
                   container = nb,
                   label=.at_tab_name,
                   expand=TRUE)
  
  stutter_tab <- ggroup(horizontal = FALSE,
                        spacing=10,
                        use.scrollwindow=FALSE,
                        container = nb,
                        label=.stutter_tab_name,
                        expand=TRUE)
  
  balance_tab <- ggroup(horizontal = FALSE,
                        spacing=10,
                        use.scrollwindow=FALSE,
                        container = nb,
                        label=.balance_tab_name,
                        expand=TRUE)
  
  concordance_tab <- ggroup(horizontal = FALSE,
                            spacing=10,
                            use.scrollwindow=FALSE,
                            container = nb,
                            label=.concordance_tab_name,
                            expand=TRUE)
  
  drop_tab <- ggroup(horizontal = FALSE,
                     spacing=10,
                     use.scrollwindow=FALSE,
                     container = nb,
                     label=.drop_tab_name,
                     expand=TRUE)
  
  mixture_tab <- ggroup(horizontal = FALSE,
                            spacing=10,
                            use.scrollwindow=FALSE,
                            container = nb,
                            label=.mixture_tab_name,
                            expand=TRUE)
  
  result_tab <- ggroup(horizontal = FALSE,
                       spacing=10,
                       use.scrollwindow=FALSE,
                       container = nb,
                       label=.result_tab_name,
                       expand=TRUE)
  
  precision_tab <- ggroup(horizontal = FALSE,
                          spacing=10,
                          use.scrollwindow=FALSE,
                          container = nb,
                          label=.precision_tab_name,
                          expand=TRUE)

  pullup_tab <- ggroup(horizontal = FALSE,
                          spacing=10,
                          use.scrollwindow=FALSE,
                          container = nb,
                          label=.pullup_tab_name,
                          expand=TRUE)
  
  # START #####################################################################
  
  glabel("", container=start_tab) # Adds some space.
  
  # STR TYPING KIT ------------------------------------------------------------
  
  about_txt <- paste("STR-validator (pronounced starvalidator) is a package ",
                     "developed for validation and process control of methods ",
                     "and instruments in a forensic genetic laboratory setting. ",
                     "This graphical user interface make it very easy to ",
                     "analyse validation data in accordance with ENFSI and SWGDAM ",
                     "guidelines.",
                     "The code has been extensively tested in order to assure correct results. ",
                     "\n\n",
                     "Created by:\n",
                     "Oskar Hansson, Department of Forensic Biology (NIPH, Norway)\n\n",
                     "General information and tutorials:\n",
                     "https://sites.google.com/site/forensicapps/strvalidator\n\n",
                     "Facebook:\n",
                     "https://www.facebook.com/pages/STR-validator/240891279451450?ref=tn_tnmn\n",
                     "https://www.facebook.com/groups/strvalidator/\n\n",
                     "Please report bugs to:\n",
                     "https://github.com/OskarHansson/strvalidator/issues\n\n",
                     "The source is hosted at GitHub:\n",
                     "https://github.com/OskarHansson/strvalidator", sep="")
  
  gtext(text=about_txt, width = NULL, height = NULL, font.attr = NULL, 
        wrap = TRUE, expand=TRUE, container = start_tab) 
  
  start_license_btn <- gbutton(text = "License", container = start_tab, expand=FALSE) 
  
  addHandlerChanged(start_license_btn, handler = function (h, ...) {
    
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
    
    gmessage(message = license_txt,
             title="License",
             icon = "info",
             parent = w) 
    
    
  } )
  
  # PROJECT MANAGER ###########################################################
  
  # Vertical main group.
  project_f1 <- ggroup(horizontal=FALSE,
                       container = project_tab,
                       spacing=10,
                       expand=TRUE)
  
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
  project_f2 <- gframe(text="Projects",
                       horizontal=TRUE,
                       spacing=10,
                       container = project_f1,
                       expand=TRUE)
  
  # Button group.
  project_g1 <- ggroup(horizontal=FALSE,
                       spacing=10,
                       container = project_f2,
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
  
  addHandlerChanged(project_open_btn, handler = function (h, ...) {
    
    # Get selected projects file name.
    val_name <- svalue(project_tbl)
    val_id <- as.numeric(project_tbl[svalue(project_tbl, index=TRUE), "Id"])
    val_prj <- .project_path_list[val_id]
    val_env <- .strvalidator_env
    
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
        
        # Move to workspace tab.
        svalue(nb)<- match(.file_tab_name, names(nb))
        
      }
      
    }
    
  } )
  
  addHandlerChanged(project_add_btn, handler = function (h, ...) {
    
    # Get selected projects file name.
    val_name <- svalue(project_tbl)
    val_id <- as.numeric(project_tbl[svalue(project_tbl, index=TRUE), "Id"])
    val_prj <- .project_path_list[val_id]
    val_env <- .strvalidator_env
    
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
  project_f3 <- gframe(text="Description",
                       horizontal=TRUE,
                       spacing=10,
                       container = project_f1,
                       expand=TRUE)
  
  # Button group.
  project_g3 <- ggroup(horizontal = FALSE, spacing = 10,
                       container = project_f3, expand = FALSE)
  
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
  project_g4 <- ggroup(horizontal=FALSE,
                       spacing=10,
                       container=project_f3,
                       expand=TRUE)
  
  # Project description window.
  proj_info_lbl <- glabel(text="Project:", anchor=c(-1 ,0),
                          container=project_g4)
  proj_info_txt <- gtext(text=NULL, height=300, expand=TRUE,
                         wrap=TRUE, container=project_g4) 
  
  # WORKSPACE #################################################################
  
  # LOADED DATASETS -----------------------------------------------------------
  
  workspace_f1 <- gframe(text = "Project",
                         markup = FALSE,
                         pos = 0,
                         horizontal=TRUE,
                         container = file_tab,
                         expand=TRUE)
  
  workspace_f1g1 <- ggroup(horizontal=FALSE,
                           container = workspace_f1,
                           expand=FALSE)
  
  ws_open_btn <- gbutton(text="Open",
                         border=TRUE,
                         container = workspace_f1g1)
  tooltip(ws_open_btn) <- "Open project"
  
  ws_save_btn <- gbutton(text="Save",
                         border=TRUE,
                         container = workspace_f1g1)
  tooltip(ws_save_btn) <- "Save project"
  
  ws_saveas_btn <- gbutton(text="Save As",
                           border=TRUE,
                           container = workspace_f1g1)
  tooltip(ws_saveas_btn) <- "Choose a location and save project"
  
  ws_import_btn <- gbutton(text="Import",
                           border=TRUE,
                           container = workspace_f1g1)
  tooltip(ws_import_btn) <- "Import data from file"
  
  ws_export_btn <- gbutton(text="Export",
                           border=TRUE,
                           container = workspace_f1g1)
  tooltip(ws_export_btn) <- "Open the export dialoge"
  
  ws_add_btn <- gbutton(text="Add",
                        border=TRUE,
                        container = workspace_f1g1)
  tooltip(ws_add_btn) <- "Merge a project with the current project"
  
  ws_refresh_btn <- gbutton(text="Refresh",
                            border=TRUE,
                            container = workspace_f1g1) 
  tooltip(ws_refresh_btn) <- "Refresh the workspace"
  
  ws_remove_btn <- gbutton(text="Delete",
                           border=TRUE,
                           container = workspace_f1g1) 
  tooltip(ws_remove_btn) <- "Delete selected object"
  
  ws_rename_btn <- gbutton(text="Rename",
                           border=TRUE,
                           container = workspace_f1g1)
  tooltip(ws_rename_btn) <- "Rename selected object"
  
  ws_view_btn <- gbutton(text="View",
                         border=TRUE,
                         container = workspace_f1g1)
  tooltip(ws_view_btn) <- "View selected object"
  
  ws_loaded_tbl <- gWidgets::gtable(items=data.frame(Object="", Size="",
                                                     stringsAsFactors=FALSE), 
                                    multiple = TRUE,
                                    chosencol = 1,
                                    expand = TRUE,
                                    container = workspace_f1) 
  
  
  addHandlerChanged(ws_rename_btn, handler = function (h, ...) {
    
    
    objectName <- svalue(ws_loaded_tbl)
    
    if(length(objectName) == 1){
      
      # Get the object to save.
      datanew <- get(objectName, envir = .strvalidator_env)
        
      # Save data.
      saveObject(name=NULL, object=datanew, suggest=objectName,
                 parent=w, remove=objectName, env=.strvalidator_env,
                 debug=debug)
      
      .refreshLoaded()
      
    } else {
      gmessage(message="Currently you can only rename one object at a time!",
               title="Error",
               icon = "error",
               parent = w) 
    }
    
  } )
  
  addHandlerChanged(ws_open_btn, handler = function (h, ...) {
    
    val_env <- .strvalidator_env
    ws_path <- gfile(text="Select a saved workspace or dataset", type="open",
                     filter = list("R files" = list(patterns = c("*.R","*.Rdata"))),
                     multi=FALSE)
    
    if(!is.na(ws_path)){
      if(file.exists(ws_path)){
        
        # Clear environment.
        remove(list=ls(envir=val_env, all.names=TRUE),
               envir=val_env, inherits=FALSE)
        
        # Load new project.
        load(file=ws_path, envir = .strvalidator_env)
        .refreshLoaded()
        .loadSavedSettings()
        
      } else {
        
        gmessage(message="The workspace file was not found",
                 title="File not found",
                 icon = "error",
                 parent = w) 
      }
    }    
    
  } )
  
  addHandlerChanged(ws_add_btn, handler = function (h, ...) {
    
    ws_path <- gfile(text="Select a saved workspace or dataset", type="open",
                     filter = list("R files" = list(patterns = c("*.R","*.Rdata"))),
                     multi=FALSE)
    
    if(!is.na(ws_path)){
      if(file.exists(ws_path)){
        
        # Add new project.
        load(file=ws_path, envir = .strvalidator_env)
        .refreshLoaded()
        .loadSavedSettings()
        
      } else {
        
        gmessage(message="The workspace file was not found",
                 title="File not found",
                 icon = "error",
                 parent = w) 
      }
    }    
    
  } )
  
  addHandlerChanged(ws_import_btn, handler = function (h, ...) {
    
    # Open GUI.
    import_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    .refreshLoaded()
    
  } )  
  
  addHandlerChanged(ws_export_btn, handler = function (h, ...) {
    
    # Open GUI.
    export_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )  
  
  addHandlerChanged(ws_refresh_btn, handler = function (h, ...) {
    
    .refreshLoaded()
    
  })
  
  addHandlerChanged(ws_view_btn, handler = function(h, ...) {
    
    # Get selected dataset name(s).
    val_obj <- svalue(ws_loaded_tbl)
    val_data <- get(val_obj, envir=.strvalidator_env)
    val_class <- class(val_data)
    
    if(debug){
      print(paste("IN:", match.call()[[1]]))
      print("Changed, ws_view_btn")
      print(val_obj)
    }
    
    if (!is.null(val_obj) && !is.na(val_obj) && length(val_obj) > 0){
      
      if("data.frame" %in% val_class){
        
        # Open GUI.
        editData_gui(env=.strvalidator_env,
                     savegui=.save_gui,
                     data=get(val_obj, envir=.strvalidator_env),
                     name=val_obj,
                     edit=FALSE, debug=debug, parent=w)
        
      } else if("ggplot" %in% val_class){
        
        # Plot object.
        print(val_data)
        
      } else {
        message(paste("Object of type", val_class, "not supported!"))
      }
      
    }
    
  } )
  
  addHandlerChanged(ws_remove_btn, handler = function(h, ...) {
    
    # Get selected dataset name(s).
    val_obj <- svalue(ws_loaded_tbl)
    
    if(debug){
      print(paste("IN:", match.call()[[1]]))
      print("Changed, ws_remove_btn")
      print("Removed:")
      print(val_obj)
    }

    if(length(val_obj) == 1){
      
      if (!is.null(val_obj) && !is.na(val_obj)){
        
        # Get active reference data frame.
        remove(list=val_obj, envir=.strvalidator_env)
        
        .refreshLoaded()
        
        
      }
      
    } else if(length(val_obj) == 0) {
      gmessage(message="No object selected!",
               title="Error",
               icon = "error",
               parent = w)
      
    } else {
      gmessage(message="Currently you can only remove one object at a time!",
               title="Error",
               icon = "error",
               parent = w) 
    }
    
  } )
  
  
  addHandlerChanged(ws_save_btn, handler = function (h, ...) {
    
    # Initiate flag.
    ok <- TRUE
    
    # Get project name if available.
    if(exists(.ws_name_variable, envir=.strvalidator_env)){
      ws_name <- get(.ws_name_variable, envir=.strvalidator_env,
                     inherits=FALSE)
    } else {
      ok <- FALSE
    }
    
    # Get project path if available.
    if(exists(.ws_path_variable, envir=.strvalidator_env)){
      ws_save_path <- get(.ws_path_variable, envir=.strvalidator_env,
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
               list=ls(envir = .strvalidator_env, all.names = TRUE),
               envir = .strvalidator_env)
          
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
  
  addHandlerChanged(ws_saveas_btn, handler = function (h, ...) {
    
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
          assign(x=.ws_name_variable, value=ws_name, envir=.strvalidator_env)
          assign(x=.ws_path_variable, value=ws_save_path, envir=.strvalidator_env)
          
          # Save settings.        
          .saveSettings()
          
          # Save project. 
          save(file=ws_full_name, 
               list=ls(envir = .strvalidator_env, all.names = TRUE),
               envir = .strvalidator_env)
          
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
  
  
  # DATASETS ------------------------------------------------------------------  
  
  workspace_f2 <- gframe(text = "Load objects from R workspace",
                         markup = FALSE,
                         pos = 0,
                         horizontal=TRUE,
                         container = file_tab,
                         expand=FALSE)
  
  workspace_f2g1 <- ggroup(horizontal=FALSE,
                           container = workspace_f2,
                           
                           expand=FALSE)
  
  glabel("", container=workspace_f2g1) # Adds some space.
  
  ws_r_refresh_btn <- gbutton(text="Refresh dropdown",
                              border=TRUE,
                              container = workspace_f2g1) 
  
  
  ws_r_load_btn <- gbutton(text="Load object",
                           border=TRUE,
                           container = workspace_f2g1) 
  
  
  glabel("", container=workspace_f2g1) # Adds some space.
  
  ws_r_drp <- gdroplist(items=c("<Select object>", 
                                listObjects(env=.strvalidator_env,
                                            obj.class=.object_classes_import)), 
                        selected = 1,
                        editable = FALSE,
                        container = workspace_f2g1) 
  
  addHandlerChanged(ws_r_refresh_btn, handler = function (h, ...) {
    
    .refreshWs()
  } )
  
  addHandlerChanged(ws_r_load_btn, handler = function(h, ...) {
    
    # Get selected dataset name.
    val_name <- svalue(ws_r_drp)
    
    if (!is.na(val_name) && !is.null(val_name)){
      
      # Load dataset.
      saveObject(name=val_name, object=get(val_name),
                 parent=w, env=.strvalidator_env, debug=debug)
      
      # Update list.
      .refreshLoaded()
      
    } 
  } )
  
  
  # STR TYPING KIT ------------------------------------------------------------
  
  # DRY LAB  ##################################################################
  
  glabel("", container=drylab_tab) # Adds some space.
  
  
  dry_grid <- glayout(container = drylab_tab)
  
  # VIEW/EDIT -----------------------------------------------------------------
  
  dry_grid[1,1] <- dry_view_btn <- gbutton(text="Edit",
                                           border=TRUE,
                                           container = dry_grid) 
  
  dry_grid[1,2] <- glabel(text="Edit or view a dataset.",
                          container=dry_grid,
                          anchor=c(-1 ,0))
  
  addHandlerChanged(dry_view_btn, handler = function(h, ...) {
    
    # Open GUI.
    editData_gui(env=.strvalidator_env, savegui=.save_gui,
                 edit=TRUE, debug=debug, parent=w)
    
  } )
  
  # MAKE KIT ------------------------------------------------------------------
  
  dry_grid[2,1] <- dry_kit_btn <- gbutton(text="Kits",
                                          border=TRUE,
                                          container = dry_grid) 
  
  dry_grid[2,2] <- glabel(text="Add new kits or edit kits file.",
                          container=dry_grid,
                          anchor=c(-1 ,0))
  
  dry_grid[3,1] <- dry_plot_kit_btn <- gbutton(text="Plot Kit",
                                               border=TRUE,
                                               container = dry_grid) 
  
  dry_grid[3,2] <- glabel(text="Plot marker ranges for kits.",
                          container=dry_grid,
                          anchor=c(-1 ,0))
  
  dry_grid[4,1] <- dry_bins_btn <- gbutton(text="Analyse Overlap",
                                           border=TRUE,
                                           container = dry_grid) 
  
  dry_grid[4,2] <- glabel(text="Compare bins overlap for kits.",
                          container=dry_grid,
                          anchor=c(-1 ,0))
  
  dry_grid[5,1] <- dry_ol_btn <- gbutton(text="Analyse OL",
                                         border=TRUE,
                                         container = dry_grid) 
  
  dry_grid[5,2] <- glabel(text="Compare risk of getting off-ladder alleles for kits.",
                          container=dry_grid,
                          anchor=c(-1 ,0))
  
  addHandlerChanged(dry_kit_btn, handler = function(h, ...) {
    
    # Open GUI.
    makeKit_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  addHandlerChanged(dry_plot_kit_btn, handler = function(h, ...) {
    
    # Open GUI.
    plotKit_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  addHandlerChanged(dry_bins_btn, handler = function(h, ...) {
    
    # Open GUI.
    calculateOverlap_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  addHandlerChanged(dry_ol_btn, handler = function(h, ...) {
    
    # Open GUI.
    calculateOL_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # EDIT  #####################################################################
  
  glabel("", container=edit_tab) # Adds some space.
  
  
  edit_grid <- glayout(container = edit_tab)
  
  # VIEW/EDIT -----------------------------------------------------------------
  
  edit_grid[1,1] <- edit_view_btn <- gbutton(text="Edit",
                                             border=TRUE,
                                             container = edit_grid) 
  
  edit_grid[1,2] <- glabel(text="Edit or view a dataset.",
                           container=edit_grid,
                           anchor=c(-1 ,0))
  
  addHandlerChanged(edit_view_btn, handler = function(h, ...) {
    
    # Open GUI.
    editData_gui(env=.strvalidator_env, savegui=.save_gui,
                 edit=TRUE, debug=debug, parent=w)
    
  } )
  
  # TRIM ----------------------------------------------------------------------
  
  edit_grid[2,1] <- edit_trim_btn <- gbutton(text="Trim",
                                             border=TRUE,
                                             container = edit_grid) 
  
  edit_grid[2,2] <- glabel(text="Trim/discard samples or columns from a dataset.",
                           container=edit_grid,
                           anchor=c(-1 ,0))
  
  
  addHandlerChanged(edit_trim_btn, handler = function(h, ...) {
    
    # Open GUI.
    trim_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # SLIM ----------------------------------------------------------------------
  
  edit_grid[3,1] <- edit_slim_btn <- gbutton(text="Slim",
                                             border=TRUE,
                                             container = edit_grid) 
  
  edit_grid[3,2] <- glabel(text="Slim a dataset to 'long' format.",
                           container=edit_grid,
                           anchor=c(-1 ,0))
  
  
  addHandlerChanged(edit_slim_btn, handler = function(h, ...) {
    
    # Open GUI.
    slim_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # FILTER --------------------------------------------------------------------
  
  edit_grid[4,1] <- edit_filter_btn <- gbutton(text="Filter",
                                               border=TRUE,
                                               container = edit_grid) 
  
  edit_grid[4,2] <- glabel(text="Filter a dataset using a reference set.",
                           container=edit_grid,
                           anchor=c(-1 ,0))
  
  
  addHandlerChanged(edit_filter_btn, handler = function(h, ...) {
    
    filterProfile_gui(env=.strvalidator_env, savegui=.save_gui, parent=w)
    
  } )
  
  # CROP ----------------------------------------------------------------------
  
  edit_grid[5,1] <- edit_crop_btn <- gbutton(text="Crop",
                                             border=TRUE,
                                             container = edit_grid) 
  
  edit_grid[5,2] <- glabel(text="Discard, or replace data.",
                           container=edit_grid,
                           anchor=c(-1 ,0))
  
  
  addHandlerChanged(edit_crop_btn, handler = function(h, ...) {
    
    # Open GUI.
    cropData_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # GUESS ---------------------------------------------------------------------
  
  edit_grid[6,1] <- edit_guess_btn <- gbutton(text="Guess",
                                              border=TRUE,
                                              container = edit_grid) 
  
  edit_grid[6,2] <- glabel(text="Guess the profile from raw DNA result.",
                           container=edit_grid,
                           anchor=c(-1 ,0))
  
  
  addHandlerChanged(edit_guess_btn, handler = function(h, ...) {
    
    guessProfile_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # DYE -----------------------------------------------------------------------
  
  edit_grid[7,1] <- edit_addDye_btn <- gbutton(text="Dye",
                                               border=TRUE,
                                               container = edit_grid) 
  
  edit_grid[7,2] <- glabel(text="Add dye information according to kit.",
                           container=edit_grid,
                           anchor=c(-1 ,0))
  
  addHandlerChanged(edit_addDye_btn, handler = function(h, ...) {
    
    # Open GUI.
    addDye_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # ADD MARKER ----------------------------------------------------------------
  
  edit_grid[8,1] <- edit_addMarker_btn <- gbutton(text="Marker",
                                                  border=TRUE,
                                                  container = edit_grid) 
  
  edit_grid[8,2] <- glabel(text="Add missing markers to dataset.",
                           container=edit_grid,
                           anchor=c(-1 ,0))
  
  addHandlerChanged(edit_addMarker_btn, handler = function(h, ...) {
    
    # Open GUI.
    addMarker_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # ADD SIZE ------------------------------------------------------------------
  
  edit_grid[9,1] <- edit_addSize_btn <- gbutton(text="Size",
                                                border=TRUE,
                                                container = edit_grid) 
  
  edit_grid[9,2] <- glabel(text="Add approximate size to alleles in a dataset.",
                           container=edit_grid,
                           anchor=c(-1 ,0))
  
  addHandlerChanged(edit_addSize_btn, handler = function(h, ...) {
    
    # Open GUI.
    addSize_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # ADD DATA -------------------------------------------------------------------
  
  edit_grid[10,1] <- edit_addData_btn <- gbutton(text="Data",
                                                 border=TRUE,
                                                 container = edit_grid) 
  
  edit_grid[10,2] <- glabel(text="Add new information to a dataset.",
                            container=edit_grid,
                            anchor=c(-1 ,0))
  
  addHandlerChanged(edit_addData_btn, handler = function(h, ...) {
    
    # Open GUI.
    addData_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # CHECK SUBSET --------------------------------------------------------------
  
  edit_grid[11,1] <- edit_check_btn <- gbutton(text="Check",
                                               border=TRUE,
                                               container = edit_grid) 
  
  edit_grid[11,2] <- glabel(text="Check the subsetting of a dataset.",
                            container=edit_grid,
                            anchor=c(-1 ,0))
  
  addHandlerChanged(edit_check_btn, handler = function(h, ...) {
    
    # Open GUI.
    checkSubset_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # COMBINE -------------------------------------------------------------------
  
  edit_grid[12,1] <- edit_combine_btn <- gbutton(text="Combine",
                                              border=TRUE,
                                              container = edit_grid) 
  
  edit_grid[12,2] <- glabel(text="Combine two datasets.",
                            container=edit_grid,
                            anchor=c(-1 ,0))
  
  addHandlerChanged(edit_combine_btn, handler = function(h, ...) {
    
    # Open GUI.
    combine_gui(env=.strvalidator_env, debug=debug, parent=w)
    
  } )

  # COLUMNS -------------------------------------------------------------------
  
  edit_grid[13,1] <- edit_columns_btn <- gbutton(text="Columns",
                                                 border=TRUE,
                                                 container = edit_grid) 
  
  edit_grid[13,2] <- glabel(text="Perform actions on columns.",
                            container=edit_grid,
                            anchor=c(-1 ,0))
  
  addHandlerChanged(edit_columns_btn, handler = function(h, ...) {
    
    # Open GUI.
    columns_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # CALCULATE HETEROZYGOUS ----------------------------------------------------
  
  edit_grid[14,1] <- edit_het_btn <- gbutton(text="Heterozygous",
                                             border=TRUE,
                                             container = edit_grid) 
  
  edit_grid[14,2] <- glabel(text="Indicate heterozygous loci for a reference dataset.",
                            container=edit_grid,
                            anchor=c(-1 ,0))
  
  addHandlerChanged(edit_het_btn, handler = function(h, ...) {
    
    # Open GUI.
    calculateHeterozygous_gui(env=.strvalidator_env, debug=debug, parent=w)
    
  } )
  
  # CALCULATE H ---------------------------------------------------------------
  
  edit_grid[15,1] <- edit_h_btn <- gbutton(text="Height",
                                           border=TRUE,
                                           container = edit_grid) 
  
  edit_grid[15,2] <- glabel(text="Calculate peak height metrics.",
                            container=edit_grid,
                            anchor=c(-1 ,0))
  
  addHandlerChanged(edit_h_btn, handler = function(h, ...) {
    
    # Open GUI.
    calculateHeight_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )

  # GENERATE EPG --------------------------------------------------------------
  
  edit_grid[16,1] <- edit_epg_btn <- gbutton(text="EPG",
                                           border=TRUE,
                                           container = edit_grid) 
  
  edit_grid[16,2] <- glabel(text="Generate EPG like plot.",
                            container=edit_grid,
                            anchor=c(-1 ,0))
  
  addHandlerChanged(edit_epg_btn, handler = function(h, ...) {
    
    # Open GUI.
    generateEPG_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )

  # AT  #######################################################################
  
  
  glabel("", container=at_tab) # Adds some space.
  
  at_grid <- glayout(container = at_tab)
  
  
  # VIEW/EDIT -----------------------------------------------------------------
  
  at_grid[1,1] <- at_view_btn <- gbutton(text="Edit", border=TRUE,
                                         container = at_grid) 
  
  at_grid[1,2] <- glabel(text="Edit or view a dataset.",
                              container=at_grid, anchor=c(-1 ,0))
  
  addHandlerChanged(at_view_btn, handler = function(h, ...) {
    
    # Open GUI.
    editData_gui(env=.strvalidator_env, savegui=.save_gui
                 , edit=TRUE, debug=debug, parent=w)
    
  } )
  
  # CALCULATE -----------------------------------------------------------------
  
  at_grid[3,1] <- at_calculate_btn <- gbutton(text="Calculate", border=TRUE,
                                              container = at_grid) 
  
  at_grid[3,2] <- glabel(text="Calculate analytical threshold (AT1, AT2, AT4).",
                              container=at_grid, anchor=c(-1 ,0))
  
  
  addHandlerChanged(at_calculate_btn, handler = function(h, ...) {
    
    # Open GUI.
    calculateAT_gui(env=.strvalidator_env, savegui=.save_gui,
                    debug=debug, parent=w)
    
  } )

  # CALCULATE -----------------------------------------------------------------
  
  at_grid[4,1] <- at6_calculate_btn <- gbutton(text="Calculate", border=TRUE,
                                              container = at_grid) 
  
  at_grid[4,2] <- glabel(text="Calculate analytical threshold (AT6).",
                         container=at_grid, anchor=c(-1 ,0))
  
  addHandlerChanged(at6_calculate_btn, handler = function(h, ...) {
    
    # Open GUI.
    calculateAT6_gui(env=.strvalidator_env, savegui=.save_gui,
                    debug=debug, parent=w)
    
  } )
  
  # PLOT AT -------------------------------------------------------------------
  
  at_grid[5,1] <- at_plot_btn <- gbutton(text="Plot", border=TRUE,
                                         container = at_grid) 
  
  at_grid[5,2] <- glabel(text="Create plots for analysed data (AT6).",
                         container=at_grid)
  
  addHandlerChanged(at_plot_btn, handler = function(h, ...) {
    
    # Open GUI.
    plotAT_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # STUTTER  ##################################################################
  
  
  glabel("", container=stutter_tab) # Adds some space.
  
  stutter_grid <- glayout(container = stutter_tab)
  
  
  # VIEW/EDIT -----------------------------------------------------------------
  
  stutter_grid[1,1] <- stutter_view_btn <- gbutton(text="Edit",
                                                   border=TRUE,
                                                   container = stutter_grid) 
  
  stutter_grid[1,2] <- glabel(text="Edit or view a dataset.",
                              container=stutter_grid,
                              anchor=c(-1 ,0))
  
  addHandlerChanged(stutter_view_btn, handler = function(h, ...) {
    
    # Open GUI.
    editData_gui(env=.strvalidator_env, savegui=.save_gui,
                 edit=TRUE, debug=debug, parent=w)
    
  } )
  
  # CALCULATE -----------------------------------------------------------------
  
  stutter_grid[3,1] <- stutter_calculate_btn <- gbutton(text="Calculate",
                                                        border=TRUE,
                                                        container = stutter_grid) 
  
  stutter_grid[3,2] <- glabel(text="Calculate stutters for a dataset.",
                              container=stutter_grid,
                              anchor=c(-1 ,0))
  
  
  addHandlerChanged(stutter_calculate_btn, handler = function(h, ...) {
    
    # Open GUI.
    calculateStutter_gui(env=.strvalidator_env, savegui=.save_gui,
                         debug=debug, parent=w)
    
  } )
  
  # PLOT STUTTER --------------------------------------------------------------
  
  stutter_grid[4,1] <- stutter_plot_btn <- gbutton(text="Plot",
                                                   border=TRUE,
                                                   container = stutter_grid) 
  
  stutter_grid[4,2] <- glabel(text="Create plots for analysed data.",
                              container=stutter_grid)
  
  addHandlerChanged(stutter_plot_btn, handler = function(h, ...) {
    
    # Open GUI.
    plotStutter_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  
  # SUMMARY TABLE -------------------------------------------------------------
  
  stutter_grid[5,1] <- stutter_table_btn <- gbutton(text="Summarize",
                                                    border=TRUE,
                                                    container = stutter_grid) 
  
  stutter_grid[5,2] <- glabel(text="Summarize stutter data in a table.",
                              container=stutter_grid)
  
  addHandlerChanged(stutter_table_btn, handler = function(h, ...) {
    
    # Open GUI.
    tableStutter_gui(env=.strvalidator_env, savegui=.save_gui,
                     debug=debug, parent=w)
    
  } )
  
  
  # BALANCE  ##################################################################
  
  glabel("", container=balance_tab) # Adds some space.
  
  balance_g1 <- glayout(container = balance_tab)
  
  # VIEW/EDIT -----------------------------------------------------------------
  
  balance_g1[1,1] <- balance_view_btn <- gbutton(text="Edit",
                                                 border=TRUE,
                                                 container = balance_g1) 
  
  balance_g1[1,2] <- glabel(text="Edit or view a dataset.",
                            container=balance_g1,
                            anchor=c(-1 ,0))
  
  addHandlerChanged(balance_view_btn, handler = function(h, ...) {
    
    # Open GUI.
    editData_gui(env=.strvalidator_env, savegui=.save_gui,
                 edit=TRUE, debug=debug, parent=w)
    
  } )
  
  # ALLELE BALANCE ============================================================
  
  balance_f2 <- gframe(text = "Intralocus and interlocus balance",
                       horizontal=FALSE, container = balance_tab) 
  
  balance_g2 <- glayout(container = balance_f2)
  
  # CALCULATE -----------------------------------------------------------------

  # FUNCTION 1.  
  balance_g2[1,1] <- balance_g2_calc_1_btn <- gbutton(text="Calculate",
                                                    border=TRUE,
                                                    container = balance_g2) 
  
  balance_g2[1,2] <- glabel(text="Calculate intra/inter locus balance for a dataset (reference required).",
                            container=balance_g2)
  
  
  addHandlerChanged(balance_g2_calc_1_btn, handler = function(h, ...) {
    
    # Open GUI.
    calculateBalance_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # FUNCTION 2.
  balance_g2[2,1] <- balance_g2_calc_2_btn <- gbutton(text="Calculate",
                                                    border=TRUE,
                                                    container = balance_g2) 
  
  balance_g2[2,2] <- glabel(text="Calculate inter locus balance for a dataset (no reference required).",
                            container=balance_g2)
  
  
  addHandlerChanged(balance_g2_calc_2_btn, handler = function(h, ...) {
    
    # Open GUI.
    calculateLb_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # PLOT ----------------------------------------------------------------------
  
  balance_g2[3,1] <- balance_g2_plot_btn <- gbutton(text="Plot",
                                                    border=TRUE,
                                                    container = balance_g2) 
  
  balance_g2[3,2] <- glabel(text="Create plots for analysed data",
                            container=balance_g2)
  
  addHandlerChanged(balance_g2_plot_btn, handler = function(h, ...) {
    
    # Open GUI.
    plotBalance_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # SUMMARY TABLE -------------------------------------------------------------
  
  balance_g2[4,1] <- balance_table_btn <- gbutton(text="Summarize",
                                                  border=TRUE,
                                                  container = balance_g2) 
  
  balance_g2[4,2] <- glabel(text="Summarize balance data in a table.",
                            container=balance_g2)
  
  addHandlerChanged(balance_table_btn, handler = function(h, ...) {
    
    # Open GUI.
    tableBalance_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # CAPILLARY BALANCE =========================================================
  
  balance_f3 <- gframe(text = "Capillary balance",
                       horizontal=FALSE, container = balance_tab) 
  
  balance_g3 <- glayout(container = balance_f3)
  
  
  # CALCULATE -----------------------------------------------------------------
  
  balance_g3[1,1] <- balance_g3_calc_btn <- gbutton(text="Calculate",
                                                    border=TRUE,
                                                    container = balance_g3) 
  
  balance_g3[1,2] <- glabel(text="Calculate capillary balance for a dataset.",
                            container=balance_g3)
  
  
  addHandlerChanged(balance_g3_calc_btn, handler = function(h, ...) {
    
    # Open GUI.
    calculateCapillary_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # PLOT ----------------------------------------------------------------------
  
  balance_g3[2,1] <- balance_g3_plot_btn <- gbutton(text="Plot",
                                                    border=TRUE,
                                                    container = balance_g3) 
  
  balance_g3[2,2] <- glabel(text="Create plots for analysed data",
                            container=balance_g3)
  
  addHandlerChanged(balance_g3_plot_btn, handler = function(h, ...) {
    
    # Open GUI.
    plotCapillary_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # SUMMARY -------------------------------------------------------------------
  
  balance_g3[3,1] <- balance_g3_tab_btn <- gbutton(text="Summarize",
                                                   border=TRUE,
                                                   container = balance_g3) 
  
  balance_g3[3,2] <- glabel(text="Create summary table for analysed data",
                            container=balance_g3)
  
  addHandlerChanged(balance_g3_tab_btn, handler = function(h, ...) {
    
    # Open GUI.
    tableCapillary_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )

  # MARKER RATIO ==============================================================
  
  balance_f4 <- gframe(text = "Marker peak height ratio",
                       horizontal=FALSE, container = balance_tab) 
  
  balance_g4 <- glayout(container = balance_f4)
  
  # CALCULATE -----------------------------------------------------------------
  
  balance_g4[1,1] <- balance_g4_calc_btn <- gbutton(text="Calculate",
                                                    border=TRUE,
                                                    container = balance_g4) 
  
  balance_g4[1,2] <- glabel(text="Calculate locus ratio for a dataset.",
                            container=balance_g4)
  
  
  addHandlerChanged(balance_g4_calc_btn, handler = function(h, ...) {
    
    # Open GUI.
    calculateRatio_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # PLOT ----------------------------------------------------------------------
  
  balance_g4[2,1] <- balance_g4_plot_btn <- gbutton(text="Plot",
                                                    border=TRUE,
                                                    container = balance_g4)
  
  balance_g4[2,2] <- glabel(text="Create plots for analysed data",
                            container=balance_g4)
  
  addHandlerChanged(balance_g4_plot_btn, handler = function(h, ...) {
    
    # Open GUI.
    plotRatio_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # CONCORDANCE  ##############################################################
  
  glabel("", container=concordance_tab) # Adds some space.
  
  conc_grid <- glayout(container = concordance_tab)
  
  # VIEW/EDIT -----------------------------------------------------------------
  
  conc_grid[1,1] <- conc_view_btn <- gbutton(text="Edit",
                                             border=TRUE,
                                             container = conc_grid) 
  
  conc_grid[1,2] <- glabel(text="Edit or view a dataset.",
                           container=conc_grid,
                           anchor=c(-1 ,0))
  
  addHandlerChanged(conc_view_btn, handler = function(h, ...) {
    
    # Open GUI.
    editData_gui(env=.strvalidator_env, savegui=.save_gui,
                 edit=TRUE, debug=debug, parent=w)
    
  } )
  
  # CALCULATE -----------------------------------------------------------------
  
  conc_grid[2,1] <- conc_calculate_btn <- gbutton(text="Calculate",
                                                  border=TRUE,
                                                  container = conc_grid) 
  
  conc_grid[2,2] <- glabel(text="Calculate concordance for multiple datasets.",
                           container=conc_grid,
                           anchor=c(-1 ,0))
  
  
  addHandlerChanged(conc_calculate_btn, handler = function(h, ...) {
    
    # Open GUI.
    calculateConcordance_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # DROPOUT  ##################################################################
  
  glabel("", container=drop_tab) # Adds some space.
  
  drop_grid <- glayout(container = drop_tab)
  
  # VIEW/EDIT -----------------------------------------------------------------
  
  drop_grid[1,1] <- drop_view_btn <- gbutton(text="Edit",
                                             border=TRUE,
                                             container = drop_grid) 
  
  drop_grid[1,2] <- glabel(text="Edit or view a dataset.",
                           container=drop_grid,
                           anchor=c(-1 ,0))
  
  addHandlerChanged(drop_view_btn, handler = function(h, ...) {
    
    # Open GUI.
    editData_gui(env=.strvalidator_env, savegui=.save_gui,
                 edit=TRUE, debug=debug, parent=w)
    
  } )
  
  # CALCULATE -----------------------------------------------------------------
  
  drop_grid[2,1] <- drop_calculate_btn <- gbutton(text="Calculate",
                                                  border=TRUE,
                                                  container = drop_grid) 
  
  drop_grid[2,2] <- glabel(text="Calculate dropouts for a dataset.",
                           container=drop_grid,
                           anchor=c(-1 ,0))
  
  
  addHandlerChanged(drop_calculate_btn, handler = function(h, ...) {
    
    # Open GUI.
    calculateDropout_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # LOGISTIC REGRESSION -------------------------------------------------------
  
  drop_grid[3,1] <- drop_model_btn <- gbutton(text="Model",
                                              border=TRUE,
                                              container = drop_grid) 
  
  drop_grid[3,2] <- glabel(text="Model dropout risk",
                           container=drop_grid)
  
  addHandlerChanged(drop_model_btn, handler = function(h, ...) {
    
    # Open GUI.
    modelDropout_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # PLOT DROPOUT --------------------------------------------------------------
  
  drop_grid[4,1] <- drop_plot_btn <- gbutton(text="Plot",
                                             border=TRUE,
                                             container = drop_grid) 
  
  drop_grid[4,2] <- glabel(text="Create plots for analysed data",
                           container=drop_grid)
  
  addHandlerChanged(drop_plot_btn, handler = function(h, ...) {
    
    # Open GUI.
    plotDropout_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  
  # SUMMARY TABLE -------------------------------------------------------------
  
  # MIXTURE  ##################################################################
  
  glabel("", container=mixture_tab) # Adds some space.
  
  mix_grid <- glayout(container = mixture_tab)
  
  # VIEW/EDIT -----------------------------------------------------------------
  
  mix_grid[1,1] <- mix_view_btn <- gbutton(text="Edit",
                                             border=TRUE,
                                             container = mix_grid) 
  
  mix_grid[1,2] <- glabel(text="Edit or view a dataset.",
                           container=mix_grid,
                           anchor=c(-1 ,0))
  
  addHandlerChanged(mix_view_btn, handler = function(h, ...) {
    
    # Open GUI.
    editData_gui(env=.strvalidator_env, savegui=.save_gui,
                 edit=TRUE, debug=debug, parent=w)
    
  } )
  
  # CALCULATE -----------------------------------------------------------------
  
  mix_grid[2,1] <- mix_calculate_btn <- gbutton(text="Calculate",
                                                  border=TRUE,
                                                  container = mix_grid) 
  
  mix_grid[2,2] <- glabel(text="Calculate mixture for a dataset.",
                           container=mix_grid,
                           anchor=c(-1 ,0))
  
  
  addHandlerChanged(mix_calculate_btn, handler = function(h, ...) {
    
    # Open GUI.
    calculateMixture_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # PLOT MIXTURE --------------------------------------------------------------
  
  # SUMMARY TABLE -------------------------------------------------------------
  
  # RESULT  ###################################################################
  
  glabel("", container=result_tab) # Adds some space.
  
  result_grid <- glayout(container = result_tab)
  
  # VIEW/EDIT -----------------------------------------------------------------
  
  result_grid[1,1] <- result_view_btn <- gbutton(text="Edit",
                                                 border=TRUE,
                                                 container = result_grid) 
  
  result_grid[1,2] <- glabel(text="Edit or view a dataset.",
                             container=result_grid,
                             anchor=c(-1 ,0))
  
  addHandlerChanged(result_view_btn, handler = function(h, ...) {
    
    # Open GUI.
    editData_gui(env=.strvalidator_env, savegui=.save_gui,
                 edit=TRUE, debug=debug, parent=w)
    
  } )
  
  # RESULT TYPE ===============================================================
  
  result_f1 <- gframe(text = "Result types",
                      horizontal=FALSE, container = result_tab) 
  
  result_g1 <- glayout(container = result_f1)
  
  
  # CALCULATE -----------------------------------------------------------------
  
  result_g1[1,1] <- result_g1_calc_btn <- gbutton(text="Calculate",
                                                  border=TRUE,
                                                  container = result_g1) 
  
  result_g1[1,2] <- glabel(text="Calculate result types for a dataset.",
                           container=result_g1,
                           anchor=c(-1 ,0))
  
  
  addHandlerChanged(result_g1_calc_btn, handler = function(h, ...) {
    
    # Open GUI.
    calculateResultType_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # PLOT RESULT TYPE ----------------------------------------------------------
  
  result_g1[2,1] <- result_g1_plot_btn <- gbutton(text="Plot",
                                                  border=TRUE,
                                                  container = result_g1) 
  
  result_g1[2,2] <- glabel(text="Create plots for analysed data",
                           container=result_g1)
  
  addHandlerChanged(result_g1_plot_btn, handler = function(h, ...) {
    
    # Open GUI.
    plotResultType_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # PEAKS =====================================================================
  
  result_f2 <- gframe(text = "Number of peaks",
                      horizontal=FALSE, container = result_tab) 
  
  result_g2 <- glayout(container = result_f2)
  
  
  # CALCULATE -----------------------------------------------------------------
  
  result_g2[1,1] <- result_g2_calc_btn <- gbutton(text="Calculate",
                                                  border=TRUE,
                                                  container = result_g2) 
  
  result_g2[1,2] <- glabel(text="Count the number of peaks in sample.",
                           container=result_g2,
                           anchor=c(-1 ,0))
  
  
  addHandlerChanged(result_g2_calc_btn, handler = function(h, ...) {
    
    # Open GUI.
    calculatePeaks_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # PLOT PEAKS ----------------------------------------------------------------
  
  result_g2[2,1] <- result_g2_plot_btn <- gbutton(text="Plot",
                                                  border=TRUE,
                                                  container = result_g2) 
  
  result_g2[2,2] <- glabel(text="Create plots for analysed data",
                           container=result_g2)
  
  addHandlerChanged(result_g2_plot_btn, handler = function(h, ...) {
    
    # Open GUI.
    plotPeaks_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  

  # PEAK HEIGHT ===============================================================
  
  result_f3 <- gframe(text = "Peak height metrics",
                      horizontal=FALSE, container = result_tab) 
  
  result_g3 <- glayout(container = result_f3)
  
  # PLOT PEAKS ----------------------------------------------------------------
  
  result_g3[1,1] <- result_g3_calc_btn <- gbutton(text="Calculate",
                                                  border=TRUE,
                                                  container = result_g3) 
  
  result_g3[1,2] <- glabel(text="Calculate average and total peak height",
                           container=result_g3)
  
  addHandlerChanged(result_g3_calc_btn, handler = function(h, ...) {
    
    # Open GUI.
    calculateHeight_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  
  # DISTRIBUTIONS =============================================================
  
  result_f4 <- gframe(text = "Distributions",
                      horizontal=FALSE, container = result_tab) 
  
  result_g4 <- glayout(container = result_f4)
  
  # PLOT PEAKS ----------------------------------------------------------------
  
  result_g4[1,1] <- result_g4_plot_btn <- gbutton(text="Plot",
                                                  border=TRUE,
                                                  container = result_g4) 
  
  result_g4[1,2] <- glabel(text="Plot distributions for analysed data",
                           container=result_g4)
  
  addHandlerChanged(result_g4_plot_btn, handler = function(h, ...) {
    
    # Open GUI.
    plotDistribution_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )

    
  # DROPIN ====================================================================
  
  result_f5 <- gframe(text = "Drop-in tools",
                      horizontal=FALSE, container = result_tab) 
  
  result_g5 <- glayout(container = result_f5)
  
  
  # CALCULATE -----------------------------------------------------------------
  
  result_g5[1,1] <- result_g5_calc_btn <- gbutton(text="Calculate",
                                                  border=TRUE,
                                                  container = result_g5) 
  
  result_g5[1,2] <- glabel(text="Find spikes in sample.",
                           container=result_g5,
                           anchor=c(-1 ,0))
  
  
  addHandlerChanged(result_g5_calc_btn, handler = function(h, ...) {
    
    # Open GUI.
    calculateSpike_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # FILTER PEAKS --------------------------------------------------------------
  
  result_g5[2,1] <- result_g5_filter_btn <- gbutton(text="Filter",
                                                  border=TRUE,
                                                  container = result_g5) 
  
  result_g5[2,2] <- glabel(text="Filter spikes from data.",
                           container=result_g5)
  
  addHandlerChanged(result_g5_filter_btn, handler = function(h, ...) {
    
    # Open GUI.
    #filterSpike_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    gmessage(message = "This function is not yet implemented", title = "No function",
             icon = "warning", parent = w)
    
    
  } )
  
  # PRECISION  ################################################################
  
  glabel("", container=precision_tab) # Adds some space.
  
  precision_grid <- glayout(container = precision_tab)
  
  # VIEW/EDIT -----------------------------------------------------------------
  
  precision_grid[1,1] <- precision_view_btn <- gbutton(text="Edit",
                                                       border=TRUE,
                                                       container = precision_grid) 
  
  precision_grid[1,2] <- glabel(text="Edit or view a dataset.",
                                container=precision_grid,
                                anchor=c(-1 ,0))
  
  addHandlerChanged(precision_view_btn, handler = function(h, ...) {
    
    # Open GUI.
    editData_gui(env=.strvalidator_env, savegui=.save_gui,
                 edit=TRUE, debug=debug, parent=w)
    
  } )
  
  # PLOT RESULT TYPE ----------------------------------------------------------
  
  precision_grid[2,1] <- precision_plot_btn <- gbutton(text="Plot",
                                                       border=TRUE,
                                                       container = precision_grid) 
  
  precision_grid[2,2] <- glabel(text="Create plots for analysed data",
                                container=precision_grid)
  
  addHandlerChanged(precision_plot_btn, handler = function(h, ...) {
    
    # Open GUI.
    plotPrecision_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # SUMMARY TABLE -------------------------------------------------------------
  
  precision_grid[3,1] <- precision_table_btn <- gbutton(text="Summarize",
                                                        border=TRUE,
                                                        container = precision_grid) 
  
  precision_grid[3,2] <- glabel(text="Summarize precision data in a table.",
                                container=precision_grid,
                                anchor=c(-1 ,0))
  
  
  addHandlerChanged(precision_table_btn, handler = function(h, ...) {
    
    # Open GUI.
    tablePrecision_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # PULLUP  ###################################################################
  
  glabel("", container=pullup_tab) # Adds some space.
  
  pull_grid <- glayout(container = pullup_tab)
  
  # VIEW/EDIT -----------------------------------------------------------------
  
  pull_grid[1,1] <- pull_view_btn <- gbutton(text="Edit",
                                             border=TRUE,
                                             container = pull_grid) 
  
  pull_grid[1,2] <- glabel(text="Edit or view a dataset.",
                           container=pull_grid,
                           anchor=c(-1 ,0))
  
  addHandlerChanged(pull_view_btn, handler = function(h, ...) {
    
    # Open GUI.
    editData_gui(env=.strvalidator_env, savegui=.save_gui,
                 edit=TRUE, debug=debug, parent=w)
    
  } )
  
  # CALCULATE -----------------------------------------------------------------
  
  pull_grid[2,1] <- pull_calculate_btn <- gbutton(text="Calculate",
                                                  border=TRUE,
                                                  container = pull_grid) 
  
  pull_grid[2,2] <- glabel(text="Calculate spectral pull-up (aka. bleed-through).",
                           container=pull_grid,
                           anchor=c(-1 ,0))
  
  addHandlerChanged(pull_calculate_btn, handler = function(h, ...) {
    
    # Open GUI.
    calculatePullup_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
  } )
  
  # PLOT PULLUP ---------------------------------------------------------------
  
  pull_grid[3,1] <- pull_plot_btn <- gbutton(text="Plot",
                                                       border=TRUE,
                                                       container = pull_grid) 
  
  pull_grid[3,2] <- glabel(text="Create plots for analysed data",
                                container=pull_grid)
  
  addHandlerChanged(pull_plot_btn, handler = function(h, ...) {
    
    # Open GUI.
    plotPullup_gui(env=.strvalidator_env, savegui=.save_gui, debug=debug, parent=w)
    
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
    
    # Check if a tab name exist and then perform tasks.
    if(length(tabName) != 0){

      if(tabName == .file_tab_name){
        
        .refreshLoaded()
        .refreshWs()
        
      }
      
      if(tabName == .project_tab_name){
        
        .updateProjectList()
        
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
      
      if(tabName == .file_tab_name){
        
        .refreshLoaded()
        .refreshWs()
        
      }
    
    } # End check.
    
  })
  
  # INTERNAL FUNCTIONS ########################################################
  
  .loadSavedSettings <- function(){
    
    # First load save flag.
    if(exists(".strvalidator_savegui", envir=.strvalidator_env, inherits = FALSE)){
      svalue(savegui_chk) <- get(".strvalidator_savegui", envir=.strvalidator_env)
    }
    
    # Then load settings if true.
    if(svalue(savegui_chk)){
      if(exists(".strvalidator_project_dir", envir=.strvalidator_env, inherits = FALSE)){
        svalue(project_fb) <- get(".strvalidator_project_dir", envir=.strvalidator_env)
      }
    }
    
    if(debug){
      print("Saved settings loaded!")
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_savegui", value=svalue(savegui_chk), envir=.strvalidator_env)
      assign(x=".strvalidator_project_dir", value=svalue(project_fb), envir=.strvalidator_env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_savegui", envir=.strvalidator_env, inherits = FALSE)){
        remove(".strvalidator_savegui", envir = .strvalidator_env)
      }
      if(exists(".strvalidator_project_dir", envir=.strvalidator_env, inherits = FALSE)){
        remove(".strvalidator_project_dir", envir = .strvalidator_env)
      }
      
      if(debug){
        print("Settings cleared!")
      }
    }
    
    if(debug){
      print("Settings saved!")
    }
    
  }
  
  .refreshWs <- function(){
    
    # Get data frames in global workspace.
    dfs <- listObjects(env=.GlobalEnv, obj.class=.object_classes_import)
    
    if(!is.null(dfs)){
      
      blockHandler(ws_r_drp)
      
      # Populate drop list.
      ws_r_drp[] <- c("<Select dataframe>", dfs)
      
      # Select first item.
      svalue(ws_r_drp, index=TRUE) <- 1 
      
      unblockHandler(ws_r_drp)
      
    }
  }
  
  .refreshLoaded <- function(){
    
    if(debug){
      print(paste("IN:", match.call()[[1]]))
    }
    
    # Get list of objects.
    dfs <- listObjects(env=.strvalidator_env, obj.class=.object_classes_view)
    
    # Get size of objects.
    dfsSize <- sapply(dfs, function(x) object.size(get(x, envir = .strvalidator_env)))
    dfsSize <- unname(dfsSize)
    dfsSize <- as.numeric(dfsSize)
    
    if(!is.null(dfs)){
      
      #blockHandler(ws_loaded_tbl) # Not working.
      
      # Populate table.
      ws_loaded_tbl[,] <- data.frame(Object=dfs, Size=dfsSize,
                                     stringsAsFactors=FALSE)
      
      #unblockHandler(ws_loaded_tbl)
      
    }
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
  
  # SHOW GUI ##################################################################
  
  # Show GUI and first tab.
  svalue(nb)<-1
  visible(w)<-TRUE
  focus(w)
  message("STR-validator graphical user interface loaded!")
  
} # END OF GUI
