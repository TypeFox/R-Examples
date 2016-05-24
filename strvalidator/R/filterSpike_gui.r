# ################################################################################
# # TODO LIST
# # TODO: ... Automatic filtering of artefact peaks...
# 
# ################################################################################
# # CHANGE LOG (last 20 changes)
# # 12.10.2015: First version.
# 
# 
# #' @title Filter Spike
# #'
# #' @description
# #' GUI wrapper for the \code{\link{filterSpike}} function.
# #'
# #' @details
# #' Simplifies the use of the \code{\link{filterSpike}} function by providing a
# #'  graphical user interface to it.
# #' 
# #' @param env environment in wich to search for data frames.
# #' @param savegui logical indicating if GUI settings should be saved in the environment.
# #' @param debug logical indicating printing debug information.
# #' @param parent widget to get focus when finished.
# #' 
# #' @export
# #' 
# #' @importFrom utils help
# #' 
# #' @return TRUE
# 
# 
# filterArtefacts_gui <- function(env=parent.frame(), savegui=NULL, debug=FALSE, parent=NULL){
#   
#   # Global variables.
#   .gData <- NULL
#   .gDataName <- NULL
#   
#   if(debug){
#     print(paste("IN:", match.call()[[1]]))
#   }
#   
#   # Main window.
#   w <- gwindow(title="Filter artefacts", visible=FALSE)
#   
#   # Runs when window is closed.
#   addHandlerDestroy(w, handler = function (h, ...) {
#     
#     # Save GUI state.
#     .saveSettings()
#     
#     # Focus on parent window.
#     if(!is.null(parent)){
#       focus(parent)
#     }
#     
#   })
#   
#   # Vertical main group.
#   gv <- ggroup(horizontal=FALSE,
#               spacing=15,
#               use.scrollwindow=FALSE,
#               container = w,
#               expand=FALSE) 
# 
#   # Help button group.
#   gh <- ggroup(container = gv, expand=FALSE, fill="both")
#   
#   savegui_chk <- gcheckbox(text="Save GUI settings", checked=FALSE, container=gh)
#   
#   addSpring(gh)
#   
#   help_btn <- gbutton(text="Help", container=gh)
#   
#   addHandlerChanged(help_btn, handler = function(h, ...) {
#     
#     # Open help page for function.
#     print(help("filterSpike_gui", help_type="html"))
#     
#   })
#   
#   # DATASET ###################################################################
#   
#   f0 <- gframe(text = "Dataset",
#                horizontal=FALSE,
#                spacing = 10,
#                container = gv) 
# 
# 
#   f0g0 <- glayout(container = f0, spacing = 1)
#   
#   f0g0[1,1] <- glabel(text="Select dataset:", container=f0g0)
#   
#   f0g0[1,2] <- f0g0_data_drp <- gdroplist(items=c("<Select dataset>",
#                                                  listObjects(env=env,
#                                                              obj.class="data.frame")),
#                                          selected = 1,
#                                          editable = FALSE,
#                                          container = f0g0)
#   
#   f0g0[1,3] <- f0g0_data_col_lbl <- glabel(text=" 0 columns",
#                                               container=f0g0)
#   
#   addHandlerChanged(f0g0_data_drp, handler = function (h, ...) {
#     
#     val_obj <- svalue(f0g0_data_drp)
#     
#     # Check if suitable.
#     ok <- is.data.frame(get(val_obj, envir=env))
#     
#     if(ok){
#       
#       # Load or change components.
#       .gData <<- get(val_obj, envir=env)
#       .gDataName <<- val_obj
#       
#       svalue(f0g0_data_col_lbl) <- paste(" ", ncol(.gData), " columns")
#       svalue(f2_name) <- paste(.gDataName, "new", sep="_")
#       
#       f1g1_col1_drp[] <- c("<Select column>", names(.gData))
#       f1g1_col2_drp[] <- c("<Select column>", names(.gData))
#       
#     } else {
#       
#       .gData <<- NULL
#       .gDataName <<- NULL
#       svalue(f0g0_data_col_lbl) <- " 0 columns"
#       svalue(f2_name) <- ""
#       
#       f1g1_col1_drp[] <- c("<Select column>")
#       f1g1_col2_drp[] <- c("<Select column>")
# 
#     }
#     
#   } )
#   
#   # COLUMNS ###################################################################
#   
#   f1 <- gframe(text="Columns", horizontal=FALSE, spacing=10, container=gv)
#   
#   f1g1 <- glayout(container=f1, spacing=1)
#   
#   f1g1[1,1] <- glabel(text="Select column 1:", container=f1g1)
#   
#   f1g1[1,2] <- f1g1_col1_drp <- gdroplist(items=c("<Select column>"),
#                                            editable = FALSE,
#                                            container = f1g1)
#   
#   f1g1[2,1] <- glabel(text="Select column 2:", container=f1g1)
#   
#   f1g1[2,2] <- f1g1_col2_drp <- gdroplist(items=c("<Select column>"),
#                                            editable = FALSE,
#                                            container = f1g1)
# 
#   addHandlerChanged(f1g1_col1_drp, handler = function (h, ...) {
#     
#     val_col <- svalue(f1g1_col1_drp)
#     
#     # Check if column exist.
#     ok <- val_col %in% names(.gData)
# 
#     # Update target column.
#     if(length(ok)>0){
#       if(ok){
#         svalue(f3g1_col_edt) <- val_col
#       } else {
#         svalue(f3g1_col_edt) <- ""
#       }
#     }
#     
#   } )
#   
#   addHandlerChanged(f1g1_col2_drp, handler = function (h, ...) {
#     
#     val_col <- svalue(f1g1_col2_drp)
#     
#     # Check if column exist.
#     ok <- val_col %in% names(.gData)
# 
#     if(length(ok)>0){
#       # Enable widgets.
#       if(ok){
#         enabled(f3g1_val_edt) <- FALSE
#       } else {
#         enabled(f3g1_val_edt) <- TRUE
#       }
#     }
#     
#   } )
#   
#   # COLUMNS ###################################################################
#   
#   f3 <- gframe(text="Options", horizontal=FALSE, spacing=10, container=gv)
#   
#   f3g1 <- glayout(container=f3, spacing=1)
#   
#   f3g1[1,1] <- glabel(text="Fixed value:", container=f3g1)
#   
#   f3g1[1,2] <- f3g1_val_edt <- gedit(text="", width=25, container=f3g1)
#   
#   f3g1[2,1] <- glabel(text="Column for new values:", container=f3g1)
#   
#   f3g1[2,2] <- f3g1_col_edt <- gedit(text="", width=25, container=f3g1)
#   
#   f3g1[3,1] <- glabel(text="Action:", container=f3g1)
#   
#   action_items <- c("&","+","*","-", "/")
#   f3g1[3,2] <- f3g1_action_drp <- gdroplist(items=action_items, selected=1,
#                                               editable=FALSE, container=f3g1)
#   
#   # NAME ######################################################################
#   
#   f2 <- gframe(text = "Save as",
#                horizontal=TRUE,
#                spacing = 5,
#                container = gv) 
#   
#   glabel(text="Save as:", container=f2)
#   f2_name <- gedit(text="", width=40, container=f2)
#   
#   # BUTTON ####################################################################
# 
#   if(debug){
#     print("BUTTON")
#   }  
#   
#   combine_btn <- gbutton(text="Execute",
#                       border=TRUE,
#                       container=gv)
#   
#   addHandlerChanged(combine_btn, handler = function(h, ...) {
#     
#     val_col1 <- svalue(f1g1_col1_drp)
#     val_col2 <- svalue(f1g1_col2_drp)
#     val_action <- svalue(f3g1_action_drp)
#     val_target <- svalue(f3g1_col_edt)
#     val_fixed <- svalue(f3g1_val_edt)
#     val_name <- svalue(f2_name)
#     
#     colOk <- val_col1 %in% names(.gData)
#     
#     if (colOk){
#       
#       datanew <- filterSpike(data=.gData, col1=val_col1, col2=val_col2,
#                          operator=val_action, fixed=val_fixed,
#                          target=val_target, debug=debug)
#       
#       # Save data.
#       saveObject(name=val_name, object=datanew, parent=w, env=env)
#       
#       if(debug){
#         print(datanew)
#         print(paste("EXIT:", match.call()[[1]]))
#       }
#       
#       # Close GUI.
#       dispose(w)
#       
#     } else {
#       
#       gmessage(message="Selected column must exist in data frame!",
#                title="Error",
#                icon = "error")      
#       
#     } 
#     
#   } )
#   
#   # INTERNAL FUNCTIONS ########################################################
#   
#   .loadSavedSettings <- function(){
#     
#     # First check status of save flag.
#     if(!is.null(savegui)){
#       svalue(savegui_chk) <- savegui
#       enabled(savegui_chk) <- FALSE
#       if(debug){
#         print("Save GUI status set!")
#       }  
#     } else {
#       # Load save flag.
#       if(exists(".strvalidator_filterSpike_gui_savegui", envir=env, inherits = FALSE)){
#         svalue(savegui_chk) <- get(".strvalidator_filterSpike_gui_savegui", envir=env)
#       }
#       if(debug){
#         print("Save GUI status loaded!")
#       }  
#     }
#     if(debug){
#       print(svalue(savegui_chk))
#     }  
#     
#     # Then load settings if true.
#     if(svalue(savegui_chk)){
#       if(exists(".strvalidator_filterSpike_gui_import_opt", envir=env, inherits = FALSE)){
#         svalue(import_opt) <- get(".strvalidator_filterSpike_gui_import_opt", envir=env)
#       }
#       if(exists(".strvalidator_filterSpike_gui_fixed", envir=env, inherits = FALSE)){
#         svalue(f3g1_val_edt) <- get(".strvalidator_filterSpike_gui_fixed", envir=env)
#       }
#       if(exists(".strvalidator_filterSpike_gui_action", envir=env, inherits = FALSE)){
#         svalue(f3g1_action_drp) <- get(".strvalidator_filterSpike_gui_action", envir=env)
#       }
#       if(debug){
#         print("Saved settings loaded!")
#       }
#     }
#     
#   }
#   
#   .saveSettings <- function(){
#     
#     # Then save settings if true.
#     if(svalue(savegui_chk)){
#       
#       assign(x=".strvalidator_filterSpike_gui_savegui", value=svalue(savegui_chk), envir=env)
#       assign(x=".strvalidator_filterSpike_gui_fixed", value=svalue(f3g1_val_edt), envir=env)
#       assign(x=".strvalidator_filterSpike_gui_action", value=svalue(f3g1_action_drp), envir=env)
#       
#     } else { # or remove all saved values if false.
#       
#       if(exists(".strvalidator_filterSpike_gui_savegui", envir=env, inherits = FALSE)){
#         remove(".strvalidator_filterSpike_gui_savegui", envir = env)
#       }
#       if(exists(".strvalidator_filterSpike_gui_fixed", envir=env, inherits = FALSE)){
#         remove(".strvalidator_filterSpike_gui_fixed", envir = env)
#       }
#       if(exists(".strvalidator_filterSpike_gui_action", envir=env, inherits = FALSE)){
#         remove(".strvalidator_filterSpike_gui_action", envir = env)
#       }
#       
#       if(debug){
#         print("Settings cleared!")
#       }
#     }
#     
#     if(debug){
#       print("Settings saved!")
#     }
#     
#   }
#   
#   # END GUI ###################################################################
# 
#   # Load GUI settings.
#   .loadSavedSettings()
#   
#   # Show GUI.
#   visible(w) <- TRUE
#   focus(w)
#   
# } # End of GUI
