################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 31.12.2015: New options 'wrap' and 'at'. 'type' replaced by 'boxplot'. 
# 29.08.2015: Added importFrom.
# 09.01.2015: Enable 'generate' after selection of new sample.
# 09.12.2014: First version.


#' @title Generate EPG
#'
#' @description
#' GUI wrapper for the \code{\link{generateEPG}} function.
#'
#' @details
#' Simplifies the use of the \code{\link{generateEPG}} function by providing a graphical 
#' user interface to it.
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
#' @importFrom utils help
#' 
#' @seealso \code{\link{generateEPG}}

generateEPG_gui <- function(env=parent.frame(), savegui=NULL, debug=FALSE, parent=NULL){
  
  # Global variables.
  .gData <- NULL
  .gPlot <- NULL
  .noSample <- "<Select sample>"
  .buttonDefault <- "Generate EPG"
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # Main window.
  w <- gwindow(title="Generate electropherogram", visible=FALSE)

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
  gv <- ggroup(horizontal=FALSE, spacing=8, use.scrollwindow=FALSE,
               container = w, expand=TRUE) 

  # Help button group.
  gh <- ggroup(container = gv, expand=FALSE, fill="both")
  
  savegui_chk <- gcheckbox(text="Save GUI settings", checked=FALSE, container=gh)
  
  addSpring(gh)
  
  help_btn <- gbutton(text="Help", container=gh)
  
  addHandlerChanged(help_btn, handler = function(h, ...) {
    
    # Open help page for function.
    print(help("generateEPG_gui", help_type="html"))
    
  })
  
  # FRAME 0 ###################################################################
  
  f0 <- gframe(text = "Dataset and kit",
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
  
  g0[1,3] <- g0_data_samples_lbl <- glabel(text=" (0 samples)", container=g0)
  
  g0[1,4] <- glabel(text=" and the kit used:", container=g0)
  
  g0[1,5] <- kit_drp <- gdroplist(items=getKit(), selected = 1,
                                  editable = FALSE, container = g0) 
  
  # Sample --------------------------------------------------------------------
  
  g0[2,1] <- glabel(text="Select sample:", container=g0)
  
  g0[2,2] <- g0_sample_drp <- gdroplist(items=.noSample, selected = 1,
                                        editable = FALSE, container = g0)
  
  # Handlers ------------------------------------------------------------------
  
  
  addHandlerChanged(g0_data_drp, handler = function (h, ...) {
    
    val_obj <- svalue(g0_data_drp)
    
    # Check if suitable.
    requiredCol <- c("Sample.Name", "Marker", "Allele")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       env=env, parent=w, debug=debug)
    
    if(ok){
      # Load or change components.
      
      # Get data.
      .gData <<- get(val_obj, envir=env)
      
      # Suggest name.
      svalue(f5_save_edt) <- paste(val_obj, "_ggplot", sep="")
      
      svalue(g0_data_samples_lbl) <- paste(" (",
                                           length(unique(.gData$Sample.Name)),
                                           " samples)", sep="")
      
      # Detect kit.
      kitIndex <- detectKit(.gData, index=TRUE)
      # Select in dropdown.
      svalue(kit_drp, index=TRUE) <- kitIndex
      
      # Populate dropdown with sample names.
      .refresh_sample_drp()
      
      # Set dataset as proposed title.
      svalue(f1_title_edt) <- paste(val_obj, " (",
                                    svalue(kit_drp), ")", sep="")
      
      # Enable buttons.
      enabled(plot_epg_btn) <- TRUE
      
    } else {
      
      # Reset components.
      .gData <<- NULL
      svalue(f5_save_edt) <- ""
      svalue(g0_data_drp, index=TRUE) <- 1
      svalue(g0_data_samples_lbl) <- " (0 samples)"
      
    }
    
  } )  
  
  addHandlerChanged(g0_sample_drp, handler = function (h, ...) {
    
    # Get selected sample name.
    val_sample <- svalue(g0_sample_drp)
    
    if(!is.null(val_sample)){
      if(val_sample != .noSample){
        # Set sample name as proposed title.
        svalue(f1_title_edt) <- paste(val_sample, " (",
                                      svalue(kit_drp), ")", sep="")
        
        # Enable buttons.
        enabled(plot_epg_btn) <- TRUE
        
      }
    }
    
  } )  

  # FRAME 1 ###################################################################
  
  f1 <- gframe(text = "Options", horizontal=FALSE, spacing = 10, container = gv)
  
  glabel(text="Plot title:", anchor=c(-1 ,0), container=f1)
  f1_title_edt <- gedit(text="", width=25, container=f1)
  
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

  # FRAME 2 ###################################################################
  
  f5 <- gframe(text = "Save as",
               horizontal=TRUE,
               spacing = 5,
               container = gv) 
  
  glabel(text="Name for result:", container=f5)
  
  f5_save_edt <- gedit(text="", container=f5)
  
  f5_save_btn <- gbutton(text = "Save as object",
                         border=TRUE,
                         container = f5) 
  
  f5_ggsave_btn <- gbutton(text = "Save as image",
                           border=TRUE,
                           container = f5) 
  
  addHandlerChanged(f5_save_btn, handler = function(h, ...) {
    
    val_name <- svalue(f5_save_edt)
    
    # Change button.
    svalue(f5_save_btn) <- "Processing..."
    enabled(f5_save_btn) <- FALSE
    
    # Save data.
    saveObject(name=val_name, object=.gPlot,
               parent=w, env=env, debug=debug)
    
    # Change button.
    svalue(f5_save_btn) <- "Object saved"
    
  } )
  
  addHandlerChanged(f5_ggsave_btn, handler = function(h, ...) {
    
    val_name <- svalue(f5_save_edt)
    
    # Save data.
    ggsave_gui(ggplot=.gPlot, name=val_name, 
               parent=w, env=env, savegui=savegui, debug=debug)
    
  } )
  
  
  # BUTTON ####################################################################
  
  
  plot_epg_btn <- gbutton(text=.buttonDefault,
                        border=TRUE,
                        container=gv)
  
  addHandlerChanged(plot_epg_btn, handler = function(h, ...) {
    
    val_name <- svalue(f5_save_edt)
    val_sample <- svalue(g0_sample_drp)
    val_kit <- svalue(kit_drp)
    val_peaks <- svalue(f1_peaks_chk)
    val_scale <- svalue(f1_scale_opt)
    val_wrap <- svalue(f1_wrap_chk)
    val_box <- svalue(f1_box_chk)
    val_collapse <- svalue(f1_collapse_chk)
    val_fix <- svalue(f1_fix_chk)
    val_ignore <- svalue(f1_ignore_chk)
    val_title <- svalue(f1_title_edt)
    val_size <- svalue(f1_size_spb)
    val_angle <- svalue(f1_angle_spb)
    val_vjust <- svalue(f1_vjust_spb)
    val_hjust <- svalue(f1_hjust_spb)
    val_expand <- svalue(f1_expand_spb)
    val_at <- svalue(f1_at_spb)
    val_data <- .gData
    
    if(!is.null(val_data)){
      
      if(val_sample != .noSample){
        # Subset selected sample.
        val_data <- val_data[val_data$Sample.Name==val_sample,]
      }
      
      # Change button.
      svalue(plot_epg_btn) <- "Processing..."
      enabled(plot_epg_btn) <- FALSE
      
      gp <- generateEPG(data = val_data, kit = val_kit, title = val_title,
                        wrap = val_wrap, boxplot = val_box, peaks = val_peaks,
                        collapse = val_collapse, silent = FALSE,
                        ignore.case = val_ignore, at = val_at,
                        scale = val_scale, limit.x = val_fix,
                        label.size = val_size, label.angle = val_angle,
                        label.vjust = val_vjust, label.hjust = val_hjust,
                        expand = val_expand, debug = debug)
        
      # Store in global variable.
      .gPlot <<- gp

      # Change button.
      svalue(plot_epg_btn) <- .buttonDefault
      enabled(plot_epg_btn) <- TRUE
      
      
    } else {
      
      message <- "A dataset must be selected. Sample is optional."
      
      gmessage(message, title="Datasets not selected",
               icon = "error",
               parent = w) 
      
    }
    
  } )

  # INTERNAL FUNCTIONS ########################################################
  
  .refresh_sample_drp <- function(){
    
    if(debug){
      print("Refresh sample dropdown")
    }
    
    # Get data frames in global workspace.
    dfs <- unique(.gData$Sample.Name)
    
    if(!is.null(dfs)){
      
      blockHandler(g0_sample_drp)
      
      # Populate drop list and select first item.
      g0_sample_drp[] <- c(.noSample, dfs)
      svalue(g0_sample_drp, index=TRUE) <- 1
      
      unblockHandler(g0_sample_drp)
      
    }
    
    if(debug){
      print("Sample dropdown refreshed!")
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
      if(exists(".strvalidator_generateEPG_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_generateEPG_gui_savegui", envir=env)
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
      if(exists(".strvalidator_generateEPG_gui_size", envir=env, inherits = FALSE)){
        svalue(f1_size_spb) <- get(".strvalidator_generateEPG_gui_size", envir=env)
      }
      if(exists(".strvalidator_generateEPG_gui_angle", envir=env, inherits = FALSE)){
        svalue(f1_angle_spb) <- get(".strvalidator_generateEPG_gui_angle", envir=env)
      }
      if(exists(".strvalidator_generateEPG_gui_vjust", envir=env, inherits = FALSE)){
        svalue(f1_vjust_spb) <- get(".strvalidator_generateEPG_gui_vjust", envir=env)
      }
      if(exists(".strvalidator_generateEPG_gui_hjust", envir=env, inherits = FALSE)){
        svalue(f1_hjust_spb) <- get(".strvalidator_generateEPG_gui_hjust", envir=env)
      }
      if(exists(".strvalidator_generateEPG_gui_expand", envir=env, inherits = FALSE)){
        svalue(f1_expand_spb) <- get(".strvalidator_generateEPG_gui_expand", envir=env)
      }
      if(exists(".strvalidator_generateEPG_gui_scales", envir=env, inherits = FALSE)){
        svalue(f1_scale_opt) <- get(".strvalidator_generateEPG_gui_scales", envir=env)
      }
      if(exists(".strvalidator_generateEPG_gui_collapse", envir=env, inherits = FALSE)){
        svalue(f1_collapse_chk) <- get(".strvalidator_generateEPG_gui_collapse", envir=env)
      }
      if(exists(".strvalidator_generateEPG_gui_fix", envir=env, inherits = FALSE)){
        svalue(f1_fix_chk) <- get(".strvalidator_generateEPG_gui_fix", envir=env)
      }
      if(exists(".strvalidator_generateEPG_gui_ignore", envir=env, inherits = FALSE)){
        svalue(f1_ignore_chk) <- get(".strvalidator_generateEPG_gui_ignore", envir=env)
      }
      if(exists(".strvalidator_generateEPG_gui_peaks", envir=env, inherits = FALSE)){
        svalue(f1_peaks_chk) <- get(".strvalidator_generateEPG_gui_peaks", envir=env)
      }
      if(exists(".strvalidator_generateEPG_gui_box", envir=env, inherits = FALSE)){
        svalue(f1_box_chk) <- get(".strvalidator_generateEPG_gui_box", envir=env)
      }
      if(exists(".strvalidator_generateEPG_gui_wrap", envir=env, inherits = FALSE)){
        svalue(f1_wrap_chk) <- get(".strvalidator_generateEPG_gui_wrap", envir=env)
      }
      if(exists(".strvalidator_generateEPG_gui_at", envir=env, inherits = FALSE)){
        svalue(f1_at_spb) <- get(".strvalidator_generateEPG_gui_at", envir=env)
      }
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_generateEPG_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_generateEPG_gui_size", value=svalue(f1_size_spb), envir=env)
      assign(x=".strvalidator_generateEPG_gui_angle", value=svalue(f1_angle_spb), envir=env)
      assign(x=".strvalidator_generateEPG_gui_vjust", value=svalue(f1_vjust_spb), envir=env)
      assign(x=".strvalidator_generateEPG_gui_hjust", value=svalue(f1_hjust_spb), envir=env)
      assign(x=".strvalidator_generateEPG_gui_expand", value=svalue(f1_expand_spb), envir=env)
      assign(x=".strvalidator_generateEPG_gui_scales", value=svalue(f1_scale_opt), envir=env)
      assign(x=".strvalidator_generateEPG_gui_collapse", value=svalue(f1_collapse_chk), envir=env)
      assign(x=".strvalidator_generateEPG_gui_fix", value=svalue(f1_fix_chk), envir=env)
      assign(x=".strvalidator_generateEPG_gui_ignore", value=svalue(f1_ignore_chk), envir=env)
      assign(x=".strvalidator_generateEPG_gui_peaks", value=svalue(f1_peaks_chk), envir=env)
      assign(x=".strvalidator_generateEPG_gui_box", value=svalue(f1_box_chk), envir=env)
      assign(x=".strvalidator_generateEPG_gui_wrap", value=svalue(f1_wrap_chk), envir=env)
      assign(x=".strvalidator_generateEPG_gui_at", value=svalue(f1_at_spb), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_generateEPG_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_generateEPG_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_generateEPG_gui_size", envir=env, inherits = FALSE)){
        remove(".strvalidator_generateEPG_gui_size", envir = env)
      }
      if(exists(".strvalidator_generateEPG_gui_angle", envir=env, inherits = FALSE)){
        remove(".strvalidator_generateEPG_gui_angle", envir = env)
      }
      if(exists(".strvalidator_generateEPG_gui_vjust", envir=env, inherits = FALSE)){
        remove(".strvalidator_generateEPG_gui_vjust", envir = env)
      }
      if(exists(".strvalidator_generateEPG_gui_hjust", envir=env, inherits = FALSE)){
        remove(".strvalidator_generateEPG_gui_hjust", envir = env)
      }
      if(exists(".strvalidator_generateEPG_gui_expand", envir=env, inherits = FALSE)){
        remove(".strvalidator_generateEPG_gui_expand", envir = env)
      }
      if(exists(".strvalidator_generateEPG_gui_scales", envir=env, inherits = FALSE)){
        remove(".strvalidator_generateEPG_gui_scales", envir = env)
      }
      if(exists(".strvalidator_generateEPG_gui_collapse", envir=env, inherits = FALSE)){
        remove(".strvalidator_generateEPG_gui_collapse", envir = env)
      }
      if(exists(".strvalidator_generateEPG_gui_fix", envir=env, inherits = FALSE)){
        remove(".strvalidator_generateEPG_gui_fix", envir = env)
      }
      if(exists(".strvalidator_generateEPG_gui_ignore", envir=env, inherits = FALSE)){
        remove(".strvalidator_generateEPG_gui_ignore", envir = env)
      }
      if(exists(".strvalidator_generateEPG_gui_peaks", envir=env, inherits = FALSE)){
        remove(".strvalidator_generateEPG_gui_peaks", envir = env)
      }
      if(exists(".strvalidator_generateEPG_gui_box", envir=env, inherits = FALSE)){
        remove(".strvalidator_generateEPG_gui_box", envir = env)
      }
      if(exists(".strvalidator_generateEPG_gui_wrap", envir=env, inherits = FALSE)){
        remove(".strvalidator_generateEPG_gui_wrap", envir = env)
      }
      if(exists(".strvalidator_generateEPG_gui_at", envir=env, inherits = FALSE)){
        remove(".strvalidator_generateEPG_gui_at", envir = env)
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
