################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 06.01.2016: Fixed theme methods not found and added more themes.
# 22.12.2015: First version.

#' @title Plot Ratio
#'
#' @description
#' GUI simplifying the creation of plots from marker ratio data.
#'
#' @details Select data to plot in the drop-down menu. 
#' Automatic plot titles can be replaced by custom titles.
#' A name for the result is automatiaclly suggested.
#' The resulting plot can be saved as either a plot object or as an image.
#' @param env environment in wich to search for data frames.
#' @param savegui logical indicating if GUI settings should be saved in the environment.
#' @param debug logical indicating printing debug information.
#' @param parent widget to get focus when finished.
#' 
#' @return TRUE
#' 
#' @export
#' 
#' @importFrom utils help
#' @importFrom stats as.formula
#' @importFrom ggplot2 ggplot aes_string facet_wrap theme_gray theme_bw
#'  theme_linedraw theme_light theme_dark theme_minimal theme_classic
#'  theme_void 
#' @importFrom graphics par
#' 
#' @seealso \url{http://docs.ggplot2.org/current/} for details on plot settings.

plotRatio_gui <- function(env=parent.frame(), savegui=NULL, debug=FALSE, parent=NULL){
  
  # Global variables.
  .gData <- NULL
  .gDataName <- NULL
  .gPlot <- NULL
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # Main window.
  w <- gwindow(title="Plot marker ratios", visible=FALSE)
  
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
    print(help("plotRatio_gui", help_type="html"))
    
  })
  
  # FRAME 0 ###################################################################
  
  f0 <- gframe(text = "Dataset",
               horizontal=TRUE,
               spacing = 5,
               container = gv) 
  
  glabel(text="Select dataset:", container=f0)

  dataset_drp <- gdroplist(items=c("<Select dataset>",
                                   listObjects(env=env,
                                               obj.class="data.frame")), 
                           selected = 1,
                           editable = FALSE,
                           container = f0) 
  
  f0_samples_lbl <- glabel(text=" (0 samples)", container=f0)
  
  addHandlerChanged(dataset_drp, handler = function (h, ...) {
    
    val_obj <- svalue(dataset_drp)

    # Check if suitable.
    requiredCol <- c("Sample.Name")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       env=env, parent=w, debug=debug)
    
    if(ok){
      # Load or change components.

      # Get data.
      .gData <<- get(val_obj, envir=env)
      
      # Suggest name.
      svalue(f5_save_edt) <- paste(val_obj, "_ggplot", sep="")
      
      svalue(f0_samples_lbl) <- paste(" (",
                                      length(unique(.gData$Sample.Name)),
                                      " samples)", sep="")
      
      # Enable buttons.
      enabled(plot_browse_btn) <- TRUE
      enabled(plot_all_btn) <- TRUE
      
    } else {
    
      # Reset components.
      .gData <<- NULL
      svalue(f5_save_edt) <- ""
      svalue(dataset_drp, index=TRUE) <- 1
      svalue(f0_samples_lbl) <- " (0 samples)"
      
    }
    
  } )  
  
  # FRAME 1 ###################################################################
  
  f1 <- gframe(text = "Options",
               horizontal=FALSE,
               spacing = 5,
               container = gv)
  
  glabel(text = "Axes scales:", container=f1)
  f1_scales_opt <- gradio(items = c("free_x", "free_y", "free"),
                        selected = 1, container=f1)

  f1_titles_chk <- gcheckbox(text="Override automatic titles.",
                             checked=FALSE, container=f1)
  
  
  addHandlerChanged(f1_titles_chk, handler = function(h, ...) {
    val <- svalue(f1_titles_chk)
    if(val){
      enabled(f1g1) <- TRUE
    } else {
      enabled(f1g1) <- FALSE
    }
  } )
  
  f1g1 <- glayout(container = f1, spacing = 1)
  enabled(f1g1) <- svalue(f1_titles_chk)
  
  f1g1[1,1] <- glabel(text="Plot title:", container=f1g1)
  f1g1[1,2] <- title_edt <- gedit(text="",
                                   width=40,
                                   container=f1g1)
  
  f1g1[2,1] <- glabel(text="X title:", container=f1g1)
  f1g1[2,2] <- x_title_edt <- gedit(text="",
                                     container=f1g1)

  f1g1[3,1] <- glabel(text="Y title:", container=f1g1)
  f1g1[3,2] <- y_title_edt <- gedit(text="",
                                     container=f1g1)

  f1g2 <- glayout(container = f1)
  f1g2[1,1] <- glabel(text="Plot theme:", anchor=c(-1 ,0), container=f1g2)
  items_theme <- c("theme_grey()","theme_bw()","theme_linedraw()",
                   "theme_light()","theme_dark()","theme_minimal()",
                   "theme_classic()","theme_void()")
  f1g2[1,2] <- f1_theme_drp <- gdroplist(items = items_theme,
                                         selected = 1,
                                         container = f1g2)
  
  # FRAME 7 ###################################################################
  
  f7 <- gframe(text = "Plot marker ratios",
               horizontal=FALSE,
               container = gv) 
  
  grid7 <- glayout(container = f7)
  
  grid7[1,1] <- plot_browse_btn <- gbutton(text="Browse",
                                           border=TRUE,
                                           container=grid7)
  tooltip(plot_browse_btn) <- "Activate the console window and use Enter key to change plot"
  
  grid7[1,2] <- plot_all_btn <- gbutton(text="Plot",
                                           border=TRUE,
                                           container=grid7) 
  tooltip(plot_all_btn) <- "Plot all data in one plot, by group if available"
  
  addHandlerChanged(plot_browse_btn, handler = function(h, ...) {
    
    enabled(plot_browse_btn) <- FALSE
    .plot(what="browse")
    enabled(plot_browse_btn) <- TRUE

  } )
  
  addHandlerChanged(plot_all_btn, handler = function(h, ...) {
    
    enabled(plot_all_btn) <- FALSE
    .plot(what="plot")
    enabled(plot_all_btn) <- TRUE
    
  } )

  # FRAME 5 ###################################################################
  
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
  
  # FUNCTIONS #################################################################
  
  
  .plot <- function(what){
    
    # Get values.
    val_data <- .gData
    val_scales <- svalue(f1_scales_opt)
    val_titles <- svalue(f1_titles_chk)
    val_title <- svalue(title_edt)
    val_xtitle <- svalue(x_title_edt)
    val_ytitle <- svalue(y_title_edt)
    val_theme <- svalue(f1_theme_drp)

    if (!is.na(val_data) && !is.null(val_data)){
      
      # Fix data for plotting.
      vecColNames <- setdiff(names(val_data), c("Sample.Name", "Group"))
      
      # Expand vectors.
      intSamples <- length(val_data$Sample.Name)
      intColumns <- length(vecColNames)
      vecSampleName <- rep(val_data$Sample.Name, intColumns)
      vecGroup <- rep(val_data$Group, intColumns)

      # Create new data.frame for plotting.
      if(is.null(vecGroup)){
        dfPlot <- data.frame(Sample.Name = vecSampleName)
      } else {
        dfPlot <- data.frame(Sample.Name = vecSampleName, Group = vecGroup)
      }
      
      # Initialise variables.
      vecMarkers <- NULL
      vecRatio <- NULL
      
      # Loop over all column names containing ratios.
      for(i in seq(along = vecColNames)){

        # Add column data to vector.
        vecRatio <- c(vecRatio, val_data[, vecColNames[i]])
        vecMarkers <- c(vecMarkers, rep(vecColNames[i], intSamples))
        
      }
      
      # Add info.
      dfPlot$Ratio <- vecRatio
      dfPlot$Marker <- vecMarkers

      # Plotting alleles for observed stutters per marker.
      if(what == "browse"){
        
        # Enable confirm. NB! Remember to disable.
        par(ask=T) 
        
        # Loop through and plot individually.      
        for(i in seq(along = vecColNames)){

          if("Group" %in% names(val_data)){

            gp <- ggplot(subset(dfPlot, Marker==vecColNames[i]))
            gp <- gp + geom_boxplot(aes_string(x="Group", y="Ratio"))

            if(val_titles){
              
              gp <- .applyPlotSettings(gp = gp, theme = val_theme, main.title = val_title,
                                       x.title = val_xtitle, y.title = val_ytitle)
              
            } else {
              
              gp <- .applyPlotSettings(gp = gp, theme = val_theme, main.title = "Marker ratio",
                                       x.title = vecColNames[i], y.title = "Ratio")
              
            }
            
          } else {

            gp <- ggplot(subset(dfPlot, Marker==vecColNames[i]))
            gp <- gp + geom_boxplot(aes_string(x="Marker", y="Ratio"))

            if(val_titles){
              
              gp <- .applyPlotSettings(gp = gp, theme = val_theme, main.title = val_title,
                                       x.title = val_xtitle, y.title = val_ytitle)
              
            } else {
              
              gp <- .applyPlotSettings(gp = gp, theme = val_theme, main.title = "Marker ratio",
                                       x.title = "Marker", y.title = "Ratio")
              
            }

          }
          
          # Show plot.          
          print(gp)
          
          # Store in global variable.
          .gPlot <<- gp

        }
        
        # Disable confirm.
        par(ask=F) 
        
        
      } else if (what == "plot") {
        
        if("Group" %in% names(val_data)){
          
          # Create plot.          
          gp <- ggplot(data = dfPlot,
                       aes_string(x="Group", y="Ratio", fill="Group"))
          gp <- gp + geom_boxplot()
          gp <- gp + facet_wrap(as.formula(paste("~ Marker")), scales = val_scales)

          if(val_titles){
            
            gp <- .applyPlotSettings(gp = gp, theme = val_theme, main.title = val_title,
                                     x.title = val_xtitle, y.title = val_ytitle)
            
          } else {
            
            gp <- .applyPlotSettings(gp = gp, theme = val_theme, main.title = "Marker ratio",
                                     x.title = NULL, y.title = "Ratio")
            
          }

        } else {

          # Create plot.          
          gp <- ggplot(data = dfPlot,
                       aes_string(x="Marker", y="Ratio"))
          gp <- gp + geom_boxplot()
          gp <- gp + facet_wrap(as.formula(paste("~ Marker")), scales = val_scales)
          
          if(val_titles){
            
            gp <- .applyPlotSettings(gp = gp, theme = val_theme, main.title = val_title,
                                     x.title = val_xtitle, y.title = val_ytitle)
            
          } else {
            
            gp <- .applyPlotSettings(gp = gp, theme = val_theme, main.title = "Marker ratio",
                                     x.title = "Marker pair", y.title = "Ratio")
            
          }

        }

        # Show plot.        
        print(gp)
        
        # Store in global variable.
        .gPlot <<- gp
        
      }
      
      # Change save button.
      svalue(f5_save_btn) <- "Save as object"
      enabled(f5_save_btn) <- TRUE
      
    } else {
      
      gmessage(message="Data frame is NULL or NA!",
               title="Error",
               icon = "error")      
      
    } 
    
  }

  # INTERNAL FUNCTIONS ########################################################

  .applyPlotSettings <- function(gp, theme = "theme_grey()",
                                 main.title = NULL, x.title = NULL, y.title = NULL){
    
    # Apply theme.
    gp <- gp + eval(parse(text = theme))

    # Apply titles.    
    gp <- gp + labs(title = main.title)
    gp <- gp + xlab(x.title)
    gp <- gp + ylab(y.title)
    
    return(gp)
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
      if(exists(".strvalidator_plotRatio_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_plotRatio_gui_savegui", envir=env)
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
      if(exists(".strvalidator_plotRatio_gui_title", envir=env, inherits = FALSE)){
        svalue(title_edt) <- get(".strvalidator_plotRatio_gui_title", envir=env)
      }
      if(exists(".strvalidator_plotRatio_gui_scales", envir=env, inherits = FALSE)){
        svalue(f1_scales_opt) <- get(".strvalidator_plotRatio_gui_scales", envir=env)
      }
      if(exists(".strvalidator_plotRatio_gui_title_chk", envir=env, inherits = FALSE)){
        svalue(f1_titles_chk) <- get(".strvalidator_plotRatio_gui_title_chk", envir=env)
      }
      if(exists(".strvalidator_plotRatio_gui_x_title", envir=env, inherits = FALSE)){
        svalue(x_title_edt) <- get(".strvalidator_plotRatio_gui_x_title", envir=env)
      }
      if(exists(".strvalidator_plotRatio_gui_y_title", envir=env, inherits = FALSE)){
        svalue(y_title_edt) <- get(".strvalidator_plotRatio_gui_y_title", envir=env)
      }
      if(exists(".strvalidator_plotRatio_gui_theme", envir=env, inherits = FALSE)){
        svalue(f1_theme_drp) <- get(".strvalidator_plotRatio_gui_theme", envir=env)
      }

      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_plotRatio_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_plotRatio_gui_scales", value=svalue(f1_scales_opt), envir=env)
      assign(x=".strvalidator_plotRatio_gui_title", value=svalue(title_edt), envir=env)
      assign(x=".strvalidator_plotRatio_gui_title_chk", value=svalue(f1_titles_chk), envir=env)
      assign(x=".strvalidator_plotRatio_gui_x_title", value=svalue(x_title_edt), envir=env)
      assign(x=".strvalidator_plotRatio_gui_y_title", value=svalue(y_title_edt), envir=env)
      assign(x=".strvalidator_plotRatio_gui_theme", value=svalue(f1_theme_drp), envir=env)

    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_plotRatio_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotRatio_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_plotRatio_gui_scales", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotRatio_gui_scales", envir = env)
      }
      if(exists(".strvalidator_plotRatio_gui_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotRatio_gui_title", envir = env)
      }
      if(exists(".strvalidator_plotRatio_gui_title_chk", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotRatio_gui_title_chk", envir = env)
      }
      if(exists(".strvalidator_plotRatio_gui_x_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotRatio_gui_x_title", envir = env)
      }
      if(exists(".strvalidator_plotRatio_gui_y_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotRatio_gui_y_title", envir = env)
      }
      if(exists(".strvalidator_plotRatio_gui_theme", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotRatio_gui_theme", envir = env)
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
