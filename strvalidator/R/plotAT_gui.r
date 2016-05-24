################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 11.11.2015: Added importFrom ggplot2.
# 29.08.2015: Added importFrom.
# 28.06.2015: Changed confidence interval level to match one-sided critical t-value.
# 01.06.2015: First version.

#' @title Plot Analytical Threshold
#'
#' @description
#' GUI simplifying the creation of plots from analytical threshold data.
#'
#' @details Select data to plot in the drop-down menu. Plot regression data
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
#' @importFrom utils help str
#' @importFrom ggplot2 stat_smooth geom_abline xlim ggplot aes_string geom_point
#' position_jitter coord_cartesian theme element_text labs xlab ylab
#' 
#' @seealso \url{http://docs.ggplot2.org/current/} for details on plot settings.

plotAT_gui <- function(env=parent.frame(), savegui=NULL, debug=FALSE, parent=NULL){
  
  # Global variables.
  .gData <- NULL
  .gDataName <- NULL
  .gPlot <- NULL
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # Main window.
  w <- gwindow(title="Plot analytical threshold", visible=FALSE)
  
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
    print(help("plotAT_gui", help_type="html"))
    
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
    requiredCol <- c("Amount", "Height", "AT6")
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
      enabled(plot_at6_btn) <- TRUE
      
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
  f1g2[1,2] <- f1_theme_drp <- gdroplist(items=c("theme_grey()","theme_bw()"),
                                         selected=1,
                                         container = f1g2)
  
  # FRAME 7 ###################################################################
  
  f7 <- gframe(text = "Plot analytical threshold data",
               horizontal=FALSE,
               container = gv) 
  
  grid7 <- glayout(container = f7)
  
  grid7[1,1] <- plot_at6_btn <- gbutton(text="Plot AT6",
                                           border=TRUE,
                                           container=grid7) 
  
  addHandlerChanged(plot_at6_btn, handler = function(h, ...) {
    
    enabled(plot_at6_btn) <- FALSE
    .plotAT(what="AT6")
    enabled(plot_at6_btn) <- TRUE

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
  
  # ADVANCED OPTIONS ##########################################################
  
  e2 <- gexpandgroup(text="Data points",
               horizontal=FALSE,
               container = f1)
  
  grid2 <- glayout(container = e2)
  
  grid2[1,1] <- glabel(text="Shape:", container=grid2)
  grid2[1,2] <- shape_spb <- gspinbutton(from=0, to=25,
                                         by=1, value=18,
                                         container=grid2)
  
  grid2[1,3] <- glabel(text="Alpha:", container=grid2)
  grid2[1,4] <- alpha_spb <- gspinbutton(from=0, to=1,
                                         by=0.01, value=1,
                                         container=grid2)

  grid2[1,5] <- glabel(text="Jitter (width):", container=grid2)
  grid2[1,6] <- jitter_txt <- gedit(text="0", width=4, container=grid2)

  # FRAME 3 ###################################################################

  e3 <- gexpandgroup(text="Axes",
                     horizontal=FALSE,
                     container = f1)
  
  grid3 <- glayout(container = e3, spacing = 1)
  
  grid3[1,1:2] <- glabel(text="Limit Y axis (min-max)", container=grid3)
  grid3[2,1] <- y_min_txt <- gedit(text="", width=5, container=grid3)
  grid3[2,2] <- y_max_txt <- gedit(text="", width=5, container=grid3)

  grid3[3,1:2] <- glabel(text="Limit X axis (min-max)", container=grid3)
  grid3[4,1] <- x_min_txt <- gedit(text="", width=5, container=grid3)
  grid3[4,2] <- x_max_txt <- gedit(text="", width=5, container=grid3)

  grid3[1,3] <- glabel(text="    ", container=grid3) # Add some space.
  
  grid3[1,4] <- glabel(text="Scales:", container=grid3)
  grid3[2:4,4] <- scales_opt <- gradio(items=c("fixed","free_x","free_y","free"),
                                      selected = 2,
                                      horizontal = FALSE,
                                      container = grid3)
  
  # FRAME 4 ###################################################################
  
  e4 <- gexpandgroup(text="X labels",
                     horizontal=FALSE,
                     container = f1)
  
  grid4 <- glayout(container = e4)
  
  grid4[1,1] <- glabel(text="Text size (pts):", container=grid4)
  grid4[1,2] <- size_txt <- gedit(text="10", width=4, container=grid4)

  grid4[1,3] <- glabel(text="Angle:", container=grid4)
  grid4[1,4] <- angle_spb <- gspinbutton (from=0, to=360, by=1,
                                         value=270,
                                         container=grid4) 

  grid4[2,1] <- glabel(text="Justification (v/h):", container=grid4)
  grid4[2,2] <- vjust_spb <- gspinbutton (from=0, to=1, by=0.1,
                                          value=0.5,
                                          container=grid4)

  grid4[2,3] <- hjust_spb <- gspinbutton (from=0, to=1, by=0.1,
                                          value=0,
                                          container=grid4)

  
  
  # FUNCTIONS #################################################################
  
  
  .plotAT <- function(what){
    
    # Get values.
    val_titles <- svalue(f1_titles_chk)
    val_title <- svalue(title_edt)
    val_xtitle <- svalue(x_title_edt)
    val_ytitle <- svalue(y_title_edt)
    val_shape <- as.numeric(svalue(shape_spb))
    val_alpha <- as.numeric(svalue(alpha_spb))
    val_jitter <- as.numeric(svalue(jitter_txt))
    val_ymin <- as.numeric(svalue(y_min_txt))
    val_ymax <- as.numeric(svalue(y_max_txt))
    val_xmin <- as.numeric(svalue(x_min_txt))
    val_xmax <- as.numeric(svalue(x_max_txt))
    val_angle <- as.numeric(svalue(angle_spb))
    val_vjust <- as.numeric(svalue(vjust_spb))
    val_hjust <- as.numeric(svalue(hjust_spb))
    val_size <- as.numeric(svalue(size_txt))
    val_scales <- svalue(scales_opt)
    val_theme <- svalue(f1_theme_drp)
    
    if(debug){
      print("ARGUMENTS:")
      print("val_title")
      print(val_title)
      print("val_xtitle")
      print(val_xtitle)
      print("val_ytitle")
      print(val_ytitle)
      print("val_shape")
      print(val_shape)
      print("val_alpha")
      print(val_alpha)
      print("val_jitter")
      print(val_jitter)
      print("val_ymin")
      print(val_ymin)
      print("val_ymax")
      print(val_ymax)
      print("val_angle")
      print(val_angle)
      print("val_vjust")
      print(val_vjust)
      print("val_hjust")
      print(val_hjust)
      print("val_size")
      print(val_size)
      print("val_scales")
      print(val_scales)
      print("str(.gData)")
      print(str(.gData))
      print("val_theme")
      print(val_theme)
    }
    
    if (!is.na(.gData) && !is.null(.gData)){
      
      # Plotting data and regression for AT6.
      if(what == "AT6"){
        
        if(val_titles){
          mainTitle <- val_title
          xTitle <- val_xtitle
          yTitle <- val_ytitle
        } else {
          if(all(is.na(.gData$Weight))){
            mainTitle <- "Linear regression"
          } else {
            mainTitle <- "Weighted linear regression"
          }
          subTitle <- paste("AT6:", round(unique(.gData$AT6),0), "(RFU)")
          xTitle <- "Amount (pg)"
          yTitle <- "Average Peak Height (RFU)"
        }

        # Create plot.
        gp <- ggplot(data=.gData, aes_string(x="Amount", y="Height"))

        # Get 
        npoints <- unique(.gData$N)
        alpha <- unique(.gData$Alpha)
        atinterc <- unique(.gData$AT6)
        
        if(all(is.na(.gData$Weight))){

          # Add regression line.
          gp <- gp + stat_smooth(aes_string(x="Amount", y="Height"),
                                 method="lm", se=TRUE, n=npoints,
                                 fullrange=TRUE, level=1-alpha*2)

          } else {

          # Add weighted regression line.
          gp <- gp + stat_smooth(aes_string(x="Amount", y="Height", weight="Weight"),
                                 method="lm", se=TRUE, n=npoints,
                                 fullrange=TRUE, level=1-alpha*2)

        }
        
        # Addthreshold line.
        gp <- gp + geom_abline(intercept=atinterc, 
                               slope=0, linetype="dotted")
        
        # Set x-axis to extend regression line.
        gp <- gp + xlim(0, max(.gData$Amount))
        
      }
      
      # Apply theme.
      gp <- gp + eval(parse(text=val_theme))
      
      # Plot settings.
      gp <- gp + geom_point(shape=val_shape, alpha=val_alpha, 
                            position=position_jitter(height = 0, width=val_jitter)) 

      # Restrict y axis.
      if(!is.na(val_ymin) && !is.na(val_ymax)){
        val_y <- c(val_ymin, val_ymax)
      } else {
        val_y <- NULL
      }
      # Restrict x axis.
      if(!is.na(val_xmin) && !is.na(val_xmax)){
        val_x <- c(val_xmin, val_xmax)
      } else {
        val_x <- NULL
      }
      # Zoom in without dropping observations.
      gp <- gp + coord_cartesian(xlim=val_x, ylim=val_y)

      if(debug){
        print(paste("Plot zoomed to xlim:", val_x, "ylim:", val_y))
      }
      
      # Titles.
      gp <- gp + theme(axis.text.x=element_text(angle=val_angle,
                                                hjust=val_hjust,
                                                vjust=val_vjust,
                                                size=val_size))
      gp <- gp + labs(title=paste(mainTitle, "\n", subTitle))
      gp <- gp + xlab(xTitle)
      gp <- gp + ylab(yTitle)

      # Restrict y axis.
      if(!is.na(val_ymin) && !is.na(val_ymax)){
        val_y <- c(val_ymin, val_ymax)
      } else {
        val_y <- NULL
      }
      # Restrict x axis.
      if(!is.na(val_xmin) && !is.na(val_xmax)){
        val_x <- c(val_xmin, val_xmax)
      } else {
        val_x <- NULL
      }
      # Zoom in without dropping observations.
      gp <- gp + coord_cartesian(xlim=val_x, ylim=val_y)
      
      if(debug){
        print(paste("Plot zoomed to xlim:", val_x, "ylim:", val_y))
      }
      
      # Show plot.        
      print(gp)
      
      # Change save button.
      svalue(f5_save_btn) <- "Save as object"
      enabled(f5_save_btn) <- TRUE
    
      # Store in global variable.
      .gPlot <<- gp
      
    } else {
      
      gmessage(message="Data frame is NULL or NA!",
               title="Error",
               icon = "error")      
      
    } 
    
  }

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
      if(exists(".strvalidator_plotAT_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_plotAT_gui_savegui", envir=env)
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
      if(exists(".strvalidator_plotAT_gui_title", envir=env, inherits = FALSE)){
        svalue(title_edt) <- get(".strvalidator_plotAT_gui_title", envir=env)
      }
      if(exists(".strvalidator_plotAT_gui_title_chk", envir=env, inherits = FALSE)){
        svalue(f1_titles_chk) <- get(".strvalidator_plotAT_gui_title_chk", envir=env)
      }
      if(exists(".strvalidator_plotAT_gui_x_title", envir=env, inherits = FALSE)){
        svalue(x_title_edt) <- get(".strvalidator_plotAT_gui_x_title", envir=env)
      }
      if(exists(".strvalidator_plotAT_gui_y_title", envir=env, inherits = FALSE)){
        svalue(y_title_edt) <- get(".strvalidator_plotAT_gui_y_title", envir=env)
      }
      if(exists(".strvalidator_plotAT_gui_points_shape", envir=env, inherits = FALSE)){
        svalue(shape_spb) <- get(".strvalidator_plotAT_gui_points_shape", envir=env)
      }
      if(exists(".strvalidator_plotAT_gui_points_alpha", envir=env, inherits = FALSE)){
        svalue(alpha_spb) <- get(".strvalidator_plotAT_gui_points_alpha", envir=env)
      }
      if(exists(".strvalidator_plotAT_gui_points_jitter", envir=env, inherits = FALSE)){
        svalue(jitter_txt) <- get(".strvalidator_plotAT_gui_points_jitter", envir=env)
      }
      if(exists(".strvalidator_plotAT_gui_axes_y_min", envir=env, inherits = FALSE)){
        svalue(y_min_txt) <- get(".strvalidator_plotAT_gui_axes_y_min", envir=env)
      }
      if(exists(".strvalidator_plotAT_gui_axes_y_max", envir=env, inherits = FALSE)){
        svalue(y_max_txt) <- get(".strvalidator_plotAT_gui_axes_y_max", envir=env)
      }
      if(exists(".strvalidator_plotAT_gui_axes_x_min", envir=env, inherits = FALSE)){
        svalue(x_min_txt) <- get(".strvalidator_plotAT_gui_axes_x_min", envir=env)
      }
      if(exists(".strvalidator_plotAT_gui_axes_x_max", envir=env, inherits = FALSE)){
        svalue(x_max_txt) <- get(".strvalidator_plotAT_gui_axes_x_max", envir=env)
      }
      if(exists(".strvalidator_plotAT_gui_axes_scales", envir=env, inherits = FALSE)){
        svalue(scales_opt) <- get(".strvalidator_plotAT_gui_axes_scales", envir=env)
      }
      if(exists(".strvalidator_plotAT_gui_xlabel_size", envir=env, inherits = FALSE)){
        svalue(size_txt) <- get(".strvalidator_plotAT_gui_xlabel_size", envir=env)
      }
      if(exists(".strvalidator_plotAT_gui_xlabel_angle", envir=env, inherits = FALSE)){
        svalue(angle_spb) <- get(".strvalidator_plotAT_gui_xlabel_angle", envir=env)
      }
      if(exists(".strvalidator_plotAT_gui_xlabel_justh", envir=env, inherits = FALSE)){
        svalue(hjust_spb) <- get(".strvalidator_plotAT_gui_xlabel_justh", envir=env)
      }
      if(exists(".strvalidator_plotAT_gui_xlabel_justv", envir=env, inherits = FALSE)){
        svalue(vjust_spb) <- get(".strvalidator_plotAT_gui_xlabel_justv", envir=env)
      }
      if(exists(".strvalidator_plotAT_gui_theme", envir=env, inherits = FALSE)){
        svalue(f1_theme_drp) <- get(".strvalidator_plotAT_gui_theme", envir=env)
      }
      
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_plotAT_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_plotAT_gui_title", value=svalue(title_edt), envir=env)
      assign(x=".strvalidator_plotAT_gui_title_chk", value=svalue(f1_titles_chk), envir=env)
      assign(x=".strvalidator_plotAT_gui_x_title", value=svalue(x_title_edt), envir=env)
      assign(x=".strvalidator_plotAT_gui_y_title", value=svalue(y_title_edt), envir=env)
      assign(x=".strvalidator_plotAT_gui_points_shape", value=svalue(shape_spb), envir=env)
      assign(x=".strvalidator_plotAT_gui_points_alpha", value=svalue(alpha_spb), envir=env)
      assign(x=".strvalidator_plotAT_gui_points_jitter", value=svalue(jitter_txt), envir=env)
      assign(x=".strvalidator_plotAT_gui_axes_y_min", value=svalue(y_min_txt), envir=env)
      assign(x=".strvalidator_plotAT_gui_axes_y_max", value=svalue(y_max_txt), envir=env)
      assign(x=".strvalidator_plotAT_gui_axes_x_min", value=svalue(x_min_txt), envir=env)
      assign(x=".strvalidator_plotAT_gui_axes_x_max", value=svalue(x_max_txt), envir=env)
      assign(x=".strvalidator_plotAT_gui_axes_scales", value=svalue(scales_opt), envir=env)
      assign(x=".strvalidator_plotAT_gui_xlabel_size", value=svalue(size_txt), envir=env)
      assign(x=".strvalidator_plotAT_gui_xlabel_angle", value=svalue(angle_spb), envir=env)
      assign(x=".strvalidator_plotAT_gui_xlabel_justh", value=svalue(hjust_spb), envir=env)
      assign(x=".strvalidator_plotAT_gui_xlabel_justv", value=svalue(vjust_spb), envir=env)
      assign(x=".strvalidator_plotAT_gui_theme", value=svalue(f1_theme_drp), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_plotAT_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotAT_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_plotAT_gui_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotAT_gui_title", envir = env)
      }
      if(exists(".strvalidator_plotAT_gui_title_chk", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotAT_gui_title_chk", envir = env)
      }
      if(exists(".strvalidator_plotAT_gui_x_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotAT_gui_x_title", envir = env)
      }
      if(exists(".strvalidator_plotAT_gui_y_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotAT_gui_y_title", envir = env)
      }
      if(exists(".strvalidator_plotAT_gui_points_shape", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotAT_gui_points_shape", envir = env)
      }
      if(exists(".strvalidator_plotAT_gui_points_alpha", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotAT_gui_points_alpha", envir = env)
      }
      if(exists(".strvalidator_plotAT_gui_points_jitter", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotAT_gui_points_jitter", envir = env)
      }
      if(exists(".strvalidator_plotAT_gui_axes_y_min", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotAT_gui_axes_y_min", envir = env)
      }
      if(exists(".strvalidator_plotAT_gui_axes_y_max", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotAT_gui_axes_y_max", envir = env)
      }
      if(exists(".strvalidator_plotAT_gui_axes_x_min", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotAT_gui_axes_x_min", envir = env)
      }
      if(exists(".strvalidator_plotAT_gui_axes_x_max", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotAT_gui_axes_x_max", envir = env)
      }
      if(exists(".strvalidator_plotAT_gui_axes_scales", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotAT_gui_axes_scales", envir = env)
      }
      if(exists(".strvalidator_plotAT_gui_xlabel_size", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotAT_gui_xlabel_size", envir = env)
      }
      if(exists(".strvalidator_plotAT_gui_xlabel_angle", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotAT_gui_xlabel_angle", envir = env)
      }
      if(exists(".strvalidator_plotAT_gui_xlabel_justh", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotAT_gui_xlabel_justh", envir = env)
      }
      if(exists(".strvalidator_plotAT_gui_xlabel_justv", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotAT_gui_xlabel_justv", envir = env)
      }
      if(exists(".strvalidator_plotAT_gui_theme", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotAT_gui_theme", envir = env)
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
