################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 11.11.2015: Added importFrom ggplot2.
# 29.08.2015: Added importFrom.
# 11.10.2014: Added 'focus', added 'parent' parameter.
# 07.08.2014: Fixed boxplot error 
#  Error: stat_boxplot requires the following missing aesthetics: x, y
# 28.06.2014: Added help button and moved save gui checkbox.
# 08.05.2014: Implemented 'checkDataset'.
# 28.02.2014: Fixed plot object not saved in '.gPlot'.
# 20.01.2014: Changed 'saveImage_gui' for 'ggsave_gui'.
# 19.11.2013: Added distribution plot.
# 28.10.2013: First version.

#' @title Plot Capillary Balance
#'
#' @description
#' GUI simplifying the creation of plots from capillary balance data.
#'
#' @details Select a dataset to plot from the drop-down meny.
#' Plot capillary balance as a dotplot, boxplot or as a distribution.
#' Automatic plot titles can be replaced by custom titles.
#' A name for the result is automatiaclly suggested.
#' The resulting plot can be saved as either a plot object or as an image.
#' @param env environment in wich to search for data frames and save result.
#' @param savegui logical indicating if GUI settings should be saved in the environment.
#' @param debug logical indicating printing debug information.
#' @param parent widget to get focus when finished.
#' 
#' @return TRUE
#' 
#' @export
#' 
#' @importFrom utils help str
#' @importFrom ggplot2 ggplot aes_string facet_grid geom_point geom_line labs
#'  geom_boxplot stat_boxplot theme geom_density coord_cartesian
#' 
#' @seealso \url{http://docs.ggplot2.org/current/} for details on plot settings.

plotCapillary_gui <- function(env=parent.frame(), savegui=NULL, debug=FALSE, parent=NULL){

  # Global variables.
  .gData <- NULL
  .gDataName <- NULL
  .gPlot <- NULL
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # Main window.
  w <- gwindow(title="Plot capillary balance", visible=FALSE)
  
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
    print(help("plotCapillary_gui", help_type="html"))
    
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
  
  f0_samples_lbl <- glabel(text=" (0 rows)", container=f0)

  addHandlerChanged(dataset_drp, handler = function (h, ...) {
    
    val_obj <- svalue(dataset_drp)
    
    # Check if suitable.
    requiredCol <- c("Instrument", "Instrument.ID", "Run",
                     "Mean.Height", "SQ", "Injection", "Capillary", "Well")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       env=env, parent=w, debug=debug)
    
    if(ok){
      
      # Load or change components.
      .gData <<- get(val_obj, envir=env)
      .gDataName <<- val_obj

      # Suggest name.
      svalue(f5_save_edt) <- paste(val_obj, "_ggplot", sep="")
      
      svalue(f0_samples_lbl) <- paste(" (",
                                      nrow(.gData),
                                      " rows)", sep="")

      # Enable buttons.
      enabled(plot_dot_btn) <- TRUE
      enabled(plot_box_btn) <- TRUE
        
    } else {

      # Reset components.
      .gData <<- NULL
      svalue(f5_save_edt) <- ""
      svalue(dataset_drp, index=TRUE) <- 1
      svalue(f0_samples_lbl) <- " (0 rows)"
      
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
      enabled(grid1) <- TRUE
    } else {
      enabled(grid1) <- FALSE
    }
  } )
  
  grid1 <- glayout(container = f1, spacing = 1)
  enabled(grid1) <- svalue(f1_titles_chk)

  grid1[1,1] <- glabel(text="Plot title:", container=grid1)
  grid1[1,2] <- title_edt <- gedit(text="",
                                   width=40,
                                   container=grid1)
  
  grid1[2,1] <- glabel(text="2nd line:", container=grid1)
  grid1[2,2] <- sub_title_edt <- gedit(text="",
                                   width=40,
                                   container=grid1)

  grid1[3,1] <- glabel(text="X title:", container=grid1)
  grid1[3,2] <- x_title_edt <- gedit(text="",
                                     container=grid1)

  grid1[4,1] <- glabel(text="Y title:", container=grid1)
  grid1[4,2] <- y_title_edt <- gedit(text="",
                                     container=grid1)

  # FRAME 7 ###################################################################
  
  f7 <- gframe(text = "Plot capillary balance data",
               horizontal=FALSE,
               container = gv) 
  
  grid7 <- glayout(container = f7)
  
  grid7[1,1] <- plot_dot_btn <- gbutton(text="Dotplot",
                                           border=TRUE,
                                           container=grid7) 
  
  grid7[1,2] <- plot_box_btn <- gbutton(text="Boxplot",
                                       border=TRUE,
                                       container=grid7) 
  
  grid7[1,3] <- plot_dst_btn <- gbutton(text="Distribution",
                                        border=TRUE,
                                        container=grid7) 
  
  addHandlerChanged(plot_dot_btn, handler = function(h, ...) {
    
    # Check if suitable for plot.
    requiredCol <- c("Instrument", "Instrument.ID", "Run", "Mean.Height",
                     "Injection", "Capillary", "Well",	"Comment")
    
    if(!all(requiredCol %in% colnames(.gData))){
      
      missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
      
      message <- paste("Additional columns required:\n",
                       paste(missingCol, collapse="\n"), sep="")
      
      gmessage(message, title="message",
               icon = "error",
               parent = w) 
      
    } else {
      
      enabled(plot_dot_btn) <- FALSE
      .plotBalance(what="dotplot")
      enabled(plot_dot_btn) <- TRUE
      
    }
      
  } )
  
  addHandlerChanged(plot_box_btn, handler = function(h, ...) {
    
    # Check if suitable for plot.
    requiredCol <- c("Instrument", "Instrument.ID", "Run", "Mean.Height",
                     "Injection", "Capillary", "Well",  "Comment")
    
    if(!all(requiredCol %in% colnames(.gData))){
      
      missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
      
      message <- paste("Additional columns required:\n",
                       paste(missingCol, collapse="\n"), sep="")
      
      gmessage(message, title="message",
               icon = "error",
               parent = w) 
      
    } else {
      
      enabled(plot_box_btn) <- FALSE
      .plotBalance(what="boxplot")
      enabled(plot_box_btn) <- TRUE
      
    }
    
  } )

  addHandlerChanged(plot_dst_btn, handler = function(h, ...) {
    
    # Check if suitable for plot.
    requiredCol <- c("Instrument", "Instrument.ID", "Run", "Mean.Height",
                     "Injection", "Capillary", "Well",  "Comment")
    
    if(!all(requiredCol %in% colnames(.gData))){
      
      missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
      
      message <- paste("Additional columns required:\n",
                       paste(missingCol, collapse="\n"), sep="")
      
      gmessage(message, title="message",
               icon = "error",
               parent = w) 
      
    } else {
      
      enabled(plot_dst_btn) <- FALSE
      .plotBalance(what="dstplot")
      enabled(plot_dst_btn) <- TRUE
      
    }
    
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
  grid2[1,2] <- e2_shape_spb <- gspinbutton(from=0, to=25,
                                         by=1, value=18,
                                         container=grid2)

  grid2[1,3] <- glabel(text="Alpha:", container=grid2)
  grid2[1,4] <- e2_alpha_spb <- gspinbutton(from=0, to=1,
                                         by=0.01, value=1,
                                         container=grid2)

  # FRAME 3 ###################################################################

  e3 <- gexpandgroup(text="Axes",
                     horizontal=FALSE,
                     container = f1)
  
  grid3 <- glayout(container = e3, spacing = 1)

  grid3[1,1:2] <- glabel(text="Limit Y axis (min-max)", container=grid3)
  grid3[2,1] <- e3_y_min_edt <- gedit(text="", width=5, container=grid3)
  grid3[2,2] <- e3_y_max_edt <- gedit(text="", width=5, container=grid3)

  # FRAME 4 ###################################################################
  
#   e4 <- gexpandgroup(text="X labels",
#                      horizontal=FALSE,
#                      container = f1)
#   
#   grid4 <- glayout(container = e4)
#   
#   grid4[1,1] <- glabel(text="Text size (pts):", container=grid4)
#   grid4[1,2] <- e4_size_edt <- gedit(text="8", width=4, container=grid4)

  # FUNCTIONS #################################################################
  
  
  .plotBalance <- function(what){
    
    # Get values.
    val_titles <- svalue(f1_titles_chk)
    val_title <- svalue(title_edt)
    val_sub_title <- svalue(sub_title_edt)
    val_x_title <- svalue(x_title_edt)
    val_y_title <- svalue(y_title_edt)
    val_shape <- as.numeric(svalue(e2_shape_spb))
    val_alpha <- as.numeric(svalue(e2_alpha_spb))
    val_ymin <- as.numeric(svalue(e3_y_min_edt))
    val_ymax <- as.numeric(svalue(e3_y_max_edt))
    #val_size <- as.numeric(svalue(e4_size_edt))
    
    if(debug){
      print("val_titles")
      print(val_titles)
      print("val_title")
      print(val_title)
      print("val_sub_title")
      print(val_sub_title)
      print("val_x_title")
      print(val_x_title)
      print("val_y_title")
      print(val_y_title)
      print("val_shape")
      print(val_shape)
      print("val_alpha")
      print(val_alpha)
      print("val_ymin")
      print(val_ymin)
      print("val_ymax")
      print(val_ymax)
#       print("val_size")
#       print(val_size)
      print("str(.gData)")
      print(str(.gData))
    }

    if (!is.na(.gData) && !is.null(.gData)){

      # Remove NA.
      if(any(is.na(.gData$Mean.Height))){
        n0 <- nrow(.gData)
        .gData <- .gData[!is.na(.gData$Mean.Height),]
        n1 <- nrow(.gData)
        print(paste(n0-n1,"NA rows removed from 'Mean.Height'."))
      }
      
      # Height must be numeric.
      if(!is.numeric(.gData$Mean.Height)){
        .gData$Mean.Height <- as.numeric(.gData$Mean.Height)
        message("'Mean.Height' not numeric, converting to numeric!")
        
      }

      # Capillary must be numeric.
      if(!is.numeric(.gData$Capillary)){
        .gData$Capillary <- as.numeric(.gData$Capillary)
        message("'Capillary' not numeric, converting to numeric!")
        
      }

      # Injection must be numeric.
      if(!is.numeric(.gData$Injection)){
        .gData$Injection <- as.numeric(.gData$Injection)
        message("'Injection' not numeric, converting to numeric!")
        
      }
      
      if(debug){
        print("Before plot: str(.gData)")
        print(str(.gData))
      }
      
      
        # Select what to plot.
        if(what == "dotplot"){
          
          if(val_titles){
            mainTitle <- val_title
            subTitle <- val_sub_title
            xTitle <- val_x_title
            yTitle <- val_y_title
          } else {
            mainTitle <- paste("Mean peak height grouped by capillary for ",
                               unique(.gData$Instrument.ID),
                               " (", unique(.gData$Instrument), ")",
                               sep="")
            subTitle <- "[dotted blue line indicate global mean, red line indicate median per capillary]"
            xTitle <- "Injection"
            yTitle <- "Mean peak height (RFU)"
          }
          
          # POINT (best for few replicates (injections))
          gp <- ggplot(.gData, aes_string(x="Injection", y="Mean.Height", group="Capillary"))
          gp <- gp + facet_grid(.~Capillary, space="fixed", scales="fixed")
          gp <- gp + geom_point(shape=val_shape, alpha=val_alpha)
          gp <- gp + geom_line()
          gp <- gp + geom_line(y=mean(.gData$Mean.Height, na.rm=TRUE), lty="dotted", color="blue")
          gp <- gp + geom_line(stat = "hline", yintercept = "median", color="red")
          gp <- gp + labs(title=paste(mainTitle, "\n", subTitle),
                          x=xTitle, y=yTitle)
          
        } else if (what == "boxplot") {
          
          if(val_titles){
            mainTitle <- val_title
            subTitle <- val_sub_title
            xTitle <- val_x_title
            yTitle <- val_y_title
          } else {
            mainTitle <- paste("Mean peak height per capillary for ",
                               unique(.gData$Instrument.ID),
                               " (", unique(.gData$Instrument), ")",
                               sep="")
            subTitle <- ""
            xTitle <- "Capillary"
            yTitle <- "Mean peak height (RFU)"
          }
          
          # BOXPLOT (best suited for many replicates (injections))
          # Note: aes allows the first two arguments to be unnamed and are assumed to be x and y (respectively)
          # Note: aes_string does not have this shortcut, and so all arguments must be named. 
          .gData$Capillary <- factor(.gData$Capillary) 
          gp <- ggplot(.gData, aes_string(x="Capillary", y="Mean.Height"))
          gp <- gp + geom_boxplot()
          gp <- gp + stat_boxplot(geom ="errorbar")
          gp <- gp + labs(title=paste(mainTitle, "\n", subTitle),
                          x=xTitle, y=yTitle)
          gp <- gp + theme(legend.position="none")
          
        } else if (what == "dstplot") {
          
          if(val_titles){
            mainTitle <- val_title
            subTitle <- val_sub_title
            xTitle <- val_x_title
            yTitle <- val_y_title
          } else {
            mainTitle <- paste("Mean peak height for ",
                               unique(.gData$Instrument.ID),
                               " (", unique(.gData$Instrument), ")",
                               sep="")
            subTitle <- ""
            xTitle <- "Mean peak height (RFU)"
            yTitle <- "Density"
          }
          
          # DISTRIBUTION PLOT.
          gp <- ggplot(.gData, aes_string(x = "Mean.Height"))
          gp <- gp + geom_density()
          gp <- gp + geom_point(y=0, aes_string(colour="Run"), shape=108, size=10, alpha=val_alpha)
          gp <- gp + labs(title=paste(mainTitle, "\n", subTitle),
                          x=xTitle, y=yTitle)
          
        }
      
        # Restrict y axis.
        if(!is.na(val_ymin) && !is.na(val_ymax)){
          # Zoom in.
          gp <- gp + coord_cartesian(ylim=c(val_ymin, val_ymax))
        }
      
        # Apply theme.
#         gp <- gp + theme(axis.text.x=element_text(size=val_size))
        
        # plot.
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
      if(exists(".strvalidator_plotCapillary_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_plotCapillary_gui_savegui", envir=env)
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
      if(exists(".strvalidator_plotCapillary_gui_title", envir=env, inherits = FALSE)){
        svalue(title_edt) <- get(".strvalidator_plotCapillary_gui_title", envir=env)
      }
      if(exists(".strvalidator_plotCapillary_gui_title_chk", envir=env, inherits = FALSE)){
        svalue(f1_titles_chk) <- get(".strvalidator_plotCapillary_gui_title_chk", envir=env)
      }
      if(exists(".strvalidator_plotCapillary_gui_sub_title", envir=env, inherits = FALSE)){
        svalue(sub_title_edt) <- get(".strvalidator_plotCapillary_gui_sub_title", envir=env)
      }
      if(exists(".strvalidator_plotCapillary_gui_x_title", envir=env, inherits = FALSE)){
        svalue(x_title_edt) <- get(".strvalidator_plotCapillary_gui_x_title", envir=env)
      }
      if(exists(".strvalidator_plotCapillary_gui_y_title", envir=env, inherits = FALSE)){
        svalue(y_title_edt) <- get(".strvalidator_plotCapillary_gui_y_title", envir=env)
      }
      if(exists(".strvalidator_plotCapillary_gui_points_shape", envir=env, inherits = FALSE)){
        svalue(e2_shape_spb) <- get(".strvalidator_plotCapillary_gui_points_shape", envir=env)
      }
      if(exists(".strvalidator_plotCapillary_gui_points_alpha", envir=env, inherits = FALSE)){
        svalue(e2_alpha_spb) <- get(".strvalidator_plotCapillary_gui_points_alpha", envir=env)
      }
      if(exists(".strvalidator_plotCapillary_gui_axes_y_min", envir=env, inherits = FALSE)){
        svalue(e3_y_min_edt) <- get(".strvalidator_plotCapillary_gui_axes_y_min", envir=env)
      }
      if(exists(".strvalidator_plotCapillary_gui_axes_y_max", envir=env, inherits = FALSE)){
        svalue(e3_y_max_edt) <- get(".strvalidator_plotCapillary_gui_axes_y_max", envir=env)
      }
#       if(exists(".strvalidator_plotCapillary_gui_xlabel_size", envir=env, inherits = FALSE)){
#         svalue(e4_size_edt) <- get(".strvalidator_plotCapillary_gui_xlabel_size", envir=env)
#       }
      
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_plotCapillary_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_plotCapillary_gui_title_chk", value=svalue(f1_titles_chk), envir=env)
      assign(x=".strvalidator_plotCapillary_gui_title", value=svalue(title_edt), envir=env)
      assign(x=".strvalidator_plotCapillary_gui_sub_title", value=svalue(sub_title_edt), envir=env)
      assign(x=".strvalidator_plotCapillary_gui_x_title", value=svalue(x_title_edt), envir=env)
      assign(x=".strvalidator_plotCapillary_gui_y_title", value=svalue(y_title_edt), envir=env)
      assign(x=".strvalidator_plotCapillary_gui_points_shape", value=svalue(e2_shape_spb), envir=env)
      assign(x=".strvalidator_plotCapillary_gui_points_alpha", value=svalue(e2_alpha_spb), envir=env)
      assign(x=".strvalidator_plotCapillary_gui_axes_y_min", value=svalue(e3_y_min_edt), envir=env)
      assign(x=".strvalidator_plotCapillary_gui_axes_y_max", value=svalue(e3_y_max_edt), envir=env)
#      assign(x=".strvalidator_plotCapillary_gui_xlabel_size", value=svalue(e4_size_edt), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_plotCapillary_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotCapillary_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_plotCapillary_gui_title_chk", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotCapillary_gui_title_chk", envir = env)
      }
      if(exists(".strvalidator_plotCapillary_gui_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotCapillary_gui_title", envir = env)
      }
      if(exists(".strvalidator_plotCapillary_gui_sub_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotCapillary_gui_sub_title", envir = env)
      }
      if(exists(".strvalidator_plotCapillary_gui_x_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotCapillary_gui_x_title", envir = env)
      }
      if(exists(".strvalidator_plotCapillary_gui_y_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotCapillary_gui_y_title", envir = env)
      }
      if(exists(".strvalidator_plotCapillary_gui_points_shape", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotCapillary_gui_points_shape", envir = env)
      }
      if(exists(".strvalidator_plotCapillary_gui_points_alpha", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotCapillary_gui_points_alpha", envir = env)
      }
      if(exists(".strvalidator_plotCapillary_gui_axes_y_min", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotCapillary_gui_axes_y_min", envir = env)
      }
      if(exists(".strvalidator_plotCapillary_gui_axes_y_max", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotCapillary_gui_axes_y_max", envir = env)
      }
#       if(exists(".strvalidator_plotCapillary_gui_xlabel_size", envir=env, inherits = FALSE)){
#         remove(".strvalidator_plotCapillary_gui_xlabel_size", envir = env)
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
  
  # Show GUI.
  visible(w) <- TRUE
  focus(w)
  
}
