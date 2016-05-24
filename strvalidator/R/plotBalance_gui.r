################################################################################
# TODO LIST
# TODO: ...

# NOTE: Column names used for calculations with data.table is declared
# in globals.R to avoid NOTES in R CMD CHECK.

################################################################################
# CHANGE LOG (last 20 changes)
# 06.01.2016: Fixed theme methods not found and added more themes.
# 04.01.2016: Fixed error object 'val_ncol' not found when val_wrap=1.
# 30.12.2015: Wrapping options changed to radio button and implemented by Dye.
# 30.12.2015: Changed default for drop sex markers to FALSE.
# 30.12.2015: Wrapped 'is.numeric' with checking that columns exist.
# 19.11.2015: Changed axes default to 'fixed' to avoid common plot error.
# 11.11.2015: Added importFrom gridExtra arrangeGrob, and ggplot2.
# 11.11.2015: Added more themes.
# 08.11.2015: Added new plot options 'Hb vs. Marker' and 'Lb vs. Marker'.
# 08.11.2015: Added options to plot all data 'Facet per marker and wrap by dye'. 
# 29.08.2015: Added importFrom.
# 17.04.2015: Included a check for plot facet error caused by all NA's.
# 14.12.2014: Updated to handle gender -> sex.marker option in getKit.
# 11.10.2014: Added 'focus', added 'parent' parameter.
# 28.06.2014: Added help button and moved save gui checkbox.
# 07.05.2014: Implemented new column 'TPH' for X-axis when plot Lb.
# 06.05.2014: Implemented 'checkDataset'.
# 05.05.2014: Fixed 'drop gender' and 'plot log' settings not saved.
# 02.05.2014: Fixed bug when 'drop gender marker=TRUE' and no defined marker.
# 25.05.2014: Fixed scale_colour_manual use first color for all in complex plot.
# 15.04.2014: Fixed now handle no observation in an entire dye channel.
# 14.04.2014: Fixed position_jitter height now fixed to zero (prev. default).
# 23.02.2014: Fixed different y max for complex plot, when supposed to be fixed.
# 23.02.2014: Implemented theme.

#' @title Plot Balance
#'
#' @description
#' GUI simplifying the creation of plots from balance data.
#'
#' @details Select a dataset to plot and the typing kit used (if not autodetected).
#' Plot heterozygote peak balance versus the average locus peak height,
#' the average profile peak height 'H', or by the difference in repeat units
#' (delta). Plot inter-locus balance versus the average locus peak height, or
#' the average profile peak height 'H'. Automatic plot titles can be replaced by
#' custom titles. Sex markers can be excluded. It is possible to plot
#' logarithmic ratios. A name for the result is automatiaclly suggested.
#' The resulting plot can be saved as either a plot object or as an image.
#' @param env environment in wich to search for data frames and save result.
#' @param savegui logical indicating if GUI settings should be saved in the environment.
#' @param debug logical indicating printing debug information.
#' @param parent widget to get focus when finished.
#' 
#' @export
#' 
#' @importFrom utils help str
#' @importFrom stats as.formula
#' @importFrom gridExtra arrangeGrob
#' @importFrom grid unit textGrob grid.newpage grid.draw
# @importFrom gtable gtable_add_grob gtable
#' @importFrom plyr rbind.fill
#' @importFrom  data.table data.table
#' @importFrom ggplot2 ggplot aes_string geom_boxplot geom_point position_jitter
#'  facet_grid facet_wrap scale_colour_manual coord_cartesian guides guide_legend
#'  theme element_text labs xlab ylab element_blank ggplotGrob theme_gray
#'  theme_bw theme_linedraw theme_light theme_dark theme_minimal theme_classic
#'  theme_void 
#' 
#' @return TRUE
#' 
#' @seealso \url{http://docs.ggplot2.org/current/} for details on plot settings.

plotBalance_gui <- function(env=parent.frame(), savegui=NULL, debug=FALSE, parent=NULL){

  # Global variables.
  .gData <- NULL
  .gDataName <- NULL
  .gPlot <- NULL
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # Main window.
  w <- gwindow(title="Plot balance", visible=FALSE)
  
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
    print(help("plotBalance_gui", help_type="html"))
    
  })
  
  # FRAME 0 ###################################################################
  
  f0 <- gframe(text = "Dataset and kit",
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

  glabel(text=" and the kit used:", container=f0)

  kit_drp <- gdroplist(items=getKit(), 
                           selected = 1,
                           editable = FALSE,
                           container = f0) 

  addHandlerChanged(dataset_drp, handler = function (h, ...) {
    
    val_obj <- svalue(dataset_drp)

    # Check if suitable.
    requiredCol <- c("Sample.Name", "Marker")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       env=env, parent=w, debug=debug)
    
    if(ok){
      # Load or change components.
      
      # Get data.
      .gData <<- get(val_obj, envir=env)
      .gDataName <<- val_obj

      # Suggest name.
      svalue(f5_save_edt) <- paste(val_obj, "_ggplot", sep="")
      
      svalue(f0_samples_lbl) <- paste(" (",
                                      length(unique(.gData$Sample.Name)),
                                      " samples)", sep="")

      # Detect kit.
      kitIndex <- detectKit(.gData, index=TRUE)
      # Select in dropdown.
      svalue(kit_drp, index=TRUE) <- kitIndex
      
      # Enable buttons.
      enabled(plot_hb_btn) <- TRUE
      enabled(plot_hb_d_btn) <- TRUE
      enabled(plot_hb_h_btn) <- TRUE
      enabled(plot_lb_btn) <- TRUE
      enabled(plot_lb_h_btn) <- TRUE
        
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
  items_theme <- c("theme_grey()","theme_bw()","theme_linedraw()",
                   "theme_light()","theme_dark()","theme_minimal()",
                   "theme_classic()","theme_void()")
  f1g2[1,2] <- f1_theme_drp <- gdroplist(items = items_theme,
                                         selected = 1,
                                         container = f1g2)

  f1_drop_chk <- gcheckbox(text="Drop sex markers", checked=FALSE, container=f1)
  
  f1_logHb_chk <- gcheckbox(text="Plot Log(balance)", checked=FALSE,
                            container=f1)
  
  f1_wrap_opt <- gradio(items = c("Do not facet or wrap",
                                  "Wrap by Dye",
                                  "Facet by Marker and wrap by Dye"),
                        horizontal = FALSE, container = f1)
  
  # FRAME 7 ###################################################################
  
  f7 <- gframe(text = "Plot Balance data",
               horizontal=FALSE,
               container = gv) 
  
  grid7 <- glayout(container = f7)
  
  grid7[1,1] <- plot_hb_btn <- gbutton(text="Hb vs. Height",
                                           border=TRUE,
                                           container=grid7) 
  
  grid7[1,2] <- plot_hb_d_btn <- gbutton(text="Hb vs. Delta",
                                       border=TRUE,
                                       container=grid7) 
  
  grid7[1,3] <- plot_hb_h_btn <- gbutton(text="Hb vs. 'H'",
                                       border=TRUE,
                                       container=grid7) 
  
  grid7[1,4] <- plot_hb_m_btn <- gbutton(text="Hb vs. Marker",
                                         border=TRUE,
                                         container=grid7) 
  
  grid7[1,5] <- plot_lb_btn <- gbutton(text="Lb vs. Height",
                                       border=TRUE,
                                       container=grid7) 

  grid7[1,6] <- plot_lb_h_btn <- gbutton(text="Lb vs. 'H'",
                                       border=TRUE,
                                       container=grid7) 
  
  grid7[1,7] <- plot_lb_m_btn <- gbutton(text="Lb vs. Marker",
                                         border=TRUE,
                                         container=grid7) 
  
  addHandlerChanged(plot_hb_btn, handler = function(h, ...) {
    
    # Check if suitable for plot.
    requiredCol <- c("Sample.Name", "Marker", "Hb", "MPH")
    
    if(!all(requiredCol %in% colnames(.gData))){
      
      missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
      
      message <- paste("Additional columns required:\n",
                       paste(missingCol, collapse="\n"), sep="")
      
      gmessage(message, title="message",
               icon = "error",
               parent = w) 
      
    } else {
      
      enabled(plot_hb_btn) <- FALSE
      .plotBalance(what="Hb")
      enabled(plot_hb_btn) <- TRUE
      
    }
      
  } )
  
  addHandlerChanged(plot_hb_d_btn, handler = function(h, ...) {
    
    # Check if suitable for plot.
    requiredCol <- c("Sample.Name", "Marker", "Delta", "Hb")
    
    if(!all(requiredCol %in% colnames(.gData))){
      
      missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
      
      message <- paste("Additional columns required:\n",
                       paste(missingCol, collapse="\n"), sep="")
      
      gmessage(message, title="message",
               icon = "error",
               parent = w) 
      
    } else {
      
      enabled(plot_hb_d_btn) <- FALSE
      .plotBalance(what="Hb_D")
      enabled(plot_hb_d_btn) <- TRUE
      
    }
    
  } )
  
  addHandlerChanged(plot_hb_h_btn, handler = function(h, ...) {
    
    # Check if suitable for plot.
    requiredCol <- c("Sample.Name", "Marker", "Hb", "H")
    
    if(!all(requiredCol %in% colnames(.gData))){
      
      missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
      
      message <- paste("Additional columns required:\n",
                       paste(missingCol, collapse="\n"), sep="")
      
      gmessage(message, title="message",
               icon = "error",
               parent = w) 
      
    } else {
      
      enabled(plot_hb_h_btn) <- FALSE
      .plotBalance(what="Hb_H")
      enabled(plot_hb_h_btn) <- TRUE
      
    }
    
  } )

  addHandlerChanged(plot_hb_m_btn, handler = function(h, ...) {
    
    # Check if suitable for plot.
    requiredCol <- c("Sample.Name", "Marker", "Hb")
    
    if(!all(requiredCol %in% colnames(.gData))){
      
      missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
      
      message <- paste("Additional columns required:\n",
                       paste(missingCol, collapse="\n"), sep="")
      
      gmessage(message, title="message",
               icon = "error",
               parent = w) 
      
    } else {
      
      enabled(plot_hb_m_btn) <- FALSE
      .plotBalance(what="Hb_M", complex = FALSE)
      enabled(plot_hb_m_btn) <- TRUE
      
    }
    
  } )
  
  addHandlerChanged(plot_lb_btn, handler = function(h, ...) {
    
    # Check if suitable for plot.
    requiredCol <- c("Sample.Name", "Marker", "Lb", "TPH")
    
    if(!all(requiredCol %in% colnames(.gData))){
      
      missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
      
      message <- paste("Additional columns required:\n",
                       paste(missingCol, collapse="\n"), sep="")
      
      gmessage(message, title="message",
               icon = "error",
               parent = w) 
      
    } else {
      
      enabled(plot_lb_btn) <- FALSE
      .plotBalance(what="Lb")
      enabled(plot_lb_btn) <- TRUE
      
    }
    
  } )

  addHandlerChanged(plot_lb_h_btn, handler = function(h, ...) {
    
    # Check if suitable for plot.
    requiredCol <- c("Sample.Name", "Marker", "Lb", "H")
    
    if(!all(requiredCol %in% colnames(.gData))){
      
      missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
      
      message <- paste("Additional columns required:\n",
                       paste(missingCol, collapse="\n"), sep="")
      
      gmessage(message, title="message",
               icon = "error",
               parent = w) 
      
    } else {
      
      enabled(plot_lb_h_btn) <- FALSE
      .plotBalance(what="Lb_H")
      enabled(plot_lb_h_btn) <- TRUE
      
    }
    
  } )

  addHandlerChanged(plot_lb_m_btn, handler = function(h, ...) {
    
    # Check if suitable for plot.
    requiredCol <- c("Sample.Name", "Marker", "Lb")
    
    if(!all(requiredCol %in% colnames(.gData))){
      
      missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
      
      message <- paste("Additional columns required:\n",
                       paste(missingCol, collapse="\n"), sep="")
      
      gmessage(message, title="message",
               icon = "error",
               parent = w) 
      
    } else {
      
      enabled(plot_lb_m_btn) <- FALSE
      .plotBalance(what="Lb_M", complex = FALSE)
      enabled(plot_lb_m_btn) <- TRUE
      
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
                                         by=0.01, value=0.60,
                                         container=grid2)

  grid2[1,5] <- glabel(text="Jitter (width):", container=grid2)
  grid2[1,6] <- e2_jitter_edt <- gedit(text="0", width=4, container=grid2)

  # FRAME 3 ###################################################################

  e3 <- gexpandgroup(text="Axes",
                     horizontal=FALSE,
                     container = f1)
  
  grid3 <- glayout(container = e3, spacing = 1)

  grid3[1,1:2] <- glabel(text="Limit Y axis (min-max)", container=grid3)
  grid3[2,1] <- e3_y_min_edt <- gedit(text="", width=5, container=grid3)
  grid3[2,2] <- e3_y_max_edt <- gedit(text="", width=5, container=grid3)

  grid3[3,1:2] <- glabel(text="Limit X axis (min-max)", container=grid3)
  grid3[4,1] <- e3_x_min_edt <- gedit(text="", width=5, container=grid3)
  grid3[4,2] <- e3_x_max_edt <- gedit(text="", width=5, container=grid3)

  grid3[1,3] <- glabel(text="    ", container=grid3) # Add some space.
  
  grid3[1,4] <- glabel(text="Scales:", container=grid3)
  grid3[2:4,4] <- e3_scales_opt <- gradio(items=c("fixed","free_x","free_y","free"),
                                      selected = 1,
                                      horizontal = FALSE,
                                      container = grid3)
  
  # FRAME 4 ###################################################################
  
  e4 <- gexpandgroup(text="X labels",
                     horizontal=FALSE,
                     container = f1)
  
  grid4 <- glayout(container = e4)
  
  grid4[1,1] <- glabel(text="Text size (pts):", container=grid4)
  grid4[1,2] <- e4_size_edt <- gedit(text="10", width=4, container=grid4)

  grid4[1,3] <- glabel(text="Angle:", container=grid4)
  grid4[1,4] <- e4_angle_spb <- gspinbutton (from=0, to=360, by=1,
                                         value=270,
                                         container=grid4) 

  grid4[2,1] <- glabel(text="Justification (v/h):", container=grid4)
  grid4[2,2] <- e4_vjust_spb <- gspinbutton (from=0, to=1, by=0.1,
                                          value=0.5,
                                          container=grid4)

  grid4[2,3] <- e4_hjust_spb <- gspinbutton (from=0, to=1, by=0.1,
                                          value=0,
                                          container=grid4)

  
  
  # FUNCTIONS #################################################################
  
  
  .plotBalance <- function(what, complex=NULL){
    
    # Get values.
    val_titles <- svalue(f1_titles_chk)
    val_title <- svalue(title_edt)
    val_xtitle <- svalue(x_title_edt)
    val_ytitle <- svalue(y_title_edt)
    val_shape <- as.numeric(svalue(e2_shape_spb))
    val_alpha <- as.numeric(svalue(e2_alpha_spb))
    val_jitter <- as.numeric(svalue(e2_jitter_edt))
    val_log <- svalue(f1_logHb_chk)
    val_ymin <- as.numeric(svalue(e3_y_min_edt))
    val_ymax <- as.numeric(svalue(e3_y_max_edt))
    val_xmin <- as.numeric(svalue(e3_x_min_edt))
    val_xmax <- as.numeric(svalue(e3_x_max_edt))
    val_angle <- as.numeric(svalue(e4_angle_spb))
    val_vjust <- as.numeric(svalue(e4_vjust_spb))
    val_hjust <- as.numeric(svalue(e4_hjust_spb))
    val_size <- as.numeric(svalue(e4_size_edt))
    val_scales <- svalue(e3_scales_opt)
    val_kit <- svalue(kit_drp)
    val_drop <- svalue(f1_drop_chk)
    val_theme <- svalue(f1_theme_drp)
    val_wrap <- svalue(f1_wrap_opt, index = TRUE)
    
    if(debug){
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
      print("val_log")
      print(val_log)
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
      print("str(.gData)")
      print(str(.gData))
      print("val_drop")
      print(val_drop)
      print("val_kit")
      print(val_kit)
      print("val_theme")
      print(val_theme)
      print("val_wrap")
      print(val_wrap)
    }
    
    # Declare variables.
    ymax <- NULL  # For complex plots.
    ymin <- NULL  # For complex plots.
    
    if (!is.na(.gData) && !is.null(.gData)){
      
      # Call functions.
      # Add color information.
      if(is.null(.gData$Dye)){
        .gData <- addColor(data=.gData, kit=val_kit, need="Dye")
        message("'Dye' is missing. Dye information added!")

      }
      # Sort by marker in kit
      .gData <- sortMarker(data=.gData,
                          kit=val_kit,
                          add.missing.levels = TRUE)

      # Drop sex markers.
      if(val_drop){

        # Get sex marker.
        sexMarkers <- getKit(val_kit, what="Sex.Marker")
        
        # Check if sexMarkers was found.
        if(length(sexMarkers) > 0){

          # Drop sex markers.
          n0 <- nrow(.gData)
          for(m in seq(along=sexMarkers)){
            .gData <- .gData[.gData$Marker != sexMarkers[m], ]
          }
          n1 <- nrow(.gData)
          message(paste(n1, " rows after removing ", n0-n1, " sex marker rows.", sep=""))
          
          # Refactor and keep order of levels.
          .gData$Marker <- factor(.gData$Marker, 
                                  levels=levels(.gData$Marker)[!levels(.gData$Marker) %in% sexMarkers])
          
        }

      }
      
      # Height must be numeric (not string).
      if("MPH" %in% names(.gData)){
        if(!is.numeric(.gData$MPH)){
          .gData$MPH <- as.numeric(as.character(.gData$MPH))
          message("'MPH' not numeric, converting to numeric.")
        }
      }
      
      # Height must be numeric (not string).
      if("TPH" %in% names(.gData)){
        if(!is.numeric(.gData$TPH)){
          .gData$TPH <- as.numeric(as.character(.gData$TPH))
          message("'TPH' not numeric, converting to numeric.")
          
        }
      }
      
      # Control complex plot.
      if(val_wrap == 1 | val_wrap == 2){
        complex <- FALSE
        message("val_wrap=1/2 overrides and set complex=FALSE")
      }
      
      # Make data frame from dataset marker levels.
      markerDye <- data.frame(Marker=levels(.gData$Marker))
      # Add colors.
      markerDye <- addColor(data=markerDye, kit=val_kit)
      # Get Marker and Dye column.
      markerDye <- markerDye[c("Marker","Dye")]
      # Extract unique elements.
      uniqueMarkerDye <- markerDye[!duplicated(markerDye),]
      # Calculate number of unique columns per dye.
      val_ncol <- unique(table(uniqueMarkerDye$Dye))

      # Check if 'simple' or 'complex' plotting:
      if(is.null(complex)){
        
        # Auto detect if complex plot.
        complex <- length(val_ncol) > 1
        
      }
      
      # Make palette.
      val_palette <- unique(getKit(val_kit, what="Color")$Color)
      val_palette <- addColor(val_palette, have="Color", need="R.Color")
      
      if(debug){
        print("Before plot: str(.gData)")
        print(str(.gData))
        print("Number of columns")
        print(val_ncol)
        print("val_palette")
        print(val_palette)
        print("levels(.gData$MPH), expect NULL")
        print(levels(.gData$MPH))
        print("levels(.gData$Hb), expect NULL")
        print(levels(.gData$Hb))
      }
      
      # Convert to log(balance).
      if(val_log){
        .gData$Hb <- log(.gData$Hb)
        .gData$Lb <- log(.gData$Lb)
      }

      # Create custom titles.
      if(val_titles){
        mainTitle <- val_title
        xTitle <- val_xtitle
        yTitle <- val_ytitle
      }
        
      # Create default titles.
      if(!val_titles){
        
        if(debug){
          print("Using default titles.")
        }
        
        if(what == "Hb"){
          
          mainTitle <- "Heterozygous balance"
          xTitle <- "Mean peak height (RFU)"
          if(val_log){
            yTitle <- "Log(Ratio)"
          } else {
            yTitle <- "Ratio"
          }
          
        } else if (what == "Hb_D") {
          
          mainTitle <- "Heterozygous balance"
          xTitle <- "Repeat difference"
          if(val_log){
            yTitle <- "Log(Ratio)"
          } else {
            yTitle <- "Ratio"
          }
          
        } else if (what == "Hb_H") {
          
          mainTitle <- "Heterozygous balance"
          xTitle <- "Average peak height 'H' (RFU)"
          if(val_log){
            yTitle <- "Log(Ratio)"
          } else {
            yTitle <- "Ratio"
          }
          
        } else if (what == "Hb_M") {
          
          mainTitle <- "Heterozygous balance"
          xTitle <- "Locus"
          if(val_log){
            yTitle <- "Log(Ratio)"
          } else {
            yTitle <- "Ratio"
          }
          
        } else if (what == "Lb") {
          
          mainTitle <- "Locus balance"
          xTitle <- "Locus peak height (RFU)"
          if(val_log){
            yTitle <- "Log(Ratio)"
          } else {
            yTitle <- "Ratio"
          }
          
        } else if (what == "Lb_H") {
          
          mainTitle <- "Locus balance"
          xTitle <- "Average peak height 'H' (RFU)"
          if(val_log){
            yTitle <- "Log(Ratio)"
          } else {
            yTitle <- "Ratio"
          }
          
        } else if (what == "Lb_M") {
          
          mainTitle <- "Locus balance"
          xTitle <- "Locus"
          if(val_log){
            yTitle <- "Log(Ratio)"
          } else {
            yTitle <- "Ratio"
          }
          
        } else {
          
          stop(paste("what =", what, "not implemented for create title!"))
          
        }
      }

      if(debug){
        print("Titles:")
        print(mainTitle)
        print(xTitle)
        print(yTitle)
      }
      
      # Construct plot differently.
      if(!complex){
        # Simple plot, equal number of markers per dye.

        if(debug){
          print("Simple plot.")
        }

        # Create data.table to check for facet error caused by all NA's.
        dt <- data.table::data.table(.gData)
        
        # Select what to plot and create default titles.
        if(what == "Hb"){
          
          gp <- ggplot(.gData, aes_string(x="MPH", y="Hb", colour="Dye"))
          val_box <- FALSE
          
          # Check for facet error caused by all NA's.
          tmp <- dt[, list(Sum=sum(Hb)), by=Marker]
          if(any(is.na(tmp$Sum))){
            message("Empty facets detected! If this leads to plot error try another scale for axes")
          }
          
        } else if (what == "Hb_D") {
          
          gp <- ggplot(.gData, aes_string(x="Delta", y="Hb", colour="Dye"))
          val_box <- FALSE
          
          # Check for facet error caused by all NA's.
          tmp <- dt[, list(Sum=sum(Hb)), by=Marker]
          if(any(is.na(tmp$Sum))){
            message("Empty facets detected! If this leads to plot error try another scale for axes")
          }
          
        } else if (what == "Hb_H") {
          
          gp <- ggplot(.gData, aes_string(x="H", y="Hb", colour="Dye"))
          val_box <- FALSE
          
          # Check for facet error caused by all NA's.
          tmp <- dt[, list(Sum=sum(Hb)), by=Marker]
          if(any(is.na(tmp$Sum))){
            message("Empty facets detected! If this leads to plot error try another scale for axes")
          }
          
        } else if (what == "Hb_M") {
          
          gp <- ggplot(.gData, aes_string(x="Marker", y="Hb", colour="Dye"))
          val_box <- TRUE
          
          # Check for facet error caused by all NA's.
          tmp <- dt[, list(Sum=sum(Hb)), by=Marker]
          if(any(is.na(tmp$Sum))){
            message("Empty facets detected! If this leads to plot error try another scale for axes")
          }
          
        } else if (what == "Lb") {
          
          gp <- ggplot(.gData, aes_string(x="TPH", y="Lb", colour="Dye"))
          val_box <- FALSE
          
          # Check for facet error caused by all NA's.
          tmp <- dt[, list(Sum=sum(Lb)), by=Marker]
          if(any(is.na(tmp$Sum))){
            message("Empty facets detected! If this leads to plot error try another scale for axes")
          }
          
        } else if (what == "Lb_H") {
          
          gp <- ggplot(.gData, aes_string(x="H", y="Lb", colour="Dye"))
          val_box <- FALSE
          
          # Check for facet error caused by all NA's.
          tmp <- dt[, list(Sum=sum(Lb)), by=Marker]
          if(any(is.na(tmp$Sum))){
            message("Empty facets detected! If this leads to plot error try another scale for axes")
          }
          
        } else if (what == "Lb_M") {
          
          gp <- ggplot(.gData, aes_string(x="Marker", y="Lb", colour="Dye"))
          val_box <- TRUE
          
          # Check for facet error caused by all NA's.
          tmp <- dt[, list(Sum=sum(Lb)), by=Marker]
          if(any(is.na(tmp$Sum))){
            message("Empty facets detected! If this leads to plot error try another scale for axes")
          }
          
        } else {
          
          stop(paste("what =", what, "not implemented for create plot!"))
          
        }
        
        if(debug){
          print("Plot created.")
        }
        
        # Apply theme.
        gp <- gp + eval(parse(text=val_theme))
        
        # Create Plot.
        if(val_box){
          
          # Create box and whisker plot.
          gp <- gp + geom_boxplot(alpha=val_alpha)
          
        } else {
          
          # Create scatter plot.
          gp <- gp + geom_point(shape=val_shape, alpha=val_alpha,
                                position=position_jitter(height = 0, width=val_jitter))
          
        }
        
        # Facet plot.
        if(val_wrap == 1){
          
          # Do nothing.
          
        } else if(val_wrap == 2){
          
          # Plot by dye per row.
          gp <- gp + facet_wrap(as.formula(paste("~", "Dye")), ncol=1,
                                drop=FALSE, scales=val_scales)

        } else if(val_wrap == 3){

          # Plot per marker one dye per row.
          gp <- gp + facet_grid("Dye ~ Marker")
          # NB! 'facet_wrap' does not seem to support strings.
          #     Use 'as.formula(paste("string1", "string2"))' as a workaround.
          gp <- gp + facet_wrap(as.formula(paste("~", "Marker")),
                                ncol=val_ncol, drop=FALSE, scales=val_scales)
          
        } else {
          
          stop("val_wrap =", val_wrap, "not implemented.")
          
        }
        
        # Add colours.
        gp <- gp + scale_colour_manual(guide=FALSE, values=val_palette)
        
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
        
        # Add titles etc.
        gp <- gp + guides(fill = guide_legend(reverse=TRUE))
        gp <- gp + theme(axis.text.x=element_text(angle=val_angle,
                                                  hjust=val_hjust,
                                                  vjust=val_vjust,
                                                  size=val_size))
        gp <- gp + labs(title=mainTitle)
        gp <- gp + xlab(xTitle)
        gp <- gp + ylab(yTitle)

        # plot.
        print(gp)
        
        # Change save button.
        svalue(f5_save_btn) <- "Save as object"
        enabled(f5_save_btn) <- TRUE
        
      } else if (complex){
        # Complex plot, unequal number of markers per dye.
        
        if(debug){
          print("Complex plot.")
        }
        
        if(val_scales %in% c("fixed","free_x")){
          # Keep Y max fixed.
          if("Hb" %in% what){
            ymax <- max(.gData$Hb, na.rm=TRUE) * 1.05
            ymin <- min(.gData$Hb, na.rm=TRUE) * 0.95
          }
          if("Lb" %in% what){
            ymax <- max(.gData$Lb, na.rm=TRUE) * 1.05
            ymin <- min(.gData$Lb, na.rm=TRUE) * 0.95
          }
        }
        
        # Get kit colors and convert to dyes.
        dyes <- unique(getKit(val_kit, what="Color")$Color)
        dyes <- addColor(dyes, have="Color", need="Dye")
        # Number of dyes.
        noDyes <- length(dyes)
        # Number of rows in table object (one for each dye + title + x title).
        noRows <- length(dyes) + 2
        
        # Create table object.
        # Note: width(1.5 for y-title, and the rest for plots)
        #       height(1.5 for plot title, equal for each plot, and 1.5 for x-title)
        g <- gtable::gtable(widths = grid::unit(c(1.5,1),c("lines","null")),
                    heights = grid::unit(c(1.5,rep(1,noDyes),1.5), c("line",rep("null",noDyes),"line")))
        
        # Add titles.        
        g <- gtable::gtable_add_grob(g, grid::textGrob(mainTitle), t=1,b=1,l=2,r=2)
        g <- gtable::gtable_add_grob(g, grid::textGrob(xTitle), t=noRows ,b=noRows ,l=2,r=2)
        g <- gtable::gtable_add_grob(g, grid::textGrob(yTitle, rot=90), t=1,b=noRows ,l=1,r=1)
        
        # Get all markers to be plotted and add dye for subsetting.
        gLevel <- data.frame(Marker=levels(.gData$Marker))
        gLevel <- addColor(gLevel, kit=val_kit)
        
        # Loop over all dyes.
        for(d in seq(along=dyes)){
          
          # Get data for current dye.
          gDataSub <- .gData[.gData$Dye == dyes[d],]

          # Get current markers/levels.
          gDyeLevel <- as.character(gLevel$Marker[gLevel$Dye==dyes[d]])
          
          # Can't handle zero rows.
          if(nrow(gDataSub)==0){
            tmp <- data.frame(Marker=gDyeLevel, Allele=NA, MPH=0, Hb=0,
                              Lb=0, H=0, Delta=0, Dye=dyes[d])
            gDataSub <- plyr::rbind.fill(gDataSub, tmp)

          }
          
          # Refactor to levels of current dye (and maintain order).
          gDataSub$Marker <- factor(gDataSub$Marker, levels=gDyeLevel)
          gDataSub$Dye <- factor(dyes[d])
          
          # Create data.table to check for facet error caused by all NA's.
          dt <- data.table::data.table(gDataSub)
          
          # Create a plot for the current subset.
          # Select what to plot.
          if(what == "Hb"){
            
            gp <- ggplot(gDataSub, aes_string(x = "MPH", y = "Hb", colour="Dye"))
            val_box <- FALSE
            
            # Check for facet error caused by all NA's.
            tmp <- dt[, list(Sum=sum(Hb)), by=Marker]
            if(any(is.na(tmp$Sum))){
              message("Empty facets detected! If this leads to plot error try another scale for axes")
            }
            
          } else if (what == "Hb_D") {
            
            gp <- ggplot(gDataSub, aes_string(x="Delta", y="Hb", colour="Dye"))
            val_box <- FALSE
            
            # Check for facet error caused by all NA's.
            tmp <- dt[, list(Sum=sum(Hb)), by=Marker]
            if(any(is.na(tmp$Sum))){
              message("Empty facets detected! If this leads to plot error try another scale for axes")
            }
            
          } else if (what == "Hb_H") {
            
            gp <- ggplot(gDataSub, aes_string(x = "H", y = "Hb", colour="Dye"))
            val_box <- FALSE
            
            # Check for facet error caused by all NA's.
            tmp <- dt[, list(Sum=sum(Hb)), by=Marker]
            if(any(is.na(tmp$Sum))){
              message("Empty facets detected! If this leads to plot error try another scale for axes")
            }
            
          } else if (what == "Hb_M") {
            
            gp <- ggplot(gDataSub, aes_string(x = "Marker", y = "Hb", colour="Dye"))
            val_box <- TRUE
            
            # Check for facet error caused by all NA's.
            tmp <- dt[, list(Sum=sum(Hb)), by=Marker]
            if(any(is.na(tmp$Sum))){
              message("Empty facets detected! If this leads to plot error try another scale for axes")
            }
            
          } else if (what == "Lb") {
            
            gp <- ggplot(gDataSub, aes_string(x = "TPH", y = "Lb", colour="Dye"))
            val_box <- FALSE
            
            # Check for facet error caused by all NA's.
            tmp <- dt[, list(Sum=sum(Lb)), by=Marker]
            if(any(is.na(tmp$Sum))){
              message("Empty facets detected! If this leads to plot error try another scale for axes")
            }
            
          } else if (what == "Lb_H") {
            
            gp <- ggplot(gDataSub, aes_string(x = "H", y = "Lb", colour="Dye"))
            val_box <- FALSE
            
            # Check for facet error caused by all NA's.
            tmp <- dt[, list(Sum=sum(Lb)), by=Marker]
            if(any(is.na(tmp$Sum))){
              message("Empty facets detected! If this leads to plot error try another scale for axes")
            }
            
          } else if (what == "Lb_M") {
            
            gp <- ggplot(gDataSub, aes_string(x = "Marker", y = "Lb", colour="Dye"))
            val_box <- TRUE
            
            # Check for facet error caused by all NA's.
            tmp <- dt[, list(Sum=sum(Lb)), by=Marker]
            if(any(is.na(tmp$Sum))){
              message("Empty facets detected! If this leads to plot error try another scale for axes")
            }
            
          } else {
            
            stop(paste("what =", what, "not implemented for create complex plot!"))
            
          }
          
          # Apply theme.
          gp <- gp + eval(parse(text=val_theme))
          
          # Create Plot.
          if(val_box){
            
            # Create box and whisker plot.
            gp <- gp + geom_boxplot(alpha=val_alpha)
            
          } else {
            
            # Create scatter plot.
            gp <- gp + geom_point(shape=val_shape, alpha = val_alpha,
                                  position = position_jitter(height = 0,
                                                             width = val_jitter))
            
          }
          
          # Wrap by marker one dye per row.
          gp <- gp + facet_grid("Dye ~ Marker", scales=val_scales, drop = FALSE) # Keep dye labels.

          # Add colours.
          gp <- gp + scale_colour_manual(guide=FALSE, values=val_palette[d], drop=FALSE)
          
          # Set margin around each plot. Note: top, right, bottom, left.
          gp <- gp + theme(plot.margin = grid::unit(c(0.25, 1.25, 0, 0), "lines"))
          
          # Restrict y axis.
          if(!is.na(val_ymin) && !is.na(val_ymax)){
            val_y <- c(val_ymin, val_ymax)
          } else {
            if(val_scales %in% c("fixed","free_x")){
              # Keep Y fixed.
              val_y <- c(ymin, ymax)
            }
            val_y <- NULL
          }
          # Restrict x axis.
          if(!is.na(val_xmin) && !is.na(val_xmax)){
            val_x <- c(val_xmin, val_xmax)
          } else {
            if(val_scales %in% c("fixed","free_x")){
              # Keep Y fixed.
              val_y <- c(ymin, ymax)
            }
            val_x <- NULL
          }
          # Zoom in without dropping observations.
          gp <- gp + coord_cartesian(xlim=val_x, ylim=val_y)
          
          if(debug){
            print(paste("Plot zoomed to xlim:", val_x, "ylim:", val_y))
          }
          
          # Remove titles, axis labels and legend.
          gp <- gp + labs(title = element_blank())
          gp <- gp + theme(axis.title.x = element_blank())
          gp <- gp + theme(axis.text.x=element_text(angle=val_angle,
                                                    hjust=val_hjust,
                                                    vjust=val_vjust,
                                                    size=val_size))
          gp <- gp + theme(axis.title.y = element_blank())
          gp <- gp + theme(legend.position="none")
          
          
          # Add plot panel to table object.  
          g <- gtable::gtable_add_grob(g,ggplotGrob(gp), t=(d+1),b=(d+1),l=2,r=2)
          
        }
        
        # Plot.
        grid::grid.newpage()
        grid::grid.draw(g)
   
        # This is step 1 in workaround to save 'complex plots':
        # Step 1: http://stackoverflow.com/a/20433318/2173340
        # Step 2: http://stackoverflow.com/a/18407452/2173340
        gp <- gridExtra::arrangeGrob(g)
        
        # Change save button.
        svalue(f5_save_btn) <- "Save as object"
        enabled(f5_save_btn) <- FALSE
        
      } # End if(complex)
      
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
      if(exists(".strvalidator_plotBalance_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_plotBalance_gui_savegui", envir=env)
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
      if(exists(".strvalidator_plotBalance_gui_title", envir=env, inherits = FALSE)){
        svalue(title_edt) <- get(".strvalidator_plotBalance_gui_title", envir=env)
      }
      if(exists(".strvalidator_plotBalance_gui_title_chk", envir=env, inherits = FALSE)){
        svalue(f1_titles_chk) <- get(".strvalidator_plotBalance_gui_title_chk", envir=env)
      }
      if(exists(".strvalidator_plotBalance_gui_x_title", envir=env, inherits = FALSE)){
        svalue(x_title_edt) <- get(".strvalidator_plotBalance_gui_x_title", envir=env)
      }
      if(exists(".strvalidator_plotBalance_gui_y_title", envir=env, inherits = FALSE)){
        svalue(y_title_edt) <- get(".strvalidator_plotBalance_gui_y_title", envir=env)
      }
      if(exists(".strvalidator_plotBalance_gui_sex", envir=env, inherits = FALSE)){
        svalue(f1_drop_chk) <- get(".strvalidator_plotBalance_gui_sex", envir=env)
      }
      if(exists(".strvalidator_plotBalance_gui_log", envir=env, inherits = FALSE)){
        svalue(f1_logHb_chk) <- get(".strvalidator_plotBalance_gui_log", envir=env)
      }
      if(exists(".strvalidator_plotBalance_gui_points_shape", envir=env, inherits = FALSE)){
        svalue(e2_shape_spb) <- get(".strvalidator_plotBalance_gui_points_shape", envir=env)
      }
      if(exists(".strvalidator_plotBalance_gui_points_alpha", envir=env, inherits = FALSE)){
        svalue(e2_alpha_spb) <- get(".strvalidator_plotBalance_gui_points_alpha", envir=env)
      }
      if(exists(".strvalidator_plotBalance_gui_points_jitter", envir=env, inherits = FALSE)){
        svalue(e2_jitter_edt) <- get(".strvalidator_plotBalance_gui_points_jitter", envir=env)
      }
      if(exists(".strvalidator_plotBalance_gui_axes_y_min", envir=env, inherits = FALSE)){
        svalue(e3_y_min_edt) <- get(".strvalidator_plotBalance_gui_axes_y_min", envir=env)
      }
      if(exists(".strvalidator_plotBalance_gui_axes_y_max", envir=env, inherits = FALSE)){
        svalue(e3_y_max_edt) <- get(".strvalidator_plotBalance_gui_axes_y_max", envir=env)
      }
      if(exists(".strvalidator_plotBalance_gui_axes_x_min", envir=env, inherits = FALSE)){
        svalue(e3_x_min_edt) <- get(".strvalidator_plotBalance_gui_axes_x_min", envir=env)
      }
      if(exists(".strvalidator_plotBalance_gui_axes_x_max", envir=env, inherits = FALSE)){
        svalue(e3_x_max_edt) <- get(".strvalidator_plotBalance_gui_axes_x_max", envir=env)
      }
      if(exists(".strvalidator_plotBalance_gui_axes_scales", envir=env, inherits = FALSE)){
        svalue(e3_scales_opt) <- get(".strvalidator_plotBalance_gui_axes_scales", envir=env)
      }
      if(exists(".strvalidator_plotBalance_gui_xlabel_size", envir=env, inherits = FALSE)){
        svalue(e4_size_edt) <- get(".strvalidator_plotBalance_gui_xlabel_size", envir=env)
      }
      if(exists(".strvalidator_plotBalance_gui_xlabel_angle", envir=env, inherits = FALSE)){
        svalue(e4_angle_spb) <- get(".strvalidator_plotBalance_gui_xlabel_angle", envir=env)
      }
      if(exists(".strvalidator_plotBalance_gui_xlabel_justh", envir=env, inherits = FALSE)){
        svalue(e4_hjust_spb) <- get(".strvalidator_plotBalance_gui_xlabel_justh", envir=env)
      }
      if(exists(".strvalidator_plotBalance_gui_xlabel_justv", envir=env, inherits = FALSE)){
        svalue(e4_vjust_spb) <- get(".strvalidator_plotBalance_gui_xlabel_justv", envir=env)
      }
      if(exists(".strvalidator_plotBalance_gui_theme", envir=env, inherits = FALSE)){
        svalue(f1_theme_drp) <- get(".strvalidator_plotBalance_gui_theme", envir=env)
      }
      if(exists(".strvalidator_plotBalance_gui_wrap", envir=env, inherits = FALSE)){
        svalue(f1_wrap_opt) <- get(".strvalidator_plotBalance_gui_wrap", envir=env)
      }

      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_plotBalance_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_plotBalance_gui_sex", value=svalue(f1_drop_chk), envir=env)
      assign(x=".strvalidator_plotBalance_gui_log", value=svalue(f1_logHb_chk), envir=env)
      assign(x=".strvalidator_plotBalance_gui_title", value=svalue(title_edt), envir=env)
      assign(x=".strvalidator_plotBalance_gui_title_chk", value=svalue(f1_titles_chk), envir=env)
      assign(x=".strvalidator_plotBalance_gui_x_title", value=svalue(x_title_edt), envir=env)
      assign(x=".strvalidator_plotBalance_gui_y_title", value=svalue(y_title_edt), envir=env)
      assign(x=".strvalidator_plotBalance_gui_points_shape", value=svalue(e2_shape_spb), envir=env)
      assign(x=".strvalidator_plotBalance_gui_points_alpha", value=svalue(e2_alpha_spb), envir=env)
      assign(x=".strvalidator_plotBalance_gui_points_jitter", value=svalue(e2_jitter_edt), envir=env)
      assign(x=".strvalidator_plotBalance_gui_axes_y_min", value=svalue(e3_y_min_edt), envir=env)
      assign(x=".strvalidator_plotBalance_gui_axes_y_max", value=svalue(e3_y_max_edt), envir=env)
      assign(x=".strvalidator_plotBalance_gui_axes_x_min", value=svalue(e3_x_min_edt), envir=env)
      assign(x=".strvalidator_plotBalance_gui_axes_x_max", value=svalue(e3_x_max_edt), envir=env)
      assign(x=".strvalidator_plotBalance_gui_axes_scales", value=svalue(e3_scales_opt), envir=env)
      assign(x=".strvalidator_plotBalance_gui_xlabel_size", value=svalue(e4_size_edt), envir=env)
      assign(x=".strvalidator_plotBalance_gui_xlabel_angle", value=svalue(e4_angle_spb), envir=env)
      assign(x=".strvalidator_plotBalance_gui_xlabel_justh", value=svalue(e4_hjust_spb), envir=env)
      assign(x=".strvalidator_plotBalance_gui_xlabel_justv", value=svalue(e4_vjust_spb), envir=env)
      assign(x=".strvalidator_plotBalance_gui_theme", value=svalue(f1_theme_drp), envir=env)
      assign(x=".strvalidator_plotBalance_gui_wrap", value=svalue(f1_wrap_opt), envir=env)

    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_plotBalance_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotBalance_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_plotBalance_gui_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotBalance_gui_title", envir = env)
      }
      if(exists(".strvalidator_plotBalance_gui_title_chk", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotBalance_gui_title_chk", envir = env)
      }
      if(exists(".strvalidator_plotBalance_gui_x_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotBalance_gui_x_title", envir = env)
      }
      if(exists(".strvalidator_plotBalance_gui_y_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotBalance_gui_y_title", envir = env)
      }
      if(exists(".strvalidator_plotBalance_gui_sex", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotBalance_gui_sex", envir = env)
      }
      if(exists(".strvalidator_plotBalance_gui_log", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotBalance_gui_log", envir = env)
      }
      if(exists(".strvalidator_plotBalance_gui_points_shape", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotBalance_gui_points_shape", envir = env)
      }
      if(exists(".strvalidator_plotBalance_gui_points_alpha", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotBalance_gui_points_alpha", envir = env)
      }
      if(exists(".strvalidator_plotBalance_gui_points_jitter", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotBalance_gui_points_jitter", envir = env)
      }
      if(exists(".strvalidator_plotBalance_gui_axes_y_min", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotBalance_gui_axes_y_min", envir = env)
      }
      if(exists(".strvalidator_plotBalance_gui_axes_y_max", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotBalance_gui_axes_y_max", envir = env)
      }
      if(exists(".strvalidator_plotBalance_gui_axes_x_min", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotBalance_gui_axes_x_min", envir = env)
      }
      if(exists(".strvalidator_plotBalance_gui_axes_x_max", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotBalance_gui_axes_x_max", envir = env)
      }
      if(exists(".strvalidator_plotBalance_gui_axes_scales", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotBalance_gui_axes_scales", envir = env)
      }
      if(exists(".strvalidator_plotBalance_gui_xlabel_size", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotBalance_gui_xlabel_size", envir = env)
      }
      if(exists(".strvalidator_plotBalance_gui_xlabel_angle", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotBalance_gui_xlabel_angle", envir = env)
      }
      if(exists(".strvalidator_plotBalance_gui_xlabel_justh", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotBalance_gui_xlabel_justh", envir = env)
      }
      if(exists(".strvalidator_plotBalance_gui_xlabel_justv", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotBalance_gui_xlabel_justv", envir = env)
      }
      if(exists(".strvalidator_plotBalance_gui_theme", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotBalance_gui_theme", envir = env)
      }
      if(exists(".strvalidator_plotBalance_gui_wrap", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotBalance_gui_wrap", envir = env)
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
