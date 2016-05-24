################################################################################
# TODO LIST
# TODO: Add nice breaks (se tutorial vWA by height), can be fixed using 'fixed_x' scale.

################################################################################
# CHANGE LOG (last 20 changes)
# 06.01.2016: Fixed theme methods not found and added more themes.
# 11.11.2015: Added importFrom grid and gridExtra arrangeGrob, and ggplot2.
# 11.11.2015: Added more themes.
# 29.08.2015: Added importFrom.
# 05.01.2015: 'Save as object' now disabled when complex plot.
# 14.12.2014: Option to drop sex markers.
# 14.12.2014: Updated to handle gender -> sex.marker option in getKit.
# 11.10.2014: Added 'focus', added 'parent' parameter.
# 28.06.2014: Added help button and moved save gui checkbox.
# 06.05.2014: Implemented 'checkDataset'.
# 05.05.2014: Fixed same color scale for all sub plots in complex plots.
# 14.04.2014: Fixed different y max for complex plot, when supposed to be fixed.
# 14.04.2014: Fixed now handle no observation in an entire dye channel.
# 14.04.2014: Fixed position_jitter height now fixed to zero (prev. default).
# 23.02.2014: Fixed different y max for complex plot, when supposed to be fixed.
# 23.02.2014: Fixed shape for 'complex' plots.
# 13.02.2014: Implemented theme.
# 20.01.2014: Implemented ggsave with workaround for complex plots.
# 30.11.2013: Added info on number of samples.
# 30.11.2013: Fixed 'complex' plot.
# 30.11.2013: Fixed 'facet_wrap' with strings.

#' @title Plot Stutter
#'
#' @description
#' GUI simplifying the creation of plots from stutter data.
#'
#' @details Select data to plot in the drop-down menu. Check that the correct
#' kit has been detected. Plot stutter data by parent allele or by peak height. 
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
#' @importFrom gridExtra arrangeGrob
# @importFrom gtable gtable_add_grob gtable gtable_filter
#' @importFrom grid unit textGrob grid.newpage grid.draw
#' @importFrom plyr rbind.fill
#' @importFrom scales pretty_breaks
#' @importFrom utils help str
#' @importFrom stats as.formula
#' @importFrom grDevices hcl
#' @importFrom ggplot2 ggplot aes_string scale_x_continuous geom_point position_jitter
#'  facet_grid facet_wrap coord_cartesian guides guide_legend theme element_text
#'  labs xlab ylab ggplotGrob scale_colour_manual element_blank theme_gray
#'  theme_bw theme_linedraw theme_light theme_dark theme_minimal theme_classic
#'  theme_void
#' 
#' @seealso \url{http://docs.ggplot2.org/current/} for details on plot settings.

plotStutter_gui <- function(env=parent.frame(), savegui=NULL, debug=FALSE, parent=NULL){
  
  # Global variables.
  .gData <- NULL
  .gDataName <- NULL
  .gPlot <- NULL
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # Main window.
  w <- gwindow(title="Plot stutter proportions", visible=FALSE)
  
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
    print(help("plotStutter_gui", help_type="html"))
    
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
    requiredCol <- c("Marker", "Allele", "HeightA", "Stutter", "Type")
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
      
      # Detect kit.
      kitIndex <- detectKit(.gData, index=TRUE)
      # Select in dropdown.
      svalue(kit_drp, index=TRUE) <- kitIndex
      
      # Enable buttons.
      enabled(plot_allele_btn) <- TRUE
      enabled(plot_height_btn) <- TRUE
      
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
  
  f1_drop_chk <- gcheckbox(text="Drop sex markers",
                           checked=TRUE,
                           container=f1)

  # FRAME 7 ###################################################################
  
  f7 <- gframe(text = "Plot stutter data",
               horizontal=FALSE,
               container = gv) 
  
  grid7 <- glayout(container = f7)
  
  grid7[1,1] <- plot_allele_btn <- gbutton(text="Ratio vs. Allele",
                                           border=TRUE,
                                           container=grid7) 
  
  grid7[1,2] <- plot_height_btn <- gbutton(text="Ratio vs. Height",
                                           border=TRUE,
                                           container=grid7) 
  
  addHandlerChanged(plot_allele_btn, handler = function(h, ...) {
    
    enabled(plot_allele_btn) <- FALSE
    .plotStutter(what="allele")
    enabled(plot_allele_btn) <- TRUE

  } )
  
  addHandlerChanged(plot_height_btn, handler = function(h, ...) {
    
    enabled(plot_height_btn) <- FALSE
    .plotStutter(what="height")
    enabled(plot_height_btn) <- TRUE
    
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
                                         by=0.01, value=0.60,
                                         container=grid2)

  grid2[1,5] <- glabel(text="Jitter (width):", container=grid2)
  grid2[1,6] <- jitter_txt <- gedit(text="0.1", width=4, container=grid2)

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
  grid4[1,2] <- size_txt <- gedit(text="8", width=4, container=grid4)

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
  
  
  .plotStutter <- function(what){
    
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
    val_kit <- svalue(kit_drp)
    val_theme <- svalue(f1_theme_drp)
    val_drop <- svalue(f1_drop_chk)
    
    # Declare variables.
    ymax <- NULL  # For complex plots.
    ymin <- NULL  # For complex plots.
    
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
      print("levels(.gData$Allele)")
      print(levels(.gData$Allele))
      print("levels(.gData$Stutter)")
      print(levels(.gData$Stutter))
      print("levels(.gData$Marker)")
      print(levels(.gData$Marker))
      print("val_theme")
      print(val_theme)
      print("val_drop")
      print(val_drop)
    }
    
    if (!is.na(.gData) && !is.null(.gData)){
      
      # Call functions.
      # Add color information.
      if(!"Dye" %in% names(.gData)){
        .gData <- addColor(data=.gData, kit=val_kit, need="Dye", debug=debug)
        message("'Dye' added to dataset!")
      }

      # Sort by marker in kit
      .gData <- sortMarker(data=.gData,
                          kit=val_kit,
                          add.missing.levels = TRUE)

      # Drop sex markers.
      if(val_drop){
        
        # Get sex markers.
        sexMarkers <- getKit(val_kit, what="Sex.Marker")
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
      
      # Create factors and round. IMPORTANT!
      .gData$Type <- factor(round(.gData$Type,2))
      
      # Sort stutter/allele factors. IMPORTANT!
      .gData$Stutter <- factor(.gData$Stutter, levels=sort(unique(as.numeric(as.character(.gData$Stutter)))))
      .gData$Allele <- factor(.gData$Allele, levels=sort(unique(as.numeric(as.character(.gData$Allele)))))
      
      # Height must be numeric (not string).
      .gData$HeightA <- as.numeric(as.character(.gData$HeightA))
      
      if(debug){
        print("BEFORE PLOTTING:")
        print("str(.gData)")
        print(str(.gData))
        print("levels(.gData$Allele)")
        print(levels(.gData$Allele))
        print("levels(.gData$Stutter)")
        print(levels(.gData$Stutter))
        print("levels(.gData$Marker)")
        print(levels(.gData$Marker))
      }

      # Check if 'simple' or 'complex' plotting:
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
      
      
      # Plotting alleles for observed stutters per marker.
      if(what == "allele"){
        
        if(val_titles){
          mainTitle <- val_title
          xTitle <- val_xtitle
          yTitle <- val_ytitle
        } else {
          mainTitle <- "Stutter ratios"
          xTitle <- "True allele"
          yTitle <- "Ratio"
        }
        
        gp <- ggplot(.gData, aes_string(x="Allele", y="Ratio", colour="Type"))
        
        
      } else if (what == "height") {
        
        if(val_titles){
          mainTitle <- val_title
          xTitle <- val_xtitle
          yTitle <- val_ytitle
        } else {
          mainTitle <- "Stutter ratios"
          xTitle <- "True allele (RFU)"
          yTitle <- "Ratio"
        }
        
        gp <- ggplot(.gData, aes_string(x="HeightA", y="Ratio", colour="Type"))
        gp <- gp + scale_x_continuous(breaks = scales::pretty_breaks())
        
      }
      
      # Apply theme.
      gp <- gp + eval(parse(text=val_theme))
      
      # Plot settings.
      gp <- gp + geom_point(shape=val_shape, alpha=val_alpha, 
                            position=position_jitter(height = 0, width=val_jitter)) 

      # Facet and keep all levels.
      gp <- gp + facet_grid("Dye ~ Marker", drop = FALSE)
      # NB! 'facet_wrap' does not seem to support strings.
      #     Use 'as.formula(paste("string1", "string2"))' as a workaround.
      gp <- gp + facet_wrap(as.formula(paste("Dye ~ Marker")), ncol=val_ncol, # Keep dye labels
                            drop=FALSE, scales=val_scales)
      #gp <- gp + facet_wrap(as.formula(paste("~ Marker")), ncol=val_ncol, # No dye labels
      #                      drop=FALSE, scales=val_scales)
      
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
      
      # Titles and legends.
      gp <- gp + guides(fill = guide_legend(reverse=TRUE))
      gp <- gp + theme(axis.text.x=element_text(angle=val_angle,
                                                hjust=val_hjust,
                                                vjust=val_vjust,
                                                size=val_size))
      gp <- gp + labs(title=mainTitle)
      gp <- gp + xlab(xTitle)
      gp <- gp + ylab(yTitle)

      # Check plot type.
      if(length(val_ncol) == 1){
        # Simple plot, equal number of markers per dye.

        if(debug){
          print(paste("Simple plot, val_ncol:",
                      paste(val_ncol, collapse=", ")))
        }
        # NB! 'facet_wrap' does not seem to support strings.
        #     Use 'as.formula(paste("string1", "string2"))' as a workaround.
        gp <- gp + facet_wrap(as.formula(paste("~ Marker")), ncol=val_ncol,
                              drop=FALSE, scales=val_scales)
        
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
      
      } else if (length(val_ncol) > 1){
        # Complex plot, unequal number of markers per dye.

        if(debug){
          print(paste("Complex plot, val_ncol:",
                      paste(val_ncol, collapse=", ")))
        }
        
        # Extract the legend from the 'simple' plot.
        guide <- gtable::gtable_filter(ggplotGrob(gp), pattern="guide")
        
        # Get y max to be able to use same scale across plots.
        ymax <- max(.gData$Ratio, na.rm=TRUE) * 1.05
        ymin <- min(.gData$Ratio, na.rm=TRUE) * 0.95
        
        if(debug){
          print("ymax:")
          print(ymax)
          print("ymin:")
          print(ymin)
        }
        
        # Get kit colors and convert to dyes.
        dyes <- unique(getKit(val_kit, what="Color")$Color)
        dyes <- addColor(dyes, have="Color", need="Dye")
        # Number of dyes.
        noDyes <- length(dyes)
        # Number of rows in table object (one per dye + title + x title).
        noRows <- length(dyes) + 2

        # Create table object.
        # Note: width(1.5 for y-title, and the rest for plots + guides)
        #       height(1.5 for plot title, equal for each plot, and 1.5 for x-title)
        g <- gtable::gtable(widths=grid::unit.c(grid::unit(1.5, "lines"),
                                          grid::unit(1, "null"),
                                          sum(guide$widths)),
                            heights = grid::unit(c(1.5,rep(1,noDyes),1.5),
                                           c("line", rep("null", noDyes), "line")))

        # Add titles.
        g <- gtable::gtable_add_grob(g, grid::textGrob(mainTitle), t=1,b=1,l=2,r=2)
        g <- gtable::gtable_add_grob(g, grid::textGrob(xTitle), t=noRows ,b=noRows ,l=2,r=2)
        g <- gtable::gtable_add_grob(g, grid::textGrob(yTitle, rot=90), t=1,b=noRows ,l=1,r=1)

        # Add the legend to the table object.
        g <- gtable::gtable_add_grob(g,guide , t=1,b=noRows,l=3,r=3)

        # Get all markers to be plotted and add dye for subsetting.
        gLevel <- data.frame(Marker=levels(.gData$Marker))
        gLevel <- addColor(gLevel, kit=val_kit)
        
        # Make palette.
        gTypeLevel <- levels(.gData$Type)
        val_palette <- .gg_color_hue(length(gTypeLevel))
        names(val_palette)  <- gTypeLevel
        
        # Loop over all dyes.
        for(d in seq(along=dyes)){
          
          # Get data for current dye.
          gDataSub <- .gData[.gData$Dye == dyes[d],]
          
          # Get current markers/levels.
          gDyeLevel <- as.character(gLevel$Marker[gLevel$Dye==dyes[d]])
          
          # Can't handle zero rows.
          if(nrow(gDataSub)==0){
            tmp <- data.frame(Marker=gDyeLevel, Allele=NA, HeightA=0, Ratio=0, Type=-1)
            gDataSub <- plyr::rbind.fill(gDataSub, tmp)
            
          }

          # Refactor to levels of current dye (and maintain order).
          gDataSub$Marker <- factor(gDataSub$Marker, levels=gDyeLevel)
          gDataSub$Dye <- factor(dyes[d])
          
          # Create a plot for the current subset.
          # Select what to plot.
          if(what == "allele"){
            # Plotting alleles for observed stutters per marker.

            # Create a plot for the current subset.
            gp <- ggplot(gDataSub, aes_string(x = "Allele", y = "Ratio", color="Type"))
            
          } else if (what == "height") {
            # Plotting true allele height for observed stutters per marker.
            
            # Create a plot for the current subset.
            gp <- ggplot(gDataSub, aes_string(x = "HeightA", y = "Ratio", color="Type"))
            gp <- gp + scale_x_continuous(breaks = scales::pretty_breaks())
            
          }
          
          # Apply theme.
          gp <- gp + eval(parse(text=val_theme))
          
          # Plot settings.
          gp <- gp + geom_point(shape=val_shape, alpha = val_alpha,
                                position = position_jitter(height = 0, width = val_jitter))

          # Add custom color palette (use same for all sub plots).
          gp <- gp + scale_colour_manual(values = val_palette)
          
          # Facet plot and keep all levels (in current dye). 
          gp <- gp + facet_grid("Dye ~ Marker", scales=val_scales, drop = FALSE) # Keep dye labels.
          #gp <- gp + facet_grid("~ Marker", scales=val_scales, drop = FALSE) # No dye labels.
          
          # Set margin around each plot. Note: top, right, bottom, left.
          gp <- gp + theme(plot.margin = grid::unit(c(0.25, 0, 0, 0), "lines"))
          
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
          g <- gtable::gtable_add_grob(g,ggplotGrob(gp),
                                       t=(d+1),b=(d+1),l=2,r=2)

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
        
      } else {
        # Not supported!
        stop(paste("Unsupported number of columns:", val_ncol))
      }
      
      
      # Store in global variable.
      .gPlot <<- gp
      
    } else {
      
      gmessage(message="Data frame is NULL or NA!",
               title="Error",
               icon = "error")      
      
    } 
    
  }

  # INTERNAL FUNCTIONS ########################################################

  # Return a number of ggplot default colors.
  .gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
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
      if(exists(".strvalidator_plotStutter_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_plotStutter_gui_savegui", envir=env)
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
      if(exists(".strvalidator_plotStutter_gui_title", envir=env, inherits = FALSE)){
        svalue(title_edt) <- get(".strvalidator_plotStutter_gui_title", envir=env)
      }
      if(exists(".strvalidator_plotStutter_gui_title_chk", envir=env, inherits = FALSE)){
        svalue(f1_titles_chk) <- get(".strvalidator_plotStutter_gui_title_chk", envir=env)
      }
      if(exists(".strvalidator_plotStutter_gui_x_title", envir=env, inherits = FALSE)){
        svalue(x_title_edt) <- get(".strvalidator_plotStutter_gui_x_title", envir=env)
      }
      if(exists(".strvalidator_plotStutter_gui_y_title", envir=env, inherits = FALSE)){
        svalue(y_title_edt) <- get(".strvalidator_plotStutter_gui_y_title", envir=env)
      }
      if(exists(".strvalidator_plotStutter_gui_points_shape", envir=env, inherits = FALSE)){
        svalue(shape_spb) <- get(".strvalidator_plotStutter_gui_points_shape", envir=env)
      }
      if(exists(".strvalidator_plotStutter_gui_points_alpha", envir=env, inherits = FALSE)){
        svalue(alpha_spb) <- get(".strvalidator_plotStutter_gui_points_alpha", envir=env)
      }
      if(exists(".strvalidator_plotStutter_gui_points_jitter", envir=env, inherits = FALSE)){
        svalue(jitter_txt) <- get(".strvalidator_plotStutter_gui_points_jitter", envir=env)
      }
      if(exists(".strvalidator_plotStutter_gui_axes_y_min", envir=env, inherits = FALSE)){
        svalue(y_min_txt) <- get(".strvalidator_plotStutter_gui_axes_y_min", envir=env)
      }
      if(exists(".strvalidator_plotStutter_gui_axes_y_max", envir=env, inherits = FALSE)){
        svalue(y_max_txt) <- get(".strvalidator_plotStutter_gui_axes_y_max", envir=env)
      }
      if(exists(".strvalidator_plotStutter_gui_axes_x_min", envir=env, inherits = FALSE)){
        svalue(x_min_txt) <- get(".strvalidator_plotStutter_gui_axes_x_min", envir=env)
      }
      if(exists(".strvalidator_plotStutter_gui_axes_x_max", envir=env, inherits = FALSE)){
        svalue(x_max_txt) <- get(".strvalidator_plotStutter_gui_axes_x_max", envir=env)
      }
      if(exists(".strvalidator_plotStutter_gui_axes_scales", envir=env, inherits = FALSE)){
        svalue(scales_opt) <- get(".strvalidator_plotStutter_gui_axes_scales", envir=env)
      }
      if(exists(".strvalidator_plotStutter_gui_xlabel_size", envir=env, inherits = FALSE)){
        svalue(size_txt) <- get(".strvalidator_plotStutter_gui_xlabel_size", envir=env)
      }
      if(exists(".strvalidator_plotStutter_gui_xlabel_angle", envir=env, inherits = FALSE)){
        svalue(angle_spb) <- get(".strvalidator_plotStutter_gui_xlabel_angle", envir=env)
      }
      if(exists(".strvalidator_plotStutter_gui_xlabel_justh", envir=env, inherits = FALSE)){
        svalue(hjust_spb) <- get(".strvalidator_plotStutter_gui_xlabel_justh", envir=env)
      }
      if(exists(".strvalidator_plotStutter_gui_xlabel_justv", envir=env, inherits = FALSE)){
        svalue(vjust_spb) <- get(".strvalidator_plotStutter_gui_xlabel_justv", envir=env)
      }
      if(exists(".strvalidator_plotStutter_gui_theme", envir=env, inherits = FALSE)){
        svalue(f1_theme_drp) <- get(".strvalidator_plotStutter_gui_theme", envir=env)
      }
      if(exists(".strvalidator_plotStutter_gui_sex", envir=env, inherits = FALSE)){
        svalue(f1_drop_chk) <- get(".strvalidator_plotStutter_gui_sex", envir=env)
      }
      
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_plotStutter_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_plotStutter_gui_title", value=svalue(title_edt), envir=env)
      assign(x=".strvalidator_plotStutter_gui_title_chk", value=svalue(f1_titles_chk), envir=env)
      assign(x=".strvalidator_plotStutter_gui_x_title", value=svalue(x_title_edt), envir=env)
      assign(x=".strvalidator_plotStutter_gui_y_title", value=svalue(y_title_edt), envir=env)
      assign(x=".strvalidator_plotStutter_gui_points_shape", value=svalue(shape_spb), envir=env)
      assign(x=".strvalidator_plotStutter_gui_points_alpha", value=svalue(alpha_spb), envir=env)
      assign(x=".strvalidator_plotStutter_gui_points_jitter", value=svalue(jitter_txt), envir=env)
      assign(x=".strvalidator_plotStutter_gui_axes_y_min", value=svalue(y_min_txt), envir=env)
      assign(x=".strvalidator_plotStutter_gui_axes_y_max", value=svalue(y_max_txt), envir=env)
      assign(x=".strvalidator_plotStutter_gui_axes_x_min", value=svalue(x_min_txt), envir=env)
      assign(x=".strvalidator_plotStutter_gui_axes_x_max", value=svalue(x_max_txt), envir=env)
      assign(x=".strvalidator_plotStutter_gui_axes_scales", value=svalue(scales_opt), envir=env)
      assign(x=".strvalidator_plotStutter_gui_xlabel_size", value=svalue(size_txt), envir=env)
      assign(x=".strvalidator_plotStutter_gui_xlabel_angle", value=svalue(angle_spb), envir=env)
      assign(x=".strvalidator_plotStutter_gui_xlabel_justh", value=svalue(hjust_spb), envir=env)
      assign(x=".strvalidator_plotStutter_gui_xlabel_justv", value=svalue(vjust_spb), envir=env)
      assign(x=".strvalidator_plotStutter_gui_theme", value=svalue(f1_theme_drp), envir=env)
      assign(x=".strvalidator_plotStutter_gui_sex", value=svalue(f1_drop_chk), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_plotStutter_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotStutter_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_plotStutter_gui_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotStutter_gui_title", envir = env)
      }
      if(exists(".strvalidator_plotStutter_gui_title_chk", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotStutter_gui_title_chk", envir = env)
      }
      if(exists(".strvalidator_plotStutter_gui_x_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotStutter_gui_x_title", envir = env)
      }
      if(exists(".strvalidator_plotStutter_gui_y_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotStutter_gui_y_title", envir = env)
      }
      if(exists(".strvalidator_plotStutter_gui_points_shape", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotStutter_gui_points_shape", envir = env)
      }
      if(exists(".strvalidator_plotStutter_gui_points_alpha", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotStutter_gui_points_alpha", envir = env)
      }
      if(exists(".strvalidator_plotStutter_gui_points_jitter", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotStutter_gui_points_jitter", envir = env)
      }
      if(exists(".strvalidator_plotStutter_gui_axes_y_min", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotStutter_gui_axes_y_min", envir = env)
      }
      if(exists(".strvalidator_plotStutter_gui_axes_y_max", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotStutter_gui_axes_y_max", envir = env)
      }
      if(exists(".strvalidator_plotStutter_gui_axes_x_min", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotStutter_gui_axes_x_min", envir = env)
      }
      if(exists(".strvalidator_plotStutter_gui_axes_x_max", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotStutter_gui_axes_x_max", envir = env)
      }
      if(exists(".strvalidator_plotStutter_gui_axes_scales", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotStutter_gui_axes_scales", envir = env)
      }
      if(exists(".strvalidator_plotStutter_gui_xlabel_size", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotStutter_gui_xlabel_size", envir = env)
      }
      if(exists(".strvalidator_plotStutter_gui_xlabel_angle", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotStutter_gui_xlabel_angle", envir = env)
      }
      if(exists(".strvalidator_plotStutter_gui_xlabel_justh", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotStutter_gui_xlabel_justh", envir = env)
      }
      if(exists(".strvalidator_plotStutter_gui_xlabel_justv", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotStutter_gui_xlabel_justv", envir = env)
      }
      if(exists(".strvalidator_plotStutter_gui_theme", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotStutter_gui_theme", envir = env)
      }
      if(exists(".strvalidator_plotStutter_gui_sex", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotStutter_gui_sex", envir = env)
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
