################################################################################
# TODO LIST
# TODO: Add colours to plot strips according to:
#       http://stackoverflow.com/questions/19440069/ggplot2-facet-wrap-strip-color-based-on-variable-in-data-set

################################################################################
# CHANGE LOG (last 20 changes)
# 06.01.2016: Fixed theme methods not found and added more themes.
# 11.11.2015: Added importFrom gridExtra arrangeGrob, and ggplot2.
# 11.11.2015: Added importFrom grid.
# 11.11.2015: Added more themes.
# 29.08.2015: Added importFrom.
# 14.12.2014: Updated to handle gender -> sex.marker option in getKit.
# 08.12.2014: First version.

#' @title Plot Pull-up
#'
#' @description
#' GUI simplifying the creation of plots from pull-up data.
#'
#' @details Select a dataset to plot and the typing kit used (if not autodetected).
#' Plot pull-up peak ratio versus the peak height of the known allele
#' Automatic plot titles can be replaced by custom titles.
#' Sex markers can be excluded. A name for the result is automatiaclly suggested.
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
#' @importFrom stats as.formula
#' @importFrom grid unit textGrob grid.newpage grid.draw
#' @importFrom gridExtra arrangeGrob
#' @importFrom ggplot2 ggplot aes_string geom_point position_jitter facet_grid
#'  facet_wrap scale_colour_manual coord_cartesian guides guide_legend theme
#'  element_text labs xlab ylab element_blank ggplotGrob theme_gray theme_bw
#'  theme_linedraw theme_light theme_dark theme_minimal theme_classic
#'  theme_void 
#' 
#' @seealso \url{http://docs.ggplot2.org/current/} for details on plot settings.

plotPullup_gui <- function(env=parent.frame(), savegui=NULL, debug=FALSE, parent=NULL){
  
  # Global variables.
  .gData <- NULL
  .gDataName <- NULL
  .gPlot <- NULL
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # Main window.
  w <- gwindow(title="Plot pull-up", visible=FALSE)
  
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
    print(help("plotPullup_gui", help_type="html"))
    
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
    requiredCol <- c("Sample.Name", "Marker", "Dye", "Allele", "Height",
                     "Size", "Data.Point", "P.Marker", "P.Dye", "P.Allele",
                     "P.Height", "P.Size", "P.Data.Point", "Delta", "Ratio")
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
      enabled(plot_height_btn) <- TRUE
      enabled(plot_allele_btn) <- TRUE
      
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
  
  f7 <- gframe(text = "Plot pull-up data",
               horizontal=FALSE,
               container = gv) 
  
  grid7 <- glayout(container = f7)
  
  grid7[1,1] <- plot_height_btn <- gbutton(text="Ratio vs. Height",
                                           border=TRUE,
                                           container=grid7) 
  
  grid7[1,2] <- plot_allele_btn <- gbutton(text="Ratio vs. Allele",
                                           border=TRUE,
                                           container=grid7) 
  
  addHandlerChanged(plot_height_btn, handler = function(h, ...) {
    
    # Check if suitable for plot.
    requiredCol <- c("Sample.Name", "Marker", "Dye", "Allele", "Height",
                     "Size", "Data.Point", "P.Marker", "P.Dye", "P.Allele",
                     "P.Height", "P.Size", "P.Data.Point", "Delta", "Ratio")
    
    if(!all(requiredCol %in% colnames(.gData))){
      
      missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
      
      message <- paste("Additional columns required:\n",
                       paste(missingCol, collapse="\n"), sep="")
      
      gmessage(message, title="message",
               icon = "error",
               parent = w) 
      
    } else {
      
      enabled(plot_height_btn) <- FALSE
      .plotPullup(what="Height")
      enabled(plot_height_btn) <- TRUE
      
    }
    
  } )
  
  addHandlerChanged(plot_allele_btn, handler = function(h, ...) {
    
    # Check if suitable for plot.
    requiredCol <- c("Sample.Name", "Marker", "Dye", "Allele", "Height",
                     "Size", "Data.Point", "P.Marker", "P.Dye", "P.Allele",
                     "P.Height", "P.Size", "P.Data.Point", "Delta", "Ratio")
    
    if(!all(requiredCol %in% colnames(.gData))){
      
      missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
      
      message <- paste("Additional columns required:\n",
                       paste(missingCol, collapse="\n"), sep="")
      
      gmessage(message, title="message",
               icon = "error",
               parent = w) 
      
    } else {
      
      enabled(plot_allele_btn) <- FALSE
      .plotPullup(what="Allele")
      enabled(plot_allele_btn) <- TRUE
      
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
                                          selected = 2,
                                          horizontal = FALSE,
                                          container = grid3)
  
  # FRAME 4 ###################################################################
  
  e4 <- gexpandgroup(text="X labels",
                     horizontal=FALSE,
                     container = f1)
  
  grid4 <- glayout(container = e4)
  
  grid4[1,1] <- glabel(text="Text size (pts):", container=grid4)
  grid4[1,2] <- e4_size_edt <- gedit(text="8", width=4, container=grid4)
  
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
  
  
  .plotPullup <- function(what){
    
    # Get values.
    val_titles <- svalue(f1_titles_chk)
    val_title <- svalue(title_edt)
    val_xtitle <- svalue(x_title_edt)
    val_ytitle <- svalue(y_title_edt)
    val_shape <- as.numeric(svalue(e2_shape_spb))
    val_alpha <- as.numeric(svalue(e2_alpha_spb))
    val_jitter <- as.numeric(svalue(e2_jitter_edt))
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
    }
    
    # Declare variables.
    ymax <- NULL  # For complex plots.
    ymin <- NULL  # For complex plots.
    
    if (!is.na(.gData) && !is.null(.gData)){
      
      # Call functions.
      # Sort by marker in kit and add Dye levels.
      .gData <- sortMarker(data=.gData,
                           kit=val_kit,
                           add.missing.levels = TRUE)
      
      # Get kit colors and convert to dyes.
      dyes <- unique(getKit(kit=val_kit, what="Color")$Color)
      dyes <- addColor(data=dyes, have="Color", need="Dye")
      
      # Factor and maintain correct order of levels.
      .gData$Dye <- factor(.gData$Dye, levels=dyes)
      .gData$P.Dye <- factor(.gData$P.Dye, levels=dyes)
      
      # Drop sex markers.
      if(val_drop){
        
        # Get sex marker.
        sexMarkers <- getKit(kit=val_kit, what="Sex.Marker")
        
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
      if(!is.numeric(.gData$Height)){
        .gData$Height <- as.numeric(as.character(.gData$Height))
        message("'Height' not numeric, converting to numeric.")
        
      }
      
      # Ratio must be numeric (not string).
      if(!is.numeric(.gData$Ratio)){
        .gData$Ratio <- as.numeric(as.character(.gData$Ratio))
        message("'Ratio' not numeric, converting to numeric.")
        
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
      
      # Make palette.
      val_palette <- unique(getKit(kit=val_kit, what="Color")$Color)
      val_palette <- addColor(data=val_palette, have="Color", need="R.Color")
      
      if(debug){
        print("Before plot: str(.gData)")
        print(str(.gData))
        print("Number of columns")
        print(val_ncol)
        print("val_palette")
        print(val_palette)
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
        
        if(what == "Height"){
          
          mainTitle <- "Pull-up ratio"
          xTitle <- "Allele peak height (RFU)"
          yTitle <- "Ratio"
          
        } else if(what == "Allele"){
          
          mainTitle <- "Pull-up ratio"
          xTitle <- "Allele designation"
          yTitle <- "Ratio"
          
        } else {
          
          stop(paste("what=", what, " not handled!"))
          
        }
        
      }
      
      if(debug){
        print("Titles:")
        print(mainTitle)
        print(xTitle)
        print(yTitle)
      }
      
      # Construct plot differently.
      if(length(val_ncol) == 1){
        # Simple plot, equal number of markers per dye.
        
        if(debug){
          print("Simple plot.")
        }
        
        # Select what to plot and create default titles.
        if(what == "Height"){
          
          gp <- ggplot(.gData, aes_string(x="Height", y="Ratio", colour="P.Dye"))
          
        } else if(what == "Allele"){
          
          gp <- ggplot(.gData, aes_string(x="Allele", y="Ratio", colour="P.Dye"))
          
        }
        
        if(debug){
          print("Plot created.")
        }
        
        # Apply theme.
        gp <- gp + eval(parse(text=val_theme))
        
        # Plot settings.
        gp <- gp + geom_point(shape=val_shape, alpha=val_alpha,
                              position=position_jitter(height = 0, width=val_jitter))
        gp <- gp + facet_grid("Dye ~ Marker")
        # NB! 'facet_wrap' does not seem to support strings.
        #     Use 'as.formula(paste("string1", "string2"))' as a workaround.
        gp <- gp + facet_wrap(as.formula(paste("~", "Marker")), ncol=val_ncol,
                             drop=FALSE, scales=val_scales)
        # Add manual scale to get colours according to 'P.Dye'.
        # NB! important to use drop=FALSE if not all are represented in data.
        gp <- gp + scale_colour_manual(guide=FALSE, values=val_palette, drop=FALSE)
        
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
          print(paste("Plot zoomed to xlim:", paste(val_x, collapse=","),
                      "ylim:", paste(val_y, collapse=",")))
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
        
      } else if (length(val_ncol) > 1){
        # Complex plot, unequal number of markers per dye.
        
        if(debug){
          print("Complex plot.")
        }
        
        if(val_scales %in% c("fixed","free_x")){
          ymax <- max(.gData$Ratio, na.rm=TRUE) * 1.05
          ymin <- min(.gData$Ratio, na.rm=TRUE) * 0.95
        }
        
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
        gLevel <- addColor(data=gLevel, kit=val_kit)
        
        # Loop over all dyes.
        for(d in seq(along=dyes)){
          
          # Get data for current dye.
          gDataSub <- .gData[.gData$Dye == dyes[d],]

          # Get current markers/levels.
          gDyeLevel <- as.character(gLevel$Marker[gLevel$Dye==dyes[d]])
          
          # Can't handle zero rows.
          if(nrow(gDataSub)==0){
            # Dummy data must have at least one value (not NA) for each property used to create plot (including facets)
            # e.g. x/y/colour
            tmp <- data.frame(Sample.Name="", Marker=gDyeLevel, Dye=dyes[d], Allele="", Height=0,
                              Size=0, Data.Point=0, P.Marker=NA, P.Dye=dyes[d], P.Allele=NA,
                              P.Height=0, P.Size=0, P.Data.Point=0, Delta=0, Ratio=0)

            # Combine with 
            gDataSub <- plyr::rbind.fill(gDataSub, tmp)
            
          }
          
          # Refactor to levels of current dye (and maintain order).
          gDataSub$Marker <- factor(gDataSub$Marker, levels=gDyeLevel)
          gDataSub$Dye <- factor(dyes[d])
          gDataSub$P.Dye <- factor(gDataSub$P.Dye, levels=dyes) # Must contain all levels.
          
          # Create a plot for the current subset.
          # Select what to plot.
          if(what == "Height"){
            
            gp <- ggplot(gDataSub, aes_string(x = "Height", y = "Ratio", colour="P.Dye"))
            
          } else if(what == "Allele"){
            
            gp <- ggplot(gDataSub, aes_string(x = "Allele", y = "Ratio", colour="P.Dye"))
            
          }
          
          # Apply theme.
          gp <- gp + eval(parse(text=val_theme))
          
          # Plot settings.
          gp <- gp + geom_point(shape=val_shape, alpha = val_alpha,
                                position = position_jitter(height = 0, width = val_jitter))
          # Add manual scale to get colours according to 'P.Dye'.
          # NB! important to use drop=FALSE if not all are represented in data.
          gp <- gp + scale_colour_manual(guide=FALSE, values=val_palette, drop=FALSE)
          gp <- gp + facet_grid("Dye ~ Marker", scales=val_scales, drop = FALSE) # Keep dye labels.
          #gp <- gp + facet_grid("~ Marker", scales=val_scales)  # No dye labels.
          
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
            print(paste("Plot zoomed to xlim:", paste(val_x, collapse=","),
                        "ylim:", paste(val_y, collapse=",")))
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
      if(exists(".strvalidator_plotPullup_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_plotPullup_gui_savegui", envir=env)
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
      if(exists(".strvalidator_plotPullup_gui_title", envir=env, inherits = FALSE)){
        svalue(title_edt) <- get(".strvalidator_plotPullup_gui_title", envir=env)
      }
      if(exists(".strvalidator_plotPullup_gui_title_chk", envir=env, inherits = FALSE)){
        svalue(f1_titles_chk) <- get(".strvalidator_plotPullup_gui_title_chk", envir=env)
      }
      if(exists(".strvalidator_plotPullup_gui_x_title", envir=env, inherits = FALSE)){
        svalue(x_title_edt) <- get(".strvalidator_plotPullup_gui_x_title", envir=env)
      }
      if(exists(".strvalidator_plotPullup_gui_y_title", envir=env, inherits = FALSE)){
        svalue(y_title_edt) <- get(".strvalidator_plotPullup_gui_y_title", envir=env)
      }
      if(exists(".strvalidator_plotPullup_gui_sex", envir=env, inherits = FALSE)){
        svalue(f1_drop_chk) <- get(".strvalidator_plotPullup_gui_sex", envir=env)
      }
      if(exists(".strvalidator_plotPullup_gui_points_shape", envir=env, inherits = FALSE)){
        svalue(e2_shape_spb) <- get(".strvalidator_plotPullup_gui_points_shape", envir=env)
      }
      if(exists(".strvalidator_plotPullup_gui_points_alpha", envir=env, inherits = FALSE)){
        svalue(e2_alpha_spb) <- get(".strvalidator_plotPullup_gui_points_alpha", envir=env)
      }
      if(exists(".strvalidator_plotPullup_gui_points_jitter", envir=env, inherits = FALSE)){
        svalue(e2_jitter_edt) <- get(".strvalidator_plotPullup_gui_points_jitter", envir=env)
      }
      if(exists(".strvalidator_plotPullup_gui_axes_y_min", envir=env, inherits = FALSE)){
        svalue(e3_y_min_edt) <- get(".strvalidator_plotPullup_gui_axes_y_min", envir=env)
      }
      if(exists(".strvalidator_plotPullup_gui_axes_y_max", envir=env, inherits = FALSE)){
        svalue(e3_y_max_edt) <- get(".strvalidator_plotPullup_gui_axes_y_max", envir=env)
      }
      if(exists(".strvalidator_plotPullup_gui_axes_x_min", envir=env, inherits = FALSE)){
        svalue(e3_x_min_edt) <- get(".strvalidator_plotPullup_gui_axes_x_min", envir=env)
      }
      if(exists(".strvalidator_plotPullup_gui_axes_x_max", envir=env, inherits = FALSE)){
        svalue(e3_x_max_edt) <- get(".strvalidator_plotPullup_gui_axes_x_max", envir=env)
      }
      if(exists(".strvalidator_plotPullup_gui_axes_scales", envir=env, inherits = FALSE)){
        svalue(e3_scales_opt) <- get(".strvalidator_plotPullup_gui_axes_scales", envir=env)
      }
      if(exists(".strvalidator_plotPullup_gui_xlabel_size", envir=env, inherits = FALSE)){
        svalue(e4_size_edt) <- get(".strvalidator_plotPullup_gui_xlabel_size", envir=env)
      }
      if(exists(".strvalidator_plotPullup_gui_xlabel_angle", envir=env, inherits = FALSE)){
        svalue(e4_angle_spb) <- get(".strvalidator_plotPullup_gui_xlabel_angle", envir=env)
      }
      if(exists(".strvalidator_plotPullup_gui_xlabel_justh", envir=env, inherits = FALSE)){
        svalue(e4_hjust_spb) <- get(".strvalidator_plotPullup_gui_xlabel_justh", envir=env)
      }
      if(exists(".strvalidator_plotPullup_gui_xlabel_justv", envir=env, inherits = FALSE)){
        svalue(e4_vjust_spb) <- get(".strvalidator_plotPullup_gui_xlabel_justv", envir=env)
      }
      if(exists(".strvalidator_plotPullup_gui_theme", envir=env, inherits = FALSE)){
        svalue(f1_theme_drp) <- get(".strvalidator_plotPullup_gui_theme", envir=env)
      }
      
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_plotPullup_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_plotPullup_gui_sex", value=svalue(f1_drop_chk), envir=env)
      assign(x=".strvalidator_plotPullup_gui_title", value=svalue(title_edt), envir=env)
      assign(x=".strvalidator_plotPullup_gui_title_chk", value=svalue(f1_titles_chk), envir=env)
      assign(x=".strvalidator_plotPullup_gui_x_title", value=svalue(x_title_edt), envir=env)
      assign(x=".strvalidator_plotPullup_gui_y_title", value=svalue(y_title_edt), envir=env)
      assign(x=".strvalidator_plotPullup_gui_points_shape", value=svalue(e2_shape_spb), envir=env)
      assign(x=".strvalidator_plotPullup_gui_points_alpha", value=svalue(e2_alpha_spb), envir=env)
      assign(x=".strvalidator_plotPullup_gui_points_jitter", value=svalue(e2_jitter_edt), envir=env)
      assign(x=".strvalidator_plotPullup_gui_axes_y_min", value=svalue(e3_y_min_edt), envir=env)
      assign(x=".strvalidator_plotPullup_gui_axes_y_max", value=svalue(e3_y_max_edt), envir=env)
      assign(x=".strvalidator_plotPullup_gui_axes_x_min", value=svalue(e3_x_min_edt), envir=env)
      assign(x=".strvalidator_plotPullup_gui_axes_x_max", value=svalue(e3_x_max_edt), envir=env)
      assign(x=".strvalidator_plotPullup_gui_axes_scales", value=svalue(e3_scales_opt), envir=env)
      assign(x=".strvalidator_plotPullup_gui_xlabel_size", value=svalue(e4_size_edt), envir=env)
      assign(x=".strvalidator_plotPullup_gui_xlabel_angle", value=svalue(e4_angle_spb), envir=env)
      assign(x=".strvalidator_plotPullup_gui_xlabel_justh", value=svalue(e4_hjust_spb), envir=env)
      assign(x=".strvalidator_plotPullup_gui_xlabel_justv", value=svalue(e4_vjust_spb), envir=env)
      assign(x=".strvalidator_plotPullup_gui_theme", value=svalue(f1_theme_drp), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_plotPullup_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPullup_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_plotPullup_gui_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPullup_gui_title", envir = env)
      }
      if(exists(".strvalidator_plotPullup_gui_title_chk", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPullup_gui_title_chk", envir = env)
      }
      if(exists(".strvalidator_plotPullup_gui_x_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPullup_gui_x_title", envir = env)
      }
      if(exists(".strvalidator_plotPullup_gui_y_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPullup_gui_y_title", envir = env)
      }
      if(exists(".strvalidator_plotPullup_gui_sex", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPullup_gui_sex", envir = env)
      }
      if(exists(".strvalidator_plotPullup_gui_points_shape", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPullup_gui_points_shape", envir = env)
      }
      if(exists(".strvalidator_plotPullup_gui_points_alpha", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPullup_gui_points_alpha", envir = env)
      }
      if(exists(".strvalidator_plotPullup_gui_points_jitter", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPullup_gui_points_jitter", envir = env)
      }
      if(exists(".strvalidator_plotPullup_gui_axes_y_min", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPullup_gui_axes_y_min", envir = env)
      }
      if(exists(".strvalidator_plotPullup_gui_axes_y_max", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPullup_gui_axes_y_max", envir = env)
      }
      if(exists(".strvalidator_plotPullup_gui_axes_x_min", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPullup_gui_axes_x_min", envir = env)
      }
      if(exists(".strvalidator_plotPullup_gui_axes_x_max", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPullup_gui_axes_x_max", envir = env)
      }
      if(exists(".strvalidator_plotPullup_gui_axes_scales", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPullup_gui_axes_scales", envir = env)
      }
      if(exists(".strvalidator_plotPullup_gui_xlabel_size", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPullup_gui_xlabel_size", envir = env)
      }
      if(exists(".strvalidator_plotPullup_gui_xlabel_angle", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPullup_gui_xlabel_angle", envir = env)
      }
      if(exists(".strvalidator_plotPullup_gui_xlabel_justh", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPullup_gui_xlabel_justh", envir = env)
      }
      if(exists(".strvalidator_plotPullup_gui_xlabel_justv", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPullup_gui_xlabel_justv", envir = env)
      }
      if(exists(".strvalidator_plotPullup_gui_theme", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPullup_gui_theme", envir = env)
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
