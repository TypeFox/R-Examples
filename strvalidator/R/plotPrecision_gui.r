################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 06.01.2016: Fixed theme methods not found and added more themes.
# 11.11.2015: Added importFrom grid unit.c, gridExtra arrangeGrob, and ggplot2.
# 11.11.2015: Added more themes.
# 29.08.2015: Added importFrom.
# 11.10.2014: Added 'focus', added 'parent' parameter.
# 12.09.2014: Filter rows with Allele=NA (Fixes issue #6).
# 28.06.2014: Added help button and moved save gui checkbox.
# 08.05.2014: Implemented 'checkDataset'.
# 23.02.2014: Fixed column check for plots.
# 06.02.2014: Implemented theme and colour.
# 06.02.2014: Implemented new dot/box plot and plot all data (not only min/max)
# 22.01.2014: Fixed bug, different y-axis max value for complex plots.
# 20.01.2014: Implemented ggsave with workaround for complex plots.
# 07.12.2013: First version.

#' @title Plot Precision
#'
#' @description
#' GUI simplifying the creation of plots from precision data.
#'
#' @details Plot precision data for size, height, or data point as dotplot or
#' boxplot. Plot per marker or all in one. Use the mean value or the allele
#' designation as x-axis labels.
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
#' @importFrom grid unit textGrob grid.newpage grid.draw unit.c
# @importFrom gtable gtable_add_grob gtable
#' @importFrom utils help str head
#' @importFrom grDevices palette
#' @importFrom stats as.formula
#' @importFrom ggplot2 ggplot geom_point aes_string geom_boxplot facet_grid
#'  facet_wrap coord_cartesian labs xlab ylab theme element_text element_blank
#'  ggplotGrob theme_gray theme_bw theme_linedraw theme_light theme_dark
#'  theme_minimal theme_classic theme_void 
#' 
#' @seealso \url{http://docs.ggplot2.org/current/} for details on plot settings.
#' 

plotPrecision_gui <- function(env=parent.frame(), savegui=NULL, debug=FALSE, parent=NULL){
  
  # Global variables.
  .gData <- NULL
  .gDataName <- NULL
  .gPlot <- NULL
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # Main window.
  w <- gwindow(title="Plot precision", visible=FALSE)
  
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
    print(help("plotPrecision_gui", help_type="html"))
    
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
  
  glabel(text=" and the kit used:", container=f0)
  
  kit_drp <- gdroplist(items=getKit(), 
                       selected = 1,
                       editable = FALSE,
                       container = f0) 
  
  addHandlerChanged(dataset_drp, handler = function (h, ...) {
    
    val_obj <- svalue(dataset_drp)
    
    # Check if suitable.
    requiredCol <- c("Marker", "Allele")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       string="OL", stringcol="Allele",
                       env=env, parent=w, debug=debug)
    
    if(ok){
      
      # Load or change components.
      .gData <<- get(val_obj, envir=env)

      # Suggest name.
      svalue(f5_save_edt) <- paste(val_obj, "_ggplot", sep="")
      
      # Detect kit.
      kitIndex <- detectKit(.gData, index=TRUE)
      # Select in dropdown.
      svalue(kit_drp, index=TRUE) <- kitIndex
      
      # Enable buttons.
      enabled(f7_size_btn) <- TRUE
      enabled(f7_height_btn) <- TRUE
      enabled(f7_data_btn) <- TRUE
      enabled(f8_size_btn) <- TRUE
      enabled(f8_height_btn) <- TRUE
      enabled(f8_data_btn) <- TRUE
          
    } else {

      # Reset components.
      .gData <<- NULL
      svalue(f5_save_edt) <- ""
      svalue(dataset_drp, index=TRUE) <- 1
      
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
  
  grid1[2,1] <- glabel(text="X title:", container=grid1)
  grid1[2,2] <- x_title_edt <- gedit(text="",
                                     container=grid1)
  
  grid1[3,1] <- glabel(text="Y title:", container=grid1)
  grid1[3,2] <- y_title_edt <- gedit(text="",
                                     container=grid1)
  
  f1_facet_chk <- gcheckbox(text="Plot per marker",
                              checked=TRUE,
                              container=f1)
  
  f1g2 <- glayout(container = f1)
  f1g2[1,1] <- glabel(text="X axis:", anchor=c(-1 ,0), container=f1g2)
  f1g2[1,2] <- f1_axis_opt <- gradio(items=c("Mean","Allele"),
                                     selected = 2,
                                     horizontal = TRUE,
                                     container = f1g2)

  f1g2[2,1] <- glabel(text="Plot theme:", anchor=c(-1 ,0), container=f1g2)
  items_theme <- c("theme_grey()","theme_bw()","theme_linedraw()",
                   "theme_light()","theme_dark()","theme_minimal()",
                   "theme_classic()","theme_void()")
  f1g2[2,2] <- f1_theme_drp <- gdroplist(items = items_theme,
                                         selected = 1,
                                         container = f1g2)
  
  
  # FRAME 7 ###################################################################
  
  f7 <- gframe(text = "Plot precision data as dotplot",
               horizontal=FALSE,
               container = gv) 
  
  grid7 <- glayout(container = f7)
  
  grid7[1,1] <- f7_size_btn <- gbutton(text="Size",
                                           border=TRUE,
                                           container=grid7) 
  
  grid7[1,2] <- f7_height_btn <- gbutton(text="Height",
                                         border=TRUE,
                                         container=grid7) 

  grid7[1,3] <- f7_data_btn <- gbutton(text="Data point",
                                         border=TRUE,
                                         container=grid7) 

  addHandlerChanged(f7_size_btn, handler = function(h, ...) {

    # Check if suitable for plot.
    requiredCol <- c("Marker", "Allele", "Size")
  
    if(!all(requiredCol %in% colnames(.gData))){
      
      missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
      
      message <- paste("Additional columns required:\n",
                       paste(missingCol, collapse="\n"), sep="")
      
      gmessage(message, title="message",
               icon = "error",
               parent = w) 
      
    } else {
      
      enabled(f7_size_btn) <- FALSE
      .plot(what="Size", how="dotplot")
      enabled(f7_size_btn) <- TRUE
      
    }
    
    # Change save button.
    svalue(f5_save_btn) <- "Save as object"
    enabled(f5_save_btn) <- TRUE
      
  } )
  
  addHandlerChanged(f7_height_btn, handler = function(h, ...) {
    
    # Check if suitable for plot.
    requiredCol <- c("Marker", "Allele", "Height")
    
    if(!all(requiredCol %in% colnames(.gData))){
      
      missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
      
      message <- paste("Additional columns required:\n",
                       paste(missingCol, collapse="\n"), sep="")
      
      gmessage(message, title="message",
               icon = "error",
               parent = w) 
      
    } else {
      
      enabled(f7_height_btn) <- FALSE
      .plot(what="Height", how="dotplot")
      enabled(f7_height_btn) <- TRUE
      
    }

    # Change save button.
    svalue(f5_save_btn) <- "Save as object"
    enabled(f5_save_btn) <- TRUE
    
  } )
  
  addHandlerChanged(f7_data_btn, handler = function(h, ...) {
    
    # Check if suitable for plot.
    requiredCol <- c("Marker", "Allele", "Data.Point")
    
    if(!all(requiredCol %in% colnames(.gData))){
      
      missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
      
      message <- paste("Additional columns required:\n",
                       paste(missingCol, collapse="\n"), sep="")
      
      gmessage(message, title="message",
               icon = "error",
               parent = w) 
      
    } else {
      
      enabled(f7_data_btn) <- FALSE
      .plot(what="Data.Point", how="dotplot")
      enabled(f7_data_btn) <- TRUE
      
    }
    
    # Change save button.
    svalue(f5_save_btn) <- "Save as object"
    enabled(f5_save_btn) <- TRUE
    
  } )
  
  # FRAME 8 ###################################################################
  
  f8 <- gframe(text = "Plot precision data as boxplot",
               horizontal=FALSE,
               container = gv) 
  
  grid8 <- glayout(container = f8)
  
  grid8[1,1] <- f8_size_btn <- gbutton(text="Size",
                                         border=TRUE,
                                         container=grid8) 
  
  grid8[1,2] <- f8_height_btn <- gbutton(text="Height",
                                           border=TRUE,
                                           container=grid8) 
  
  grid8[1,3] <- f8_data_btn <- gbutton(text="Data point",
                                         border=TRUE,
                                         container=grid8) 
  
  addHandlerChanged(f8_size_btn, handler = function(h, ...) {

    # Check if suitable for plot.
    requiredCol <- c("Marker", "Allele", "Size")
    
    if(!all(requiredCol %in% colnames(.gData))){
      
      missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
      
      message <- paste("Additional columns required:\n",
                       paste(missingCol, collapse="\n"), sep="")
      
      gmessage(message, title="message",
               icon = "error",
               parent = w) 
      
    } else {
      
      enabled(f8_size_btn) <- FALSE
      .plot(what="Size", how="boxplot")
      enabled(f8_size_btn) <- TRUE
      
    }

    # Change save button.
    svalue(f5_save_btn) <- "Save as object"
    enabled(f5_save_btn) <- TRUE
    
  } )
  
  addHandlerChanged(f8_height_btn, handler = function(h, ...) {

    # Check if suitable for plot.
    requiredCol <- c("Marker", "Allele", "Height")
    
    if(!all(requiredCol %in% colnames(.gData))){
      
      missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
      
      message <- paste("Additional columns required:\n",
                       paste(missingCol, collapse="\n"), sep="")
      
      gmessage(message, title="message",
               icon = "error",
               parent = w) 
      
    } else {
      
      enabled(f8_height_btn) <- FALSE
      .plot(what="Height", how="boxplot")
      enabled(f8_height_btn) <- TRUE
    
    }
    
    # Change save button.
    svalue(f5_save_btn) <- "Save as object"
    enabled(f5_save_btn) <- TRUE
    
  } )
  
  addHandlerChanged(f8_data_btn, handler = function(h, ...) {
    
    # Check if suitable for plot.
    requiredCol <- c("Marker", "Allele", "Data.Point")
    
    if(!all(requiredCol %in% colnames(.gData))){
      
      missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
      
      message <- paste("Additional columns required:\n",
                       paste(missingCol, collapse="\n"), sep="")
      
      gmessage(message, title="message",
               icon = "error",
               parent = w) 
      
    } else {
      
      enabled(f8_data_btn) <- FALSE
      .plot(what="Data.Point", how="boxplot")
      enabled(f8_data_btn) <- TRUE
    
    }
    
    # Change save button.
    svalue(f5_save_btn) <- "Save as object"
    enabled(f5_save_btn) <- TRUE
    
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
  
  grid2[1,5] <- glabel(text="Colour:", container=grid2)
  grid2[1,6] <- colour_drp <- gdroplist(items=c("white",palette()),
                                    selected=2,
                                    container = grid2)
  
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
  
  
  .plot <- function(what, how){
    
    # Get values.
    val_titles <- svalue(f1_titles_chk)
    val_title <- svalue(title_edt)
    val_xtitle <- svalue(x_title_edt)
    val_ytitle <- svalue(y_title_edt)
    val_shape <- as.numeric(svalue(shape_spb))
    val_alpha <- as.numeric(svalue(alpha_spb))
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
    val_data <- .gData
    val_facet <- svalue(f1_facet_chk)
    val_colour <- svalue(colour_drp)
    val_axis <- svalue(f1_axis_opt)
    val_theme <- svalue(f1_theme_drp)
    
    if(debug){
      print("ARGUMENTS:")
      print("what")
      print(what)
      print("how")
      print(how)
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
      print("colour")
      print(val_colour)
      print("val_axis")
      print(val_axis)
      print("val_theme")
      print(val_theme)
    }
    
    if(is.factor(val_data$Marker)){
      
    }
    
    if (!is.null(val_data) && !is.na(val_data)){
      
      if(debug){
        print("BEFORE PLOTTING:")
        print("str(val_data)")
        print(str(val_data))
        print("levels(val_data$Allele)")
        print(levels(val_data$Allele))
        print("levels(val_data$Marker)")
        print(levels(val_data$Marker))
      }
      
      # Calculate mean and deviation ------------------------------------------
      
      what.mean <- paste(what, "Mean", sep=".")
      
      # Add new column.
      val_data[[what.mean]] <- NA
      val_data[[what]] <- as.numeric(val_data[[what]])
      
      # Remove NA's
      if(any(is.na(val_data$Allele))){
        tmp1 <- nrow(val_data)
        val_data <- val_data[!is.na(val_data$Allele), ]
        tmp2 <- nrow(val_data)
        message(paste(tmp1-tmp2, "NA rows removed!"))
      }
      
      # Get all markers.
      marker <- unique(val_data$Marker)
      
      # Loop over all markers.
      for(m in seq(along=marker)){
        
        # Select current marker.
        selMarker <- val_data$Marker == marker[m]
        
        # Get unique alleles.
        allele <- unique(val_data[selMarker, ]$Allele)
        
        # Loop over all alleles.        
        for(a in seq(along=allele)){
          
          # Select current allele.
          selAllele <- val_data$Allele == allele[a]
          
          # Combine selection.
          selection <- selMarker & selAllele
          
          # Calculate mean and save in dataframe.
          val_data[selection, ][[what.mean]] <- mean(val_data[selection, ][[what]], na.rm=TRUE) 
          
        }
      }
      
      # Calculate deviation.
      dev <- val_data[[what]] - val_data[[what.mean]]
      val_data <- data.frame(Marker=val_data$Marker,
                             Allele=val_data$Allele,
                             Mean=val_data[[what.mean]],
                             Deviation=dev,
                             stringsAsFactors=FALSE)
      
      # Make sorted allele factors (use low values for X/Y).
      numericAlleles <- unique(val_data$Allele)
      numericAlleles <- gsub("X", "0.0001", numericAlleles, ignore.case = TRUE)
      numericAlleles <- gsub("Y", "0.0002", numericAlleles, ignore.case = TRUE)
      orderedAlleles <- order(as.numeric(numericAlleles))
      val_data$Allele <- factor(val_data$Allele,levels=unique(val_data$Allele)[orderedAlleles])
      
      # TODO: NB! although plotting min/max values are very tidy,
      #       it may not be very informative. Include as an option?
      #         # Calculate deviation
      #         neg <- val_data$Size.Min - val_data$Size.Mean
      #         pos <- val_data$Size.Max - val_data$Size.Mean
      #         val_data <- data.frame(Marker=rep(val_data$Marker, 2),
      #                                Mean=val_data$Size.Mean,
      #                                Value=c(neg, pos),
      #                                Deviation=rep(c("Min","Max"), each=nrow(val_data)),
      #                                stringsAsFactors=FALSE)
      
      if(debug){
        print("val_data after calculating deviation:")
        print(str(val_data))
        print(head(val_data))
      }
      
      # End calculate mean and deviation --------------------------------------
      
      # Call functions.
      # Add color information.
      if(!"Dye" %in% names(val_data)){
        val_data <- addColor(data=val_data, kit=val_kit,
                             need="Dye", debug=debug)
        message("'Dye' added to dataset!")
      }
      
      # Sort by marker in kit
      val_data <- sortMarker(data=val_data,
                             kit=val_kit,
                             add.missing.levels = TRUE)
      
      
      if(debug){
        print("AFTER SORT MARKERS:")
        print("str(val_data)")
        print(str(val_data))
        print("levels(val_data$Allele)")
        print(levels(val_data$Allele))
        print("levels(val_data$Marker)")
        print(levels(val_data$Marker))
      }
      
      # Create titles..
      if(val_titles){
        mainTitle <- val_title
        xTitle <- val_xtitle
        yTitle <- val_ytitle
      } else {
        
        if(what == "Size"){
          
          mainTitle <- "Allele size range for allelic ladders"

          yTitle <- "Deviation from mean (bp)"
          
          if(val_axis == "Mean"){
            xTitle <- "Mean size in basepair (bp)"
          } else if(val_axis == "Allele"){
            xTitle <- "Allele"
          } else {
            warning(paste("val_axis=", val_axis, "not implemented!"))
          }
          
        } else if(what == "Height"){
          
          mainTitle <- "Allele height range for allelic ladders"
          
          yTitle <- "Deviation from mean (RFU)"
          
          if(val_axis == "Mean"){
            xTitle <- "Mean height in relative fluorescent units (RFU)"
          } else if(val_axis == "Allele"){
            xTitle <- "Allele"
          } else {
            warning(paste("val_axis=", val_axis, "not implemented!"))
          }
          
        } else if(what == "Data.Point"){
          
          mainTitle <- "Allele data point range for allelic ladders"
          
          yTitle <- "Deviation from mean"
          
          if(val_axis == "Mean"){
            xTitle <- "Mean scan number in data points"
          } else if(val_axis == "Allele"){
            xTitle <- "Allele"
          } else {
            warning(paste("val_axis=", val_axis, "not implemented!"))
          }
          
        } else {
          
          warning(paste("what=", what, "not implemented!")) 
          
        }
        
      }

      # TODO: NB! although plotting min/max values are very tidy,
      #       it may not be very informative. Include as an option?
#       # Create plot.
#       gp <- ggplot(val_data, aes_string(x="Mean", y="Value", color="Deviation"),
#                    shape=val_shape, alpha=val_alpha)
#       gp <- gp + geom_point()
      
      # Create plot.
      if(how == "dotplot"){
        
        # Create dotplot.
        gp <- ggplot(val_data)
        gp <- gp + geom_point(aes_string(x=val_axis, y="Deviation"),
                              alpha=val_alpha, shape=val_shape, colour=val_colour)
        # gp <- gp + facet_wrap(~Marker) # TODO: is this needed?
        
      } else if(how == "boxplot"){
        
        # Create boxplot (per allele).
        gp <- ggplot(val_data)
        gp <- gp + geom_boxplot( aes_string(x=val_axis, y="Deviation"),
                                 alpha=val_alpha, shape=val_shape, fill=val_colour)
        
      } else {
        
        warning(paste("how=", how, "not implemented!")) 
        
      }
      
      # Apply theme.
      gp <- gp + eval(parse(text=val_theme))
      
      # Default is simple plot (length(val_ncol)==1).
      val_ncol <- 0
      
      # Facet plot.
      if(val_facet){

        # Check if 'simple' or 'complex' plotting:
        # Get Marker and Dye column.
        markerDye <- val_data[c("Marker","Dye")]
        # Extract unique elements.
        uniqueMarkerDye <- markerDye[!duplicated(markerDye),]
        # Calculate number of unique columns per dye.
        val_ncol <- unique(table(uniqueMarkerDye$Dye))
        
        # Facet plot.        
        gp <- gp + facet_grid("Dye ~ Marker")
        # NB! 'facet_wrap' does not seem to support strings.
        #     Use 'as.formula(paste("string1", "string2"))' as a workaround.
        gp <- gp + facet_wrap(as.formula(paste("~ Marker")), ncol=val_ncol,
                              drop=FALSE, scales=val_scales)
        
      }
      
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
      
      # Titles and legends.
      gp <- gp + labs(title=mainTitle)
      gp <- gp + xlab(xTitle)
      gp <- gp + ylab(yTitle)
      gp <- gp + theme(axis.text.x=element_text(angle=val_angle,
                                                hjust=val_hjust,
                                                vjust=val_vjust,
                                                size=val_size))
      
      # Check plot type.
      if(length(val_ncol) == 1){
        # Simple plot, equal number of markers per dye.
        
        if(debug){
          print(paste("Simple plot, val_ncol:",
                      paste(val_ncol, collapse=", ")))
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
        
        # With guide:
        # Extract the legend from the 'simple' plot.
        #guide <- gtable::gtable_filter(ggplotGrob(gp), pattern="guide")
        
        # Get y max/min to be able to use same scale across plots.
        yMax <- max(val_data$Deviation, na.rm=TRUE)
        yMin <- min(val_data$Deviation, na.rm=TRUE)
        
        if(debug){
          print("yMax/yMin:")
          print(yMax)
          print(yMin)
        }
      
        # Get kit colors and convert to dyes.
        dyes <- unique(getKit(val_kit, what="Color")$Color)
        dyes <- addColor(dyes, have="Color", need="Dye")
        # Number of dyes.
        noDyes <- length(dyes)
        # Number of rows in table object (one per dye + title + x title).
        noRows <- length(dyes) + 2
        
        # Create table object.
        # With guide:
        # Note: width(1.5 for y-title, and the rest for plots + guides)
        #       height(1.5 for plot title, equal for each plot, and 1.5 for x-title)
#         g <- gtable::gtable(widths=grid::unit.c(grid::unit(1.5, "lines"),
#                                                 grid::unit(1, "null"),
#                                                 sum(guide$widths)),
#                             heights = grid::unit(c(1.5,rep(1,noDyes),1.5),
#                                                  c("line", rep("null", noDyes), "line")))
        
        # Without guide:
        # Note: width(1.5 for y-title, and the rest for plots + margin)
        #       height(1.5 for plot title, equal for each plot, and 1.5 for x-title)
        g <- gtable::gtable(widths=grid::unit.c(grid::unit(1.5, "lines"),
                                                grid::unit(1, "null"),grid::unit(1.5, "lines")),
                            heights = grid::unit(c(1.5,rep(1,noDyes),1.5),
                                                 c("line", rep("null", noDyes), "line")))
        # Add titles.
        g <- gtable::gtable_add_grob(g, grid::textGrob(mainTitle), t=1,b=1,l=2,r=2)
        g <- gtable::gtable_add_grob(g, grid::textGrob(xTitle), t=noRows ,b=noRows ,l=2,r=2)
        g <- gtable::gtable_add_grob(g, grid::textGrob(yTitle, rot=90), t=1,b=noRows ,l=1,r=1)
        
        # With guide:
        # Add the legend to the table object.
        #g <- gtable::gtable_add_grob(g,guide , t=1,b=noRows,l=3,r=3)
      
        if(debug){
          print("Complex plot area created. Looping over dyes.")
        }
        
        # Loop over all dyes.
        for(d in seq(along=dyes)){
          
          # Create a plot for the current subset.
          if(how == "dotplot"){

            # Create a plot for the current subset.
            gp <- ggplot(subset(val_data, val_data$Dye == dyes[d]))
            gp <- gp + geom_point(aes_string(x = val_axis, y = "Deviation"),
                                  alpha = val_alpha, shape=val_shape, colour=val_colour)
            
          } else if(how == "boxplot"){
            
            # Create a plot for the current subset.
            gp <- ggplot(subset(val_data, val_data$Dye == dyes[d])) 
            gp <- gp + geom_boxplot(aes_string(x = val_axis, y = "Deviation"),
                                    alpha = val_alpha, shape=val_shape, fill=val_colour)
            
          } else {
            
            warning(paste("how=", how, "not implemented!")) 
            
          }
          
          if(debug){
            print(paste("Complex plot base pane", d, "created."))
          }
          
          # Plot settings.
          gp <- gp + facet_grid("Dye ~ Marker", scales=val_scales)
          
          # Apply theme.
          gp <- gp + eval(parse(text=val_theme))
          
          # Set margin around each plot. Note: top, right, bottom, left.
          gp <- gp + theme(plot.margin = grid::unit(c(0.25, 0, 0, 0), "lines"))

          # Restrict y axis.
          if(!is.na(val_ymin) && !is.na(val_ymax)){
            val_y <- c(val_ymin, val_ymax)
          } else {
            # Make scales work on multiple plots.
            if(val_scales != "free" && val_scales != "free_y"){
              val_y <- c(yMin, yMax)
            } else {
              val_y <- NULL
            }
          }
          # Restrict x axis.
          if(!is.na(val_xmin) && !is.na(val_xmax)){
            val_x <- c(val_xmin, val_xmax)
          } else {
            val_x <- NULL
          }
          # Zoom in without dropping observations.
          gp <- gp + coord_cartesian(xlim=val_x, ylim=val_y)

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
      if(exists(".strvalidator_plotPrecision_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_plotPrecision_gui_savegui", envir=env)
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
      if(exists(".strvalidator_plotPrecision_gui_title", envir=env, inherits = FALSE)){
        svalue(title_edt) <- get(".strvalidator_plotPrecision_gui_title", envir=env)
      }
      if(exists(".strvalidator_plotPrecision_gui_title_chk", envir=env, inherits = FALSE)){
        svalue(f1_titles_chk) <- get(".strvalidator_plotPrecision_gui_title_chk", envir=env)
      }
      if(exists(".strvalidator_plotPrecision_gui_x_title", envir=env, inherits = FALSE)){
        svalue(x_title_edt) <- get(".strvalidator_plotPrecision_gui_x_title", envir=env)
      }
      if(exists(".strvalidator_plotPrecision_gui_y_title", envir=env, inherits = FALSE)){
        svalue(y_title_edt) <- get(".strvalidator_plotPrecision_gui_y_title", envir=env)
      }
      if(exists(".strvalidator_plotPrecision_gui_points_shape", envir=env, inherits = FALSE)){
        svalue(shape_spb) <- get(".strvalidator_plotPrecision_gui_points_shape", envir=env)
      }
      if(exists(".strvalidator_plotPrecision_gui_points_alpha", envir=env, inherits = FALSE)){
        svalue(alpha_spb) <- get(".strvalidator_plotPrecision_gui_points_alpha", envir=env)
      }
      if(exists(".strvalidator_plotPrecision_gui_points_colour", envir=env, inherits = FALSE)){
        svalue(colour_drp) <- get(".strvalidator_plotPrecision_gui_points_colour", envir=env)
      }
      if(exists(".strvalidator_plotPrecision_gui_axes_y_min", envir=env, inherits = FALSE)){
        svalue(y_min_txt) <- get(".strvalidator_plotPrecision_gui_axes_y_min", envir=env)
      }
      if(exists(".strvalidator_plotPrecision_gui_axes_y_max", envir=env, inherits = FALSE)){
        svalue(y_max_txt) <- get(".strvalidator_plotPrecision_gui_axes_y_max", envir=env)
      }
      if(exists(".strvalidator_plotPrecision_gui_axes_x_min", envir=env, inherits = FALSE)){
        svalue(x_min_txt) <- get(".strvalidator_plotPrecision_gui_axes_x_min", envir=env)
      }
      if(exists(".strvalidator_plotPrecision_gui_axes_x_max", envir=env, inherits = FALSE)){
        svalue(x_max_txt) <- get(".strvalidator_plotPrecision_gui_axes_x_max", envir=env)
      }
      if(exists(".strvalidator_plotPrecision_gui_axes_scales", envir=env, inherits = FALSE)){
        svalue(scales_opt) <- get(".strvalidator_plotPrecision_gui_axes_scales", envir=env)
      }
      if(exists(".strvalidator_plotPrecision_gui_xlabel_size", envir=env, inherits = FALSE)){
        svalue(size_txt) <- get(".strvalidator_plotPrecision_gui_xlabel_size", envir=env)
      }
      if(exists(".strvalidator_plotPrecision_gui_xlabel_angle", envir=env, inherits = FALSE)){
        svalue(angle_spb) <- get(".strvalidator_plotPrecision_gui_xlabel_angle", envir=env)
      }
      if(exists(".strvalidator_plotPrecision_gui_xlabel_justh", envir=env, inherits = FALSE)){
        svalue(hjust_spb) <- get(".strvalidator_plotPrecision_gui_xlabel_justh", envir=env)
      }
      if(exists(".strvalidator_plotPrecision_gui_xlabel_justv", envir=env, inherits = FALSE)){
        svalue(vjust_spb) <- get(".strvalidator_plotPrecision_gui_xlabel_justv", envir=env)
      }
      if(exists(".strvalidator_plotPrecision_gui_facet", envir=env, inherits = FALSE)){
        svalue(f1_facet_chk) <- get(".strvalidator_plotPrecision_gui_facet", envir=env)
      }
      if(exists(".strvalidator_plotPrecision_gui_theme", envir=env, inherits = FALSE)){
        svalue(f1_theme_drp) <- get(".strvalidator_plotPrecision_gui_theme", envir=env)
      }
      
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_plotPrecision_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_plotPrecision_gui_title", value=svalue(title_edt), envir=env)
      assign(x=".strvalidator_plotPrecision_gui_title_chk", value=svalue(f1_titles_chk), envir=env)
      assign(x=".strvalidator_plotPrecision_gui_x_title", value=svalue(x_title_edt), envir=env)
      assign(x=".strvalidator_plotPrecision_gui_y_title", value=svalue(y_title_edt), envir=env)
      assign(x=".strvalidator_plotPrecision_gui_points_shape", value=svalue(shape_spb), envir=env)
      assign(x=".strvalidator_plotPrecision_gui_points_alpha", value=svalue(alpha_spb), envir=env)
      assign(x=".strvalidator_plotPrecision_gui_points_colour", value=svalue(colour_drp), envir=env)
      assign(x=".strvalidator_plotPrecision_gui_axes_y_min", value=svalue(y_min_txt), envir=env)
      assign(x=".strvalidator_plotPrecision_gui_axes_y_max", value=svalue(y_max_txt), envir=env)
      assign(x=".strvalidator_plotPrecision_gui_axes_x_min", value=svalue(x_min_txt), envir=env)
      assign(x=".strvalidator_plotPrecision_gui_axes_x_max", value=svalue(x_max_txt), envir=env)
      assign(x=".strvalidator_plotPrecision_gui_axes_scales", value=svalue(scales_opt), envir=env)
      assign(x=".strvalidator_plotPrecision_gui_xlabel_size", value=svalue(size_txt), envir=env)
      assign(x=".strvalidator_plotPrecision_gui_xlabel_angle", value=svalue(angle_spb), envir=env)
      assign(x=".strvalidator_plotPrecision_gui_xlabel_justh", value=svalue(hjust_spb), envir=env)
      assign(x=".strvalidator_plotPrecision_gui_xlabel_justv", value=svalue(vjust_spb), envir=env)
      assign(x=".strvalidator_plotPrecision_gui_facet", value=svalue(f1_facet_chk), envir=env)
      assign(x=".strvalidator_plotPrecision_gui_theme", value=svalue(f1_theme_drp), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_plotPrecision_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPrecision_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_plotPrecision_gui_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPrecision_gui_title", envir = env)
      }
      if(exists(".strvalidator_plotPrecision_gui_title_chk", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPrecision_gui_title_chk", envir = env)
      }
      if(exists(".strvalidator_plotPrecision_gui_x_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPrecision_gui_x_title", envir = env)
      }
      if(exists(".strvalidator_plotPrecision_gui_y_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPrecision_gui_y_title", envir = env)
      }
      if(exists(".strvalidator_plotPrecision_gui_points_shape", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPrecision_gui_points_shape", envir = env)
      }
      if(exists(".strvalidator_plotPrecision_gui_points_alpha", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPrecision_gui_points_alpha", envir = env)
      }
      if(exists(".strvalidator_plotPrecision_gui_points_colour", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPrecision_gui_points_colour", envir = env)
      }
      if(exists(".strvalidator_plotPrecision_gui_axes_y_min", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPrecision_gui_axes_y_min", envir = env)
      }
      if(exists(".strvalidator_plotPrecision_gui_axes_y_max", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPrecision_gui_axes_y_max", envir = env)
      }
      if(exists(".strvalidator_plotPrecision_gui_axes_x_min", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPrecision_gui_axes_x_min", envir = env)
      }
      if(exists(".strvalidator_plotPrecision_gui_axes_x_max", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPrecision_gui_axes_x_max", envir = env)
      }
      if(exists(".strvalidator_plotPrecision_gui_axes_scales", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPrecision_gui_axes_scales", envir = env)
      }
      if(exists(".strvalidator_plotPrecision_gui_xlabel_size", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPrecision_gui_xlabel_size", envir = env)
      }
      if(exists(".strvalidator_plotPrecision_gui_xlabel_angle", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPrecision_gui_xlabel_angle", envir = env)
      }
      if(exists(".strvalidator_plotPrecision_gui_xlabel_justh", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPrecision_gui_xlabel_justh", envir = env)
      }
      if(exists(".strvalidator_plotPrecision_gui_xlabel_justv", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPrecision_gui_xlabel_justv", envir = env)
      }
      if(exists(".strvalidator_plotPrecision_gui_facet", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPrecision_gui_facet", envir = env)
      }
      if(exists(".strvalidator_plotPrecision_gui_theme", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPrecision_gui_theme", envir = env)
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

