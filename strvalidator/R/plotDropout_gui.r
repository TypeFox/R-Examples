################################################################################
# TODO LIST
# TODO: Number of decimals on x axis as an option.
# TODO: Custom colors.
# TODO: Just one plot button, and a dropdown to select column to sort by. Ascend/descen+numeric/character
# TODO: ...NOT FINISHED!!!
#      (need to change if/when a preferred dropout method has beed decided)

################################################################################
# CHANGE LOG (last 20 changes)
# 11.11.2015: Added importFrom ggplot2.
# 29.08.2015: Added importFrom.
# 16.05.2015: Fixed issue#10 colors hardcoded as ESX17 for dotplot.
# 11.10.2014: Added 'focus', added 'parent' parameter.
# 28.06.2014: Added help button and moved save gui checkbox.
# 08.05.2014: Implemented 'checkDataset'.
# 15.04.2014: Fixed position_jitter height now fixed to zero (prev. default).
# 17.02.2014: Fixed NA in title for ecdp.
# 17.02.2014: Fixed heatmap by 'H' loosing samples with equal 'H'.
# 20.01.2014: Changed 'saveImage_gui' for 'ggsave_gui'.
# 05.11.2013: Fixed not possible to limit both y/x axes.
# 04.11.2013: Added edcf plot.
# 01.11.2013: Added 'override titles' option.
# 23.10.2013: Added save as image.
# 20.10.2013: Added plot by sample name. Fixed x-label font size not changing.
# 18.09.2013: Updated to support new 'addColor' function, replacing 'addDye'.
# 17.09.2013: Added missing plot by concentration.
# 18.07.2013: Check before overwrite object.
# 15.07.2013: Save as ggplot object to workspace instead of image.
# 15.07.2013: Added save GUI settings.

#' @title Plot Drop-out Events
#'
#' @description
#' GUI simplifying the creation of plots from dropout data.
#'
#' @details Plot dropout data as heatmap arranged by, average peak height, 
#' amount, concentration, or sample name. It is also possible to plot the
#' empirical cumulative distribution (ecdp) of the peak heights of surviving heterozygote
#' alleles (with dropout of the parter allele), or a dotplot of all dropout events.
#' The peak height of homozygote alleles can be included in the ecdp.
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
#' @importFrom scales pretty_breaks
#' @importFrom utils help str
#' @importFrom grDevices rgb
#' @importFrom ggplot2 ggplot aes_string geom_tile scale_fill_manual guides
#'  guide_legend theme element_text labs ylab xlab scale_y_discrete scale_x_discrete
#'  stat_ecdf scale_colour_discrete scale_x_continuous scale_y_continuous
#'  coord_cartesian geom_point position_jitter scale_colour_manual
#' 
#' @seealso \url{http://docs.ggplot2.org/current/} for details on plot settings.

plotDropout_gui <- function(env=parent.frame(), savegui=NULL, debug=FALSE, parent=NULL){
  
  # Global variables.
  .gData <- NULL
  .gDataColumns <- NULL
  .gPlot <- NULL
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
    #print(head(data))
  }
  
  # Main window.  
  w <- gwindow(title="Plot dropout data", visible=FALSE)
  
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
    print(help("plotDropout_gui", help_type="html"))
    
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
    requiredCol <- c("Sample.Name", "Marker", "Allele", "Height",
                     "Dropout", "Rfu", "Heterozygous")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       env=env, parent=w, debug=debug)
    
    if(ok){
      
      # Load or change components.
      .gData <<- get(val_obj, envir=env)
      .gDataColumns <<- names(.gData)
      
      # Suggest name.
      svalue(f5_save_edt) <- paste(val_obj, "_ggplot", sep="")
      # Detect kit.
      kitIndex <- detectKit(.gData, index=TRUE)
      # Select in dropdown.
      svalue(kit_drp, index=TRUE) <- kitIndex
      
      # Enable plot buttons.
      enabled(f7_plot_h_btn) <- TRUE
      enabled(f7_plot_amount_btn) <- TRUE
      enabled(f7_plot_conc_btn) <- TRUE
      enabled(f7_plot_sample_btn) <- TRUE
      enabled(f8_plot_ecdf_btn) <- TRUE
      enabled(f8_plot_dot_btn) <- TRUE
      
    } else {
      
      # Reset components.
      .gData <<- NULL
      .gDataColumns <<- NULL
      svalue(f5_save_edt) <- ""
      
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
  grid1[1,2] <- f1_title_edt <- gedit(text="",
                                   width=40,
                                   container=grid1)
  
  grid1[2,1] <- glabel(text="X title:", container=grid1)
  grid1[2,2] <- f1_xtitle_edt <- gedit(text="",
                                     container=grid1)

  grid1[3,1] <- glabel(text="Y title:", container=grid1)
  grid1[3,2] <- f1_ytitle_edt <- gedit(text="",
                                     container=grid1)

  # FRAME 7 ###################################################################
  
  f7 <- gframe(text = "Plot heatmap by",
               horizontal=TRUE,
               container = gv) 
  
  f7_plot_h_btn <- gbutton(text="Average peak height",
                        border=TRUE,
                        container=f7) 

  f7_plot_amount_btn <- gbutton(text="Amount",
                             border=TRUE,
                             container=f7) 

  f7_plot_conc_btn <- gbutton(text="Concentration",
                           border=TRUE,
                           container=f7) 

  f7_plot_sample_btn <- gbutton(text="Sample",
                             border=TRUE,
                             container=f7) 
  
  addHandlerChanged(f7_plot_h_btn, handler = function(h, ...) {

    # Check if suitable for plot.
    requiredCol <- c("Sample.Name", "Marker", "Dropout", "H")
    
    if(!all(requiredCol %in% colnames(.gData))){
      
      missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]

      message <- paste("Additional columns required:\n",
                       paste(missingCol, collapse="\n"), sep="")
      
      gmessage(message, title="message",
               icon = "error",
               parent = w) 
      
    } else {

      enabled(f7_plot_h_btn) <- FALSE
      .plotStutter(what="heat_h")
      enabled(f7_plot_h_btn) <- TRUE
      
    }
    
    # Change save button.
    svalue(f5_save_btn) <- "Save as object"
    enabled(f5_save_btn) <- TRUE
    
  } )

  addHandlerChanged(f7_plot_amount_btn, handler = function(h, ...) {
    
    # Check if suitable for plot.
    requiredCol <- c("Sample.Name", "Marker", "Dropout", "Amount")
    
    if(!all(requiredCol %in% colnames(.gData))){
      
      missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]

      message <- paste("Additional columns required:\n",
                       paste(missingCol, collapse="\n"), sep="")
      
      gmessage(message, title="message",
               icon = "error",
               parent = w) 
      
    } else {
      
      enabled(f7_plot_amount_btn) <- FALSE
      .plotStutter(what="heat_amount")
      enabled(f7_plot_amount_btn) <- TRUE
      
    }
    
    # Change save button.
    svalue(f5_save_btn) <- "Save as object"
    enabled(f5_save_btn) <- TRUE
    
  } )
  
  addHandlerChanged(f7_plot_conc_btn, handler = function(h, ...) {
    
    # Check if suitable for plot.
    requiredCol <- c("Sample.Name", "Marker", "Dropout", "Concentration")
    
    if(!all(requiredCol %in% colnames(.gData))){
      
      missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]

      message <- paste("Additional columns required:\n",
                       paste(missingCol, collapse="\n"), sep="")
      
      gmessage(message, title="message",
               icon = "error",
               parent = w) 
      
    } else {
      
      enabled(f7_plot_conc_btn) <- FALSE
      .plotStutter(what="heat_conc")
      enabled(f7_plot_conc_btn) <- TRUE
      
    }
    
    # Change save button.
    svalue(f5_save_btn) <- "Save as object"
    enabled(f5_save_btn) <- TRUE
    
  } )
  
  addHandlerChanged(f7_plot_sample_btn, handler = function(h, ...) {
    
    # Check if suitable for plot.
    requiredCol <- c("Sample.Name", "Marker", "Dropout")
    
    if(!all(requiredCol %in% colnames(.gData))){
      
      missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
      
      message <- paste("Additional columns required:\n",
                       paste(missingCol, collapse="\n"), sep="")
      
      gmessage(message, title="message",
               icon = "error",
               parent = w) 
      
    } else {
      
      enabled(f7_plot_sample_btn) <- FALSE
      .plotStutter(what="sample")
      enabled(f7_plot_sample_btn) <- TRUE
      
    }
    
    # Change save button.
    svalue(f5_save_btn) <- "Save as object"
    enabled(f5_save_btn) <- TRUE
    
  } )

  # FRAME 8 ###################################################################
  
  f8 <- gframe(text = "Other plots",
               horizontal=TRUE,
               container = gv) 
  
  f8_plot_ecdf_btn <- gbutton(text="ecdp",
                              border=TRUE,
                              container=f8)
  
  f8_hom_chk <- gcheckbox(text="Plot homozygous peaks.",
                          checked=FALSE,
                          container=f8)
  
  f8_plot_dot_btn <- gbutton(text="Dotplot",
                              border=TRUE,
                              container=f8)
  
  addHandlerChanged(f8_plot_ecdf_btn, handler = function(h, ...) {
    
    # Check if suitable for plot.
    requiredCol <- c("Sample.Name", "Marker", "Dropout", "Height", "Heterozygous")
    
    if(!all(requiredCol %in% colnames(.gData))){
      
      missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
      
      message <- paste("Additional columns required:\n",
                       paste(missingCol, collapse="\n"), sep="")
      
      gmessage(message, title="message",
               icon = "error",
               parent = w) 
      
    } else {
      
      enabled(f8_plot_ecdf_btn) <- FALSE
      .plotStutter(what="ecdf")
      enabled(f8_plot_ecdf_btn) <- TRUE
      
    }
    
    # Change save button.
    svalue(f5_save_btn) <- "Save as object"
    enabled(f5_save_btn) <- TRUE
    
  } )
  
  addHandlerChanged(f8_plot_dot_btn, handler = function(h, ...) {
    
    # Check if suitable for plot.
    requiredCol <- c("Sample.Name", "Marker", "Dropout", "Height", "Heterozygous")
    
    if(!all(requiredCol %in% colnames(.gData))){
      
      missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
      
      message <- paste("Additional columns required:\n",
                       paste(missingCol, collapse="\n"), sep="")
      
      gmessage(message, title="message",
               icon = "error",
               parent = w) 
      
    } else {
      
      enabled(f8_plot_dot_btn) <- FALSE
      .plotStutter(what="dot")
      enabled(f8_plot_dot_btn) <- TRUE
      
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
  
  # FRAME 3 ###################################################################
  
  e3 <- gexpandgroup(text="Axes (applies to continous axes)",
                     horizontal=FALSE,
                     container = f1)
  
  grid3 <- glayout(container = e3, spacing = 1)
  
  grid3[1,1:2] <- glabel(text="Limit Y axis (min-max)", container=grid3)
  grid3[2,1] <- e3_y_min_edt <- gedit(text="", width=5, container=grid3)
  grid3[2,2] <- e3_y_max_edt <- gedit(text="", width=5, container=grid3)
  
  grid3[3,1:2] <- glabel(text="Limit X axis (min-max)", container=grid3)
  grid3[4,1] <- e3_x_min_edt <- gedit(text="", width=5, container=grid3)
  grid3[4,2] <- e3_x_max_edt <- gedit(text="", width=5, container=grid3)
  
  # FRAME 4 ###################################################################
  
  e4 <- gexpandgroup(text="X labels",
                     horizontal=FALSE,
                     container = f1)
  
  grid4 <- glayout(container = e4)
  
  grid4[1,1] <- glabel(text="Text size (pts):", container=grid4)
  grid4[1,2] <- e4_size_txt <- gedit(text="10", width=4, container=grid4)

  grid4[1,3] <- glabel(text="Angle:", container=grid4)
  grid4[1,4] <- e4_angle_spb <- gspinbutton (from=0, to=360, by=1,
                                         value=270,
                                         container=grid4) 

  grid4[2,1] <- glabel(text="Justification (v/h):", container=grid4)
  grid4[2,2] <- e4_vjust_spb <- gspinbutton (from=0, to=1, by=0.1,
                                          value=0.3,
                                          container=grid4)

  grid4[2,3] <- e4_hjust_spb <- gspinbutton (from=0, to=1, by=0.1,
                                          value=0,
                                          container=grid4)

  # FUNCTIONS #################################################################
  
  
  .plotStutter <- function(what){
    
    # Get values.
    val_titles <- svalue(f1_titles_chk)
    val_title <- svalue(f1_title_edt)
    val_xtitle <- svalue(f1_xtitle_edt)
    val_ytitle <- svalue(f1_ytitle_edt)
    val_angle <- as.numeric(svalue(e4_angle_spb))
    val_vjust <- as.numeric(svalue(e4_vjust_spb))
    val_hjust <- as.numeric(svalue(e4_hjust_spb))
    val_size <- as.numeric(svalue(e4_size_txt))
    val_kit <- svalue(kit_drp)
    val_hom <- svalue(f8_hom_chk)
    val_ymin <- as.numeric(svalue(e3_y_min_edt))
    val_ymax <- as.numeric(svalue(e3_y_max_edt))
    val_xmin <- as.numeric(svalue(e3_x_min_edt))
    val_xmax <- as.numeric(svalue(e3_x_max_edt))
    
    if(debug){
      print("val_title")
      print(val_title)
      print("val_xtitle")
      print(val_xtitle)
      print("val_ytitle")
      print(val_ytitle)
      print("val_angle")
      print(val_angle)
      print("val_vjust")
      print(val_vjust)
      print("val_hjust")
      print(val_hjust)
      print("val_size")
      print(val_size)
      print("val_hom")
      print(val_hom)
      print("str(.gData)")
      print(str(.gData))
    }
    
    
    if (!is.na(.gData) && !is.null(.gData)){
      
      
      # Call functions.
      
      # Color information.
      if(is.null(.gData$Dye)){
        .gData <- addColor(data=.gData, kit=val_kit, need="Dye")
      }

      # Sort by marker in kit
      .gData <- sortMarker(data=.gData,
                          kit=val_kit,
                          add.missing.levels = TRUE)
      
      
      if(debug){
        print("Before plot: str(.gData)")
        print(str(.gData))
      }

      # Create custom titles.
      if(val_titles){
        mainTitle <- val_title
        xTitle <- val_xtitle
        yTitle <- val_ytitle
      } 
        
      # Select what to plot and create default titles.
      if(what == "heat_h"){

        # Create default titles.
        if(!val_titles){
          mainTitle <- "Allele and locus dropout"
          xTitle <- "Average peak height 'H' (RFU)"
          yTitle <- "Marker"
        }

        # Sort according to H.
        if(!is.numeric(.gData$H)){
          .gData$H <- as.numeric(.gData$H)
          message("'H' converted to numeric.")
        }
        .gData <- .gData[order(.gData$H),]
        
        # Add H to sample name.
        .gData$Sample.Name <- paste(.gData$H, " (", .gData$Sample.Name, ")", sep="")
        
        # Create factors.
        .gData$Dropout <- factor(.gData$Dropout, levels=c(0,1,2))
        .gData$Sample.Name <- factor(.gData$Sample.Name,
                                     levels=unique(.gData$Sample.Name))
        
        # Create x labels.
        xlabels <- .gData[!duplicated(.gData[, c("Sample.Name", "H")]), ]$H
        xlabels <- round(as.double(xlabels), digits=3)
        
        # Define colours.
        col <- c(rgb(0,0.737,0), rgb(1,0.526,1), rgb(0.526,0,0.526))
        
        # Create plot.
        gp <- ggplot(.gData, aes_string(x = "Sample.Name", y = "Marker", fill = "Dropout"))
        gp <- gp + geom_tile(colour = "white") #OK
        gp <- gp + scale_fill_manual(values=col, name="Dropout", breaks=c("0", "1", "2"),
                                     labels=c("none", "allele", "locus"))
        gp <- gp + guides(fill = guide_legend(reverse=TRUE)) # OK
        gp <- gp + theme(axis.text.x=element_text(angle=val_angle,
                                                  hjust = val_hjust,
                                                  vjust = val_vjust,
                                                  size = val_size))
        
        gp <- gp + labs(title=mainTitle)
        gp <- gp + ylab(yTitle)
        gp <- gp + xlab(xTitle)
        
        # Reverse y-axis and relabel x-ticks.
        gp <- gp + scale_y_discrete(limits = rev(levels(.gData$Marker))) + 
          scale_x_discrete(labels=formatC(xlabels, 0, format = "f")) +
          theme(axis.text.x=element_text(family="sans", face="bold", size=val_size))
        
      } else if (what == "heat_amount") {
        
        # Create default titles.
        if(!val_titles){
          mainTitle <- "Allele and locus dropout"
          xTitle <- "Amount amplified DNA (ng)"
          yTitle <- "Marker"
        }
        
        # Sort according to average amount of DNA
        if(!is.numeric(.gData$Amount)){
          .gData$Amount <- as.numeric(.gData$Amount)
          message("'Amount' converted to numeric.")
        }
        .gData <- .gData[order(.gData$Amount),]
        
        # Add amount to sample name.
        .gData$Sample.Name <- paste(.gData$Amount, " (", .gData$Sample.Name, ")", sep="")
        
        # Create factors.
        .gData$Dropout <- factor(.gData$Dropout, levels=c(0,1,2))
        .gData$Sample.Name <- factor(.gData$Sample.Name,
                                     levels=unique(.gData$Sample.Name))
        
        # Create x labels.
        xlabels <- .gData[!duplicated(.gData[, c("Sample.Name", "Amount")]), ]$Amount
        xlabels <- round(as.double(xlabels), digits=3)
        
        # Define colours.
        col <- c(rgb(0,0.737,0), rgb(1,0.526,1), rgb(0.526,0,0.526))

        # Create plot.
        gp <- ggplot(.gData, aes_string(x = "Sample.Name", y = "Marker", fill = "Dropout"))
        gp <- gp + geom_tile(colour = "white") #OK
        gp <- gp + scale_fill_manual(values=col, name="Dropout", breaks=c("0", "1", "2"),
                                               labels=c("none", "allele", "locus"))
        gp <- gp + guides(fill = guide_legend(reverse=TRUE)) # OK
        gp <- gp + theme(axis.text.x=element_text(angle=val_angle,
                                                  hjust = val_hjust,
                                                  vjust = val_vjust,
                                                  size = val_size))
        gp <- gp + labs(title=mainTitle)
        gp <- gp + ylab(yTitle)
        gp <- gp + xlab(xTitle)

        # Reverse y-axis and relabel x-ticks.
        gp <- gp + scale_y_discrete(limits = rev(levels(.gData$Marker))) + 
          scale_x_discrete(labels=formatC(xlabels, 3, format = "f")) +
          theme(axis.text.x=element_text(family="sans", face="bold", size=val_size))
        
        
      } else if (what == "heat_conc") {
        # Sort according to concentration of DNA.

        # Create default titles.
        if(!val_titles){
          mainTitle <- "Allele and locus dropout"
          xTitle <- "Concentration (ng/uL)"
          yTitle <- "Marker"
        }
        
        # Sort according to concentration.
        if(!is.numeric(.gData$Concentration)){
          .gData$Concentration <- as.numeric(.gData$Concentration)
          message("'Concentration' converted to numeric.")
        }
        .gData <- .gData[order(.gData$Concentration),]
        
        # Add concentration to sample name.
        .gData$Sample.Name <- paste(.gData$Concentration, " (", .gData$Sample.Name, ")", sep="")
        
        # Create factors.
        .gData$Dropout <- factor(.gData$Dropout, levels=c(0,1,2))
        .gData$Sample.Name <- factor(.gData$Sample.Name,
                                     levels=unique(.gData$Sample.Name))
        
        # Create x labels.
        xlabels <- .gData[!duplicated(.gData[, c("Sample.Name", "Concentration")]), ]$Concentration
        xlabels <- round(as.double(xlabels), digits=4)
        
        # Define colours.
        col <- c(rgb(0,0.737,0), rgb(1,0.526,1), rgb(0.526,0,0.526))
        
        # Create plot.
        gp <- ggplot(.gData, aes_string(x = "Sample.Name", y = "Marker", fill = "Dropout"))
        gp <- gp + geom_tile(colour = "white") #OK
        gp <- gp + scale_fill_manual(values=col, name="Dropout", breaks=c("0", "1", "2"),
                                     labels=c("none", "allele", "locus"))
        gp <- gp + guides(fill = guide_legend(reverse=TRUE)) # OK
        gp <- gp + theme(axis.text.x=element_text(angle=val_angle,
                                                  hjust = val_hjust,
                                                  vjust = val_vjust,
                                                  size = val_size))
        gp <- gp + labs(title=mainTitle)
        gp <- gp + ylab(yTitle)
        gp <- gp + xlab(xTitle)
        
        # Reverse y-axis and relabel x-ticks.
        gp <- gp + scale_y_discrete(limits = rev(levels(.gData$Marker))) + 
          scale_x_discrete(labels=formatC(xlabels, 4, format = "f")) +
          theme(axis.text.x=element_text(family="sans", face="bold", size=val_size))

      } else if (what == "sample") {
        # Sort according to sample name.
        
        # Create default titles.
        if(!val_titles){
          mainTitle <- "Allele and locus dropout"
          xTitle <- "Sample name"
          yTitle <- "Marker"
        }
        
        # Sort according to sample name.
        .gData <- .gData[order(.gData$Sample.Name),]
        
        # Create factors.
        .gData$Dropout <- factor(.gData$Dropout)
        
        # Create x labels.
        xlabels <- .gData[!duplicated(.gData[, "Sample.Name"]), ]$Sample.Name
        
        # Define colours.
        col <- c(rgb(0,0.737,0), rgb(1,0.526,1), rgb(0.526,0,0.526))
        
        # Create plot.
        gp <- ggplot(.gData, aes_string(x = "Sample.Name", y = "Marker", fill = "Dropout"))
        gp <- gp + geom_tile(colour = "white") #OK
        gp <- gp + scale_fill_manual(values=col, name="Dropout", breaks=c("0", "1", "2"),
                                     labels=c("none", "allele", "locus"))
        gp <- gp + guides(fill = guide_legend(reverse=TRUE)) # OK
        gp <- gp + theme(axis.text.x=element_text(angle=val_angle,
                                                  hjust = val_hjust,
                                                  vjust = val_vjust,
                                                  size = val_size))
        gp <- gp + labs(title=mainTitle)
        gp <- gp + ylab(yTitle)
        gp <- gp + xlab(xTitle)
        
        # Reverse y-axis and relabel x-ticks.
        gp <- gp + scale_y_discrete(limits = rev(levels(.gData$Marker))) + 
          scale_x_discrete(labels=xlabels) +
          theme(axis.text.x=element_text(family="sans", face="bold", size=val_size))
        
      } else if (what == "ecdf") {
        # Plot empirical cumulative distribution.

        # Remove NA in dropout col.
        # NB! THIS HAS TO BE CHANGED WHEN A DROPOUT MODEL HAS BEEN SELECTED!
        n0 <- nrow(.gData)
        .gData <- .gData[!is.na(.gData$Dropout),]
        n1 <- nrow(.gData)
        message(paste("Analyse ", n1,
                      " rows (",n0-n1,
                      " rows with NA in Dropout removed.",
                      sep=""))
        
        if(val_hom){
          
          # Remove locus dropouts.
          # NB! THIS HAS TO BE CHANGED WHEN A DROPOUT MODEL HAS BEEN SELECTED!
          n0 <- nrow(.gData)
          .gData <- .gData[!is.na(.gData$Height),]
          n1 <- nrow(.gData)
          message(paste("Analyse ", n1,
                        " rows (",n0-n1,
                        " NA rows i.e. locus dropout, removed from column 'Height').",
                        sep=""))
          
          # Remove locus dropout=2.
          # NB! THIS HAS TO BE CHANGED WHEN A DROPOUT MODEL HAS BEEN SELECTED!
          n0 <- nrow(.gData)
          .gData <- .gData[.gData$Dropout!=2,]
          n1 <- nrow(.gData)
          message(paste("Analyse ", n1,
                        " rows (",n0-n1,
                        " rows with 2 in Dropout removed.",
                        sep=""))

          # Remove heterozygous loci without dropout.
          n0 <- nrow(.gData)
          .gData <- .gData[!(.gData$Heterozygous==1 & .gData$Dropout==0),]
          n1 <- nrow(.gData)
          message(paste("Analyse ", n1,
                        " rows (",n0-n1,
                        " heterozygous loci without dropout removed.",
                        sep=""))
          
        } else {
          
          # Remove non-dropouts.
          n0 <- nrow(.gData)
          .gData <- .gData[.gData$Dropout==1,]
          n1 <- nrow(.gData)
          message(paste("Analyse ", n1,
                        " rows (",n0-n1,
                        " non-dropouts removed).",
                        sep=""))

        }

        # Create plot.
        if(val_hom){
          
          # Create default titles.
          if(!val_titles){
            mainTitle <- paste("Empirical cumulative distribution for",
                               sum(.gData$Dropout==1) ,
                               "heterozygous alleles (with dropout) and",
                               sum(.gData$Dropout==0) , "homozygous peaks")
            xTitle <- "Peak height (RFU)"
            yTitle <- "Cumulative probability"
          }

          # NB! Convert numeric to character (continous to discrete).
          # To avoid Error: Continuous value supplied to discrete scale.
          .gData$Heterozygous <- as.character(.gData$Heterozygous)

          # With homozygous data and heterozygous dropout data.
          gp <- ggplot(data=.gData, aes_string(x="Height", color="Heterozygous"))
          gp <- gp + stat_ecdf(data=subset(.gData, .gData$Heterozygous=="0"))
          gp <- gp + stat_ecdf(data=subset(.gData, .gData$Heterozygous=="1"))

          # Add legend.
          gp <- gp + scale_colour_discrete(name="Alleles",
                                           breaks=c("0", "1"),
                                           labels=c("Homozygous", "Heterozygous"))
          
        } else {

          # Create default titles.
          if(!val_titles){
            mainTitle <- paste("Empirical cumulative distribution for",
                               sum(.gData$Dropout==1),
                               "heterozygous alleles (with dropout of the sister allele)")
            xTitle <- "Peak height of surviving allele (RFU)"
            yTitle <- "Cumulative probability"
          }

          # With heterozygous dropout data.
          gp <- ggplot(.gData) + stat_ecdf(aes_string(x="Height"))
          
        }
        # TODO: Add optional threshold line.
        #Fn(t) = #{xi <= t}/n = 1/n sum(i=1,n) Indicator(xi <= t).
        # x = rfu, Fn(t) = probability
        # Or bootstrap for "confidence" interval..
        
        # Add titles and settings.        
        gp <- gp + theme(axis.text.x=element_text(angle=val_angle,
                                                  vjust = val_vjust,
                                                  size = val_size))
        gp <- gp + labs(title=mainTitle)
        gp <- gp + ylab(yTitle)
        gp <- gp + xlab(xTitle)
        gp <- gp + scale_x_continuous(breaks = scales::pretty_breaks())
        gp <- gp + scale_y_continuous(breaks = seq(0, 1, 0.1))
        
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
        
        
      } else if (what == "dot") {
        # Plot dropouts per locus.
        
        # NA heights.
        # NB! THIS HAS TO BE CHANGED WHEN A DROPOUT MODEL HAS BEEN SELECTED!
        n0 <- nrow(.gData)
        .gData <- .gData[!is.na(.gData$Height),]
        n1 <- nrow(.gData)
        message(paste("Analyse ", n1,
                      " rows (",n0-n1,
                      " NA rows removed from column 'Height').",
                      sep=""))
        
        # NA Dropouts.
        # NB! THIS HAS TO BE CHANGED WHEN A DROPOUT MODEL HAS BEEN SELECTED!
        n0 <- nrow(.gData)
        .gData <- .gData[!is.na(.gData$Dropout),]
        n1 <- nrow(.gData)
        message(paste("Analyse ", n1,
                      " rows (",n0-n1,
                      " NA rows removed from column 'Dropout').",
                      sep=""))

        # Remove non-dropouts.
        n0 <- nrow(.gData)
        .gData <- .gData[.gData$Dropout==1,]
        n1 <- nrow(.gData)
        message(paste("Analyse ", n1,
                      " rows (",n0-n1,
                      " non-dropouts removed).",
                      sep=""))
          
        # Create default titles.
        if(!val_titles){
          mainTitle <- paste(nrow(.gData),
                             "heterozygous alleles with dropout of the sister allele")
          xTitle <- "Locus"
          yTitle <- "Peak height of surviving allele (RFU)"
        }
        
        # Create plot.
        plotColor <- getKit(kit=val_kit, what="Color")
        plotColor <- unique(plotColor$Color)
        plotColor <- addColor(plotColor, need="R.Color", have="Color")

        # Create plot.
        gp <- ggplot(data=.gData, aes_string(x="Marker", y="Height"))
        
        # NB! This colour is only a grouping variable, NOT plot color.
        gp <- gp + geom_point(data=.gData, mapping = aes_string(colour = "Dye"),
                              position = position_jitter(height = 0, width = 0.2)) 

        # Specify colour values must be strings, NOT factors!
        # NB! The plot colours are specified as here as strings.
        # NB! Custom colours work on DATA AS SORTED FACTOR + COLOR CHARACTER.
        gp <- gp + scale_colour_manual(guide=FALSE, values=as.character(plotColor), drop=FALSE)

        # Add titles and settings.        
        gp <- gp + theme(axis.text.x=element_text(angle=val_angle,
                                                  vjust = val_vjust,
                                                  size = val_size))
        gp <- gp + labs(title=mainTitle)
        gp <- gp + ylab(yTitle)
        gp <- gp + xlab(xTitle)
        #gp <- gp + scale_y_continuous(breaks = seq(0, 1, 0.1))

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
        
      } else if (what == "heat_mx") {
        
#         if(!val_titles){
#           mainTitle <- "Allele and locus dropout"
#           xTitle <- "Mixture proportion (Mx)"
#           yTitle <- "Marker"
#         }
#         
#         .gData <- .gData[order(.gData$Sample.Name),]
#         
#         .gData$Dropout <- factor(.gData$Dropout)
#         
        # Mx Data:
        # .gData <- newdata[order(newdata$Ratio),]
        # .gData <- newdata[order(newdata$Proportion),]
      
        # Mx data:
        # .gData$Sample.Name<-paste(.gData$Ratio, " (", .gData$Sample.Name, ")", sep="")
        # .gData$Sample.Name<-paste(.gData$Proportion, " (", .gData$Sample.Name, ")", sep="")
      
        # Mx data:
        # .gData <- .gData [order(.gData$Ratio),]
        # .gData <- .gData [order(.gData$Proportion),]
      
        # Mx data SGM Plus.
        # .gData<-addColor()
        # .gData<-sortMarker(.gData,"SGM Plus")
      
        # Mx Data:
        # xlabels<-.gData[!duplicated(.gData[, c("Sample.Name", "Ratio")]), ]$Ratio
        # xlabels<-.gData[!duplicated(.gData[, c("Sample.Name", "Proportion")]), ]$Proportion
      
        # Mx data:
        # hm.title <- "Heatmap: allele and locus dropout for 'F' SGM Plus (3500)"
        # hm.xlab <- "Proportion"
        # Mx data:
      
        #gp <- gp + scale_y_discrete(limits = rev(levels(.gData$Marker))) + 
        #  scale_x_discrete(labels=formatC(xlabels, 4, format = "f")) +
        #  theme(axis.text.x=element_text(angle=-90, hjust = 0, vjust = 0.4, size = 10))

      }

      # Draw plot.    
      print(gp)
      
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
      if(exists(".strvalidator_plotDropout_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_plotDropout_gui_savegui", envir=env)
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
      if(exists(".strvalidator_plotDropout_gui_title", envir=env, inherits = FALSE)){
        svalue(f1_title_edt) <- get(".strvalidator_plotDropout_gui_title", envir=env)
      }
      if(exists(".strvalidator_plotDropout_gui_title_chk", envir=env, inherits = FALSE)){
        svalue(f1_titles_chk) <- get(".strvalidator_plotDropout_gui_title_chk", envir=env)
      }
      if(exists(".strvalidator_plotDropout_gui_x_title", envir=env, inherits = FALSE)){
        svalue(f1_xtitle_edt) <- get(".strvalidator_plotDropout_gui_x_title", envir=env)
      }
      if(exists(".strvalidator_plotDropout_gui_y_title", envir=env, inherits = FALSE)){
        svalue(f1_ytitle_edt) <- get(".strvalidator_plotDropout_gui_y_title", envir=env)
      }
#       if(exists(".strvalidator_plotDropout_gui_points_shape", envir=env, inherits = FALSE)){
#         svalue(shape_txt) <- get(".strvalidator_plotDropout_gui_points_shape", envir=env)
#       }
#       if(exists(".strvalidator_plotDropout_gui_points_alpha", envir=env, inherits = FALSE)){
#         svalue(alpha_txt) <- get(".strvalidator_plotDropout_gui_points_alpha", envir=env)
#       }
#       if(exists(".strvalidator_plotDropout_gui_points_jitterh", envir=env, inherits = FALSE)){
#         svalue(jitterh_txt) <- get(".strvalidator_plotDropout_gui_points_jitterh", envir=env)
#       }
#       if(exists(".strvalidator_plotDropout_gui_points_jitterv", envir=env, inherits = FALSE)){
#         svalue(jitterv_txt) <- get(".strvalidator_plotDropout_gui_points_jitterv", envir=env)
#       }
      if(exists(".strvalidator_plotDropout_gui_axes_y_min", envir=env, inherits = FALSE)){
        svalue(e3_y_min_edt) <- get(".strvalidator_plotDropout_gui_axes_y_min", envir=env)
      }
      if(exists(".strvalidator_plotDropout_gui_axes_y_max", envir=env, inherits = FALSE)){
        svalue(e3_y_max_edt) <- get(".strvalidator_plotDropout_gui_axes_y_max", envir=env)
      }
      if(exists(".strvalidator_plotDropout_gui_axes_x_min", envir=env, inherits = FALSE)){
        svalue(e3_x_min_edt) <- get(".strvalidator_plotDropout_gui_axes_x_min", envir=env)
      }
      if(exists(".strvalidator_plotDropout_gui_axes_x_max", envir=env, inherits = FALSE)){
        svalue(e3_x_max_edt) <- get(".strvalidator_plotDropout_gui_axes_x_max", envir=env)
      }
      if(exists(".strvalidator_plotDropout_gui_xlabel_size", envir=env, inherits = FALSE)){
        svalue(e4_size_txt) <- get(".strvalidator_plotDropout_gui_xlabel_size", envir=env)
      }
      if(exists(".strvalidator_plotDropout_gui_xlabel_angle", envir=env, inherits = FALSE)){
        svalue(e4_angle_spb) <- get(".strvalidator_plotDropout_gui_xlabel_angle", envir=env)
      }
      if(exists(".strvalidator_plotDropout_gui_xlabel_justh", envir=env, inherits = FALSE)){
        svalue(e4_hjust_spb) <- get(".strvalidator_plotDropout_gui_xlabel_justh", envir=env)
      }
      if(exists(".strvalidator_plotDropout_gui_xlabel_justv", envir=env, inherits = FALSE)){
        svalue(e4_vjust_spb) <- get(".strvalidator_plotDropout_gui_xlabel_justv", envir=env)
      }
      if(exists(".strvalidator_plotDropout_gui_hom", envir=env, inherits = FALSE)){
        svalue(f8_hom_chk) <- get(".strvalidator_plotDropout_gui_hom", envir=env)
      }
      
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_plotDropout_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_plotDropout_gui_title", value=svalue(f1_title_edt), envir=env)
      assign(x=".strvalidator_plotDropout_gui_title_chk", value=svalue(f1_titles_chk), envir=env)
      assign(x=".strvalidator_plotDropout_gui_x_title", value=svalue(f1_xtitle_edt), envir=env)
      assign(x=".strvalidator_plotDropout_gui_y_title", value=svalue(f1_ytitle_edt), envir=env)
#       assign(x=".strvalidator_plotDropout_gui_points_plot", value=svalue(e2_plotpoints_chk), envir=env)
#       assign(x=".strvalidator_plotDropout_gui_points_shape", value=svalue(shape_txt), envir=env)
#       assign(x=".strvalidator_plotDropout_gui_points_alpha", value=svalue(alpha_txt), envir=env)
#       assign(x=".strvalidator_plotDropout_gui_points_jitterh", value=svalue(jitterh_txt), envir=env)
#       assign(x=".strvalidator_plotDropout_gui_points_jitterv", value=svalue(jitterv_txt), envir=env)
      assign(x=".strvalidator_plotDropout_gui_axes_y_min", value=svalue(e3_y_min_edt), envir=env)
      assign(x=".strvalidator_plotDropout_gui_axes_y_max", value=svalue(e3_y_max_edt), envir=env)
      assign(x=".strvalidator_plotDropout_gui_axes_x_min", value=svalue(e3_x_min_edt), envir=env)
      assign(x=".strvalidator_plotDropout_gui_axes_x_max", value=svalue(e3_x_max_edt), envir=env)
      assign(x=".strvalidator_plotDropout_gui_xlabel_size", value=svalue(e4_size_txt), envir=env)
      assign(x=".strvalidator_plotDropout_gui_xlabel_angle", value=svalue(e4_angle_spb), envir=env)
      assign(x=".strvalidator_plotDropout_gui_xlabel_justh", value=svalue(e4_hjust_spb), envir=env)
      assign(x=".strvalidator_plotDropout_gui_xlabel_justv", value=svalue(e4_vjust_spb), envir=env)
      assign(x=".strvalidator_plotDropout_gui_hom", value=svalue(f8_hom_chk), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_plotDropout_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotDropout_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_plotDropout_gui_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotDropout_gui_title", envir = env)
      }
      if(exists(".strvalidator_plotDropout_gui_title_chk", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotDropout_gui_title_chk", envir = env)
      }
      if(exists(".strvalidator_plotDropout_gui_x_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotDropout_gui_x_title", envir = env)
      }
      if(exists(".strvalidator_plotDropout_gui_y_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotDropout_gui_y_title", envir = env)
      }
#       if(exists(".strvalidator_plotDropout_gui_points_plot", envir=env, inherits = FALSE)){
#         remove(".strvalidator_plotDropout_gui_points_plot", envir = env)
#       }
#       if(exists(".strvalidator_plotDropout_gui_points_shape", envir=env, inherits = FALSE)){
#         remove(".strvalidator_plotDropout_gui_points_shape", envir = env)
#       }
#       if(exists(".strvalidator_plotDropout_gui_points_alpha", envir=env, inherits = FALSE)){
#         remove(".strvalidator_plotDropout_gui_points_alpha", envir = env)
#       }
#       if(exists(".strvalidator_plotDropout_gui_points_jitterh", envir=env, inherits = FALSE)){
#         remove(".strvalidator_plotDropout_gui_points_jitterh", envir = env)
#       }
#       if(exists(".strvalidator_plotDropout_gui_points_jitterv", envir=env, inherits = FALSE)){
#         remove(".strvalidator_plotDropout_gui_points_jitterv", envir = env)
#       }
      if(exists(".strvalidator_plotDropout_gui_axes_y_min", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotDropout_gui_axes_y_min", envir = env)
      }
      if(exists(".strvalidator_plotDropout_gui_axes_y_max", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotDropout_gui_axes_y_max", envir = env)
      }
      if(exists(".strvalidator_plotDropout_gui_axes_x_min", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotDropout_gui_axes_x_min", envir = env)
      }
      if(exists(".strvalidator_plotDropout_gui_axes_x_max", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotDropout_gui_axes_x_max", envir = env)
      }
      if(exists(".strvalidator_plotDropout_gui_xlabel_size", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotDropout_gui_xlabel_size", envir = env)
      }
      if(exists(".strvalidator_plotDropout_gui_xlabel_angle", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotDropout_gui_xlabel_angle", envir = env)
      }
      if(exists(".strvalidator_plotDropout_gui_xlabel_justh", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotDropout_gui_xlabel_justh", envir = env)
      }
      if(exists(".strvalidator_plotDropout_gui_xlabel_justv", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotDropout_gui_xlabel_justv", envir = env)
      }
      if(exists(".strvalidator_plotDropout_gui_hom", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotDropout_gui_hom", envir = env)
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
