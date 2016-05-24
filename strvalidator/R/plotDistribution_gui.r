################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 06.01.2016: Fixed theme methods not found and added more themes.
# 11.11.2015: Added importFrom ggplot2.
# 11.11.2015: Added more themes.
# 07.10.2015: NA's now removed prior to plotting, and from number of observations.
# 29.08.2015: Added importFrom.
# 11.06.2015: Fixed title for histogram plot.
# 09.06.2015: Fixed 'overlay boxplot' not saved.
# 04.05.2015: Added 'Histogram'.
# 11.10.2014: Added 'focus', added 'parent' parameter.
# 28.06.2014: Added help button and moved save gui checkbox.
# 11.05.2014: Fixed boxplot bug, box not drawn.
# 08.05.2014: Implemented 'checkDataset'.
# 23.02.2014: No column required.
# 23.02.2014: Conversion of 'Height', 'Size', and 'Data.Point' to numeric.
# 23.02.2014: Fixed boxplot width bug.
# 09.02.2014: Automatically try to guess x title from chosen column.
# 07.02.2014: First version.

#' @title Plot Distribution
#'
#' @description
#' GUI simplifying the creation of distribution plots.
#'
#' @details Plot the distribution of data as cumulative distribution function,
#' probability density function, or count. First select a dataset, then select
#' a group (if any), finally select a column to plot the distribution of.
#' It is possible to overlay a boxplot.
#' Various smoothing kernels and bandwidths can be specified.
#' More info on kernels and bandwidth: \url{http://www.inside-r.org/r-doc/stats/density}
#' Automatic plot titles can be replaced by custom titles.
#' A name for the result is automatiaclly suggested.
#' The resulting plot can be saved as either a plot object or as an image.
#' @param env environment in wich to search for data frames and save result.
#' @param savegui logical indicating if GUI settings should be saved in the environment.
#' @param debug logical indicating printing debug information.
#' @param parent widget to get focus when finished.
#' 
#' @export
#' 
#' @importFrom utils help str head
#' @importFrom ggplot2 qplot ggplot aes_string stat_ecdf geom_density ggplot_build
#'  geom_boxplot geom_segment geom_point labs theme_gray theme_bw
#'  theme_linedraw theme_light theme_dark theme_minimal theme_classic
#'  theme_void 
#' 
#' @return TRUE

plotDistribution_gui <- function(env=parent.frame(), savegui=NULL, debug=FALSE, parent=NULL){
  
  # Global variables.
  .gData <- NULL
  .gDataName <- NULL
  .gPlot <- NULL
  .palette <- c("Set1","Set2","Set3","Accent","Dark2",
                "Paired","Pastel1", "Pastel2")
  # Qualitative palette, do not imply magnitude differences between legend
  # classes, and hues are used to create the primary visual differences 
  # between classes. Qualitative schemes are best suited to representing
  # nominal or categorical data.
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # Main window.
  w <- gwindow(title="Plot distributions", visible=FALSE)
  
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
    print(help("plotDistribution_gui", help_type="html"))
    
  })
  
  # FRAME 0 ###################################################################
  
  f0 <- gframe(text = "Dataset",
               horizontal=TRUE,
               spacing = 5,
               container = gv) 
  
  f0g0 <- glayout(container = f0)
  
  f0g0[1,1] <- glabel(text="Select dataset:", container=f0g0)
  
  f0g0[1,2] <- dataset_drp <- gdroplist(items=c("<Select dataset>",
                                                listObjects(env=env,
                                                            obj.class="data.frame")),
                                        selected = 1,
                                        editable = FALSE,
                                        container = f0g0) 
  
  f0g0[1,3] <- f0_samples_lbl <- glabel(text=" (0 rows)", container=f0g0)
  
  f0g0[2,1] <- glabel(text="Select group:", container=f0g0)
  f0g0[2,2] <- f0_group_drp <- gdroplist(items="<Select group>",
                                         selected = 1, container=f0g0)
  f0g0[2,3] <- f0_rows_lbl <- glabel(text=" (0 rows)", container=f0g0)
  
  f0g0[3,1] <- glabel(text="Select column:", container=f0g0)
  f0g0[3,2] <- f0_column_drp <- gdroplist(items="<Select column>",
                                          selected = 1, container=f0g0)
  
  
  addHandlerChanged(dataset_drp, handler = function (h, ...) {
    
    val_obj <- svalue(dataset_drp)
    
    # Check if suitable.
    requiredCol <- NULL
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       env=env, parent=w, debug=debug)
    
    if(ok){
      
      # Load or change components.
      .gData <<- get(val_obj, envir=env)
      .gDataName <<- val_obj
      
      # Refresh column in drop lists.
      .refresh_column_drp()
      
      # Suggest name.
      svalue(f5_save_edt) <- paste(val_obj, "_ggplot", sep="")
      
      # Get number of observations.
      svalue(f0_samples_lbl) <- paste(" (", nrow(.gData), " rows)", sep="")
      
      # Get number of observations in subset.
      val <- svalue(f0_group_drp)
      if(val %in% names(.gData)){
        rows <- nrow(.gData[.gData$Group==val, ])
        svalue(f0_rows_lbl) <- paste(" (", rows, " rows)", sep="")
      }
      
      # Enable buttons.
      enabled(f7_ecdf_btn) <- TRUE
      enabled(f7_pdf_btn) <- TRUE
      enabled(f7_qplot_btn) <- TRUE
      
    } else {
      
      # Reset components.
      .gData <<- NULL
      svalue(f5_save_edt) <- ""
      svalue(f0_samples_lbl) <- " (0 rows)"
      
    }    
    
  } )  
  
  addHandlerChanged(f0_group_drp, handler = function (h, ...) {
    
    val <- svalue(f0_group_drp)
    rows <- nrow(.gData[.gData$Group==val, ])
    
    # Update number of observations.
    svalue(f0_rows_lbl) <- paste(" (", rows, " rows)", sep="")
    
  } )  
  
  addHandlerChanged(f0_column_drp, handler = function (h, ...) {
    
    # Enable buttons.
    enabled(f7_ecdf_btn) <- TRUE
    enabled(f7_pdf_btn) <- TRUE
    enabled(f7_qplot_btn) <- TRUE
    
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
  f1g1[1,2] <- f1_title_edt <- gedit(text="",
                                     width=40,
                                     container=f1g1)
  
  f1g1[2,1] <- glabel(text="X title:", container=f1g1)
  f1g1[2,2] <- f1_xtitle_edt <- gedit(text="",
                                      container=f1g1)
  
  f1g1[3,1] <- glabel(text="Y title:", container=f1g1)
  f1g1[3,2] <- f1_ytitle_edt <- gedit(text="",
                                      container=f1g1)
  
  f1g2 <- glayout(container = f1, spacing = 1)
  f1g2[1,1] <- glabel(text="Plot theme:", anchor=c(-1 ,0), container=f1g2)
  items_theme <- c("theme_grey()","theme_bw()","theme_linedraw()",
                   "theme_light()","theme_dark()","theme_minimal()",
                   "theme_classic()","theme_void()")
  f1g2[1,2] <- f1_theme_drp <- gdroplist(items = items_theme,
                                         selected = 1,
                                         container = f1g2)
  
  f1e1 <- gexpandgroup(text = "Boxplot", horizontal=FALSE,
                       spacing = 5, container = f1) 
  
  f1_box_chk <- gcheckbox(text="Overlay boxplot",
                          checked=TRUE,
                          container=f1e1)
  
  f1g3 <- glayout(container = f1e1, spacing = 1)
  f1g3[1,1] <- glabel(text="Width of boxplot:", container=f1g3)
  f1g3[1,2] <- f1_width_spn <- gspinbutton(from=0, to=1, by=0.01, value=0.25,
                                           container=f1g3)
  
  f1e2 <- gexpandgroup(text = "Distribution function", horizontal=FALSE,
                       spacing = 5, container = f1) 
  
  f1g4 <- glayout(container = f1e2, spacing = 1)
  
  f1_kernel <- c("gaussian", "rectangular", "triangular", "epanechnikov", "biweight", "cosine","optcosine") 
  f1g4[1,1] <- glabel(text="Smoothing kernel:", container=f1g4)
  f1g4[1,2] <- f1_kernel_drp <- gdroplist(items=f1_kernel,
                                          selected = 1, container=f1g4)
  
  f1_adjust <- c(4,2,1,0.5,0.25) 
  f1g4[2,1] <- glabel(text="Adjust bandwidth:", container=f1g4)
  f1g4[2,2] <- f1_adjustbw_cbo <- gcombobox(items=f1_adjust,
                                            selected = 3, editable = TRUE,
                                            container=f1g4)
  
  f1e3 <- gexpandgroup(text = "Histogram", horizontal=FALSE,
                       spacing = 5, container = f1) 
  
  f1g5 <- glayout(container = f1e3, spacing = 1)
  
  f1g5[1,1] <- glabel(text="Adjust binwidth:", container=f1g5)
  f1g5[1,2] <- f1_binwidth_edt <- gedit(text="1", container=f1g5)
  
  
  # FRAME 7 ###################################################################
  
  f7 <- gframe(text = "Plot distribution",
               horizontal=FALSE,
               container = gv) 
  
  f7g7 <- glayout(container = f7)
  
  f7g7[1,1] <- f7_ecdf_btn <- gbutton(text="CDF", border=TRUE, container=f7g7) 
  
  addHandlerChanged(f7_ecdf_btn, handler = function(h, ...) {
    
    val_column <- svalue(f0_column_drp)
    
    if(val_column == "<Select column>"){
      
      gmessage(message="A data column must be specified!",
               title="Error",
               icon = "error")      
      
    } else {
      
      enabled(f7_ecdf_btn) <- FALSE
      .plot(how="cdf")
      enabled(f7_ecdf_btn) <- TRUE
      
    }
    
  } )
  
  f7g7[1,2] <- f7_pdf_btn <- gbutton(text="PDF", border=TRUE, container=f7g7) 
  
  addHandlerChanged(f7_pdf_btn, handler = function(h, ...) {
    
    val_column <- svalue(f0_column_drp)
    
    if(val_column == "<Select column>"){
      
      gmessage(message="A data column must be specified!",
               title="Error",
               icon = "error")      
      
    } else {
      
      enabled(f7_pdf_btn) <- FALSE
      .plot(how="pdf")
      enabled(f7_pdf_btn) <- TRUE
      
    }
    
  } )
  
  f7g7[1,3] <- f7_qplot_btn <- gbutton(text="Histogram", border=TRUE, container=f7g7) 
  
  addHandlerChanged(f7_qplot_btn, handler = function(h, ...) {
    
    val_column <- svalue(f0_column_drp)
    
    if(val_column == "<Select column>"){
      
      gmessage(message="A data column must be specified!",
               title="Error",
               icon = "error")      
      
    } else {
      
      enabled(f7_qplot_btn) <- FALSE
      .plot(how="qplot")
      enabled(f7_qplot_btn) <- TRUE
      
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
  
  # FUNCTIONS #################################################################
  
  .plot <- function(how){
    
    # Get values.
    val_data <- .gData
    val_titles <- svalue(f1_titles_chk)
    val_title <- svalue(f1_title_edt)
    val_x_title <- svalue(f1_xtitle_edt)
    val_y_title <- svalue(f1_ytitle_edt)
    val_theme <- svalue(f1_theme_drp)
    val_group <- svalue(f0_group_drp)
    val_column <- svalue(f0_column_drp)
    val_kernel <- svalue(f1_kernel_drp)
    val_adjustbw <- as.numeric(svalue(f1_adjustbw_cbo))
    val_boxplot <- svalue(f1_box_chk)
    val_width <- svalue(f1_width_spn)
    val_binwidth <- as.numeric(svalue(f1_binwidth_edt))
    
    if(debug){
      print("val_titles")
      print(val_titles)
      print("val_title")
      print(val_title)
      print("val_x_title")
      print(val_x_title)
      print("val_y_title")
      print(val_y_title)
      print("val_kernel")
      print(val_kernel)
      print("val_column")
      print(val_column)
      print("str(val_data)")
      print(str(val_data))
    }
    
    # Check if data.
    if (!is.na(val_data) && !is.null(val_data)){
      
      if(debug){
        print("Before plot: str(val_data)")
        print(str(val_data))
        print(head(val_data))
      }
      
      # Get number of observations.
      nb <- nrow(val_data)
      
      # Get data for selected group.
      if("Group" %in% names(val_data)){
        
        # Store nb of observations.
        nb0 <- nb
        
        # Subset according to group.
        val_data <- val_data[val_data$Group==val_group, ]
        
        # Update number of observations.
        nb <- nrow(val_data)
        
        # Show message.
        message(paste("Subset group = '", val_group,
                      "'', removed ", nb0-nb, " rows.", sep=""))
        
      }
      
      # Convert to numeric.
      if(!is.numeric(val_data[, val_column])){
        val_data[, val_column] <- as.numeric(val_data[, val_column])
        message(paste(val_column, " converted to numeric."))
      }
      
      if(debug){
        print("After subsetting (val_data)")
        print(str(val_data))
        print(head(val_data))
      }
      
      # Remove NA's
      if(any(is.na(val_data[, val_column]))){
        
        # Store nb of observations.
        nb0 <- nb
        
        # Update number of observations.
        nb <- nrow(val_data[!is.na(val_data[val_column]),])
        
        # Show message.
        message(paste("Removed ", nb0-nb, " NA rows.", sep=""))
        
        if(debug){
          print("After subsetting (val_data)")
          print(str(val_data))
          print(head(val_data))
        }
        
      }

      # Create titles.
      if(val_titles){
        
        if(debug){
          print("Custom titles")
        }
        
        mainTitle <- val_title
        xTitle <- val_x_title
        yTitle <- val_y_title
        
      } else {
        
        if(debug){
          print("Default titles")
        }
        
        # Different titles.
        if(how == "cdf"){
          
          mainTitle <- paste("Cumulative density function (",
                             nb, " observations)", sep="")
          
          yTitle <- "Proportion"
          
        } else if(how=="pdf"){
          
          mainTitle <- paste("Probability density function (",
                             nb, " observations)", sep="")
          
          yTitle <- "Proportion"
          
        } else if(how=="qplot"){
          
          mainTitle <- paste("Histogram (",
                             nb, " observations)", sep="")
          
          yTitle <- "Count"
          
        } else {
          
          warning(paste("how=", how, "not implemented for titles!"))
          
        }
        
        # Different X axis depending on chosen column.
        if(val_column == "Height"){
          
          xTitle <- "Peak height (RFU)"
          
        } else if(val_column == "Size"){
          
          xTitle <- "Fragment size (bp)"
          
        } else if(val_column == "Data.Point"){
          
          xTitle <- "Data point"
          
        } else {
          
          xTitle <- val_column
          
        }
        
      }
      
      # Create plots.
      if(how == "cdf"){
        
        if(debug){
          print("Create cdf plot")
        }
        
        # ECDP
        gp <- ggplot(data=val_data, aes_string(x=val_column))
        gp <- gp + stat_ecdf()
        
      } else if(how == "pdf"){
        
        if(debug){
          print("Create pdf plot")
        }
        
        gp <- ggplot(data=val_data, aes_string(x=val_column))
        # More info on kernels and bandwidth: http://www.inside-r.org/r-doc/stats/density
        gp <- gp + geom_density(aes_string(x=val_column), kernel=val_kernel, adjust=val_adjustbw)
        
      } else if(how == "qplot"){
        
        # Get data as vector.
        dv <- val_data[, val_column]
        
        if(debug){
          print("Create Histogram")
          str(head(dv))
        }
        
        # Create plot.
        gp <- qplot(x=dv, binwidth=val_binwidth)
        
      } else {
        
        warning(paste("how=", how, "not implemented for plots!"))
        
      }
      
      if(debug){
        print("Plot created")
      }
      
      # Overlay boxplot.
      if(val_boxplot){
        
        if(debug){
          print("Overlay boxplot")
        }
        
        # Extract information from plot:
        gb <- ggplot_build(gp)
        ywidth <- max(gb$data[[1]]$y, na.rm=TRUE) * (val_width / 2)
        ymean <- max(gb$data[[1]]$y, na.rm=TRUE) / 2
        
        # Create a normal boxplot.
        gbox <- ggplot(data=val_data, aes_string(x=1, y=val_column))
        gbox <- gbox + geom_boxplot()
        
        # Extract information from boxplot.
        gb <- ggplot_build(gbox)
        xmax <- gb$data[[1]]$ymax
        xmin <- gb$data[[1]]$ymin
        left <- gb$data[[1]]$lower
        middle  <- gb$data[[1]]$middle
        right <- gb$data[[1]]$upper
        dots <- unlist(gb$data[[1]]$outliers)
        
        val_box <- data.frame(xmin=xmin, xmax = xmax,
                              ymin=ymean-ywidth, ymax=ymean+ywidth, ymean=ymean,
                              left=left, middle=middle, right=right)
        
        
        if(debug){
          print("val_box")
          print(val_box)
          print("dots")
          print(dots)
        }
        
        # Manually overlay a boxplot:
        # Add box.
        # Should work...
        #        gp <- gp + geom_polygon(data=val_box, aes_string(x = c("left","left","right","right"),
        #                                                         y = c("ymin","ymax","ymax","ymin")),
        #                                color=1, alpha=0)
        #        gp <- gp + geom_rect(data=val_box, aes_string(xmin = "left", xmax="right",
        #                                                      ymin = "ymin", ymax="ymax"),
        #                             color=1, alpha=0)
        # Add top.
        gp <- gp + geom_segment(data=val_box, aes_string(x="left", y="ymax",
                                                         xend="right", yend="ymax"))
        
        # Add bottom.
        gp <- gp + geom_segment(data=val_box, aes_string(x="left", y="ymin",
                                                         xend="right", yend="ymin"))
        
        # Add left.
        gp <- gp + geom_segment(data=val_box, aes_string(x="left", y="ymin",
                                                         xend="left", yend="ymax"))
        
        # Add right.
        gp <- gp + geom_segment(data=val_box, aes_string(x="right", y="ymin",
                                                         xend="right", yend="ymax"))
        
        # Add median.
        gp <- gp + geom_segment(data=val_box, aes_string(x="middle", y="ymin",
                                                         xend="middle", yend="ymax"))
        # Add whiskers.
        gp <- gp + geom_segment(data=val_box, aes_string(x="xmin", y="ymean",
                                                         xend="left", yend="ymean"))
        gp <- gp + geom_segment(data=val_box, aes_string(x="xmax", y="ymean",
                                                         xend="right", yend="ymean"))
        # Add outliers.
        out <- data.frame(x = dots, y = rep(ymean, length(dots)))
        gp <- gp + geom_point(data=out, aes_string(x = "x", y = "y"))
        
        if(debug){
          print("Boxplot created")
        }
        
      } # End if boxplot.
      
      # Add titles.
      gp <- gp + labs(title=mainTitle, x=xTitle, y=yTitle, fill=NULL)
      
      # Apply theme.
      gp <- gp + eval(parse(text=val_theme))
      
      # plot.
      print(gp)
      
      # Store in global variable.
      .gPlot <<- gp
      
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
  
  .refresh_column_drp <- function(){
    
    if(debug){
      print("Refresh group and column dropdown")
    }
    
    # Get data frames in global workspace.
    groups <- unique(.gData$Group)
    columns <- names(.gData)
    
    if(!is.null(groups)){
      
      blockHandler(f0_group_drp)
      
      # Populate drop list.
      f0_group_drp[] <- c("<Select group>", groups)
      
      unblockHandler(f0_group_drp)
      
      if(debug){
        print("Group dropdown refreshed!")
      }
      
    }
    
    if(!is.null(columns)){
      
      blockHandler(f0_column_drp)
      
      # Populate drop list.
      f0_column_drp[] <- c("<Select column>", columns)
      
      unblockHandler(f0_column_drp)
      
      if(debug){
        print("Column dropdown refreshed!")
      }
      
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
      if(exists(".strvalidator_plotDistribution_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_plotDistribution_gui_savegui", envir=env)
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
      if(exists(".strvalidator_plotDistribution_gui_title", envir=env, inherits = FALSE)){
        svalue(f1_title_edt) <- get(".strvalidator_plotDistribution_gui_title", envir=env)
      }
      if(exists(".strvalidator_plotDistribution_gui_title_chk", envir=env, inherits = FALSE)){
        svalue(f1_titles_chk) <- get(".strvalidator_plotDistribution_gui_title_chk", envir=env)
      }
      if(exists(".strvalidator_plotDistribution_gui_x_title", envir=env, inherits = FALSE)){
        svalue(f1_xtitle_edt) <- get(".strvalidator_plotDistribution_gui_x_title", envir=env)
      }
      if(exists(".strvalidator_plotDistribution_gui_y_title", envir=env, inherits = FALSE)){
        svalue(f1_ytitle_edt) <- get(".strvalidator_plotDistribution_gui_y_title", envir=env)
      }
      if(exists(".strvalidator_plotDistribution_gui_box", envir=env, inherits = FALSE)){
        svalue(f1_box_chk) <- get(".strvalidator_plotDistribution_gui_box", envir=env)
      }
      if(exists(".strvalidator_plotDistribution_gui_kernel", envir=env, inherits = FALSE)){
        svalue(f1_kernel_drp) <- get(".strvalidator_plotDistribution_gui_kernel", envir=env)
      }
      if(exists(".strvalidator_plotDistribution_gui_theme", envir=env, inherits = FALSE)){
        svalue(f1_theme_drp) <- get(".strvalidator_plotDistribution_gui_theme", envir=env)
      }
      if(exists(".strvalidator_plotDistribution_gui_width", envir=env, inherits = FALSE)){
        svalue(f1_width_spn) <- get(".strvalidator_plotDistribution_gui_width", envir=env)
      }
      if(exists(".strvalidator_plotDistribution_gui_binwidth", envir=env, inherits = FALSE)){
        svalue(f1_binwidth_edt) <- get(".strvalidator_plotDistribution_gui_binwidth", envir=env)
      }
      
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_plotDistribution_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_plotDistribution_gui_title_chk", value=svalue(f1_titles_chk), envir=env)
      assign(x=".strvalidator_plotDistribution_gui_title", value=svalue(f1_title_edt), envir=env)
      assign(x=".strvalidator_plotDistribution_gui_x_title", value=svalue(f1_xtitle_edt), envir=env)
      assign(x=".strvalidator_plotDistribution_gui_y_title", value=svalue(f1_ytitle_edt), envir=env)
      assign(x=".strvalidator_plotDistribution_gui_box", value=svalue(f1_box_chk), envir=env)
      assign(x=".strvalidator_plotDistribution_gui_kernel", value=svalue(f1_kernel_drp), envir=env)
      assign(x=".strvalidator_plotDistribution_gui_theme", value=svalue(f1_theme_drp), envir=env)
      assign(x=".strvalidator_plotDistribution_gui_width", value=svalue(f1_width_spn), envir=env)
      assign(x=".strvalidator_plotDistribution_gui_binwidth", value=svalue(f1_binwidth_edt), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_plotDistribution_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotDistribution_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_plotDistribution_gui_title_chk", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotDistribution_gui_title_chk", envir = env)
      }
      if(exists(".strvalidator_plotDistribution_gui_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotDistribution_gui_title", envir = env)
      }
      if(exists(".strvalidator_plotDistribution_gui_x_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotDistribution_gui_x_title", envir = env)
      }
      if(exists(".strvalidator_plotDistribution_gui_y_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotDistribution_gui_y_title", envir = env)
      }
      if(exists(".strvalidator_plotDistribution_gui_box", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotDistribution_gui_box", envir = env)
      }
      if(exists(".strvalidator_plotDistribution_gui_kernel", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotDistribution_gui_kernel", envir = env)
      }
      if(exists(".strvalidator_plotDistribution_gui_theme", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotDistribution_gui_theme", envir = env)
      }
      if(exists(".strvalidator_plotDistribution_gui_width", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotDistribution_gui_width", envir = env)
      }
      if(exists(".strvalidator_plotDistribution_gui_binwidth", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotDistribution_gui_binwidth", envir = env)
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
