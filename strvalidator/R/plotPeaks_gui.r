################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (last 20 changes)
# 11.11.2015: Added importFrom ggplot2.
# 29.08.2015: Added importFrom.
# 11.10.2014: Added 'focus', added 'parent' parameter.
# 28.06.2014: Added help button and moved save gui checkbox.
# 08.05.2014: Implemented 'checkDataset'.
# 20.01.2014: Changed 'saveImage_gui' for 'ggsave_gui'.
# 12.01.2014: First version.

#' @title Plot Peaks
#'
#' @description
#' GUI simplifying the creation of plots from result type data.
#'
#' @details Plot result type data. It is possible to customise titles and font
#' size. Data can be plotted as as frequency or proportion. The values can be
#' printed on the plot with custom number of decimals. There are several 
#' colour palettes to chose from. 
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
#' @importFrom plyr count
#' @importFrom utils help str
#' @importFrom ggplot2 ggplot aes_string theme_grey geom_bar scale_fill_brewer
#'  labs geom_text theme
#' 

plotPeaks_gui <- function(env=parent.frame(), savegui=NULL, debug=FALSE, parent=NULL){

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
  w <- gwindow(title="Plot peaks", visible=FALSE)
  
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
    print(help("plotPeaks_gui", help_type="html"))
    
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
    requiredCol <- c("Sample.Name", "Peaks", "Group", "Id")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       env=env, parent=w, debug=debug)
    
    if(ok){
      
      # Load or change components.
      .gData <<- get(val_obj, envir=env)
      .gDataName <<- val_obj

      # Suggest name.
      svalue(f5_save_edt) <- paste(val_obj, "_ggplot", sep="")
      
      svalue(f0_samples_lbl) <- paste(" (",
                                      length(unique(.gData$Id)),
                                      " samples)", sep="")

      # Enable buttons.
      enabled(plot_btn) <- TRUE
        
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
               spacing = 10,
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

  f1_prop_chk <- gcheckbox(text="Plot proportion",
                              checked=TRUE,
                              container=f1)
  
  grid2 <- glayout(container = f1, spacing = 1)
  grid2[1,1] <- glabel(text="Base font size (pts):", container=grid2)
  grid2[1,2] <- f1_base_size_edt <- gedit(text="18", width=4, container=grid2)
  
  grid3 <- glayout(container = f1, spacing = 1)
  grid3[1,1] <- glabel(text="Colour palette:", container=grid3)
  grid3[1,2] <- f1_palette_drp <- gdroplist(items=.palette,
                                            selected = 1,
                                            editable = FALSE,
                                            container = grid3)

  grid4 <- glayout(container = f1, spacing = 1)
  grid4[1,1] <- f1_print_chk <- gcheckbox(text="Print values as bar labels", checked=TRUE, container=grid4)
  grid4[2,1] <- glabel(text="Number of decimals for bar labels:", container=grid4)
  grid4[2,2] <- f1_decimal_spb <- gspinbutton(from=0, to=9, by=1, value=4,
                                              container=grid4)
  grid4[3,1] <- glabel(text="Font size for bar labels (pts):", container=grid4)
  grid4[3,2] <- f1_lab_size_edt <- gedit(text="5", width=4, container=grid4)
  
  # FRAME 7 ###################################################################
  
  plot_btn <- gbutton(text="Plot", border=TRUE, container=gv) 
  
  addHandlerChanged(plot_btn, handler = function(h, ...) {
    
      enabled(plot_btn) <- FALSE
      .plotBalance()
      enabled(plot_btn) <- TRUE

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
  
  .plotBalance <- function(){
    
    # Get values.
    val_titles <- svalue(f1_titles_chk)
    val_title <- svalue(f1_title_edt)
    val_x_title <- svalue(f1_xtitle_edt)
    val_y_title <- svalue(f1_ytitle_edt)
    val_base_size <- as.numeric(svalue(f1_base_size_edt))
    val_lab_size <- as.numeric(svalue(f1_lab_size_edt))
    val_palette <- svalue(f1_palette_drp)
    val_decimals <- as.numeric(svalue(f1_decimal_spb))
    val_print <- svalue(f1_print_chk)
    val_prop <- svalue(f1_prop_chk)
    
    if(debug){
      print("val_titles")
      print(val_titles)
      print("val_title")
      print(val_title)
      print("val_x_title")
      print(val_x_title)
      print("val_y_title")
      print(val_y_title)
      print("val_base_size")
      print(val_base_size)
      print("val_lab_size")
      print(val_lab_size)
      print("val_palette")
      print(val_palette)
      print("val_decimals")
      print(val_decimals)
      print("val_print")
      print(val_print)
      print("val_prop")
      print(val_prop)
      print("str(.gData)")
      print(str(.gData))
    }

    # Check if data.
    if (!is.na(.gData) && !is.null(.gData)){

      if(debug){
        print("Before plot: str(.gData)")
        print(str(.gData))
      }

      # Prepare data.
      # Get one row from each sample for plotting.
      .gData <- .gData[!duplicated(.gData[,'Id']),]
      
      # Create titles.
      if(val_titles){
        mainTitle <- val_title
        xTitle <- val_x_title
        yTitle <- val_y_title
      } else {
        numberOfSamples <- nrow(.gData)
        mainTitle <- paste("Analysis of peaks from",
                           numberOfSamples, "samples")
        xTitle <- "Group"
        if(val_prop){
          yTitle <- "Proportion"
        } else {
          yTitle <- "Count"
        }
      }

      # Count samples per group.
      .gData <- plyr::count(.gData, vars="Group")
      #Calculate frequencies.
      if(val_prop){
        .gData$freq <- .gData$freq / sum(.gData$freq)
      }
      .gData$lab <- round(.gData$freq, val_decimals)
      
      # Create plot.
      gp <- ggplot(.gData, aes_string(x = "Group", y="freq", fill = "Group"))
      gp <- gp + theme_grey(base_size = val_base_size) 
      gp <- gp + geom_bar(stat = "identity", position = "stack")
      
      # Add color.
      gp <- gp + scale_fill_brewer(palette = val_palette) # NB! only 9 colors.

      # Add titles.
      gp <- gp + labs(title=mainTitle, x=xTitle, y=yTitle, fill=NULL)
      
      # Print value labels on bars.
      if(val_print){
        gp <- gp + geom_text(aes_string(x="Group", y="freq",
                                 ymax="freq", label="lab", 
                                 hjust=0.5, vjust=0), size=val_lab_size)
      }

      # Remove legend.
      gp <- gp + theme(legend.position="none")

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
      if(exists(".strvalidator_plotPeaks_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_plotPeaks_gui_savegui", envir=env)
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
      if(exists(".strvalidator_plotPeaks_gui_title", envir=env, inherits = FALSE)){
        svalue(f1_title_edt) <- get(".strvalidator_plotPeaks_gui_title", envir=env)
      }
      if(exists(".strvalidator_plotPeaks_gui_title_chk", envir=env, inherits = FALSE)){
        svalue(f1_titles_chk) <- get(".strvalidator_plotPeaks_gui_title_chk", envir=env)
      }
      if(exists(".strvalidator_plotPeaks_gui_x_title", envir=env, inherits = FALSE)){
        svalue(f1_xtitle_edt) <- get(".strvalidator_plotPeaks_gui_x_title", envir=env)
      }
      if(exists(".strvalidator_plotPeaks_gui_y_title", envir=env, inherits = FALSE)){
        svalue(f1_ytitle_edt) <- get(".strvalidator_plotPeaks_gui_y_title", envir=env)
      }
      if(exists(".strvalidator_plotPeaks_gui_base_size", envir=env, inherits = FALSE)){
        svalue(f1_base_size_edt) <- get(".strvalidator_plotPeaks_gui_base_size", envir=env)
      }
      if(exists(".strvalidator_plotPeaks_gui_label_size", envir=env, inherits = FALSE)){
        svalue(f1_lab_size_edt) <- get(".strvalidator_plotPeaks_gui_label_size", envir=env)
      }
      if(exists(".strvalidator_plotPeaks_gui_print", envir=env, inherits = FALSE)){
        svalue(f1_print_chk) <- get(".strvalidator_plotPeaks_gui_print", envir=env)
      }
      if(exists(".strvalidator_plotPeaks_gui_print", envir=env, inherits = FALSE)){
        svalue(f1_prop_chk) <- get(".strvalidator_plotPeaks_gui_prop", envir=env)
      }
      if(exists(".strvalidator_plotPeaks_gui_palette", envir=env, inherits = FALSE)){
        svalue(f1_palette_drp) <- get(".strvalidator_plotPeaks_gui_palette", envir=env)
      }
      if(exists(".strvalidator_plotPeaks_gui_decimal", envir=env, inherits = FALSE)){
        svalue(f1_decimal_spb) <- get(".strvalidator_plotPeaks_gui_decimal", envir=env)
      }
            
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_plotPeaks_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_plotPeaks_gui_title_chk", value=svalue(f1_titles_chk), envir=env)
      assign(x=".strvalidator_plotPeaks_gui_title", value=svalue(f1_title_edt), envir=env)
      assign(x=".strvalidator_plotPeaks_gui_x_title", value=svalue(f1_xtitle_edt), envir=env)
      assign(x=".strvalidator_plotPeaks_gui_y_title", value=svalue(f1_ytitle_edt), envir=env)
      assign(x=".strvalidator_plotPeaks_gui_base_size", value=svalue(f1_base_size_edt), envir=env)
      assign(x=".strvalidator_plotPeaks_gui_label_size", value=svalue(f1_lab_size_edt), envir=env)
      assign(x=".strvalidator_plotPeaks_gui_print", value=svalue(f1_print_chk), envir=env)
      assign(x=".strvalidator_plotPeaks_gui_prop", value=svalue(f1_prop_chk), envir=env)
      assign(x=".strvalidator_plotPeaks_gui_palette", value=svalue(f1_palette_drp), envir=env)
      assign(x=".strvalidator_plotPeaks_gui_decimal", value=svalue(f1_decimal_spb), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_plotPeaks_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPeaks_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_plotPeaks_gui_title_chk", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPeaks_gui_title_chk", envir = env)
      }
      if(exists(".strvalidator_plotPeaks_gui_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPeaks_gui_title", envir = env)
      }
      if(exists(".strvalidator_plotPeaks_gui_x_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPeaks_gui_x_title", envir = env)
      }
      if(exists(".strvalidator_plotPeaks_gui_y_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPeaks_gui_y_title", envir = env)
      }
      if(exists(".strvalidator_plotPeaks_gui_base_size", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPeaks_gui_base_size", envir = env)
      }
      if(exists(".strvalidator_plotPeaks_gui_label_size", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPeaks_gui_label_size", envir = env)
      }
      if(exists(".strvalidator_plotPeaks_gui_print", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPeaks_gui_print", envir = env)
      }
      if(exists(".strvalidator_plotPeaks_gui_prop", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPeaks_gui_prop", envir = env)
      }
      if(exists(".strvalidator_plotPeaks_gui_palette", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPeaks_gui_palette", envir = env)
      }
      if(exists(".strvalidator_plotPeaks_gui_decimal", envir=env, inherits = FALSE)){
        remove(".strvalidator_plotPeaks_gui_decimal", envir = env)
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
