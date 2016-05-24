################################################################################
# TODO LIST
# TODO: ...

# NOTE: Column names used for calculations with data.table is declared
# in globals.R to avoid NOTES in R CMD CHECK.

################################################################################
# CHANGE LOG (last 20 changes)
# 11.11.2015: Added importFrom ggplot2.
# 21.10.2015: Added attributes.
# 28.08.2015: Added importFrom
# 18.08.2015: Changed label for AT4 option.
# 04.05.2015: First version.

#' @title Calculate Analytical Threshold
#'
#' @description
#' GUI wrapper for the \code{\link{blockAT}} and \code{\link{calculateAT}} function.
#'
#' @details
#' Simplifies the use of the \code{\link{calculateAT}} and
#'  \code{\link{calculateAT}} function by providing a graphical user interface.
#'  In addition there are integrated control functions.
#' 
#' @param env environment in wich to search for data frames and save result.
#' @param savegui logical indicating if GUI settings should be saved in the environment.
#' @param debug logical indicating printing debug information.
#' @param parent widget to get focus when finished.
#' 
#' @return TRUE
#' 
#' @export
#' 
#' @importFrom utils help head str
#' @importFrom graphics title
#' @importFrom ggplot2 ggtitle scale_shape_discrete ggplot facet_wrap geom_point
#'  aes_string scale_colour_manual geom_rect
#' 
#' @seealso \code{\link{calculateAT}}, \code{\link{blockAT}},
#'  \code{\link{checkSubset}}


calculateAT_gui <- function(env=parent.frame(), savegui=NULL,
                            debug=FALSE, parent=NULL){
  
  # Global variables.
  .gData <- NULL
  .gSamples <- NULL
  .gDataPrep <- NULL
  .gPlot <- NULL
  .gRef <- NULL
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # WINDOW ####################################################################
  
  if(debug){
    print("WINDOW")
  }  
  
  # Main window.
  w <- gwindow(title="Calculate analytical threshold", visible=FALSE)
  
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
               spacing=5,
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
    print(help("calculateAT_gui", help_type="html"))
    
  })
  
  # FRAME 0 ###################################################################
  
  if(debug){
    print("FRAME 0")
  }  
  
  f0 <- gframe(text = "Datasets",
               horizontal=TRUE,
               spacing = 5,
               container = gv) 
  
  g0 <- glayout(container = f0, spacing = 1)
  
  # Dataset -------------------------------------------------------------------
  
  g0[1,1] <- glabel(text="Select dataset:", container=g0)
  
  dfs <- c("<Select a dataset>", listObjects(env=env, obj.class="data.frame"))
  
  g0[1,2] <- g0_data_drp <- gdroplist(items=dfs, 
                                      selected = 1,
                                      editable = FALSE,
                                      container = g0)
  g0[1,3] <- g0_data_samples_lbl <- glabel(text=" 0 samples", container=g0)
  
  addHandlerChanged(g0_data_drp, handler = function (h, ...) {
    
    val_obj <- svalue(g0_data_drp)
    
    # Check if suitable.
    requiredCol <- c("Dye.Sample.Peak", "Sample.File.Name", "Marker", "Allele",
                     "Marker", "Height", "Data.Point")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       slim=TRUE, slimcol=c("Allele", "Height","Data.Point"),
                       env=env, parent=w, debug=debug)
    
    if(ok){
      # Load or change components.
      
      # get dataset.
      .gData <<- get(val_obj, envir=env)
      svalue(g0_data_samples_lbl) <- paste(length(unique(.gData$Sample.File.Name)),
                                           "samples.")
      .refresh_sample_drp()
      .gDataPrep <- NULL # Erase any previously prepared data.
      
      # Suggest a name for result.
      svalue(f4_save1_edt) <- paste(val_obj, "_at", sep="")
      svalue(f4_save2_edt) <- paste(val_obj, "_rank", sep="")
      
    } else {
      
      # Reset components.
      .gData <<- NULL
      svalue(g0_data_drp, index=TRUE) <- 1
      svalue(g0_data_samples_lbl) <- " 0 samples"
      svalue(f4_save1_edt) <- ""
      svalue(f4_save2_edt) <- ""
      svalue(f4_save3_edt) <- ""
      .refresh_sample_drp()
      .gDataPrep <- NULL # Erase any previously prepared data.
      
    }
    
  } )  
  
  # Reference -----------------------------------------------------------------
  
  g0[2,1] <- glabel(text="Select reference dataset:", container=g0)
  
  # NB! dfs defined in previous section.
  g0[2,2] <- g0_ref_drp <- gdroplist(items=dfs, 
                                     selected = 1,
                                     editable = FALSE,
                                     container = g0)
  
  g0[2,3] <- g0_ref_samples_lbl <- glabel(text=" 0 references", container=g0)
  
  addHandlerChanged(g0_ref_drp, handler = function (h, ...) {
    
    val_obj <- svalue(g0_ref_drp)
    
    # Check if suitable.
    requiredCol <- c("Sample.Name", "Marker", "Allele")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       slim=TRUE, slimcol="Allele",
                       env=env, parent=w, debug=debug)
    
    if(ok){
      # Load or change components.
      
      .gRef <<- get(val_obj, envir=env)
      svalue(g0_ref_samples_lbl) <- paste(length(unique(.gRef$Sample.Name)),
                                          "samples.")
      
    } else {
      
      # Reset components.
      .gRef <<- NULL
      svalue(g0_ref_drp, index=TRUE) <- 1
      svalue(g0_ref_samples_lbl) <- " 0 references"
      
    }
    
  } )  
  
  # CHECK ---------------------------------------------------------------------
  
  if(debug){
    print("CHECK")
  }  
  
  g0[3,2] <- g0_check_btn <- gbutton(text="Check subsetting",
                                     border=TRUE,
                                     container=g0)
  
  addHandlerChanged(g0_check_btn, handler = function(h, ...) {
    
    # Get values.
    val_data <- .gData
    val_ref <- .gRef
    val_ignore <- svalue(f1_ignore_chk)
    val_word <- svalue(f1_word_chk)
    
    if (!is.null(.gData) || !is.null(.gRef)){
      
      chksubset_w <- gwindow(title = "Check subsetting",
                             visible = FALSE, name=title,
                             width = NULL, height= NULL, parent=w,
                             handler = NULL, action = NULL)
      
      chksubset_txt <- checkSubset(data=val_data,
                                   ref=val_ref,
                                   console=FALSE,
                                   ignore.case=val_ignore,
                                   word=val_word)
      
      gtext (text = chksubset_txt, width = NULL, height = 300, font.attr = NULL, 
             wrap = FALSE, container = chksubset_w)
      
      visible(chksubset_w) <- TRUE
      
    } else {
      
      gmessage(message="Data frame is NULL!\n\n
               Make sure to select a dataset and a reference set",
               title="Error",
               icon = "error")      
      
    } 
    
  } )
  
  # FRAME 1 ###################################################################
  
  if(debug){
    print("FRAME 1")
  }  
  
  f1 <- gframe(text = "Options",
               horizontal=FALSE,
               spacing = 10,
               container = gv) 
  
  f1_ignore_chk <- gcheckbox(text="Ignore case",
                             checked=TRUE,
                             container=f1)
  
  f1_word_chk <- gcheckbox(text="Add word boundaries",
                           checked = FALSE,
                           container = f1)
  
  # LAYOUT --------------------------------------------------------------------
  
  f1g1 <- glayout(container = f1, spacing = 1)
  
  f1g1[1,1] <- f1_block_h_chk <- gcheckbox(text="Block high peaks",
                                           checked = TRUE, container = f1g1)
  f1g1[1,2] <- glabel(text="Block all peaks above (RFU): ", anchor=c(-1, 0), container=f1g1)
  
  f1g1[1,3] <- f1_block_h_edt <- gedit(text="200", width=6, container=f1g1)
  
  f1g1[2,1] <- f1_block_chk <- gcheckbox(text="Block sample alleles",
                                         checked=TRUE, container=f1g1)
  f1g1[3,1] <- f1_block_d_chk <- gcheckbox(text="Block sample alleles per dye channel",
                                           checked=TRUE, container=f1g1)
  f1g1[2,2] <- glabel(text="Range (data points) around known alleles:", anchor=c(-1, 0), container=f1g1)
  f1g1[2,3] <- f1_block_spb <- gspinbutton(from=0, to=100, by=10, value=50, container=f1g1)
  
  f1g1[4,1] <- f1_block_ils_chk <- gcheckbox(text="Block ILS peaks",
                                             checked=TRUE, container=f1g1)
  f1g1[4,2] <- glabel(text="Range (data points) around known peak: ", anchor=c(-1, 0), container=f1g1)
  f1g1[4,3] <- f1_block_ils_spb <- gspinbutton(from=0, to=100, by=20, value=10, container=f1g1)
  
  # LAYOUT --------------------------------------------------------------------
  
  f1g2 <- glayout(container = f1, spacing = 1)
  
  f1g2[1,1] <- glabel(text="Confidence level 'k' (AT1): ", container=f1g2)
  f1g2[1,2] <- f1_k_spb <- gspinbutton(from=0, to=100, by=1, value=3, container=f1g2)
  
  
  f1g2[2,1] <- glabel(text="Percentile rank threshold (AT2): ", container=f1g2)
  f1g2[2,2] <- f1_t_spb <- gspinbutton(from=0, to=1, by=0.01, value=0.99, container=f1g2)
  
  f1g2[3,1] <- glabel(text="Upper confidence 'alpha' (AT4): ", container=f1g2)
  f1g2[3,2] <- f1_a_spb <- gspinbutton(from=0, to=1, by=0.01, value=0.01, container=f1g2)
  
  # Handlers ------------------------------------------------------------------
  
  addHandlerChanged(f1_block_h_chk, handler = function(h, ...) {
    
    # Update otions.
    .refresh_options()
    
  })
  
  addHandlerChanged(f1_block_chk, handler = function(h, ...) {
    
    # Update otions.
    .refresh_options()
    
  })
  
  addHandlerChanged(f1_block_ils_chk, handler = function(h, ...) {
    
    # Update otions.
    .refresh_options()
    
  })
  
  
  
  # FRAME 3 ###################################################################
  
  if(debug){
    print("FRAME 3")
  }  
  
  f3 <- gframe(text = "Prepare data and check blocking",
               horizontal=TRUE,
               spacing = 5,
               container = gv) 
  
  block_btn <- gbutton(text="Prepare and block", border=TRUE,
                       container=f3)
  
  f3_sample_drp <- gdroplist(items="<Select sample>", selected=1,
                             editable=FALSE, container=f3)
  
  save_btn <- gbutton(text="Save plot", border=TRUE, container=f3)
  
  addHandlerChanged(block_btn, handler = function(h, ...) {
    
    # Get values.
    val_data <- .gData
    val_ref <- .gRef
    val_block_h <- svalue(f1_block_h_chk)
    val_block <- svalue(f1_block_chk)
    val_block_d <- svalue(f1_block_d_chk)
    val_block_ils <- svalue(f1_block_ils_chk)
    val_height <- as.numeric(svalue(f1_block_h_edt))
    val_range <- svalue(f1_block_spb)
    val_range_ils <- svalue(f1_block_ils_spb)
    val_ignore <- svalue(f1_ignore_chk)
    val_word <- svalue(f1_word_chk)
    
    if(debug){
      print("Read Values:")
      print("val_data")
      print(head(val_data))
      print("val_block_h")
      print(val_block_h)
      print("val_block")
      print(val_block)
      print("val_range")
      print(val_range)
      print("val_block_d")
      print(val_block_d)
      print("val_block_ils")
      print(val_block_ils)
      print("val_range_ils")
      print(val_range_ils)
      print("val_ignore")
      print(val_ignore)
      print("val_word")
      print(val_word)
    }
    
    # Prepare data.
    .gDataPrep <<- blockAT(data=val_data, ref=val_ref,
                           block.height=val_block_h,
                           height=val_height,
                           block.sample=val_block,
                           per.dye=val_block_d,
                           range.sample=val_range,
                           block.ils=val_block_ils,
                           range.ils=val_range_ils,
                           ignore.case=val_ignore,
                           word=val_word,
                           debug=debug)
    
    # Unselect sample.
    svalue(f3_sample_drp, index=TRUE) <- 1
    
  })
  
  addHandlerChanged(f3_sample_drp, handler = function(h, ...) {
    
    # Get values.
    val_sample <- svalue(f3_sample_drp)
    
    if(!is.null(.gDataPrep) & !is.null(val_sample)){
      
      # Get values.
      val_block_h <- svalue(f1_block_h_chk)
      val_block <- svalue(f1_block_chk)
      val_block_d <- svalue(f1_block_d_chk)
      val_block_ils <- svalue(f1_block_ils_chk)
      val_range <- svalue(f1_block_spb)
      val_range_ils <- svalue(f1_block_ils_spb)
      
      if(val_sample %in% unique(.gDataPrep$Sample.File.Name)){
        
        # Must come after 'val_sample'.
        val_data <- subset(.gDataPrep, Sample.File.Name==val_sample)
        
        if(debug){
          print("Read Values:")
          print("val_data")
          print(head(val_data))
          print("val_sample")
          print(val_sample)
          print("val_block_h")
          print(val_block_h)
          print("val_block")
          print(val_block)
          print("val_range")
          print(val_range)
          print("val_block_d")
          print(val_block_d)
          print("val_block_ils")
          print(val_block_ils)
          print("val_range_ils")
          print(val_range_ils)
        }
        
        # Get all dyes.
        dyes <- as.character(unique(val_data$Dye))
        colorsKit <- unique(getKit("ESX17", what="Color")$Color)
        dyesKit <- addColor(colorsKit, have="Color", need="Dye")
        dyeILS <- setdiff(dyes, dyesKit)
        
        # Refactor and keep order of levels.
        val_data$Dye <- factor(val_data$Dye, levels=unique(val_data$Dye))
        
        # Create plot.
        gp <- ggplot(data=val_data)
        gp <- gp + ggtitle(paste("Blocked data for", val_sample))
        gp <- gp + facet_wrap( ~ Dye, ncol=1, scales = "fixed", drop = FALSE)
        if(val_block_h){
          # Change shape, color, and legend.
          gp <- gp + geom_point(aes_string(x="Data.Point", y="Height",
                                           colour="H.Block", shape="H.Block"))
          gp <- gp + scale_shape_discrete(name  ="Peaks",
                                          breaks=c(FALSE, TRUE),
                                          labels=c("Included", "Excluded"))
          gp <- gp + scale_colour_manual(values=c("black", "red"),
                                         name="Peaks",
                                         breaks=c(FALSE, TRUE),
                                         labels=c("Included", "Excluded"))
        } else {
          # Use default color and shape.
          gp <- gp + geom_point(aes_string(x="Data.Point", y="Height"))
        }
        
        if(debug){
          print("Plot created!")
        }
        
        if(val_block_ils){
          
          # ILS blocking data frame for plot:
          dfIls <- val_data[val_data$ILS == TRUE, ]
          ilsDye <- unique(dfIls$Dye)
          dpBlock <- dfIls$Data.Point
          dyeBlock <- rep(unique(val_data$Dye), each=length(dpBlock))
          dpBlock <- rep(dpBlock, length(unique(val_data$Dye)))
          dfBlock <- data.frame(Dye=dyeBlock, Data.Point=dpBlock,
                                Xmin=dpBlock - val_range_ils,
                                Xmax=dpBlock + val_range_ils)
          
          if(nrow(dfBlock) > 0){
            
            # Add blocking range to plot.
            gp <- gp + geom_rect(data = dfBlock,
                                 aes_string(ymin = -Inf, ymax = Inf,
                                            xmin = "Xmin", xmax = "Xmax"),
                                 alpha = 0.2,
                                 fill = addColor(ilsDye, have="Dye", need="Color"))
            
          }
          
          if(debug){
            print(str(dfBlock))
            print("ILS blocked!")
          }
          
        }
        
        if(val_block){
          # Sample blocking data frame for plot:
          dfSample <- val_data[val_data$S.Block==TRUE, ]
          
          if(val_block_d){
            
            # Loop over dyes and add block ranges.
            for(d in seq(along=dyesKit)){
              # Get data points for selected sample.
              dpBlock <- dfSample$Data.Point[dfSample$Dye==dyesKit[d]]
              dyeBlock <- rep(dyesKit[d], length(dpBlock))
              dfBlock <- data.frame(Dye=dyeBlock, Data.Point=dpBlock,
                                    Xmin=dpBlock - val_range,
                                    Xmax=dpBlock + val_range)
              
              if(nrow(dfBlock) > 0){
                
                # Add blocking range to plot.
                gp <- gp + geom_rect(data = dfBlock,
                                     aes_string(ymin = -Inf, ymax = Inf,
                                                xmin = "Xmin", xmax = "Xmax"),
                                     alpha = 0.2, fill = colorsKit[d])
                
              }
              
            }
            
            if(debug){
              print(str(dfBlock))
              print("Sample blocked per dye!")
            }
            
          } else {
            
            dpBlock <- dfSample$Data.Point
            dyeBlock <- rep(dyesKit, each=length(dpBlock))
            dpBlock <- rep(dpBlock, length(dyesKit))
            dfBlock <- data.frame(Dye=dyeBlock, Data.Point=dpBlock,
                                  Xmin=dpBlock - val_range,
                                  Xmax=dpBlock + val_range)
            
            if(nrow(dfBlock) > 0){
              
              # Add blocking range to plot.
              gp <- gp + geom_rect(data = dfBlock,
                                   aes_string(ymin = -Inf, ymax = Inf,
                                              xmin = "Xmin", xmax = "Xmax"),
                                   alpha = 0.2, fill="red")
              
            }
            
            if(debug){
              print(str(dfBlock))
              print("Sample blocked across dyes!")
            }
            
          }
          
        }
        
        # Show plot.
        print(gp)
        
        # Save plot object.
        .gPlot <<- gp
        
      } # End 'sample exist' if.
      
    } # End 'data exist' if.
    
  })
  
  addHandlerChanged(save_btn, handler = function(h, ...) {
    
    # Get sample name.
    val_name <- svalue(f3_sample_drp) 
    
    if(!is.null(.gPlot)){
      
      # Save data.
      ggsave_gui(ggplot=.gPlot, name=val_name, parent=w, env=env,
                 savegui=savegui, debug=debug)
      
    } else {
      
      message <- "Click 'Prepare and block' and select a sample before saving."
      gmessage(message, title="No plot!",
               icon = "info", parent = w) 
      
    }
    
  })
  
  
  # FRAME 4 ###################################################################
  
  if(debug){
    print("FRAME 4")
  }  
  
  f4 <- gframe(text = "Save as",
               horizontal=FALSE,
               spacing = 5,
               container = gv) 
  
  glabel(text="Name for result:", anchor=c(-1, 0), container=f4)
  
  f4_save1_edt <- gedit(text="", container=f4)
  
  glabel(text="Name for percentile rank list:", anchor=c(-1, 0), container=f4)
  
  f4_save2_edt <- gedit(text="", container=f4)
  
  # BUTTON ####################################################################
  
  if(debug){
    print("BUTTON")
  }  
  
  calculate_btn <- gbutton(text="Calculate",
                           border=TRUE,
                           container=gv)
  
  addHandlerChanged(calculate_btn, handler = function(h, ...) {
    
    # Get values.
    if(is.null(.gDataPrep)){
      val_data <- .gData
    } else {
      val_data <- .gDataPrep
    }
    val_ref <- .gRef
    val_ignore <- svalue(f1_ignore_chk)
    val_word <- svalue(f1_word_chk)
    val_block_h <- svalue(f1_block_h_chk)
    val_block <- svalue(f1_block_chk)
    val_block_d <- svalue(f1_block_d_chk)
    val_block_ils <- svalue(f1_block_ils_chk)
    val_height <- as.numeric(svalue(f1_block_h_edt))
    val_range <- svalue(f1_block_spb)
    val_range_ils <- svalue(f1_block_ils_spb)
    val_k <- svalue(f1_k_spb)
    val_t <- svalue(f1_t_spb)
    val_a <- svalue(f1_a_spb)
    val_name1 <- svalue(f4_save1_edt)
    val_name2 <- svalue(f4_save2_edt)
    
    if(debug){
      print("Read Values:")
      print("val_data")
      print(head(val_data))
      print("val_ref")
      print(head(val_ref))
      print("val_ignore")
      print(val_ignore)
      print("val_word")
      print(val_word)
      print("val_block_h")
      print(val_block_h)
      print("val_height")
      print(val_height)
      print("val_block")
      print(val_block)
      print("val_range")
      print(val_range)
      print("val_block_d")
      print(val_block_d)
      print("val_block_ils")
      print(val_block_ils)
      print("val_range_ils")
      print(val_range_ils)
      print("val_k")
      print(val_k)
      print("val_t")
      print(val_t)
      print("val_a")
      print(val_a)
      print("val_name1")
      print(val_name1)
      print("val_name2")
      print(val_name2)
    }
    
    # Check if data.
    if(!is.null(val_data)){
      
      # If sample blocking is active a reference must be selected.
      if((val_block & !is.null(.gRef)) | (!val_block)){
        
        
        # Change button.
        svalue(calculate_btn) <- "Processing..."
        enabled(calculate_btn) <- FALSE
        
        datanew <- calculateAT(data=val_data,
                               ref=val_ref,
                               block.height=val_block_h,
                               height=val_height,
                               block.sample=val_block,
                               per.dye=val_block_d,
                               range.sample=val_range,
                               block.ils=val_block_ils,
                               range.ils=val_range_ils,
                               k=val_k,
                               rank.t=val_t,
                               alpha=val_a,
                               ignore.case=val_ignore,
                               word=val_word,
                               debug=debug)
        

        # Add attributes.
        attr(datanew[[1]], which="calculateAT_gui, data") <- svalue(g0_data_drp)
        attr(datanew[[1]], which="calculateAT_gui, ref") <- svalue(g0_ref_drp)
        attr(datanew[[1]], which="calculateAT_gui, k") <- val_k
        attr(datanew[[1]], which="calculateAT_gui, rank.t") <- val_t
        attr(datanew[[1]], which="calculateAT_gui, alpha") <- val_a
        attr(datanew[[1]], which="calculateAT_gui, block.height") <- val_block_h
        attr(datanew[[1]], which="calculateAT_gui, height") <- val_height
        attr(datanew[[1]], which="calculateAT_gui, block") <- val_block
        attr(datanew[[1]], which="calculateAT_gui, range.sample") <- val_range
        attr(datanew[[1]], which="calculateAT_gui, block.ils") <- val_block_ils
        attr(datanew[[1]], which="calculateAT_gui, range.ils") <- val_range_ils
        attr(datanew[[1]], which="calculateAT_gui, per.dye") <- val_block_d
        attr(datanew[[1]], which="calculateAT_gui, ignore.case") <- val_ignore
        attr(datanew[[1]], which="calculateAT_gui, word") <- val_word
        
        attr(datanew[[2]], which="calculateAT_gui, data") <- svalue(g0_data_drp)
        attr(datanew[[2]], which="calculateAT_gui, ref") <- svalue(g0_ref_drp)
        attr(datanew[[2]], which="calculateAT_gui, rank.t") <- val_t
        attr(datanew[[2]], which="calculateAT_gui, block.height") <- val_block_h
        attr(datanew[[2]], which="calculateAT_gui, height") <- val_height
        attr(datanew[[2]], which="calculateAT_gui, block") <- val_block
        attr(datanew[[2]], which="calculateAT_gui, range.sample") <- val_range
        attr(datanew[[2]], which="calculateAT_gui, block.ils") <- val_block_ils
        attr(datanew[[2]], which="calculateAT_gui, range.ils") <- val_range_ils
        attr(datanew[[2]], which="calculateAT_gui, per.dye") <- val_block_d
        attr(datanew[[2]], which="calculateAT_gui, ignore.case") <- val_ignore
        attr(datanew[[2]], which="calculateAT_gui, word") <- val_word
        
        # Save data.
        saveObject(name=val_name1, object=datanew[[1]], parent=w, env=env)
        saveObject(name=val_name2, object=datanew[[2]], parent=w, env=env)
        
        
        if(debug){
          print(str(datanew))
          print(head(datanew))
          print(paste("EXIT:", match.call()[[1]]))
        }
        
        # Close GUI.
        dispose(w)
        
      } else {
        
        message <- "If 'Block sample alleles' is active a reference dataset must be selected."
        
        gmessage(message, title="Datasets not selected",
                 icon = "error",
                 parent = w) 
        
      }
      
    } else {
      
      message <- "A dataset must be selected."
      
      gmessage(message, title="Datasets not selected",
               icon = "error",
               parent = w) 
      
    }
    
  } )
  
  # INTERNAL FUNCTIONS ########################################################
  
  .refresh_sample_drp <- function(){
    
    if(debug){
      print("Refresh sample dropdown")
    }
    
    # Get data frames in global workspace.
    samples <- unique(.gData$Sample.File.Name)
    
    if(!is.null(samples)){
      
      # Populate drop list.
      f3_sample_drp[] <- c("<Select sample>", samples)
      
    } else {
      
      # Populate drop list.
      f3_sample_drp[] <- c("<Select sample>")
      
    }
    
    if(debug){
      print("Samples dropdown refreshed!")
    }
  }
  
  .refresh_options <- function(){
    
    if(debug){
      print("Refresh options")
    }
    
    val_block_h <- svalue(f1_block_h_chk)
    val_block <- svalue(f1_block_chk)
    val_block_d <- svalue(f1_block_d_chk)
    val_block_ils <- svalue(f1_block_ils_chk)
    
    # Update dependent widgets.
    if(val_block_h){
      enabled(f1_block_h_edt) <- TRUE
    } else {
      enabled(f1_block_h_edt) <- FALSE
    }
    
    # Update dependent widgets.
    if(val_block){
      enabled(f1_block_d_chk) <- TRUE
      enabled(f1_block_spb) <- TRUE
    } else {
      enabled(f1_block_d_chk) <- FALSE
      enabled(f1_block_spb) <- FALSE
    }
    
    # Update dependent widgets.
    if(val_block_ils){
      enabled(f1_block_ils_spb) <- TRUE
    } else {
      enabled(f1_block_ils_spb) <- FALSE
    }

    if(debug){
      print("Options refreshed!")
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
      if(exists(".strvalidator_calculateAT_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_calculateAT_gui_savegui", envir=env)
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
      if(exists(".strvalidator_calculateAT_gui_block_h", envir=env, inherits = FALSE)){
        svalue(f1_block_h_chk) <- get(".strvalidator_calculateAT_gui_block_h", envir=env)
      }
      if(exists(".strvalidator_calculateAT_gui_block", envir=env, inherits = FALSE)){
        svalue(f1_block_chk) <- get(".strvalidator_calculateAT_gui_block", envir=env)
      }
      if(exists(".strvalidator_calculateAT_gui_block_ils", envir=env, inherits = FALSE)){
        svalue(f1_block_ils_chk) <- get(".strvalidator_calculateAT_gui_block_ils", envir=env)
      }
      if(exists(".strvalidator_calculateAT_gui_dye", envir=env, inherits = FALSE)){
        svalue(f1_block_d_chk) <- get(".strvalidator_calculateAT_gui_dye", envir=env)
      }
      if(exists(".strvalidator_calculateAT_gui_height", envir=env, inherits = FALSE)){
        svalue(f1_block_h_edt) <- get(".strvalidator_calculateAT_gui_height", envir=env)
      }
      if(exists(".strvalidator_calculateAT_gui_range", envir=env, inherits = FALSE)){
        svalue(f1_block_spb) <- get(".strvalidator_calculateAT_gui_range", envir=env)
      }
      if(exists(".strvalidator_calculateAT_gui_range_ils", envir=env, inherits = FALSE)){
        svalue(f1_block_ils_spb) <- get(".strvalidator_calculateAT_gui_range_ils", envir=env)
      }
      if(exists(".strvalidator_calculateAT_gui_k", envir=env, inherits = FALSE)){
        svalue(f1_k_spb) <- get(".strvalidator_calculateAT_gui_k", envir=env)
      }
      if(exists(".strvalidator_calculateAT_gui_t", envir=env, inherits = FALSE)){
        svalue(f1_t_spb) <- get(".strvalidator_calculateAT_gui_t", envir=env)
      }
      if(exists(".strvalidator_calculateAT_gui_a", envir=env, inherits = FALSE)){
        svalue(f1_a_spb) <- get(".strvalidator_calculateAT_gui_a", envir=env)
      }
      if(exists(".strvalidator_calculateAT_gui_ignore", envir=env, inherits = FALSE)){
        svalue(f1_ignore_chk) <- get(".strvalidator_calculateAT_gui_ignore", envir=env)
      }
      if(exists(".strvalidator_calculateAT_gui_word", envir=env, inherits = FALSE)){
        svalue(f1_word_chk) <- get(".strvalidator_calculateAT_gui_word", envir=env)
      }
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_calculateAT_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_calculateAT_gui_block_h", value=svalue(f1_block_h_chk), envir=env)
      assign(x=".strvalidator_calculateAT_gui_block", value=svalue(f1_block_chk), envir=env)
      assign(x=".strvalidator_calculateAT_gui_block_ils", value=svalue(f1_block_ils_chk), envir=env)
      assign(x=".strvalidator_calculateAT_gui_dye", value=svalue(f1_block_d_chk), envir=env)
      assign(x=".strvalidator_calculateAT_gui_height", value=svalue(f1_block_h_edt), envir=env)
      assign(x=".strvalidator_calculateAT_gui_range", value=svalue(f1_block_spb), envir=env)
      assign(x=".strvalidator_calculateAT_gui_range_ils", value=svalue(f1_block_ils_spb), envir=env)
      assign(x=".strvalidator_calculateAT_gui_k", value=svalue(f1_k_spb), envir=env)
      assign(x=".strvalidator_calculateAT_gui_t", value=svalue(f1_t_spb), envir=env)
      assign(x=".strvalidator_calculateAT_gui_a", value=svalue(f1_a_spb), envir=env)
      assign(x=".strvalidator_calculateAT_gui_ignore", value=svalue(f1_ignore_chk), envir=env)
      assign(x=".strvalidator_calculateAT_gui_word", value=svalue(f1_word_chk), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_calculateAT_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateAT_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_calculateAT_gui_block_h", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateAT_gui_block_h", envir = env)
      }
      if(exists(".strvalidator_calculateAT_gui_block", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateAT_gui_block", envir = env)
      }
      if(exists(".strvalidator_calculateAT_gui_block_ils", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateAT_gui_block_ils", envir = env)
      }
      if(exists(".strvalidator_calculateAT_gui_dye", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateAT_gui_dye", envir = env)
      }
      if(exists(".strvalidator_calculateAT_gui_height", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateAT_gui_height", envir = env)
      }
      if(exists(".strvalidator_calculateAT_gui_range", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateAT_gui_range", envir = env)
      }
      if(exists(".strvalidator_calculateAT_gui_range_ils", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateAT_gui_range_ils", envir = env)
      }
      if(exists(".strvalidator_calculateAT_gui_k", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateAT_gui_k", envir = env)
      }
      if(exists(".strvalidator_calculateAT_gui_t", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateAT_gui_t", envir = env)
      }
      if(exists(".strvalidator_calculateAT_gui_a", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateAT_gui_a", envir = env)
      }
      if(exists(".strvalidator_calculateAT_gui_ignore", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateAT_gui_ignore", envir = env)
      }
      if(exists(".strvalidator_calculateAT_gui_word", envir=env, inherits = FALSE)){
        remove(".strvalidator_calculateAT_gui_word", envir = env)
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
  
  # Update otions.
  .refresh_options()
  
  # Show GUI.
  visible(w) <- TRUE
  focus(w)
  
}
