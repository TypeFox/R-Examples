################################################################################
# TODO LIST
# TODO: implement log(P(D)) scale. Quite tricky...
# TODO: parameter perLocus? Not priority since easy to make subset of dataframe.

################################################################################
# CHANGE LOG (last 20 changes)
# 11.11.2015: Added importFrom ggplot2.
# 29.08.2015: Added importFrom.
# 19.08.2015: Added more information to the description.
# 18.06.2015: Rounded printed probabilities to three decimals.
# 20.03.2015: Rounded printed conservative drop-out threshold to integer.
# 05.01.2015: Changed check of suggested package ResourceSelection in accordance
#             with Writing R extensions v 3.2.1 section 1.1.3.1.
# 14.12.2014: Added option to use average peak height 'H'.
# 14.12.2014: Updated to handle gender -> sex.marker option in getKit.
# 11.10.2014: Added 'focus', added 'parent' parameter.
# 28.06.2014: Added help button and moved save gui checkbox.
# 28.06.2014: Changed notation on plot to be more correct.
# 08.05.2014: Implemented 'checkDataset'.
# 18.02.2014: Implemented conserative T estimate.
# 18.02.2014: Removed erroneously implemented prediction interval for T.
# 27.01.2014: Fixed bug not checking required columns upon selection of dataset.
# 20.01.2014: Changed 'saveImage_gui' for 'ggsave_gui'.
# 16.01.2014: Changed according to new column names in 'calculateDropout' res.
# 13.11.2013: Implemented 'Hosmer-Lemeshow test'.
# 13.11.2013: Nicer code for multiple input choices to model + bug in limit x/y.
# 06.11.2013: Fixed prediction interval for log.

#' @title Model And Plot Drop-out Events
#'
#' @description
#' Model the probability of drop-out and plot graphs.
#'
#' @details
#' \code{\link{calculateDropout}} score drop-out events relative to a user
#' defined LDT in four different ways:
#' (1) by reference to the low molecular weight allele (Method1),
#' (2) by reference to the high molecular weight allele (Method2),
#' (3) by reference to a random allele (MethodX), and
#' (4) by reference to the locus (MethodL).
#' Options 1-3 are recommended by the DNA commision (see reference),
#' while option 4 is included for experimental purposes.
#' Options 1-3 may discard many dropout events while option 4 catches all
#' drop-out events. On the other hand options 1-3 can score events below
#' the LDT, while option 4 cannot, making accurate predictions possible
#' below the LDT.
#' 
#' Method X/1/2 records the peak height of the partner allele to be used as
#' the explanatory variable in the logistic regression. The locus method L also
#' do this when there has been a drop-out, if not the the mean peak height for
#' the locus is used. Peak heights for the locus method are stored in a
#' separate column.
#' 
#' Using the scored drop-out events and the peak heights of the surviving
#' alleles the probability of drop-out can be modelled by logistic regression
#' as described in Appendix B the first reference.
#' P(dropout|H) = B0 + B1*H, where 'H' is the peak height or log(peak height).
#' This produces a plot with the predicted probabilities for a range of peak heights.
#' There are options to print the model parameters, mark the stochastic
#' threshold at a specified probability of drop-out, include the underlying
#' observations, and to calculate a specified prediction interval.
#' A conservative estimate of the stochastic threshold can be calculated
#' from the prediction interval: the risk of observing a drop-out probability
#' greater than the specified threshold limit, at the conservative peak height,
#' is less than a specified value (e.g. 1-0.95=0.05). By default the gender
#' marker is excluded from the dataset used for modelling, and the peak height
#' is used as explanatory variable. The average peak height 'H' can be used
#' instead of the allele/locus peak height. Optionally, the logarithm of the peak
#' height can be used. To evaluate the goodness of fit for the logistic
#' regression the Hosmer-Lemeshow test is used. A value below 0.05 indicates
#' a poor fit. Alternatives to the logistic regression method are discussed in
#' reference 4 and 5. [13,32].
#' 
#' Explanation of the result:
#' Dropout - all alleles are scored according to the limit of detection threshold (LDT).
#' This is the observations and is not used for modelling.
#' Rfu - peak height of the surviving allele.
#' MethodX - a random reference allele is selected and drop-out is scored in
#' relation to the the partner allele.
#' Method1 - the low molecular weight allele is selected and drop-out is
#' scored if the high molecular weight allele is missing.
#' Method2 - the high molecular weight allele is selected and drop-out is
#' scored if the low molecular weight allele is missing.
#' MethodL - drop-out is scored per locus i.e. drop-out if any allele is missing.
#' MethodL.Ph - peak height of the surviving allele if one allele has droped out,
#' or the average peak height if no drop-out.
#' 
#' @param env environment in wich to search for data frames and save result.
#' @param savegui logical indicating if GUI settings should be saved in the environment.
#' @param debug logical indicating printing debug information.
#' @param parent widget to get focus when finished.
#' 
#' @references
#' Peter Gill et.al.,
#'  DNA commission of the International Society of Forensic Genetics:
#'  Recommendations on the evaluation of STR typing results that may
#'  include drop-out and/or drop-in using probabilistic methods,
#'  Forensic Science International: Genetics, Volume 6, Issue 6, December 2012,
#'  Pages 679-688, ISSN 1872-4973, 10.1016/j.fsigen.2012.06.002.
#' \url{http://www.sciencedirect.com/science/article/pii/S1872497312001354}
#' @references
#' Peter Gill, Roberto Puch-Solis, James Curran,
#'  The low-template-DNA (stochastic) threshold-Its determination relative to
#'  risk analysis for national DNA databases,
#'  Forensic Science International: Genetics, Volume 3, Issue 2, March 2009,
#'  Pages 104-111, ISSN 1872-4973, 10.1016/j.fsigen.2008.11.009.
#' \url{http://www.sciencedirect.com/science/article/pii/S1872497308001798}
#' @references
#' H. DW Jr., S. Lemeshow, Applied Logistic Regression, John Wiley & Sons, 2004.
#' @references
#' A.A. Westen, L.J.W. Grol, J. Harteveld, A.S. Matai, P. de Knijff, T. Sijen,
#'  Assessment of the stochastic threshold, back- and forward stutter filters
#'  and low template techniques for NGM,
#'  Forensic Science International: Genetetics, Volume 6, Issue 6 December 2012,
#'  Pages 708-715, ISSN 1872-4973, 10.1016/j.fsigen.2012.05.001.
#'  \url{http://www.sciencedirect.com/science/article/pii/S1872497312001172}
#' @references
#' R. Puch-Solis, A.J. Kirkham, P. Gill, J. Read, S. Watson, D. Drew,
#'  Practical determination of the low template DNA threshold,
#'  Forensic Science International: Genetetics, Volume 5, Issue 5, November 2011,
#'  Pages 422-427, ISSN 1872-4973, 10.1016/j.fsigen.2010.09.001.
#'  \url{http://www.sciencedirect.com/science/article/pii/S1872497310001511}
#' 
#' @return TRUE
#' 
#' @export
#' 
#' @importFrom utils help str head
#' @importFrom grDevices palette
#' @importFrom stats glm binomial fitted predict plogis qnorm
#' @importFrom ggplot2 ggplot aes_string geom_line geom_point position_jitter
#'  geom_ribbon geom_segment geom_text coord_cartesian theme element_text labs
#'  xlab ylab
#'  
#' @seealso \code{\link{calculateDropout}}, \code{\link{plotDropout_gui}}

modelDropout_gui <- function(env=parent.frame(), savegui=NULL, debug=FALSE, parent=NULL){

  # Global variables.
  .gData <- NULL
  .gPlot <- NULL
  
  if(debug){
    print(paste("IN:", match.call()[[1]]))
  }
  
  # Main window.
  w <- gwindow(title="Plot dropout prediction", visible=FALSE)
  
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
    print(help("modelDropout_gui", help_type="html"))
    
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
  
  glabel(text=" and the kit used:", container=f0)
  
  kit_drp <- gdroplist(items=getKit(), 
                       selected = 1,
                       editable = FALSE,
                       container = f0) 

  addHandlerChanged(dataset_drp, handler = function (h, ...) {
    
    val_obj <- svalue(dataset_drp)
    
    # Check if suitable.
    requiredCol <- c("Height")
    ok <- checkDataset(name=val_obj, reqcol=requiredCol,
                       env=env, parent=w, debug=debug)
    
    if(ok){
      
      # Load or change components.
      .gData <<- get(val_obj, envir=env)
      
      ph_range <- NA
      if("Heterozygous" %in% names(.gData)){
        # Make sure numeric, then find range for heterozygotes.
        ph_range <- range(as.numeric(.gData$Height[.gData$Heterozygous==1 & .gData$Dropout!=2]),
                          na.rm=TRUE)
      } else {
        # Make sure numeric, then find min and max.
        ph_range <- range(as.numeric(.gData$Height), na.rm=TRUE)
      }
      svalue(f1g2_low_lbl) <- ph_range[1]
      svalue(f1g2_high_lbl) <- ph_range[2]
      
      # Suggest name.
      svalue(f5_save_edt) <- paste(val_obj, "_ggplot", sep="")
      
      # Detect kit.
      kitIndex <- detectKit(.gData, index=TRUE)
      # Select in dropdown.
      svalue(kit_drp, index=TRUE) <- kitIndex
      
      # Check additional required columns and enable/disable plot button.
      .checkColumns()
        
    } else {

      # Reset components.
      .gData <<- NULL
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
      enabled(f1g1) <- TRUE
    } else {
      enabled(f1g1) <- FALSE
    }
  } )
  
  f1g1 <- glayout(container = f1, spacing = 1)
  enabled(f1g1) <- svalue(f1_titles_chk)

  # Legends
  f1g1[1,1] <- glabel(text="Plot title:", container=f1g1)
  f1g1[1,2] <- f1_title_edt <- gedit(text="",
                                     width=60,
                                     container=f1g1)
  
  f1g1[2,1] <- glabel(text="X title:", container=f1g1)
  f1g1[2,2] <- f1_x_title_edt <- gedit(text="",
                                       width=60,
                                       container=f1g1)
  
  f1g1[3,1] <- glabel(text="Y title:", container=f1g1)
  f1g1[3,2] <- f1_y_title_edt <- gedit(text="",
                                       width=60,
                                       container=f1g1)
  

  # Group 2.
  f1g2 <- ggroup(horizontal = TRUE, spacing = 5,  container = f1)
  glabel(text="Dataset peak height range:", container = f1g2)
  f1g2_low_lbl <- glabel(text="", width = 6, container = f1g2)
  glabel(text="-", container = f1g2)
  f1g2_high_lbl <- glabel(text="", width = 6, container = f1g2)
  glabel(text=" RFU", container = f1g2)
  

  # Other options.
  log_model <- gcheckbox(text="Log (Height)", checked=FALSE, container=f1)
  
  f1_sex_chk <- gcheckbox(text="Exclude sex markers",
                               checked = TRUE,
                               container = f1)
  
  glabel(text=paste("NB! Currently, the recommended methods are the first three options.\n",
                    "The fourth alternative has not been evaluated by the DNA Commission.",
                    "\nSee details for more information."),
         anchor=c(-1 ,0), container=f1)
  
  glabel(text="Model drop-out from scoring method:", anchor=c(-1 ,0), container=f1)
  f1_column_opt <- gradio(items=c("Relative a random allele and peak height of surviving allele",
                                  "Relative the low molecular weight allele and peak height of surviving allele",
                                  "Relative the high molecular weight allele and peak height of surviving allele",
                                  "Relative the locus and peak height of surviving allele, or mean locus peak height"),
                          selected = 2,
                          horizontal = FALSE,
                          container = f1)

  f1_h_chk <- gcheckbox(text="Use average peak height 'H' instead of allele/locus peak hight",
                                 checked = FALSE,
                                 container = f1)
  
  f1_printmodel_chk <- gcheckbox(text="Print model",
                                 checked = FALSE,
                                 container = f1)
  
  addHandlerChanged(f1_column_opt, handler = function(h, ...) {
    
    .checkColumns()
    
  } )
  
  addHandlerChanged(f1_h_chk, handler = function(h, ...) {
    
    .checkColumns()
    
  } )
  
  # FRAME 7 ###################################################################
  
  f7 <- gframe(text = "Plot drop-out data",
               horizontal=FALSE,
               container = gv) 
  
  f7g1 <- glayout(container = f7)
  
  f7g1[1,1] <- f7_plot_drop_btn <- gbutton(text="Plot predicted drop-out probability",
                                           border=TRUE,
                                           container=f7g1) 
  
  
  addHandlerChanged(f7_plot_drop_btn, handler = function(h, ...) {
    
    if(!is.null(.gData)){
      enabled(f7_plot_drop_btn) <- FALSE
      svalue(f7_plot_drop_btn) <- "Processing..."
      .plotDrop()
      svalue(f7_plot_drop_btn) <- "Plot predicted drop-out probability"
      enabled(f7_plot_drop_btn) <- TRUE
    } else {
      message <- paste("Select a drop-out dataset")
      
      gmessage(message, title="Could not find dataset",
               icon = "error",
               parent = w) 
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

  # EXPAND 1 ##################################################################

  e1 <- gexpandgroup(text="Drop-out prediction and threshold",
                     horizontal=FALSE,
                     container = f1)
  
  # FRAME 1 -------------------------------------------------------------------
  # DROPOUT THRESHOLD
  
  e1f1 <- gframe(text = "", horizontal = FALSE, container = e1) 

  # Group 2.
  e1f1g2 <- ggroup(horizontal = TRUE, spacing = 5,  container = e1f1)
  
  e1f1_threshold_chk <- gcheckbox(text="Mark threshold @ P(D):",
                                               checked = TRUE,
                                               container = e1f1g2)

  e1f1_risk_spn <- gspinbutton (from=0, to=1, by=0.001,
                                                 value=0.05,
                                                 container=e1f1g2)
  
  # Group 3.
  e1f1g3 <- ggroup(horizontal = TRUE, spacing = 5,  container = e1f1)
  
  e1_linetypes <- c("blank", "solid", "dashed", "dotted", "dotdash","longdash","twodash")
  
  glabel("Line type", container = e1f1g3)
  
  e1f1_t_linetype_drp <- gdroplist(items=e1_linetypes,
                                     selected=2,
                                     container = e1f1g3)
  
  glabel("Line colour", container = e1f1g3) 
  
  e1f1_t_linecolor_drp <- gdroplist(items=palette(),
                                selected=2,
                                container = e1f1g3)
  
  # Group 4.
  e1f1g4 <- ggroup(horizontal = TRUE, spacing = 5,  container = e1f1)
  
  e1f1_print_chk <- gcheckbox(text="Print threshold value",
                                               checked = TRUE,
                                               container = e1f1g4)
  
  # FRAME 2 -------------------------------------------------------------------
  # PREDICTION INTERVAL
  
  e1f2 <- gframe(text = "", horizontal = FALSE, container = e1) 

  # Group 1.
  e1f2g1 <- ggroup(horizontal = TRUE, spacing = 5,  container = e1f2)
  
  
  glabel(text="Prediction interval:", container = e1f2g1)
  
  e1f2_conf_spn <- gspinbutton (from=0, to=1, by=0.001,
                                                 value=0.950,
                                                 container = e1f2g1)
  
  # Group 2.
  e1f2g2 <- ggroup(horizontal = TRUE, spacing = 5,  container = e1f2)
  
  e1f2_print_interval_chk <- gcheckbox(text="Print conservative T value",
                                           checked = TRUE,
                                           container = e1f2g2)
  
  # Group 3.
  e1f2g3 <- ggroup(horizontal = TRUE, spacing = 5,  container = e1f2)
  
  e1f2_mark_interval_chk <- gcheckbox(text="Draw prediction interval:",
                                                checked = TRUE,
                                                container = e1f2g3)
  
  glabel("Alpha", container = e1f2g3) 
  
  e1f2_interval_spb <- gspinbutton (from=0, to=1, by=0.01,
                                                     value=0.25,
                                                     container=e1f2g3)
  
  glabel("Fill colour", container = e1f2g3) 
  e1f2_interval_drp <- gdroplist(items=palette(),
                                             selected=2,
                                             container = e1f2g3)
  
  # EXPAND 2 ##################################################################
  
  e2 <- gexpandgroup(text="Data points",
                     horizontal=FALSE,
                     container = f1)
  
  e2f1 <- gframe(text = "", horizontal = FALSE, container = e2) 
  
  e2g1 <- glayout(container = e2f1)
  
  e2g1[1,1] <- e2g1_plotpoints_chk <- gcheckbox(text="Plot data points",
                                               checked = TRUE,
                                               container = e2g1)
  e2g1[1,2] <- glabel(text="Shape:", container=e2g1)
  e2g1[1,3] <- e2g1_shape_spb <- gspinbutton(from=0, to=25,
                                         by=1, value=18,
                                         container=e2g1)
    
  e2g1[1,4] <- glabel(text="Alpha:", container=e2g1)
  e2g1[1,5] <- e2g1_alpha_spb <- gspinbutton(from=0, to=1,
                                         by=0.01, value=0.60,
                                         container=e2g1)
  
  e2g1[1,6] <- glabel(text="Jitter (h/v):", container=e2g1)
  e2g1[1,7] <- e2g1_jitterh_edt <- gedit(text="0", width=4, container=e2g1)
  e2g1[1,8] <- e2g1_jitterv_edt <- gedit(text="0", width=4, container=e2g1)
  
  # EXPAND 3 ##################################################################
  
  e3 <- gexpandgroup(text="Axes",
                     horizontal=FALSE,
                     container = f1)
  
  e3f1 <- gframe(text = "", horizontal = FALSE, container = e3) 

  glabel(text="NB! Must provide both min and max value.",
         anchor=c(-1 ,0), container=e3f1)

  e3g1 <- glayout(container = e3f1, spacing = 1)
  e3g1[1,1:2] <- glabel(text="Limit Y axis (min-max)", container=e3g1)
  e3g1[2,1] <- e3g1_y_min_edt <- gedit(text="", width=5, container=e3g1)
  e3g1[2,2] <- e3g1_y_max_edt <- gedit(text="", width=5, container=e3g1)
  
  e3g1[3,1:2] <- glabel(text="Limit X axis (min-max)", container=e3g1)
  e3g1[4,1] <- e3g1_x_min_edt <- gedit(text="", width=5, container=e3g1)
  e3g1[4,2] <- e3g1_x_max_edt <- gedit(text="", width=5, container=e3g1)
  
  # FRAME 4 ###################################################################
  
  e4 <- gexpandgroup(text="X labels",
                     horizontal=FALSE,
                     container = f1)
  
  e4f1 <- gframe(text = "", horizontal = FALSE, container = e4) 

  e4g1 <- glayout(container = e4f1)
  
  e4g1[1,1] <- glabel(text="Text size (pts):", container=e4g1)
  e4g1[1,2] <- e4g1_size_edt <- gedit(text="8", width=4, container=e4g1)
  
  e4g1[1,3] <- glabel(text="Angle:", container=e4g1)
  e4g1[1,4] <- e4g1_angle_spb <- gspinbutton (from=0, to=360, by=1,
                                          value=0,
                                          container=e4g1) 
  
  e4g1[2,1] <- glabel(text="Justification (v/h):", container=e4g1)
  e4g1[2,2] <- e4g1_vjust_spb <- gspinbutton (from=0, to=1, by=0.1,
                                          value=0.5,
                                          container=e4g1)
  
  e4g1[2,3] <- e4g1_hjust_spb <- gspinbutton (from=0, to=1, by=0.1,
                                          value=0.5,
                                          container=e4g1)
  
  
  # FUNCTIONS #################################################################
  
  .plotDrop <- function(){
    
    logModel <- svalue(log_model)
    
    # Get values.
    val_p_dropout <- svalue(e1f1_risk_spn)
    val_predint <- svalue(e1f2_conf_spn) 
    val_predline <- svalue(e1f1_t_linetype_drp)
    val_predcol <- svalue(e1f1_t_linecolor_drp)
    
    val_titles <- svalue(f1_titles_chk)
    val_title <- svalue(f1_title_edt)
    val_xtitle <- svalue(f1_x_title_edt)
    val_ytitle <- svalue(f1_y_title_edt)
    val_column <- svalue(f1_column_opt, index=TRUE)
    val_sex <- svalue(f1_sex_chk)
    val_shape <- as.numeric(svalue(e2g1_shape_spb))
    val_alpha <- as.numeric(svalue(e2g1_alpha_spb))
    val_jitterh <- as.numeric(svalue(e2g1_jitterh_edt))
    val_jitterv <- as.numeric(svalue(e2g1_jitterv_edt))
    val_xmin <- as.numeric(svalue(e3g1_x_min_edt))
    val_xmax <- as.numeric(svalue(e3g1_x_max_edt))
    val_ymin <- as.numeric(svalue(e3g1_y_min_edt))
    val_ymax <- as.numeric(svalue(e3g1_y_max_edt))
    val_angle <- as.numeric(svalue(e4g1_angle_spb))
    val_vjust <- as.numeric(svalue(e4g1_vjust_spb))
    val_hjust <- as.numeric(svalue(e4g1_hjust_spb))
    val_size <- as.numeric(svalue(e4g1_size_edt))
    val_model <- svalue(f1_printmodel_chk)
    val_points <- svalue(e2g1_plotpoints_chk)
    val_threshold <- svalue(e1f1_threshold_chk)
    val_threshold_print <- svalue(e1f1_print_chk)
    val_prediction_interval <- svalue(e1f2_mark_interval_chk)
    val_prediction_print <- svalue(e1f2_print_interval_chk)
    val_interval_col <- svalue(e1f2_interval_drp)
    val_interval_alpha <- svalue(e1f2_interval_spb)
    val_h <- svalue(f1_h_chk)
    
    # Calculate values.
    val_pi_alpha <- 1 - val_predint
    
    if(debug){
      print("val_title")
      print(val_title)
      print("val_xtitle")
      print(val_xtitle)
      print("val_ytitle")
      print(val_ytitle)
      print("val_column")
      print(val_column)
      print("val_sex")
      print(val_sex)
      print("val_shape")
      print(val_shape)
      print("val_alpha")
      print(val_alpha)
      print("val_xmin")
      print(val_xmin)
      print("val_xmax")
      print(val_xmax)
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
      print("Before cleaning:")
      print(str(.gData))
    }

    # MODEL ###################################################################

    # Make copy of selected data frame (to allow re-plotting).
    obsData <- .gData

    # Get data for the selected analysis.
    if(val_column == 1){

      obsData$Dep <- .gData$MethodX
      obsData$Exp <- .gData$Height
      
    } else if(val_column == 2){
      
      obsData$Dep <- .gData$Method1
      obsData$Exp <- .gData$Height
      
    } else if(val_column == 3){
      
      obsData$Dep <- .gData$Method2
      obsData$Exp <- .gData$Height
      
    } else if(val_column == 4){
      
      obsData$Dep <- .gData$MethodL
      obsData$Exp <- .gData$MethodL.Ph
      
    }
    
    if(val_h){
      # Use average peak height.
      obsData$Exp <- .gData$H
    }
    
    # Clean -------------------------------------------------------------------

    message("Model drop-out for dataset with:")
    message(paste(nrow(obsData), " rows.", sep=""))
    
    # Remove homozygous loci
    if("Heterozygous" %in% names(obsData)){
      n0 <- nrow(obsData)
      obsData <- obsData[obsData$Heterozygous == 1, ]
      n1 <- nrow(obsData)
      message(paste(n1, " rows after removing ", n0-n1, " homozygous rows.", sep=""))
    }
    
    # Remove locus droput.
    if("Dep" %in% names(obsData)){
      n0 <- nrow(obsData)
      obsData <- obsData[obsData$Dep != 2, ]
      n1 <- nrow(obsData)
      message(paste(n1, " rows after removing ", n0-n1, " locus drop-out rows.", sep=""))
    }

    # Remove sex markers.
    if(val_sex){
      n0 <- nrow(obsData)
      sexMarkers <- getKit(kit=svalue(kit_drp), what="Sex.Marker")
      for(m in seq(along=sexMarkers)){
        obsData <- obsData[obsData$Marker != sexMarkers[m], ]
      }
      n1 <- nrow(obsData)
      message(paste(n1, " rows after removing ", n0-n1, " sex marker rows.", sep=""))
    }

    # Remove NA Explanatory.
    if(any(is.na(obsData$Exp))){
      n0 <- nrow(obsData)
      obsData <- obsData[!is.na(obsData$Exp), ]
      n1 <- nrow(obsData)
      message(paste(n1, " rows after removing ", n0-n1, " NA rows in explanatory column.", sep=""))
    }
  
    # Remove NA Dependent.
    if(any(is.na(obsData$Dep))){
      n0 <- nrow(obsData)
      obsData <- obsData[!is.na(obsData$Dep), ]
      n1 <- nrow(obsData)
      message(paste(n1, " rows after removing ", n0-n1, " NA rows in dependent column.", sep=""))
    }

    message(paste(nrow(obsData), " rows in total for analysis.", sep=""))
    
    if(debug){
      print("After cleaning:")
      print(str(obsData))
      print("NA in Exp/Dep:")
      print(any(is.na(obsData$Exp)))
      print(any(is.na(obsData$Dep)))
    }
    
    # Model -------------------------------------------------------------------

    # Build prediction range for smoother curve.
    val_pred_xmin <- min(obsData$Exp)
    val_pred_xmax <- max(obsData$Exp)
    xplot <- seq(val_pred_xmin, val_pred_xmax)
    predRange <- data.frame(Exp=xplot)
    
    # Create data for modelling.
    modData <- obsData
    # Convert to log values.
    if(logModel){
      modData$Exp <- log(obsData$Exp)
      predRange$Exp <- log(predRange$Exp)
    }

    # Perform logistic regression on the selected column.
    dropoutModel <- glm(Dep~Exp, family=binomial("logit"), data=modData)
    sumfit <- summary(dropoutModel)
    
    # Calculate model score.
    hosOk <- FALSE
    if(requireNamespace("ResourceSelection", quietly = TRUE)){
      #p-value <0.05 rejects the model.
      hos <- ResourceSelection::hoslem.test(dropoutModel$y, fitted(dropoutModel)) 
      hosOk <- TRUE
    }
        
    # Titles.
    if(val_titles){
      mainTitle <- val_title
      xTitle <- val_xtitle
      yTitle <- val_ytitle
    } else {
      if(val_h){
        mainTitle <- "Drop-out probability as a function of average peak height"
        xTitle <- "Average peak height 'H', (RFU)"
      } else {
        mainTitle <- "Drop-out probability as a function of present-allele height"
        xTitle <- "Peak height, (RFU)"
      }
      yTitle <- "Drop-out probability, P(D)"
    }
    
    # Extract model parameters.
    b0 <- sumfit$coefficients[1]
    b1 <- sumfit$coefficients[2]
    
    if(debug){
      print("Model summary:")
      print(sumfit)
    }
    
    
    if(debug){
      print("b0")
      print(b0)
      print("b1")
      print(b1)
    }
    
    # Calculate probabilities for the prediction range.
    ypred <-  predict(dropoutModel, predRange, type="link", se.fit=TRUE)

    # Calculate the prediction interval.
    ylower <- plogis(ypred$fit - qnorm(1-val_pi_alpha/2)*ypred$se)  # Lower confidence limit.
    yupper <- plogis(ypred$fit + qnorm(1-val_pi_alpha/2)*ypred$se)  # Upper confidence limit.

    # Calculate conservative prediction curve.
    yconservative <- plogis(ypred$fit + qnorm(1-val_pi_alpha)*ypred$se)

    # Calculate y values for plot.
    yplot <- plogis(ypred$fit)

    # Create legend text.
    legendModel <- paste("Model parameters: \u03B20=", round(b0,3),
                         ", \u03B21=", round(b1,3), sep="")
    
    if(hosOk){
      # Add Hosmer-Lemeshow test.
      legendModel <- paste(legendModel,
                           "\nHosmer-Lemeshow test: p = ", round(hos$p.value, 4),
                           sep="")
    } else {
      message("Package 'ResourceSelection' is required for Hosmer-Lemeshow test.")
    }
    
    # Save prediction in a dataframe.
    predictionDf <- data.frame(Exp=xplot, Prob=yplot, yupper=yupper, ylower=ylower)

    if(debug){
      print("predictionDf:")
      print(head(predictionDf, 100))
      print(str(predictionDf))
    }
    
    # Calculate dropout threshold T.
    if(logModel){
      drop_py <- log(val_p_dropout) - log(1 - val_p_dropout)
      t_dropout <- exp((drop_py - b0) / b1)
    }  else {
      t_dropout <- (log(val_p_dropout / (1 - val_p_dropout)) - b0) / b1
    }

    if(debug){
      print(paste("t_dropout =", t_dropout))
    }
    
    if(!logModel){
      if(t_dropout < 0){
        if(debug){
          print("t_dropout < 0 -> NA")
        }
        t_dropout <- NA # Can't handle negative values.
      }
    }

    # Calculate conservative threshold at P(D).
    t_dropout_cons <- xplot[min(which(yconservative < val_p_dropout))]
    # rfu1 <- predrange[min(which( yupper <val_p_dropout))] #too conservative?
    
    # PLOT ####################################################################
    
    if (!is.na(predictionDf) && !is.null(predictionDf)){
      
      # Plotting global dropout probability.
      gp <- ggplot(data = predictionDf, 
                   aes_string(y = "Prob", x="Exp")) + geom_line() 

      # Plot observed data points (heterozygotes).
      if(val_points){
        
          gp <- gp + geom_point(data=obsData, aes_string(x="Exp", y="Dep"),
                                shape=val_shape, alpha=val_alpha, 
                                position=position_jitter(width=val_jitterh,
                                                         height=val_jitterv)) 
      }
      
      # Prediction interval.
      if(val_prediction_interval){
        
        gp <- gp + geom_ribbon(data = predictionDf,
                               aes_string(y = "Prob", ymin = "ylower", ymax = "yupper"),
                               fill = val_interval_col,
                               alpha = val_interval_alpha) 
        
      }
 
      # Dropout threshold.      
      if(val_threshold){
 
        # Initiate.
        thresholdLegend <- ""
        
        # Create threshold label.
        thresholdLegend <- paste("P(dropout|T=",
                                 round(t_dropout, 0),
                                 ")=",
                                 round(val_p_dropout, 3),
                                 sep="")

        # Add prediction interval.
        if(val_prediction_print){
          thresholdLegend <- paste(thresholdLegend,
                                   "\n P(dropout>", round(val_p_dropout, 3),
                                   "|T=",
                                   round(t_dropout_cons,0),
                                   ")<",
                                   val_pi_alpha * 100, "%",
                                   sep="")
          
        }

        # Make data frame.
        if(is.na(t_dropout)){
          t_height <- 0
        } else {
          t_height <- t_dropout
        }
        thresholdLabel <- data.frame(Exp = t_height,
                                     Prob = val_p_dropout,
                                     label = thresholdLegend)
        
        if(debug){
          print("thresholdLabel")
          print(thresholdLabel)
        }
        
        if(!is.na(t_dropout)){

          if(debug){
            print("Mark threshold")
          }
          
          # Horizontal threshold line.
          if(!is.na(val_xmin) && !is.na(val_xmax)){
            xtemp <- val_xmin
          } else {
            xtemp <- 0
          }
          
          # Add horizontal threshold line.
          gp <- gp + geom_segment(data=thresholdLabel,
                                  aes_string(x = xtemp, y = "Prob",
                                             xend = "Exp",
                                             yend = "Prob"),
                                  color = val_predcol,
                                  linetype=val_predline)

          if(debug){
            print("Horizontal line added")
          }
          
          # Vertical threshold line.
          if(!is.na(val_ymin) && !is.na(val_ymax)){
            ytemp <- val_ymin
          } else {
            ytemp <- 0
          }
          
          # Add vertical threshold line.
          gp <- gp + geom_segment(data=thresholdLabel,
                                  aes_string(x = "Exp", y = ytemp, 
                                             xend = "Exp",
                                             yend = "Prob"),
                                  color = val_predcol,
                                  linetype=val_predline)
          
          if(debug){
            print("Vertical line added")
          }
          
        }

        # Print threshold label.
        if(val_threshold_print){
          gp <- gp + geom_text(data = thresholdLabel,
                             aes_string(x = "Exp", y = "Prob", label = "label"),
                             hjust=0, vjust=0)
          if(debug){
            print("Threshold printed")
          }
        }
        
      }
 
      # Print dropout model.
      if(val_model){
        
        if(debug){
          print("Print model")
        }
        
        # Create data frame.
        modelLabel <- data.frame(Exp = t_dropout,
                                 Prob = val_p_dropout,
                                 label = legendModel,
                                 xmax= val_pred_xmax)
        # Add model text.
        gp <- gp + geom_text(data=modelLabel, aes_string(x = Inf, y = Inf, label = "label"),
                        hjust=1, vjust=1)
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
      if(debug){
        print(paste("Zoom plot xmin/xmax,ymin/ymax:",
                    paste(val_x,collapse="/"),
                    ",",
                    paste(val_y,collapse="/")))
        print(str(val_x))
        print(str(val_y))
      }
      # Zoom in without dropping observations.
      gp <- gp + coord_cartesian(xlim=val_x, ylim=val_y)

      
      if(debug){
        print("Apply theme and labels")
      }
      
      # Apply theme.
      gp <- gp + theme(axis.text.x=element_text(angle=val_angle,
                                                hjust=val_hjust,
                                                vjust=val_vjust,
                                                size=val_size))
 
      # Add titles and labels.
      gp <- gp + labs(title=mainTitle)
      gp <- gp + xlab(xTitle)
      gp <- gp + ylab(yTitle)
      
      if(debug){
        print("Plot")
      }
      
      # Plot.
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
  
  .checkColumns <- function(){
    
    val_col <- svalue(f1_column_opt, index=TRUE)
    val_h <- svalue(f1_h_chk)
    requiredCol <- NULL
    missingCol <- NULL
    
    
    if(!is.null(.gData)){
      # Enable button.
      enabled(f7_plot_drop_btn) <- TRUE
      svalue(f7_plot_drop_btn) <- "Plot predicted drop-out probability"

      # Check available modelling columns.
      if(val_h){
        requiredCol <- c("H")
        if(!all(requiredCol %in% colnames(.gData))){
          missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
        }
      }
        
      # Check available modelling columns and enable/select.
      if(val_col == 1){
        
        requiredCol <- c("MethodX")
        if(!all(requiredCol %in% colnames(.gData))){
          missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
        }
        
      } else if(val_col == 2){
        
        requiredCol <- c("Method1")
        if(!all(requiredCol %in% colnames(.gData))){
          missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
        }
        
      } else if(val_col == 3){
        
        requiredCol <- c("Method2")
        if(!all(requiredCol %in% colnames(.gData))){
          missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
        }
        
      } else if(val_col == 4){
        
        requiredCol <- c("MethodL", "MethodL.Ph")
        if(!all(requiredCol %in% colnames(.gData))){
          missingCol <- requiredCol[!requiredCol %in% colnames(.gData)]
        }
        
      } else {
        
        message <- paste("Selection not supported!")
        
        gmessage(message, title="Error",
                 icon = "error",
                 parent = w) 
        
        # Disable button.
        enabled(f7_plot_drop_btn) <- FALSE
      }
      
      if(!is.null(missingCol)){
        
        message <- paste("Dataset is ok for drop-out analysis.\n",
                        "However, additional columns are required for this analysis:\n",
                         paste(missingCol, collapse="\n"),
                         "\n\nPlease try modelling using another scoring method.",
                         sep="")
        
        gmessage(message, title="message",
                 icon = "info",
                 parent = w)
        
        # Disable button.
        enabled(f7_plot_drop_btn) <- FALSE
        
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
      if(exists(".strvalidator_modelDropout_gui_savegui", envir=env, inherits = FALSE)){
        svalue(savegui_chk) <- get(".strvalidator_modelDropout_gui_savegui", envir=env)
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
      if(exists(".strvalidator_modelDropout_gui_title", envir=env, inherits = FALSE)){
        svalue(f1_title_edt) <- get(".strvalidator_modelDropout_gui_title", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_title_chk", envir=env, inherits = FALSE)){
        svalue(f1_titles_chk) <- get(".strvalidator_modelDropout_gui_title_chk", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_x_title", envir=env, inherits = FALSE)){
        svalue(f1_x_title_edt) <- get(".strvalidator_modelDropout_gui_x_title", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_y_title", envir=env, inherits = FALSE)){
        svalue(f1_y_title_edt) <- get(".strvalidator_modelDropout_gui_y_title", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_sex", envir=env, inherits = FALSE)){
        svalue(f1_sex_chk) <- get(".strvalidator_modelDropout_gui_sex", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_column", envir=env, inherits = FALSE)){
        svalue(f1_column_opt) <- get(".strvalidator_modelDropout_gui_column", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_print_model", envir=env, inherits = FALSE)){
        svalue(f1_printmodel_chk) <- get(".strvalidator_modelDropout_gui_print_model", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_mark_threshold", envir=env, inherits = FALSE)){
        svalue(e1f1_threshold_chk) <- get(".strvalidator_modelDropout_gui_mark_threshold", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_risk", envir=env, inherits = FALSE)){
        svalue(e1f1_risk_spn) <- get(".strvalidator_modelDropout_gui_risk", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_print_threshold", envir=env, inherits = FALSE)){
        svalue(e1f1_print_chk) <- get(".strvalidator_modelDropout_gui_print_threshold", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_t_line", envir=env, inherits = FALSE)){
        svalue(e1f1_t_linetype_drp) <- get(".strvalidator_modelDropout_gui_t_line", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_t_color", envir=env, inherits = FALSE)){
        svalue(e1f1_t_linecolor_drp) <- get(".strvalidator_modelDropout_gui_t_color", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_print_interval", envir=env, inherits = FALSE)){
        svalue(e1f2_print_interval_chk) <- get(".strvalidator_modelDropout_gui_print_interval", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_mark_interval", envir=env, inherits = FALSE)){
        svalue(e1f2_mark_interval_chk) <- get(".strvalidator_modelDropout_gui_mark_interval", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_interval_alpha", envir=env, inherits = FALSE)){
        svalue(e1f2_interval_spb) <- get(".strvalidator_modelDropout_gui_interval_alpha", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_interval_color", envir=env, inherits = FALSE)){
        svalue(e1f2_interval_drp) <- get(".strvalidator_modelDropout_gui_interval_color", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_points_plot", envir=env, inherits = FALSE)){
        svalue(e2g1_plotpoints_chk) <- get(".strvalidator_modelDropout_gui_points_plot", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_points_shape", envir=env, inherits = FALSE)){
        svalue(e2g1_shape_spb) <- get(".strvalidator_modelDropout_gui_points_shape", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_points_alpha", envir=env, inherits = FALSE)){
        svalue(e2g1_alpha_spb) <- get(".strvalidator_modelDropout_gui_points_alpha", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_points_jitterh", envir=env, inherits = FALSE)){
        svalue(e2g1_jitterh_edt) <- get(".strvalidator_modelDropout_gui_points_jitterh", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_points_jitterv", envir=env, inherits = FALSE)){
        svalue(e2g1_jitterv_edt) <- get(".strvalidator_modelDropout_gui_points_jitterv", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_axes_y_min", envir=env, inherits = FALSE)){
        svalue(e3g1_y_min_edt) <- get(".strvalidator_modelDropout_gui_axes_y_min", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_axes_y_max", envir=env, inherits = FALSE)){
        svalue(e3g1_y_max_edt) <- get(".strvalidator_modelDropout_gui_axes_y_max", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_axes_x_min", envir=env, inherits = FALSE)){
        svalue(e3g1_x_min_edt) <- get(".strvalidator_modelDropout_gui_axes_x_min", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_axes_x_max", envir=env, inherits = FALSE)){
        svalue(e3g1_x_max_edt) <- get(".strvalidator_modelDropout_gui_axes_x_max", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_xlabel_size", envir=env, inherits = FALSE)){
        svalue(e4g1_size_edt) <- get(".strvalidator_modelDropout_gui_xlabel_size", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_xlabel_angle", envir=env, inherits = FALSE)){
        svalue(e4g1_angle_spb) <- get(".strvalidator_modelDropout_gui_xlabel_angle", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_xlabel_justh", envir=env, inherits = FALSE)){
        svalue(e4g1_hjust_spb) <- get(".strvalidator_modelDropout_gui_xlabel_justh", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_xlabel_justv", envir=env, inherits = FALSE)){
        svalue(e4g1_vjust_spb) <- get(".strvalidator_modelDropout_gui_xlabel_justv", envir=env)
      }
      if(exists(".strvalidator_modelDropout_gui_h", envir=env, inherits = FALSE)){
        svalue(f1_h_chk) <- get(".strvalidator_modelDropout_gui_h", envir=env)
      }
      
      if(debug){
        print("Saved settings loaded!")
      }
    }
    
  }
  
  .saveSettings <- function(){
    
    # Then save settings if true.
    if(svalue(savegui_chk)){
      
      assign(x=".strvalidator_modelDropout_gui_savegui", value=svalue(savegui_chk), envir=env)
      assign(x=".strvalidator_modelDropout_gui_title", value=svalue(f1_title_edt), envir=env)
      assign(x=".strvalidator_modelDropout_gui_title_chk", value=svalue(f1_titles_chk), envir=env)
      assign(x=".strvalidator_modelDropout_gui_x_title", value=svalue(f1_x_title_edt), envir=env)
      assign(x=".strvalidator_modelDropout_gui_y_title", value=svalue(f1_y_title_edt), envir=env)
      assign(x=".strvalidator_modelDropout_gui_sex", value=svalue(f1_sex_chk), envir=env)
      assign(x=".strvalidator_modelDropout_gui_column", value=svalue(f1_column_opt), envir=env)
      assign(x=".strvalidator_modelDropout_gui_print_model", value=svalue(f1_printmodel_chk), envir=env)
      assign(x=".strvalidator_modelDropout_gui_mark_threshold", value=svalue(e1f1_threshold_chk), envir=env)
      assign(x=".strvalidator_modelDropout_gui_print_threshold", value=svalue(e1f1_print_chk), envir=env)
      assign(x=".strvalidator_modelDropout_gui_risk", value=svalue(e1f1_risk_spn), envir=env)
      assign(x=".strvalidator_modelDropout_gui_t_line", value=svalue(e1f1_t_linetype_drp), envir=env)
      assign(x=".strvalidator_modelDropout_gui_t_color", value=svalue(e1f1_t_linecolor_drp), envir=env)
      assign(x=".strvalidator_modelDropout_gui_print_interval", value=svalue(e1f2_print_interval_chk), envir=env)
      assign(x=".strvalidator_modelDropout_gui_mark_interval", value=svalue(e1f2_mark_interval_chk), envir=env)
      assign(x=".strvalidator_modelDropout_gui_interval_alpha", value=svalue(e1f2_interval_spb), envir=env)
      assign(x=".strvalidator_modelDropout_gui_interval_color", value=svalue(e1f2_interval_drp), envir=env)
      assign(x=".strvalidator_modelDropout_gui_points_plot", value=svalue(e2g1_plotpoints_chk), envir=env)
      assign(x=".strvalidator_modelDropout_gui_points_shape", value=svalue(e2g1_shape_spb), envir=env)
      assign(x=".strvalidator_modelDropout_gui_points_alpha", value=svalue(e2g1_alpha_spb), envir=env)
      assign(x=".strvalidator_modelDropout_gui_points_jitterh", value=svalue(e2g1_jitterh_edt), envir=env)
      assign(x=".strvalidator_modelDropout_gui_points_jitterv", value=svalue(e2g1_jitterv_edt), envir=env)
      assign(x=".strvalidator_modelDropout_gui_axes_y_min", value=svalue(e3g1_y_min_edt), envir=env)
      assign(x=".strvalidator_modelDropout_gui_axes_y_max", value=svalue(e3g1_y_max_edt), envir=env)
      assign(x=".strvalidator_modelDropout_gui_axes_x_min", value=svalue(e3g1_x_min_edt), envir=env)
      assign(x=".strvalidator_modelDropout_gui_axes_x_max", value=svalue(e3g1_x_max_edt), envir=env)
      assign(x=".strvalidator_modelDropout_gui_xlabel_size", value=svalue(e4g1_size_edt), envir=env)
      assign(x=".strvalidator_modelDropout_gui_xlabel_angle", value=svalue(e4g1_angle_spb), envir=env)
      assign(x=".strvalidator_modelDropout_gui_xlabel_justh", value=svalue(e4g1_hjust_spb), envir=env)
      assign(x=".strvalidator_modelDropout_gui_xlabel_justv", value=svalue(e4g1_vjust_spb), envir=env)
      assign(x=".strvalidator_modelDropout_gui_h", value=svalue(f1_h_chk), envir=env)
      
    } else { # or remove all saved values if false.
      
      if(exists(".strvalidator_modelDropout_gui_savegui", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_savegui", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_title", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_title_chk", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_title_chk", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_x_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_x_title", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_y_title", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_y_title", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_sex", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_sex", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_column", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_column", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_print_model", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_print_model", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_mark_threshold", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_mark_threshold", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_print_threshold", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_print_threshold", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_risk", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_risk", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_t_line", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_t_line", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_t_color", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_t_color", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_print_interval", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_print_interval", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_mark_interval", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_mark_interval", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_interval_alpha", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_interval_alpha", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_interval_color", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_interval_color", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_points_plot", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_points_plot", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_points_shape", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_points_shape", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_points_alpha", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_points_alpha", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_points_jitterh", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_points_jitterh", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_points_jitterv", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_points_jitterv", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_axes_y_min", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_axes_y_min", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_axes_y_max", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_axes_y_max", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_axes_x_min", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_axes_x_min", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_axes_x_max", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_axes_x_max", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_xlabel_size", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_xlabel_size", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_xlabel_angle", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_xlabel_angle", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_xlabel_justh", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_xlabel_justh", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_xlabel_justv", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_xlabel_justv", envir = env)
      }
      if(exists(".strvalidator_modelDropout_gui_h", envir=env, inherits = FALSE)){
        remove(".strvalidator_modelDropout_gui_h", envir = env)
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
