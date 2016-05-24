#' Scatmat Tour Plotting
#' Plots the Scatmat Tour in tab g8
#'
#' Sets up the interface for the scatterplot matrix tour
#'
#' @keywords internal
#' @author Bei Huang\email{beihuang@@iastate.edu}, Di Cook \email{dicook@@iastate.edu}, and Hadley Wickham \email{hadley@@rice.edu}
# =============================== Gui_scatmat ================================
.interface_scatmat = function(g8,data, w){

  # =============== Function: update_tour_scatmat ==================
  tour <- NULL
  tour_anim <- NULL
  update_tour_scatmat <- function(...) {
    tour <<- .create_mat_tour(data,
      var_selected = svalue(Variables_scatmat),
      cat_selected = svalue(Class_scatmat), 
      projdim_selected = svalue(Projections_scatmat),
      tour_type = svalue(TourType_scatmat),
      guided_type = svalue(GuidedType_scatmat),
      lambda = svalue(LambdaValue_scatmat),
      aps = svalue(sl_scatmat)
    )
    tour_anim <<- with(tour, new_tour(data, tour_path))

    tour$display$init(tour$data)
  # tour$display$render_frame()

    TRUE
  }
  # --------------------- End of update_tour_scatmat ----------------
  
  # ================= Function: draw_frame_scatmat ==================
  draw_frame_scatmat <- function(...) {
    # if there's no tour, don't draw anything
    if (is.null(tour)) return(FALSE)

    tour_step <- tour_anim(svalue(sl_scatmat) / 33)
    if (is.null(tour_step$proj)) return(FALSE)

    if (find_platform()$os == "win") {
      tour$display$render_frame()
    } else {
      tour$display$render_transition()
    }
    with(tour_step, tour$display$render_data(tour$data, proj, target))
    Sys.sleep(1/33)

    TRUE
  }
  # ---------------------- End of draw_frame_scatmat -----------------
  
  num <- sapply(data, is.numeric)
  # ==================Controls==========================
  vbox_scatmat <- glayout(container = g8)
    # Variable selection column
    vbox_scatmat[1, 1, anchor = c(-1, 0)] <- "Variable Selection"
    vbox_scatmat[2, 1] <- Variables_scatmat <- gcheckboxgroup(names(data[num]),
      checked = TRUE, horizontal = FALSE)
    tooltip(Variables_scatmat) <- "Select variables to display in the nD Tour."

    vbox_scatmat[3, 1, anchor = c(-1, 0)] <- "Class Selection"
    vbox_scatmat [4, 1, anchor = c(-1, 0)] <- Class_scatmat <- gtable(names(data)[!num], 
      multiple = TRUE)
    tooltip(Class_scatmat) <- "Select a class variable to classify the points."

    # Tour selection column
    vbox_scatmat[1, 2, anchor=c(-1, 0)] <- "Tour Type"
    tour_types <- c("Grand", "Little", "Local", "Guided")
    vbox_scatmat[2, 2] <- TourType_scatmat <- gradio(tour_types)
    tooltip(TourType_scatmat) <- "Select a nD Tour type."


  #Guided indices selection
  vbox_scatmat[3, 2, anchor=c(-1, 0)] <- "Guided indices"
  IntIndex <-c("holes","cmass","lda_pp","pda_pp")
  vbox_scatmat [4, 2, anchor=c(-1,-1)] <-  GuidedType_scatmat <- gdroplist(IntIndex)
  tooltip(GuidedType_scatmat) <- "Select an index type for guided tour."

  # Lambda selection
  vbox_scatmat[3, 3, anchor=c(-1, 0)] <-"Lambda"
  vbox_scatmat[4, 3] <- LambdaValue_scatmat <- gslider(from=0, to = 1, by = 0.01,value=0.02)
  tooltip(LambdaValue_scatmat) <- "Select lambda's value to calculate pda index."

    # dimension control
    vbox_scatmat[1, 3, anchor = c(-1, 0)] <- "Projection Dimension"
    projections <- c(2:length(data[num]))
    vbox_scatmat[2, 3, anchor = c(-1, 0)] <- Projections_scatmat <- gradio(projections)
    tooltip(Projections_scatmat) <- "Select pojection dimension for the nD tour."	


    # speed and pause
    vbox_scatmat[5,1, anchor = c(-1, 0)] <- "Speed"
    vbox_scatmat[6,1, expand = TRUE] <- sl_scatmat <- gslider(from = 0, to = 5, by = 0.1, value = 1)
    tooltip(sl_scatmat) <- "Drag to set the speed of the nD Tour."

    vbox_scatmat[6, 2] <- chk_pause_scatmat <- gcheckbox("Pause",
      handler = function(h, ...) pause_scatmat(svalue(h$obj)))
      tooltip(chk_pause_scatmat) <- "Click here to pause or continue the nD Tour."

    # buttons control
    anim_id <- NULL
    pause_scatmat <- function(paused) {
      svalue(chk_pause_scatmat) <- paused
      if (paused) {
        gtkIdleRemove(anim_id)
        anim_id <- NULL
      } else {
        if (!is.null(anim_id)) gtkIdleRemove(anim_id)
        anim_id <<- gIdleAdd(draw_frame_scatmat)
      }
    }
    buttonGroup_scatmat <- ggroup(horizontal = FALSE, container = vbox_scatmat)

    # addSpace(buttonGroup,10)
    button1_scatmat<- gbutton("Apply", container = buttonGroup_scatmat, handler = function(...){
      print("apply from gui_scatmat")
      pause_scatmat(FALSE)
      update_tour_scatmat()
    })
    tooltip(button1_scatmat) <- "Click here to update the options."

    # addSpace(buttonGroup,10)
    button2_scatmat<- gbutton("Quit",container = buttonGroup_scatmat, handler = function(...) {
      pause_scatmat(TRUE)
      dispose(w)
    })
    tooltip(button2_scatmat) <- "Click here to close this window."                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         

    vbox_scatmat[5:6, 3, anchor = c(0, 1)] <- buttonGroup_scatmat

}
# ============================ End of Gui_scatmat ==============
