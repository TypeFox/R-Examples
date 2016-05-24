#' Density Tour Plotting
#' Plots the Density Tour in tab g2
#'
#' Sets up the interface for the density tour
#'
#' @keywords internal
#' @author Bei Huang\email{beihuang@@iastate.edu}, Di Cook \email{dicook@@iastate.edu}, and Hadley Wickham \email{hadley@@rice.edu}
# =============================== Gui_density ==========================================
.interface_density = function(g2,data, w){

  # =============== Function: update_tour_density ===============
  tour <- NULL
  tour_anim <- NULL
  update_tour_density <- function(...) {
    tour <<- .create_1d_tour(data,
      var_selected = svalue(Variables_density),
      cat_selected = svalue(Class_density), 
      method_selected = svalue(MethodType),
      center_selected = svalue(CenterType),
      tour_type = svalue(TourType_density),
      guided_type = svalue(GuidedType_density),
      lambda = svalue(LambdaValue_density),
      aps = svalue(sl_density)
    )
    tour_anim <<- with(tour, new_tour(data, tour_path, start=basis_random(ncol(data), 1)))

    tour$display$init(tour$data)
    tour$display$render_frame()

    TRUE
  }
  # ----------------- End of update_tour_density -----------------


  # =============== Function: draw_frame_density ==================
  draw_frame_density <- function(...) {
    # if there's no tour, don't draw anything
    if (is.null(tour)) return(FALSE)

    tour_step <- tour_anim(svalue(sl_density) / 33)
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
  # -------------------- End of draw_frame_density -----------------

  num <- sapply(data, is.numeric)

  # ==================Controls==========================
  vbox_density <- glayout(container = g2)

  # Variable selection column
  vbox_density[1, 1, anchor = c(-1, 0)] <- "Variable Selection"
  vbox_density[2, 1] <- Variables_density <- gcheckboxgroup(names(data[num]),
    checked = TRUE, horizontal = FALSE)
  tooltip(Variables_density) <- "Select variables to display in the 1D Tour."

  vbox_density[3, 1, anchor = c(-1, 0)] <- "Class Selection"
  vbox_density[4, 1, anchor = c(-1, 0)] <- Class_density <- gtable(names(data)[!num], 
    multiple = TRUE)
  tooltip(Class_density) <- "Select a class variable to classify the data."

  # Tour selection column
  vbox_density[1, 2, anchor=c(-1, 0)] <- "Tour Type"
  tour_types <- c("Grand", "Little", "Local", "Guided")
  vbox_density[2, 2] <- TourType_density <- gradio(tour_types)
  tooltip(TourType_density) <- "Select a 1D Tour type."

  vbox_density[3, 2, anchor=c(-1, 0)] <- "Guided indices"
  IntIndex <-c("holes","cmass","lda_pp","pda_pp")
  vbox_density[4, 2, anchor=c(-1,-1)] <-  GuidedType_density <- gdroplist(IntIndex)
  tooltip(GuidedType_density) <- "Select an index type for guided tour."

  vbox_density[3,3, anchor=c(-1, 0)] <-"Lambda"
  vbox_density[4,3] <- LambdaValue_density <- gslider(from=0, to = 1, by = 0.01,value=0.02)
  tooltip(LambdaValue_density) <- "Select lambda's value to calculate pda index."

  # speed and pause
  vbox_density[5,1, anchor = c(-1, 0)] <- "Speed"
  vbox_density[6,1, expand = TRUE] <- sl_density <- gslider(from = 0, to = 5, by = 0.1, value = 1)
  tooltip(sl_density) <- "Drag to set the speed of the 1D Tour."

  vbox_density[5, 3] <- chk_pause_density <- gcheckbox("Pause",
    handler = function(h, ...) pause_density(svalue(h$obj)))
  tooltip(chk_pause_density) <- "Click here to pause or continue the 1D Tour."

  # method control
  vbox_density[1, 3, anchor = c(-1, 0)] <- "Method Type"
  method_types <- c("density","hist","ash")
  vbox_density[2, 3, anchor = c(-1, 0)] <- MethodType <- gradio(method_types)
  tooltip(MethodType) <- "Select a display method for the 1D tour."

  # center control
  vbox_density[5,2, anchor=c(-1,0)] <- "Center or Not"
  center_types <- c("TRUE", "FALSE")
  vbox_density[6,2, anchor=c(-1,0)] <- CenterType <- gradio(center_types)
  tooltip(CenterType) <- "Choose to center the display or let it wander."

  # buttons control
  anim_id <- NULL
  pause_density <- function(paused) {
    svalue(chk_pause_density) <- paused
    if (paused) {
      gtkIdleRemove(anim_id)
      anim_id <<- NULL
    } else {
      if (!is.null(anim_id)) gtkIdleRemove(anim_id)
      anim_id <<- gIdleAdd(draw_frame_density)
    }
  }

  buttonGroup_density <- ggroup(horizontal = FALSE, container = vbox_density)

  # addSpace(buttonGroup_density,10)
  button1_density <-gbutton("Apply", container = buttonGroup_density, handler = function(...) {
    print("apply from gui_density")
    opar <- par(mfrow = c(1,1))
    pause_density(FALSE)
    update_tour_density()
  })
  tooltip(button1_density ) <- "Click here to update the options."

  # addSpace(buttonGroup,10)
  button2_density <-gbutton("Quit",container = buttonGroup_density, handler = function(...) {
    pause_density(TRUE)
    dispose(w)
  })
  tooltip(button2_density ) <- "Click here to close this window."

  vbox_density[6, 3, anchor = c(0, 1)] <- buttonGroup_density

}
# --------------------------- End of Gui_density -----------------------------------
