#' Stereo Tour Plotting
#' Plots the Stereo Tour in tab g6
#'
#' Sets up the interface for the stereo tour
#'
#' @keywords internal
#' @author Bei Huang\email{beihuang@@iastate.edu}, Di Cook \email{dicook@@iastate.edu}, and Hadley Wickham \email{hadley@@rice.edu} 
# =============================== Gui_stereo ================================
.interface_stereo = function(g6, data, w){

  # ================= Function: update_tour_stereo ==================
  tour <- NULL
  tour_anim <- NULL
  update_tour_stereo <- function(...) {
    tour <<- .create_stereo_tour(data,
      var_selected = svalue(Variables_stereo),
      cat_selected = svalue(Class_stereo),
      tour_type = svalue(TourType_stereo),
      guided_type = svalue(GuidedType_stereo),
      lambda = svalue(LambdaValue_stereo),
      aps = svalue(sl_stereo)
    )
    tour_anim <<- with(tour, new_tour(data, tour_path))

    tour$display$init(tour$data)
    tour$display$render_frame()

    TRUE
  }
  # -------------------- End of update_tour_stereo -----------------
  
  # ================= Function: draw_frame_stereo ==================
  draw_frame_stereo <- function(...) {
    # if there's no tour, don't draw anything
    if (is.null(tour)) return(FALSE)

    tour_step <- tour_anim(svalue(sl_stereo) / 33)
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
  # -------------------- End of draw_frame_stereo -----------------

  num <- sapply(data, is.numeric)
# ==================Controls==========================
vbox_stereo <- glayout(container = g6)
  # Variable selection column
  vbox_stereo[1, 1, anchor = c(-1, 0)] <- "Variable Selection"
  vbox_stereo[2, 1] <- Variables_stereo <- gcheckboxgroup(names(data[num]),
    checked = TRUE, horizontal = FALSE)
  tooltip(Variables_stereo) <- "Select variables to display in the 3D Tour."

  vbox_stereo[3, 1, anchor = c(-1, 0)] <- "Class Selection"
  vbox_stereo[4, 1, anchor = c(-1, 0)] <- Class_stereo <- gtable(names(data)[!num], 
    multiple = TRUE)
  tooltip(Class_stereo) <- "Select a class variable to classify the points."

  # Tour selection column
  vbox_stereo[1, 2, anchor=c(-1, 0)] <- "Tour Type"
  tour_types <- c("Grand", "Little", "Local","Guided")
  vbox_stereo[2, 2] <- TourType_stereo <- gradio(tour_types)
  tooltip(TourType_stereo) <- "Select a 3D Tour type."

  vbox_stereo[3, 2, anchor=c(-1, 0)] <- "Guided indices"
  IntIndex <-c("holes","cmass","lda_pp","pda_pp")
  vbox_stereo[4, 2, anchor=c(-1,-1)] <-  GuidedType_stereo <- gdroplist(IntIndex)
  tooltip(GuidedType_stereo) <- "Select an index type for guided tour."

  vbox_stereo[3,3, anchor=c(-1, 0)] <-"Lambda"
  vbox_stereo[4,3] <- LambdaValue_stereo <- gslider(from=0, to = 1, by = 0.01,value=0.02)
  tooltip(LambdaValue_stereo) <- "Select lambda's value to calculate pda index."

  # speed and pause
  vbox_stereo[5,1, anchor = c(-1, 0)] <- "Speed"
  vbox_stereo[6,1, expand = TRUE] <- sl_stereo <- gslider(from = 0, to = 5, by = 0.1, value = 1)
  tooltip(sl_stereo) <- "Drag to set the speed of the 3D Tour."

  vbox_stereo[6, 2] <- chk_pause_stereo <- gcheckbox("Pause",
    handler = function(h, ...) pause_stereo(svalue(h$obj)))
  tooltip(chk_pause_stereo) <- "Click here to pause or continue the 3D Tour."

  # buttons control
  anim_id <- NULL
  pause_stereo <- function(paused) {
    svalue(chk_pause_stereo) <- paused
    if (paused) {
      gtkIdleRemove(anim_id)
      anim_id <- NULL
    } else {
      if (!is.null(anim_id)) gtkIdleRemove(anim_id)
      anim_id <<- gIdleAdd(draw_frame_stereo)
    }
  }
  buttonGroup_stereo <- ggroup(horizontal = FALSE, container = vbox_stereo)

  # addSpace(buttonGroup,10)
  button1_stereo<- gbutton("Apply", container = buttonGroup_stereo, handler = function(...) {
    print("apply from gui_stereo")
    opar <- par(mfrow = c(1,1))
    pause_stereo(FALSE)
    update_tour_stereo()
  })
  tooltip(button1_stereo) <- "Click here to update the options."

  # addSpace(buttonGroup,10)
  button2_stereo<- gbutton("Quit",container = buttonGroup_stereo, handler = function(...) {
    pause_stereo(TRUE)
    dispose(w)
  })
  tooltip(button2_stereo) <- "Click here to close this window."

  vbox_stereo[5:6, 3, anchor = c(0, -1)] <- buttonGroup_stereo
}
# ----------------------------- End of Gui_stereo ------------------------------
