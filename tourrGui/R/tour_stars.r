#' Stars Tour Plotting
#' Plots the Stars Tour in tab g4
#'
#' Sets up the interface for the stars tour
#'
#' @keywords internal
#' @author Bei Huang\email{beihuang@@iastate.edu}, Di Cook \email{dicook@@iastate.edu}, and Hadley Wickham \email{hadley@@rice.edu} 
# =============================== Gui_stars==============================
.interface_stars = function(g4, data, w){
  # =============== Function: update_tour_stars ==================
  tour <- NULL
  tour_anim <- NULL
  update_tour_stars <- function(...) {
    tour <<- .create_stars_tour(data,
      var_selected = svalue(Variables_stars),
      cat_selected = svalue(Class_stars), 
      dim_selected = svalue(Dimensions_stars),
      nstar_selected = svalue(Star_stars),
      tour_type = svalue(TourType_stars),
      guided_type = svalue(GuidedType_stars),
      lambda = svalue(LambdaValue_stars),
      aps = svalue(sl_stars)
    )
    tour_anim <<- with(tour, new_tour(data, tour_path))

    tour$display$init(tour$data)
    tour$display$render_frame()

    TRUE
  }
  # --------------------- End of update_tour_stars ----------------
  
  # ================= Function: draw_frame_stars ==================
  draw_frame_stars <- function(...) {
    # if there's no tour, don't draw anything
    if (is.null(tour)) return(FALSE)

    tour_step <- tour_anim(svalue(sl_stars) / 33)
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
  # -------------------- End of draw_frame_stars -----------------
  
  num <- sapply(data, is.numeric)
  # ================== Controls ==========================
  vbox_stars <- glayout(container = g4)

  # Variable selection column
  vbox_stars[1, 1, anchor = c(-1, 0)] <- "Variable Selection"
  vbox_stars[2, 1] <- Variables_stars <- gcheckboxgroup(names(data[num]),
    checked = TRUE, horizontal = FALSE)
  tooltip(Variables_stars) <- "Select variables to display in the nD Tour."

  vbox_stars[3, 1, anchor = c(-1, 0)] <- "Class Selection"
  vbox_stars[4, 1, anchor = c(-1, 0)] <- Class_stars <- gtable(names(data)[!num], 
    multiple = TRUE)
  tooltip(Class_stars) <- "Select a class variable to classify the data."

  # Tour selection column
  vbox_stars[1, 2, anchor=c(-1, 0)] <- "Tour Type"
  tour_types <- c("Grand", "Little","Local", "Guided")
  vbox_stars[2, 2] <- TourType_stars <- gradio(tour_types)
  tooltip(TourType_stars) <- "Select a nD Tour type."

  #Guided indices selection
  vbox_stars[3, 2, anchor=c(-1, 0)] <- "Guided indices"
  IntIndex <-c("holes","cmass","lda_pp","pda_pp")
  vbox_stars[4, 2, anchor=c(-1,-1)] <-  GuidedType_stars <- gdroplist(IntIndex)
  tooltip(GuidedType_stars) <- "Select an index type for guided tour."

  # Lambda selection
  vbox_stars[3, 3, anchor=c(-1, 0)] <-"Lambda"
  vbox_stars[4, 3] <- LambdaValue_stars <- gslider(from=0, to = 1, by = 0.01,value=0.02)
  tooltip(LambdaValue_stars) <- "Select lambda's value to calculate pda index."

  # dimension control
  vbox_stars[1, 3, anchor = c(-1, 0)] <- "Choose Dimension"
  dimensions <- c(3:length(data[num]))
  vbox_stars[2, 3, anchor = c(-1, 0)] <- Dimensions_stars <- gradio(dimensions)
  tooltip(Dimensions_stars) <- "Select dimension number for displaying the nD Tour."

  vbox_stars[5, 1, anchor = c(-1, 0)] <- "Choose Star Number"
  vbox_stars[6, 1, expand = TRUE] <- Star_stars <- gslider(from = 3, to = nrow(data), by = 1, value = 4 )
  tooltip(Star_stars) <- "Drag to choose the face number."

  # speed and pause
  vbox_stars[5,2, anchor = c(-1, 0)] <- "Speed"
  vbox_stars[6,2, expand = TRUE] <- sl_stars <- gslider(from = 0, to = 5, by = 0.1, value = 1)
  tooltip(sl_stars) <- "Drag to set the speed of the nD Tour."

  vbox_stars[5, 3] <- chk_pause_stars <- gcheckbox("Pause",
    handler = function(h, ...) pause_stars(svalue(h$obj)))
  tooltip(chk_pause_stars) <- "Click here to pause or continue the nD Tour."

  # buttons control
  anim_id <- NULL
  pause_stars <- function(paused) {
    svalue(chk_pause_stars) <- paused
    if (paused) {
      gtkIdleRemove(anim_id)
      anim_id <- NULL
    } else {
      if (!is.null(anim_id)) gtkIdleRemove(anim_id)
      anim_id <<- gIdleAdd(draw_frame_stars)
    }
  }
  buttonGroup_stars <- ggroup(horizontal = FALSE, container = vbox_stars)

  # addSpace(buttonGroup,10)
  button1_stars<- gbutton("Apply", container = buttonGroup_stars, handler = function(...){
    print("apply from gui_stars")
    opar <- par(mfrow = c(1,1))
    pause_stars(FALSE)
    update_tour_stars()
  })
  tooltip(button1_stars) <- "Click here to update the options."

  # addSpace(buttonGroup,10)
  button2_stars<- gbutton("Quit",container = buttonGroup_stars, handler = function(...) {
    pause_stars(TRUE)
    dispose(w)
  })
  tooltip(button2_stars) <- "Click here to update the options."

  vbox_stars[6,3, anchor = c(0, 1)] <- buttonGroup_stars
  
}
# ------------------------- End of Gui_stars ----------------------------
