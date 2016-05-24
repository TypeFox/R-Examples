#' Scatterplot Tour Plotting
#' Plots the scatterplot Tour in tab g1
#'
#' Sets up the interface for the xy tour
#'
#' @keywords internal
#' @author Bei Huang\email{beihuang@@iastate.edu}, Di Cook \email{dicook@@iastate.edu}, and Hadley Wickham \email{hadley@@rice.edu} 
# =============================== Gui_xy ================================
.interface_xy = function (g1,data,w) {

  # ================= Function: update_tour_xy =======================
  tour <- NULL
  tour_anim <- NULL
  update_tour_xy <- function(...) {
    tour <<- .create_xy_tour(data,
    var_selected = svalue(Variables_xy),
    cat_selected = svalue(Class_xy),
    axes_location = svalue(dl_xy),
    tour_type = svalue(TourType_xy),
    guided_type = svalue(GuidedType_xy),
    lambda = svalue(LambdaValue_xy),
    aps = svalue(sl_xy)
    )
    tour_anim <<- with(tour, new_tour(data, tour_path, start=basis_random(ncol(data), 2)))

    tour$display$init(tour$data)
    tour$display$render_frame()
    TRUE
  }
  # ---------------------- End of update_tour_xy ---------------------


  # ================= Function: draw_frame_xy =======================
  draw_frame_xy <- function(...) {
    # if there's no tour, don't draw anything
    if (is.null(tour)) return(FALSE)

    tour_step <- tour_anim(svalue(sl_xy) / 33)
    if (is.null(tour_step$proj)) return(FALSE)

    if (os == "win") {
      tour$display$render_frame()
    } else {
      tour$display$render_transition()
    }
    with(tour_step, tour$display$render_data(tour$data, proj, target))
    Sys.sleep(1/33)

    TRUE
  }
  # ---------------------- End of draw_frame_xy ---------------------

  os <- find_platform()$os
  num <- sapply(data, is.numeric)
  # ================== Controls for tour_xy ==========================
  vbox_xy1 <- glayout(container = g1)
  # Variable selection column
  vbox_xy1[1, 1, anchor = c(-1, 0)] <- "Variable Selection"
  vbox_xy1[2, 1] <- Variables_xy <- gcheckboxgroup(names(data[num]),
    checked = TRUE, horizontal = FALSE)
  tooltip(Variables_xy) <- "Select variables to display in the 2D Tour."

  vbox_xy1[3, 1, anchor = c(-1, 0)] <- "Class Selection"
  vbox_xy1[4, 1, anchor = c(-1, 0)] <- Class_xy <- gtable(names(data)[!num],
    multiple = TRUE)
  tooltip(Class_xy) <- "Select a class variable to color the points."

  # Tour selection column
  vbox_xy1[1, 2, anchor=c(-1, 0)] <- "Tour Type"
  tour_types <- c("Grand", "Little", "Local", "Guided")
  vbox_xy1[2, 2] <- TourType_xy <- gradio(tour_types)
  tooltip(TourType_xy) <- "Select a 2D Tour type."

  vbox_xy1[3, 2, anchor=c(-1, 0)] <- "Guided indices"
  IntIndex <-c("holes","cmass","lda_pp","pda_pp")
  vbox_xy1[4, 2, anchor=c(-1,-1)] <-  GuidedType_xy <- gdroplist(IntIndex)
  tooltip(GuidedType_xy) <- "Select an index type for guided tour."

  vbox_xy1[3,3, anchor=c(-1, 0)] <-"Lambda"
  vbox_xy1[4,3] <- LambdaValue_xy <- gslider(from=0, to = 1, by = 0.01,value=0.02)
  #svalue(LambdaValue) <- 0.02
  tooltip(LambdaValue_xy) <- "Select lambda's value to calculate pda index."


  # Speed and pause
  vbox_xy1[5,1, anchor = c(-1, 0)] <- "Speed"
  vbox_xy1[6,1, expand = TRUE] <- sl_xy <- gslider(from = 0, to = 5, by = 0.1, value = 1)
  tooltip(sl_xy) <- "Drag to set the speed of the 2D Tour."

  vbox_xy1[6, 2] <- chk_pause_xy <- gcheckbox("Pause",
    handler = function(h, ...) pause_xy(svalue(h$obj)))
  tooltip(chk_pause_xy) <- "Click here to pause or continue the 2D Tour."

  # Axes control
  vbox_xy1[1,3, anchor=c(-1,0)] <- "Axes Locations"
  locations <- c("center", "bottomleft", "off")
  vbox_xy1[2,3, anchor=c(-1,0)] <- dl_xy <- gradio(locations)
  tooltip(dl_xy) <- "Select a location for the 2D Tour axes."


  # Buttons control
  anim_id <- NULL
  pause_xy <- function(paused) {
    svalue(chk_pause_xy) <- paused
    if (paused) {
      gtkIdleRemove(anim_id)
      anim_id <<- NULL
    } else {
      if (!is.null(anim_id)) gtkIdleRemove(anim_id)
      anim_id <<- gIdleAdd(draw_frame_xy)
    }
  }

  buttonGroup_xy <- ggroup(horizontal = FALSE, container = vbox_xy1)

  # Apply button & handler
  button1_xy<- gbutton("Apply", container = buttonGroup_xy, handler = function(...) {

    print("apply from gui_xy")
    if(is.null(anim_id))
    	pause_xy(FALSE)
    opar <- par(mfrow = c(1,1))
    update_tour_xy()
  })
  tooltip(button1_xy) <- "Click here to update the options."

  button2_xy<- gbutton("Quit",container = buttonGroup_xy, handler = function(...) {
    pause_xy(TRUE)
    dispose(w)
  })
  tooltip(button2_xy) <- "Click here to update the options."

  message1 <- gbutton("Help",container = buttonGroup_xy, handler = function(...) {
gmessage("GUI_xy allows user to control a dynamic plot by using a checkbox, a ratiobox, a table, a slider and some bottons. And it could easily be extended.
It's much more convenient for users to just click on this simple GUI instead of trying to figure out how to write the proper auguments for their desirable graphics.",
title="gui_help",icon="info")
  })

  vbox_xy1[5:6, 3, anchor = c(0, 1)] <- buttonGroup_xy

}
