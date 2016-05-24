#' Faces Tour Plotting
#' Plots the Faces Tour in tab g3, internal function used by gui_tour()
#'
#' Sets up the gui for the faces tour
#'
#' @keywords internal
#' @author Bei Huang\email{beihuang@@iastate.edu}, Di Cook \email{dicook@@iastate.edu}, and Hadley Wickham \email{hadley@@rice.edu} 
# =============================== Gui_faces ====================================
.interface_faces = function(g3, data, w){

  # =============== Function: update_tour_faces ==================
  tour <- NULL
  tour_anim <- NULL
  update_tour_faces <- function(...) {
    tour <<- .create_face_tour(data,
      var_selected = svalue(Variables_faces),
      # VarIndex = svalue(Variables_faces, index = TRUE),
      cat_selected = svalue(Class_faces), 
      dim_selected = svalue(Dimensions_faces),
      nface_selected = svalue(Faces_faces),
      tour_type = svalue(TourType_faces),
      guided_type = svalue(GuidedType_faces),
      lambda = svalue(LambdaValue_faces),
      aps = svalue(sl_faces)
    )
    tour_anim <<- with(tour, new_tour(data, tour_path))

    tour$display$init(tour$data)
    tour$display$render_frame()

    TRUE
  }
  # --------------------- End of update_tour_xy ------------------
  
  # ================= Function: draw_frame =======================
  draw_frame_faces <- function(...) {
    # if there's no tour, don't draw anything
    if (is.null(tour)) return(FALSE)

    tour_step <- tour_anim(svalue(sl_faces) / 33)
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
  # ---------------------- End of draw_frame ---------------------
  
  num <- sapply(data, is.numeric)
  # ==================Controls==========================
  vbox_faces <- glayout(container = g3)
    # Variable selection column
    vbox_faces[1, 1, anchor = c(-1, 0)] <- "Variable Selection"
    vbox_faces[2, 1] <- Variables_faces <- gcheckboxgroup(names(data[num]),
      checked = TRUE, horizontal = FALSE)
    tooltip(Variables_faces ) <- "Select variables to display in the nD Tour."

    vbox_faces[3, 1, anchor = c(-1, 0)] <- "Class Selection"
    vbox_faces[4, 1, anchor = c(-1, 0)] <- Class_faces <- gtable(names(data)[!num], 
    multiple = TRUE)
    tooltip(Class_faces) <- "Select a class variable to classify the data."

    # Tour selection column
    vbox_faces[1, 2, anchor=c(-1, 0)] <- "Tour Type"
    tour_types <- c("Grand", "Little", "Local", "Guided")
    vbox_faces[2, 2] <- TourType_faces <- gradio(tour_types)
    tooltip(TourType_faces) <- "Select a nD Tour type."

  #Guided indices selection
    vbox_faces[3, 2, anchor=c(-1, 0)] <- "Guided indices"
    IntIndex <-c("holes","cmass","lda_pp","pda_pp")
    vbox_faces[4, 2, anchor=c(-1,-1)] <-  GuidedType_faces <- gdroplist(IntIndex)
    tooltip(GuidedType_faces) <- "Select an index type for guided tour."

  # Lambda selection
    vbox_faces[3, 3, anchor=c(-1, 0)] <-"Lambda"
    vbox_faces[4, 3] <- LambdaValue_faces <- gslider(from=0, to = 1, by = 0.01,value=0.02)
    tooltip(LambdaValue_faces) <- "Select lambda's value to calculate pda index."
   
    # speed and pause
    vbox_faces[5,1, anchor = c(-1, 0)] <- "Speed"
    vbox_faces[6,1, expand = TRUE] <- sl_faces <- gslider(from = 0, to = 5, by = 0.1, value = 1)
    tooltip(sl_faces) <- "Drag to set the speed of the nD Tour."

    vbox_faces[5, 3] <- chk_pause_faces <- gcheckbox("Pause",
      handler = function(h, ...) pause_faces(svalue(h$obj)))
    tooltip(chk_pause_faces) <- "Click here to pause or continue the nD Tour."

    vbox_faces[5, 2, anchor = c(-1, 0)] <- "Choose Face Number"
    vbox_faces[6, 2, expand = TRUE] <- Faces_faces <- gslider(from = 2, to = nrow(data), by = 1, value = 4 )
    tooltip(Faces_faces) <- "Drag to choose the face number."

    # dimension control
    vbox_faces[1, 3, anchor = c(-1, 0)] <- "Choose Dimension"
    dimensions <- c(2:length(data[num]))
    vbox_faces[2, 3, anchor = c(-1, 0)] <- Dimensions_faces <- gradio(dimensions)
    tooltip(Dimensions_faces) <- "Select dimension number n for displaying the nD Tour."

    # buttons control
    anim_id <- NULL
    pause_faces <- function(paused) {
      svalue(chk_pause_faces) <- paused
      if (paused) {
        gtkIdleRemove(anim_id)
        anim_id <<- NULL
      } else {
        if (!is.null(anim_id)) gtkIdleRemove(anim_id)
        anim_id <<- gIdleAdd(draw_frame_faces)
      }
    }
    buttonGroup_faces <- ggroup(horizontal = FALSE, container = vbox_faces)

    # addSpace(buttonGroup,10)
    button1_faces<- gbutton("Apply", container = buttonGroup_faces, handler = function(...){
      print("apply from gui_faces")
      pause_faces(FALSE)
      update_tour_faces()
    })
    tooltip(button1_faces) <- "Click here to update the options."

    # addSpace(buttonGroup,10)
    button2_faces<- gbutton("Quit",container = buttonGroup_faces, handler = function(...) {
      pause_faces(TRUE)
      dispose(w)
    })
    tooltip(button2_faces) <- "Click here to close this window."                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         

    vbox_faces[6, 3, anchor = c(0, 1)] <- buttonGroup_faces
}
# -------------------------- End of Gui_faces ----------------------------------
