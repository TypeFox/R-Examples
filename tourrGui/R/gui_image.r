##' Image Tour GUI                                   
##' Displays an Image Tour GUI                       
##'
##' This GUI allows users to control an image tour plot by simply moving and clicking their mouses. 
##' The Tour Type radio buttons contains four different tour types. They are the Grand Tour, Little Tour, Local Tour and Guided Tour. We can 
##' only choose one type a time. For the Guided Tour, we need to choose an index from the droplist to specify which particular search type is desired. 
##' The default index would be holes. 
##' The Speed slider can control the speed of the 1D tour. Simply dragging the mouse along the slider, changes the speed from slow to fast.
##' The Pause check box allow users to pause the dynamic 1D tour and have a close examination on the details.
##' The Apply button allows users to update the 1D tour, when it doesn't automatically update.
##' The Quit button allows users to close thie GUI window.
##' The Help button provides information about the tour and also what this GUI can do.
##' Tooltips will pop up when the mouse is moved over the GUI, which give hints about the functionality of the different GUI elements.
##' 
##' @param data a 3d array, the first two dimensions are locations on a grid, and the 3rd dimension gives the observations to be mixed with the tour.defaults to ozone dataset
##' @param ... other arguments passed on to \code{\link{animate}} and \code{\link{display_xy}}
##' @author Bei Huang\email{beihuang@@iastate.edu}, Di Cook \email{dicook@@iastate.edu}, and Hadley Wickham \email{hadley@@rice.edu}
##' @keywords display_image
##' @references Bei Huang, Dianne Cook, Hadley Wickham (2012).
##'   tourrGui: A gWidgets GUI for the Tour to Explore High-Dimensional
##'   Data Using Low-Dimensional Projections. Journal of Statistical
##'   Software, 49(6), 1-12. \url{http://www.jstatsoft.org/v49/i06/}.
##' @export
##' @examples
##' \dontrun{gui_image(ozone)}
gui_image <- function(data = ozone, ...) {
  require(tourr)
  require(gWidgets)
  require(RGtk2)
  options("guiToolkit"="RGtk2")

  os <- find_platform()$os
  num <- sapply(data, is.numeric)
  
  tour <- NULL
  tour_anim <- NULL
  update_tour <- function(...) {
    tour <<- .create_image_tour(data,
      tour_type = svalue(TourType),
      guided_type = svalue(GuidedType),
      aps = svalue(sl)
    )
    tour_anim <<- with(tour, new_tour(data, tour_path))
    
    tour$display$init(tour$data)
    tour$display$render_frame()
    
    TRUE
  }
  
  draw_frame <- function(...) {
    # if there's no tour, don't draw anything
    if (is.null(tour)) return(TRUE)  

    tour_step <- tour_anim(svalue(sl) / 33)
    if (os == "win") {
      tour$display$render_frame()
    } else {
      tour$display$render_transition()      
    }
    with(tour_step, tour$display$render_data(tour$data, proj, target))
    Sys.sleep(1/33)
    
    TRUE
  }
  
  
  # ==================Controls==========================
  w <- gwindow("2D Tour plot example", visible = FALSE)
  vbox <- glayout(container = w)

  # Tour selection column
  vbox[1, 1, anchor=c(-1, 0)] <- "Tour Type"
  tour_types <- c("Grand", "Little","Local","Guided")
  vbox[2, 1] <- TourType <- gradio(tour_types)
  tooltip(TourType) <- "Select a 1D Tour type."

  vbox[1, 2, anchor=c(-1, 0)] <- "Guided indices"
  IntIndex <-c("holes","cmass")
  vbox[2, 2, anchor=c(-1,-1)] <-  GuidedType <- gdroplist(IntIndex)
  tooltip(GuidedType) <- "Select an index type for guided tour."

  # speed and pause
  vbox[3,1, anchor = c(-1, 0)] <- "Speed"
  vbox[4,1, expand = TRUE] <- sl <- gslider(from = 0, to = 5, by = 0.1, value = 1)
  tooltip(sl) <- "Drag to set the speed of the 1D Tour."
  
  vbox[3,2] <- chk_pause <- gcheckbox("Pause", 
    handler = function(h, ...) pause(svalue(h$obj)))
  tooltip(chk_pause) <- "Click here to pause or continue the 1D Tour."

  # buttons control
  anim_id <- NULL
  pause <- function(paused) {
    svalue(chk_pause) <- paused
    if (paused) {
      gtkIdleRemove(anim_id)
      anim_id <<- NULL
    } else {
      if (!is.null(anim_id)) gtkIdleRemove(anim_id)
      anim_id <<- gIdleAdd(draw_frame)
    }
  }
  buttonGroup <- ggroup(horizontal = FALSE, container = vbox)  
  
  # addSpace(buttonGroup,10)
  button1<- gbutton("Apply", container = buttonGroup, handler = update_tour)
  tooltip(button1) <- "Click here to update the options."
  
  # addSpace(buttonGroup,10)
  button2<- gbutton("Quit",container = buttonGroup, handler = function(...) {
    pause(TRUE)
    dispose(w)
  })
  tooltip(button2) <- "Click here to close this window."

  # addSpace(buttonGroup,10)
  message1<-gbutton("Help",container = buttonGroup, handler = function(...) {
gmessage("The tour is a movie of low dimensional projections of high dimensional data. The projections are usually 1-, 2-, or 3-dimensional. They are used to expose interesting features of the high-dimensional data, such as outliers, clusters, and nonlinear dependencies.

When the projection dimension is 2, the data is usually shown as a scatterplot. Densities or histograms are used to display 1-dimensional projections. Projections of 3 or higher dimensions can be shown as stereo, parallel coordinates, scatterplot matrices or icons.

There are several different types of tours: grand, guided, little, and local. The grand tour generates a random path, while the guided uses an index on interest such as holes, central mass, lda or pda to guide the choice of projections to particular structure. The little tour moves between existing variables, only covering a subset of all the space. The local tour contrains the choice of projection to be those near the current view.

The GUI allows user to control the tour by checkboxes for the variable selection, slider for the speed, and toggle boxes for pause.",
title="gui_help",icon="info")
  })
tooltip(message1) <- "Click here for help."

  vbox[4, 2, anchor = c(0, 1)] <- buttonGroup
  
  # If on a mac, open a Cairo device, if there's not already one open
  # The cairo device has a much better refresh rate than Quartz
  if (find_platform()$os == "mac" && names(dev.cur()) != "Cairo") {
    require(Cairo)
    CairoX11()
  }
  
  update_tour()
  pause(FALSE)
  visible(w) <- TRUE
  
  invisible()
}

##' Image Tour Plotting
##' Plots the Image Tour
##'
##' @keywords internal
##' @author Bei Huang\email{beihuang@@iastate.edu}, Di Cook \email{dicook@@iastate.edu}, and Hadley Wickham \email{hadley@@rice.edu} 
.create_image_tour <- function(data, tour_type, guided_type, aps) {
 if (aps > 9999) {
   gmessage("Please quit", icon = "warning")
   return()
 }

  xs <- dim(data)[1]
  ys <- dim(data)[2]
  zs <- dim(data)[3]
  dim(data) <- c(xs * ys, zs)

  display <- display_image(xs, ys)
  # Work out which type of tour to use
  tour <- switch(tour_type,
    "Grand" = grand_tour(1), 
    "Little" = little_tour(), 
    "Local" = local_tour(),
    "Guided" = switch(guided_type, "holes"=guided_tour(holes), 
				"cmass"=guided_tour(cmass))
  )
      
  list(
    data = rescale(data),
    tour_path = tour,
    display = display,
    aps = aps
  )
}
