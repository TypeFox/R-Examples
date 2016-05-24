##' Stereo Tour GUI                                   
##' Displays an Stereo Tour GUI                       
##'
##' This GUI allows users to control the Stereo tour by simply moving and clicking their mouses.
##' The Variable Selection checkboxes contains all the numeric variables, and at least three of them need to be checked to make the display work.
##' All the categorical variables go to the Class Seclection box. We should select the class variable by double clicking the variable names. 
##' The Tour Type radio buttons contains four different tour types. They are the Grand Tour, Little Tour, Local Tour and Guided Tour. We can 
##' only choose one type a time. For the Guided Tour, we need to choose an index from the droplist to specify which particular search type is desired. 
##' The default index would be holes. For tour type Guided(lda_pp) and Guided(pda_pp), we also need to specify class variable first, and the Guided(pda_pp) 
##' is also controlled by another parameter, lambda. Lambda ranges from 0 to 1, with default at 0.02. A value of 0 will make the tour operate like Guided(lda_pp). 
##' For high-dimensional data a value closer to 1 would be advised.
##' The Speed slider can control the speed of the 3D tour. Simply dragging the mouse along the slider, changes the speed from slow to fast.
##' The Pause check box allow users to pause the dynamic 3D tour and have a close examination on the details.
##' The Apply button allows users to update the 3D tour, when it doesn't automatically update.
##' The Quit button allows users to close thie GUI window.
##' The Help button provides information about the tour and also what this GUI can do.
##' Tooltips will pop up when the mouse is moved over the GUI, which give hints about the functionality of the different GUI elements.
##' 
##' @param data matrix, or data frame containing numeric columns, defaults to flea dataset
##' @param ... other arguments passed on to \code{\link{animate}} and \code{\link{display_xy}}
##' @author Bei Huang\email{beihuang@@iastate.edu}, Di Cook \email{dicook@@iastate.edu}, and Hadley Wickham \email{hadley@@rice.edu} 
##' @keywords display_stereo
##' @references Bei Huang, Dianne Cook, Hadley Wickham (2012).
##'   tourrGui: A gWidgets GUI for the Tour to Explore High-Dimensional
##'   Data Using Low-Dimensional Projections. Journal of Statistical
##'   Software, 49(6), 1-12. \url{http://www.jstatsoft.org/v49/i06/}.
##' @export
##' @examples
##' \dontrun{gui_stereo(flea)}
gui_stereo <- function(data = flea, ...) {
  require(tourr) 
  require(gWidgets)
  require(RGtk2)
  options("guiToolkit"="RGtk2")

  os <- find_platform()$os
  num <- sapply(data, is.numeric)
  
  tour <- NULL
  tour_anim <- NULL
  update_tour <- function(...) {
    tour <<- .create_stereo_tour(data,
      var_selected = svalue(Variables),
      cat_selected = svalue(Class), 
      tour_type = svalue(TourType),
      guided_type = svalue(GuidedType),
      lambda = svalue(LambdaValue),
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

  # Variable selection column
  vbox[1, 1, anchor = c(-1, 0)] <- "Variable Selection"
  vbox[2, 1] <- Variables <- gcheckboxgroup(names(data[num]), 
    checked = TRUE, horizontal = FALSE)
  tooltip(Variables) <- "Select variables to display in the 3D Tour."

  vbox[3, 1, anchor = c(-1, 0)] <- "Class Selection"
  vbox[4, 1, anchor = c(-1, 0)] <- Class <- gtable(names(data)[!num], 
    multiple = TRUE)
  tooltip(Class) <- "Select a class variable to classify the points."

  # Tour selection column
  vbox[1, 2, anchor=c(-1, 0)] <- "Tour Type"
  tour_types <- c("Grand", "Little", "Local","Guided")
  vbox[2, 2] <- TourType <- gradio(tour_types)
  tooltip(TourType) <- "Select a 3D Tour type."

  vbox[3, 2, anchor=c(-1, 0)] <- "Guided indices"
  IntIndex <-c("holes","cmass","lda_pp","pda_pp")
  vbox[4, 2, anchor=c(-1,-1)] <-  GuidedType <- gdroplist(IntIndex)
  tooltip(GuidedType) <- "Select an index type for guided tour."

  vbox[3,3, anchor=c(-1, 0)] <-"Lambda"
  vbox[4,3] <- LambdaValue <- gslider(from=0, to = 1, by = 0.01,value=0.02)
  #svalue(LambdaValue) <- 0.02
  tooltip(LambdaValue) <- "Select lambda's value to calculate pda index."

  # speed and pause
  vbox[5,1, anchor = c(-1, 0)] <- "Speed"
  vbox[6,1, expand = TRUE] <- sl <- gslider(from = 0, to = 5, by = 0.1, value = 1)
  tooltip(sl) <- "Drag to set the speed of the 3D Tour."
  
  vbox[6, 2] <- chk_pause <- gcheckbox("Pause", 
    handler = function(h, ...) pause(svalue(h$obj)))
  tooltip(chk_pause) <- "Click here to pause or continue the 3D Tour."

  # buttons control
  anim_id <- NULL
  pause <- function(paused) {
    svalue(chk_pause) <- paused
    if (paused) {
      gtkIdleRemove(anim_id)
      anim_id <- NULL
    } else {
      if (!is.null(anim_id)) gtkIdleRemove(anim_id)
      anim_id <<- gIdleAdd(draw_frame)
    }
  }
  buttonGroup <- ggroup(horizontal = FALSE, container = vbox)  
  
  # addSpace(buttonGroup,10)
  button1<- gbutton("Apply", container = buttonGroup, handler = function(...) {
    pause(FALSE)
    update_tour()
  })
  tooltip(button1) <- "Click here to update the options."

  # addSpace(buttonGroup,10)
  button2<- gbutton("Quit", container = buttonGroup, handler = function(...) {
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

  vbox[5:6, 3, anchor = c(0, -1)] <- buttonGroup
  
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

##' Stereo Tour Plotting
##' Plots the Stereo Tour
##'
##' Initializes the stero tour
##'
##' @keywords internal
##' @author Bei Huang\email{beihuang@@iastate.edu}, Di Cook \email{dicook@@iastate.edu}, and Hadley Wickham \email{hadley@@rice.edu} 
.create_stereo_tour <- function(data, var_selected,cat_selected, tour_type, guided_type, lambda, aps) {
  if (length(var_selected) < 3) {
    gmessage("Please select at least three variables", icon = "warning")
    return()
  }

  blue = rgb(0, 0.91, 0.89)
  red = rgb(0.98, 0.052, 0)
  display <- display_stereo(blue = rgb(0, 0.91, 0.89), red = rgb(0.98, 0.052, 0))


  # Work out which type of tour to use
  tour <- switch(tour_type,
    "Grand" = grand_tour(3),
    "Little" = little_tour(3),
    "Guided" = switch(guided_type, "holes"=guided_tour(holes), 
				"cmass"=guided_tour(cmass),
				"lda_pp" = guided_tour(lda_pp(data[,cat_selected])),
				"pda_pp" = guided_tour(pda_pp(data[,cat_selected],lambda))),
    "Local" = local_tour(3)
  )
  
  sel <- data[var_selected]
  # Sphere the data if we're using a guided tour
  if (length(grep(tour_type, "Guided")) > 0) {
    sel <- sphere(sel)
  }
     
  list(
    data = rescale(sel),
    tour_path = tour,
    display = display,
    aps = aps
  )
}
