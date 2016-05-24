##' Scatterplot Tour GUI
##' A graphical user interface enabling interactive control of a scatterplot tour.
##'
##' This GUI allows users to control the scatterplot tour by simply moving and clicking their mouses.
##' The Variable Selection checkboxes contains all the numeric variables, and at least three of them need to be checked to make the display work.
##' All the categorical variables go to the Class Selection box. We should select the class variable by double clicking the variable names. 
##' If users don't specify the class variable, the selected numeric variables will be considered as one class, and all points will be black.
##' After users specify the class variable, the points will be considered as different classes according to this 
##' categorical variable, and these will appear as rainbow colors in the scatterplot tour.
##' The Tour Type radio buttons contains four different tour types. They are the Grand Tour, Little Tour, Local Tour and Guided Tour. We can 
##' only choose one type a time. For the Guided Tour, we need to choose an index from the droplist to specify which particular search type is desired. 
##' The default index is holes. For tour type Guided(lda_pp) and Guided(pda_pp), we also need to specify class variable first. The Guided(pda_pp) 
##' is also controlled by another parameter, lambda. Lambda ranges from 0 to 1, with default at 0.02. A value of 0 will make the tour operate like Guided(lda_pp). 
##' For very high-dimensional data a value closer to 1 would be advised.
##' The Axes Locations column contains three types. Users can specify where tour axes will be displayed. The choices are center, bottomleft and off.
##' The Speed slider can control the speed of the 2D tour. Simply dragging the mouse along the slider, changes the speed from slow to fast.
##' The Pause check box allow users to pause the dynamic 2D tour and have a close examination on the details.
##' The Apply button allows users to update the 2D tour, when it doesn't automatically update.
##' The Quit button allows users to close thie GUI window.
##' The Help button provides information about the tour and also what this GUI can do.
##' Tooltips will pop up when the mouse is moved over the GUI, which give hints about the functionality of the different GUI elements.
##'
##' @param data matrix, or data frame containing numeric columns, defaults to flea dataset
##' @param ... other arguments passed on to \code{\link{animate}} and \code{\link{display_xy}}
##' @author Bei Huang\email{beihuang@@iastate.edu}, Di Cook \email{dicook@@iastate.edu}, and Hadley Wickham \email{hadley@@rice.edu}
##' @keywords display_xy
##' @references Bei Huang, Dianne Cook, Hadley Wickham (2012).
##'   tourrGui: A gWidgets GUI for the Tour to Explore High-Dimensional
##'   Data Using Low-Dimensional Projections. Journal of Statistical
##'   Software, 49(6), 1-12. \url{http://www.jstatsoft.org/v49/i06/}.
##' @export
##' @examples
##' \dontrun{gui_xy(flea)}
gui_xy <- function(data = flea, ...) {
  require(tourr)
  require("colorspace")
  require("gWidgets")
  require("RGtk2")
  options("guiToolkit"="RGtk2")
 

  os <- find_platform()$os
  num <- sapply(data, is.numeric)
  
  tour <- NULL
  tour_anim <- NULL
  update_tour <- function(...) {
    tour <<- .create_xy_tour(data,
      var_selected = svalue(Variables), 
      cat_selected = svalue(Class), 
      axes_location = svalue(dl),
      tour_type = svalue(TourType),
      guided_type = svalue(GuidedType),
      lambda = svalue(LambdaValue),
      aps = svalue(sl)
    )

    tour_anim <<- with(tour, new_tour(data, tour_path, start=basis_random(ncol(data), 2)))

    tour$display$init(tour$data)
    tour$display$render_frame()    

    TRUE
  }
  
  draw_frame <- function(...) {
    # if there's no tour, don't draw anything
    if (is.null(tour)) return(FALSE)  

    tour_step <- tour_anim(svalue(sl) / 33)
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
  
  
  # ==================Controls==========================
  w <- gwindow("2D Tour plot example", visible = FALSE)
  vbox <- glayout(container = w)

  # Variable selection column
  # Divide all the variable names into two parts according to numeric variables and categorical variables.
  # The numeric variable names go under Variable Selection column, 
  # and the categorical variable names go under Class Selection column.

  vbox[1, 1, anchor = c(-1, 0)] <- "Variable Selection"
  vbox[2, 1] <- Variables <- gcheckboxgroup(names(data[num]), 
    checked = TRUE, horizontal = FALSE)
  tooltip(Variables) <- "Select variables to display in the 2D Tour."

  vbox[3, 1, anchor = c(-1, 0)] <- "Class Selection"
  vbox[4, 1, anchor = c(-1, 0)] <- Class <- gtable(names(data)[!num], 
    multiple = TRUE)
  tooltip(Class) <- "Select a class variable to color the points."


  # Tour selection column
  # The are seven kinds of tour type which allow users to switch willingfully.
 
  vbox[1, 2, anchor=c(-1, 0)] <- "Tour Type"
  tour_types <- c("Grand", "Little", "Local", "Guided")
  vbox[2, 2] <- TourType <- gradio(tour_types)
  tooltip(TourType) <- "Select a 2D Tour type."

  # Guided Tour index selection column
  # Indices including "holes", "cmass", "lda", "pda".

  vbox[3, 2, anchor=c(-1, 0)] <- "Guided indices"
  IntIndex <- c("holes", "cmass", "lda_pp", "pda_pp")
  vbox[4, 2, anchor=c(-1,-1)] <- GuidedType <- gdroplist(IntIndex)
  tooltip(GuidedType) <- "Select an index type for guided tour."
  
  # Lambda selection column
  # Lambda's range is from 0 to 1.

  vbox[3,3, anchor=c(-1, 0)] <-"Lambda"
  vbox[4,3] <- LambdaValue <- gslider(from=0, to = 1, by = 0.01,value=0.02)
  #svalue(LambdaValue) <- 0.02
  tooltip(LambdaValue) <- "Select lambda's value to calculate pda index."

  # speed and pause
  # This slider can control the speed of the 2D tour, which ranged from 0 to 5.

  vbox[5,1, anchor = c(-1, 0)] <- "Speed"
  vbox[6,1, expand = TRUE] <- sl <- gslider(from = 0, to = 5, by = 0.1, value = 1)
  tooltip(sl) <- "Drag to set the speed of the 2D Tour."
 
  # Pause box allow users to pause the dynamic 2D tour and have a close examination on the details.
  vbox[6, 2] <- chk_pause <- gcheckbox("Pause", 
    handler = function(h, ...) pause(svalue(h$obj)))
  tooltip(chk_pause) <- "Click here to pause or continue the 2D Tour."

  # axes control
  # There are three kinds of axes locations which allow users to switch willingfully.
  vbox[1,3, anchor=c(-1,0)] <- "Axes Locations"
  locations <- c("center", "bottomleft", "off")
  vbox[2,3, anchor=c(-1,0)] <- dl <- gradio(locations)
  tooltip(dl) <- "Select a location for the 2D Tour axes."

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
  button1<- gbutton("Apply", container = buttonGroup, handler = function(...) {
    pause(FALSE)
    update_tour()
  })
  tooltip(button1) <- "Click here to update the options."
 
  # addSpace(buttonGroup,10)
  button2<-gbutton("Quit", container = buttonGroup, handler = function(...) {
    pause(TRUE)
    dispose(w)
  })
  tooltip(button2) <- "Click here to close this window."
 

  # addSpace(buttonGroup,10)
  message1<-gbutton("Help", container = buttonGroup, handler = function(...) {
gmessage("The tour is a movie of low dimensional projections of high dimensional data. The projections are usually 1-, 2-, or 3-dimensional. They are used to expose interesting features of the high-dimensional data, such as outliers, clusters, and nonlinear dependencies.

When the projection dimension is 2, the data is usually shown as a scatterplot. Densities or histograms are used to display 1-dimensional projections. Projections of 3 or higher dimensions can be shown as stereo, parallel coordinates, scatterplot matrices or icons.

There are several different types of tours: grand, guided, little, and local. The grand tour generates a random path, while the guided uses an index on interest such as holes, central mass, lda or pda to guide the choice of projections to particular structure. The little tour moves between existing variables, only covering a subset of all the space. The local tour contrains the choice of projection to be those near the current view.

The GUI allows user to control the tour by checkboxes for the variable selection, slider for the speed, and toggle boxes for pause.",
title="gui_help",icon="info")
  })
tooltip(message1) <- "Click here for help."


  vbox[5:6, 3, anchor = c(0, 1)] <- buttonGroup
  
  # If on a mac, open a Cairo device, if there's not already one open
  # The cairo device has a much better refresh rate than Quartz
  if (find_platform()$os == "mac" && names(dev.cur()) != "Cairo") {
    require(Cairo)
    CairoX11()
  } else if (length(dev.list()) == 0) {
    # Open new display if necessary
    dev.new()
    # Turn off display list to maximise speed
    dev.control(displaylist = "inhibit")
  }
  
  update_tour()
  pause(FALSE)
  visible(w) <- TRUE
  
  invisible()
}

##' Scatterplot Tour Plotting
##' Plots the scatterplot Tour
##'
##' Creates tour for gui_xy
##'
##' @keywords internal
##' @author Bei Huang\email{beihuang@@iastate.edu}, Di Cook \email{dicook@@iastate.edu}, and Hadley Wickham \email{hadley@@rice.edu} 
.create_xy_tour <- function(data, var_selected, cat_selected, axes_location, tour_type, guided_type, lambda, aps) {
  if (length(var_selected) < 3) {
    gmessage("Please select at least three variables", icon = "warning")
    return()
  }
  
  # Work out point colours
  cat <- data[cat_selected]
  if (length(cat_selected) > 0) {
    # collapse to single variable if multiple selected
    int <- interaction(cat, drop = TRUE)
    pal <- rainbow_hcl(length(levels(int)))
    col <- pal[as.numeric(int)]
  } else {
    col <- "black"
  }
  
  display <- display_xy(axes = axes_location, center = TRUE, col = col)

  # Work out which type of tour to use
  tour <- switch(tour_type,
    "Grand" = grand_tour(), 
    "Little" = little_tour(), 
    "Guided" = switch(guided_type, "holes" = guided_tour(holes), 
				"cmass" = guided_tour(cmass),
				"lda_pp" = guided_tour(lda_pp(data[,cat_selected])),
				"pda_pp" = guided_tour(pda_pp(data[,cat_selected],lambda))),
    # "Local" = local_tour()
    "Local" = local_tour(basis_init(length(var_selected), 2))
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

