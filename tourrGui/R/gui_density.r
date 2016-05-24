##' Density Tour GUI                                   
##' Displays a Density Tour GUI                       
##'
##' This GUI allows users to control the density tour by simply moving and clicking their mouses.
##' The Variable Selection checkboxes contains all the numeric variables, and at least two of them need to be checked to make the display work.
##' All the categorical variables go to the Class Selection box. We should select the class variable by double clicking the variable names. 
##' Color isn't implemented with the density tour yet. 
##' The Tour Type radio buttons contains four different tour types. They are the Grand Tour, Little Tour, Local Tour and Guided Tour. We can 
##' only choose one type a time. For the Guided Tour, we need to choose an index from the droplist to specify which particular search type is desired. 
##' The default index would be holes. For tour type Guided(lda_pp) and Guided(pda_pp), we also need to specify class variable first, and the Guided(pda_pp) 
##' is also controlled by another parameter, lambda. Lambda ranges from 0 to 1, with default at 0.02. A value of 0 will make the tour operate like Guided(lda_pp). 
##' For high-dimensional data a value closer to 1 would be advised.
##' The Method Type radio buttons contains three different display methods. They are histogram, density plot and ash plot. 
##' The distribution of data projected into 1d can be displayed correspondingly as a histogram, kernel density estimate and average shifted histogram.
##' The Axes Locations column contains two choices, TRUE and FALSE. TRUE means the tour will center at the middle of x-axes, 
##' FALSE means the tour will wander to the left and right. The default value is TRUE.
##' The Speed slider can control the speed of the 1D tour. Simply dragging the mouse along the slider, changes the speed from slow to fast.
##' The Pause check box allow users to pause the dynamic 1D tour and have a close examination on the details.
##' The Apply button allows users to update the 1D tour, when it doesn't automatically update.
##' The Quit button allows users to close thie GUI window.
##' The Help button provides information about the tour and also what this GUI can do.
##' Tooltips will pop up when the mouse is moved over the GUI, which give hints about the functionality of the different GUI elements.
##' 
##' @param data matrix, or data frame containing numeric columns, defaults to flea dataset
##' @param ... other arguments passed on to \code{\link{animate}} and \code{\link{display_xy}}
##' @author Bei Huang\email{beihuang@@iastate.edu}, Di Cook \email{dicook@@iastate.edu}, and Hadley Wickham \email{hadley@@rice.edu} 
##' @keywords display_density
##' @references Bei Huang, Dianne Cook, Hadley Wickham (2012).
##'   tourrGui: A gWidgets GUI for the Tour to Explore High-Dimensional
##'   Data Using Low-Dimensional Projections. Journal of Statistical
##'   Software, 49(6), 1-12. \url{http://www.jstatsoft.org/v49/i06/}.
##' @export
##' @examples
##' \dontrun{gui_density(flea)}
gui_density <- function(data = flea, ...) {
  require(tourr)
  require(gWidgets)
  require(RGtk2)
  options("guiToolkit"="RGtk2")
  require(ash)

  
  os <- find_platform()$os
  num <- sapply(data, is.numeric)
  
  tour <- NULL
  tour_anim <- NULL
  update_tour <- function(...) {
    tour <<- .create_1d_tour(data,
      var_selected = svalue(Variables), 
      cat_selected = svalue(Class), 
      method_selected = svalue(MethodType),
      center_selected = svalue(CenterType),
      tour_type = svalue(TourType),
      guided_type = svalue(GuidedType),
      lambda = svalue(LambdaValue),
      aps = svalue(sl)
    )
    tour_anim <<- with(tour, new_tour(data, tour_path, start=basis_random(ncol(data), 1)))
    
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
  vbox[1, 1, anchor = c(-1, 0)] <- "Variable Selection"
  vbox[2, 1] <- Variables <- gcheckboxgroup(names(data[num]), 
    checked = TRUE, horizontal = FALSE)
  tooltip(Variables) <- "Select variables to display in the 1D Tour."
  
  vbox[3, 1, anchor = c(-1, 0)] <- "Class Selection"
  vbox[4, 1, anchor = c(-1, 0)] <- Class <- gtable(names(data)[!num], 
    multiple = TRUE)
  tooltip(Class) <- "Select a class variable to classify the data."


  # Tour selection column
  vbox[1, 2, anchor=c(-1, 0)] <- "Tour Type"
  tour_types <- c("Grand", "Little", "Local", "Guided")
  vbox[2, 2] <- TourType <- gradio(tour_types)
  tooltip(TourType) <- "Select a 1D Tour type."

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
  tooltip(sl) <- "Drag to set the speed of the 1D Tour."
  
  vbox[5, 3] <- chk_pause <- gcheckbox("Pause", 
    handler = function(h, ...) pause(svalue(h$obj)))
  tooltip(chk_pause) <- "Click here to pause or continue the 1D Tour."

  # method control
  vbox[1, 3, anchor = c(-1, 0)] <- "Method Type"
  method_types <- c("density","hist","ash")
  vbox[2, 3, anchor = c(-1, 0)] <- MethodType <- gradio(method_types)
  tooltip(MethodType) <- "Select a display method for the 1D tour."
    
  # center control
  vbox[5,2, anchor=c(-1,0)] <- "Axes Locations"
  center_types <- c("TRUE", "FALSE")
  vbox[6,2, anchor=c(-1,0)] <- CenterType <- gradio(center_types)
  tooltip(CenterType) <- "Choose to center the display or let it wander."

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
  button2<- gbutton("Quit",container = buttonGroup, handler = function(...) {
    pause(TRUE)
    dispose(w)
  })
  tooltip(button2) <- "Click here to close this window."

  # addSpace(buttonGroup,10)
  message1_den<-gbutton("Help",container = buttonGroup, handler = function(...) {
gmessage("The tour is a movie of low dimensional projections of high dimensional data. The projections are usually 1-, 2-, or 3-dimensional. They are used to expose interesting features of the high-dimensional data, such as outliers, clusters, and nonlinear dependencies.

When the projection dimension is 2, the data is usually shown as a scatterplot. Densities or histograms are used to display 1-dimensional projections. Projections of 3 or higher dimensions can be shown as stereo, parallel coordinates, scatterplot matrices or icons.

There are several different types of tours: grand, guided, little, and local. The grand tour generates a random path, while the guided uses an index on interest such as holes, central mass, lda or pda to guide the choice of projections to particular structure. The little tour moves between existing variables, only covering a subset of all the space. The local tour contrains the choice of projection to be those near the current view.

The GUI allows user to control the tour by checkboxes for the variable selection, slider for the speed, and toggle boxes for pause.",
title="gui_help",icon="info")
  })
tooltip(message1_den) <- "Click here for help."

  vbox[6, 3, anchor = c(0, 1)] <- buttonGroup
  
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

##' For generating 1D projections
##'
##' The 1D projection is used for the density tour
##'
##' @keywords internal
##' @author Bei Huang\email{beihuang@@iastate.edu}, Di Cook \email{dicook@@iastate.edu},and Hadley Wickham \email{hadley@@rice.edu} 
.create_1d_tour <- function(data, var_selected, cat_selected, method_selected, center_selected, tour_type, guided_type, lambda, aps) {
  if (length(var_selected) < 2) {
    gmessage("Please select at least two variables", icon = "warning")
    return()
  }
   
  display <- display_dist(method = method_selected, center = center_selected ,col=col) 

  # Work out which type of tour to use
  tour <- switch(tour_type,
    "Grand" = grand_tour(1), 
    "Little" = little_tour(), 
    "Guided" = switch(guided_type,
               "holes"=guided_tour(holes, 1), 
	       "cmass"=guided_tour(cmass, 1),
	       "lda_pp" = guided_tour(lda_pp(data[,cat_selected]), 1),
	       "pda_pp" = guided_tour(pda_pp(data[,cat_selected],lambda)), 1),
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
