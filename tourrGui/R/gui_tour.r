##' Tour GUI                                   
##' Displays all types of Tour GUI in different tabs                     
##'
##' Combines all of the tour gui's into one, putting a separate tab for each.
##' 
##' @param data matrix, or data frame containing numeric columns, defaults to flea dataset
##' @param ... other arguments passed on to \code{\link{animate}} and \code{\link{display_xy}}
##' @author Bei Huang\email{beihuang@@iastate.edu}, Di Cook \email{dicook@@iastate.edu}, and Hadley Wickham \email{hadley@@rice.edu} 
##' @keywords display
##' @references Bei Huang, Dianne Cook, Hadley Wickham (2012).
##'   tourrGui: A gWidgets GUI for the Tour to Explore High-Dimensional
##'   Data Using Low-Dimensional Projections. Journal of Statistical
##'   Software, 49(6), 1-12. \url{http://www.jstatsoft.org/v49/i06/}.
##' @export
##' @examples
##' \dontrun{gui_tour(flea)}
gui_tour<- function(data = flea, ...) {
  require(tourr)
  require(colorspace)
  require(gWidgets)
  require(RGtk2)
  options("guiToolkit"="RGtk2")
  require(ash)
  require(TeachingDemos)  
  

  os <- find_platform()$os
  num <- sapply(data, is.numeric)
  tour <- NULL
  tour_anim <- NULL

  w <- gwindow("2D Tour plot example", visible = TRUE)
  mw = gnotebook(container = w, closebuttons = TRUE)  
  g8 = ggroup(container = mw, horizontal = FALSE,label="gui_scatmat") 
  g7 = ggroup(container = mw, horizontal = FALSE,label="gui_pcp")
  g6 = ggroup(container = mw, horizontal = FALSE,label="gui_stereo")
  g5 = ggroup(container = mw, horizontal = FALSE,label="gui_andrews")
  g4 = ggroup(container = mw, horizontal = FALSE,label="gui_stars")
  g3 = ggroup(container = mw, horizontal = FALSE,label="gui_faces")
  g2 = ggroup(container = mw, horizontal = FALSE,label="gui_density")
  g1 = ggroup(container = mw, horizontal = FALSE,label="gui_xy") 
    
  
  .interface_xy(g1,data,w)
  .interface_density(g2,data,w)
  .interface_faces(g3,data,w)
  .interface_stars(g4,data,w)
  .interface_andrews(g5,data,w)
  .interface_stereo(g6,data,w)
  .interface_pcp(g7,data,w)
  .interface_scatmat(g8,data,w)
  
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
  
  
  # pause_xy(FALSE)
  # pause_density(FALSE)
  visible(w) <- TRUE
  invisible()
}
