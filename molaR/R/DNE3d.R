#' Plot results of a DNE analysis of a surface 
#' 
#' plotting function
#' 
#' @param DNE_File An object that stores the output of the DNE
#' function
#' @param setRange User-defined range for plotting color scheme, see
#' Details
#' @param logColors Logical that log transforms the color scheme
#' @param edgeMask Logical that colors edge faces black to indicate their
#' lack of contribution to the total Dirichlet normal energy
#' @param outlierMask Logical that colors outlier faces dark gray to
#' indicate their lack of contribution to the Dirichlet normal energy
#' @param legend Logical indicating whether or not a legend
#' shold be displayed
#' @param legendScale numeric value setting the relative size of the legend
#' similar in function to cex
#' @param leftOffset numeric value between -1 and 1 setting the degree of
#' offset for the plotted surface to the left. Larger values set further to right. 
#' @param fieldofview Passes an argument to par3d changing the field of
#' view in degrees of the resulting rgl 
#'
#' @details This function creates a heat map on the mesh surface
#' corresponding to the Dirichlet normal energy of each face calculated by
#' the DNE function. Hottest colors represent highest normal energy
#' values
#'
#' Dirichlet normal energies for the faces of a mesh surface tend to be
#' positively skewed, with a small proportion of the faces contributing
#' much of the total energy for the surface. When logColors is enabled the
#' function colorizes based on the log transformed Dirichlet normal
#' energies, allowing for finer resolution between faces near the mode of
#' the energy per face distribution. Disabling logColors will display the
#' untransformed Dirichlet normal energies.
#'
#' The legend will update to reflect the other arguments chosen by the
#' user. Colors currently display in the legend in bins, however the colors
#' used in the displayed mesh surface are on a continuum. Ideally, the
#' legend should reflect a continuous stretch of color from the lowest
#' calculated Dirichlet normal energy to the highest. Future versions will
#' adjust the legend to this more intuitive display.
#'
#' By default, the function sets the lowest Dirichlet normal energy
#' calculated among all faces to a cool color and the highest normal energy
#' calculated among all faces to red, and then colors the remaining faces
#' on a continuous color spectrum between these two end points using
#' either absolute or log transformed Dirichlet normal energy values
#' (depending on the status of logColors). Since the scale is relative to the
#' energies of the input surface, visual comparisons cannot directly be
#' made between  multiple plots of different surfaces. The setRange
#' argument allows users to define the minimum and maximum of the
#' plotting color scheme and use it in multiple plots. This enables the
#' direct comparison of different surfaces to one another with red equal to
#' the user-defined maximum and a cool color equal to the user-defined
#' minimum. The user should choose reasonable bounds for the
#' maximum and minimum that are near the maximum and minimum
#' Dirichlet normal energies calculated for their surfaces. setRange will
#' not accept negative values.
#' 
#' The leftOffset value sets how far to the left the surface will appear, intended
#' to help avoid overlap with the legend. Defaults to 0.75.
#' 
#' legendScale sets the relative size of the scale in the same way cex works
#' 
#' fieldofview is set to a default of 0, which is an isometric projection.
#' Increasing it alters the degree of parallax in the perspective view, up to
#' a maximum of 179 degrees.
#'
#'
#' @import
#' rgl grDevices graphics utils
#' 
#' @export
#' DNE3d

DNE3d <- function (DNE_File, setRange = c(0, 0), logColors = TRUE, edgeMask = TRUE, 
                       outlierMask = TRUE, legend = TRUE, legendScale = 1, leftOffset = 1,
                       fieldofview = 0) 
{
  if(length(DNE_File$Boundary_Values)==1){edgeMask <- FALSE}
  plyFile <- DNE_File$plyFile
  DNEs <- DNE_File$Face_Values$Dirichlet_Energy_Densities*DNE_File$Face_Values$Face_Areas
  if (setRange[1] == 0 && setRange[2] == 0 && logColors == FALSE) {
    color_range <- DNEs * (1/max(DNEs))
    color_range <- (color_range * -1 + 1) * 0.575
    DNE_colors <- hsv(color_range)
    legend_labels <- round(seq(max(DNEs), min(DNEs), l = 10), digits = 4)
    scaled <- FALSE
  }
  if (setRange[1] == 0 && setRange[2] == 0 && logColors == TRUE) {
    zeroes <- which(DNEs == 0)
    DNEs[zeroes] <- 1e-06
    logDNEs <- log(DNEs)
    logDNEs <- logDNEs + abs(min(logDNEs))
    color_range <- logDNEs * (1/max(logDNEs))
    color_range <- (color_range * -1 + 1) * 0.575
    DNE_colors <- hsv(color_range)
    legend_labels <- round(exp(seq(log(max(DNEs)), log(min(DNEs)), l = 10)), digits = 4)
    scaled <- FALSE
  }
  if (setRange[1] != 0 | setRange[2] != 0 && logColors == FALSE) {
    setMax <- max(setRange)
    setMin <- min(setRange)
    if (setMin < 0) {stop("Negative values not accepted for face energy range")}
    if (setMax < max(DNEs)) {warning("setRange max is less than highest calculated face energy")}
    if (setMin > min(DNEs)) {
      warning("setRange min is greater than lowest calculated face energy")
    }
    color_range <- DNEs * (1/setMax)
    color_range <- (color_range * -1 + 1) * 0.575
    DNE_colors <- hsv(color_range)
    legend_labels <- round(seq(setMax, setMin, l = 10), digits = 4)
    scaled <- TRUE
  }
  if (setRange[1] != 0 | setRange[2] != 0 && logColors == TRUE) {
    setMax <- max(setRange)
    setMin <- min(setRange)
    if (setMin < 0) {stop("Negative values not accepted for face energy range")}
    if (setMax < max(DNEs)) {warning("setRange max is less than highest calculated face energy")}
    if (setMin > min(DNEs)) {
      warning("setRange min is greater than lowest calculated face energy")
    }
    zeroes <- which(DNEs == 0)
    DNEs[zeroes] <- 1e-06
    logDNEs <- log(DNEs)
    Top <- log(setMax) + abs(min(logDNEs))
    logDNEs <- logDNEs + abs(min(logDNEs))
    if (max(logDNEs) > setMax) {
      Adjustment <- max(logDNEs)/Top
      color_range <- logDNEs * (1/max(logDNEs)) * Adjustment
      color_range <- (color_range * -1 + 1) * 0.575
      DNE_colors <- hsv(color_range)
    }
    if (max(logDNEs) <= setMax) {
      color_range <- logDNEs * (1/setMax)
      color_range <- (color_range * -1 + 1) * 0.675
      DNE_colors <- hsv(color_range)
    }
    if (setMin == 0) {setMin <- 1e-06}
    legend_labels <- round(exp(seq(log(setMax), log(setMin), l = 10)), digits = 4)
    scaled <- TRUE
  }
  if (edgeMask == TRUE) {
    edges <- as.numeric(rownames(DNE_File$Boundary_Values))
    DNE_colors[edges] <- "#000000"
  }
  if (outlierMask == TRUE) {
    outliers <- as.numeric(rownames(DNE_File$Outliers))
    DNE_colors[outliers] <- "#505050"
  }
  open3d()
  par3d(windowRect = c(100, 100, 900, 900))
  colormatrix <- rep(DNE_colors, 3)
  colormatrix <- matrix(colormatrix, nrow = 3, byrow = TRUE)
  shade3d(plyFile, color = colormatrix, shininess = 110)
  if (legend == TRUE) {
    if(legendScale <= 0){stop("legendScale must be a positive number")}
    if(legendScale > 1.25){
      warning("legendScale greater than 1.25 will restrict legend visibility")
    }
    molaR_bgplot(DNE_Legend(DNELabels = rev(legend_labels), scaled = scaled,
                       edgeMask = edgeMask, outlierMask = outlierMask, logColors = logColors,
                       size = legendScale))
  }
  if (leftOffset > 1) {warning("Left offset greater than 1 may restrict mesh visibility")}
  if (leftOffset < -1) {warning("Left offset less than -1 may restrict mesh visibility")}
  rgl.viewpoint(fov = fieldofview)
  ZView <- par3d("observer")[3]
  XView <- leftOffset * ZView *0.05
  observer3d(XView, 0, ZView)
}