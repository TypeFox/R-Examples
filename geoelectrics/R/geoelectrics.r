# 3D-Visualization of Geoelectric Resistivity Measurement Profiles
# Anja Kleebaum

#' @import lattice
#' @import rgl
#' @import fields
#' @import methods
#' @importFrom grid grid.text
#' @importFrom grDevices colorRamp colorRampPalette rgb
#' @importFrom graphics axis legend plot points
#' @importFrom stats lm
#' @importFrom utils read.table

###---Settings---####
pointsize <- 10
colors <- c("blue", "green", "yellow", "orange", "red", "purple")
logBase <- 10