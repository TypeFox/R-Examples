#' Plot 3D and 2D areas of a mesh used to calculate relief index
#'
#' A function that plots a three-dimensional model of the mesh
#' surface and includes a footprint of the two-dimensional area for
#' visual comparison. 
#' 
#' @param RFI_Output An object that stores the output of the RFI
#' function
#' @param displacement Moves the surface footprint some
#' proportion of the height of the mesh. 0 is no
#' displacement. Expects a value, negative values displace
#' the footprint downward. 
#' @param SurfaceColor changes the color of the 3D surface mesh
#' @param FootColor changes color of the 2D surface footprint
#' @param FootPts logical indicating whether to plot the
#' flattened points of the footprint from the original ply file
#' @param FootPtsColor color for the plotted footprint points
#' @param Opacity adjusts the opacity of the 3D mesh
#' surface
#' @param legend Logical indicating whether or not to include a
#' legend of the colors chosen to represent the 3D surface and
#' footprint
#' @param legendScale cex style numeric relative scaling factor for the legend
#' @param leftOffset how numeric between -1 and 1 for which to offset the surface
#' relative to the legend.
#' @param fieldofview Passes an argument to par3d changing the
#' field of view in degrees of the resulting rgl window. 
#' 
#' @details This function can help to visualize the three-dimensional and two
#' dimensional areas that are used in calculating the relief index of a surface by
#' displaying both at the same time. The RFI function must be performed first.
#' 
#' Opacity can be adjusted in a range from fully opaque (1) to fully
#' transparent (0) in order to help visualize the footprint. The vertical placement of
#' the footprint along the Z axis can be altered with displace depending on how the
#' user wishes to view the surface, or on the original mesh orientation.
#' 
#' fieldofview is set to a default of 0, which is an isometric projection. Increasing it
#' alters the degree of parallax in the perspective view, up to a maximum of 179
#' degrees.
#'
#' @import
#' rgl
#' 
#' @export
#' RFI3d



RFI3d <- function (RFI_Output, displacement = -1.9, SurfaceColor = "gray",
                   FootColor = "red", FootPts = FALSE, FootPtsColor = "black",
                   Opacity = 1, legend = F, legendScale = 1, leftOffset = 0,
                   fieldofview = 0) 
{
  plyFile <- RFI_Output$plyFile
  Vertices <- plyFile$vb
  x <- Vertices[1, ] - mean(Vertices[1, ])
  y <- Vertices[2, ] - mean(Vertices[2, ])
  z <- Vertices[3, ] - mean(Vertices[3, ])
  n <- rep(1, length(Vertices))
  Shifted <- as.matrix(rbind(x, y, z, n))
  ShiftedPly <- plyFile
  ShiftedPly$vb <- Shifted
  FootColor = FootColor
  open3d()
  par3d(windowRect = c(100, 100, 900, 900))
  shade3d(ShiftedPly, color = SurfaceColor, alpha = Opacity, shininess=110)
  if (legend == T) {
    if(legendScale <= 0){stop("legendScale must be a positive number")}
    if(legendScale > 1.25){
      warning("legendScale greater than 1.25 will restrict legend visibility")
    }
    molaR_bgplot(RFI_Legend(surfCol = SurfaceColor, footCol = FootColor,
                                   size = legendScale, opac = Opacity))
  }
  FootprintPts <- RFI_Output$Flattened_Pts
  MeshHeight <- abs(max(plyFile$vb[3, ]) - min(plyFile$vb[3, ]))
  displaceDist <- displacement*0.5*MeshHeight
  zpts <- FootprintPts[,3] + displaceDist
  xyz <- cbind(FootprintPts[,1:2], zpts)
  if(FootPts==TRUE){
    points3d(xyz, color=FootPtsColor)
  }
  FootprintVertices <- t(cbind(xyz, rep(1, length(zpts))))
  triangles <- t(RFI_Output$Footprint_Triangles)
  Footprint <- list(vb=FootprintVertices, it=triangles, primitivetype="triangle", material=NULL)
  class(Footprint) <- c("mesh3d", "shape3d")
  shade3d(Footprint, color=FootColor)
  rgl.viewpoint(fov = fieldofview)
  if (leftOffset > 1) {warning("Left offset greater than 1 may restrict mesh visibility")}
  if (leftOffset < -1) {warning("Left offset less than -1 may restrict mesh visibility")}
  rgl.viewpoint(fov = fieldofview)
  ZView <- par3d("observer")[3]
  XView <- leftOffset * ZView *0.055
  observer3d(XView, 0, ZView)
}