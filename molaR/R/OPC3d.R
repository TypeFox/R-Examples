#' Plot results of OPC analysis of a surface
#'
#' A function that produces a three-dimensional rendering of face
#' orientation on a surface. The OPC function will identify the
#' orientations of mesh faces and assign them to patches. It must be
#' performed prior to using the OPC3d function. 
#'
#' @param OPC_Output_Object An object that stores the output of
#' the OPC function
#' @param binColors Allows the user to change the colors filled in for
#' each directional bin 
#' @param patchOutline logical whether or not to outline the patches
#' @param outlineColor parameter designating which color to outline the patches in
#' @param maskDiscard logical indicating whether to discard the unused patches
#' @param legend Logical indicating whether or not a legend should
#' be displayed
#' @param legendScale cex style scaling factor for the legend
#' @param legendTextCol parameter designating color for the legend text
#' @param legendLineCol parameter designating the color for the legend lines
#' @param leftOffset numeric parameters disginating how far to offset the surface
#' @param fieldofview Passes an argument to par3d changing the
#' field of view in dregrees of the resulting rgl window
#' 
#' @details This function will assign a uniform color to all faces on the mesh
#' surface that share one of the 8 orientations identified by the OPC function. The
#' function returns a colored shade3d of the mesh so that patches can be visually
#' inspected. Future versions will include the option to black out patches not
#' included in the orientation patch count.
#' 
#' Several legend plotting options are availble including customizing the line and
#' text colors using color names with legendTextCol and legendLineCol, both default
#' to black. legendScale works like cex for setting the size of the relative size
#' of the legend.
#' 
#' leftOffset will determine how far the plotted surface is moved to the left to
#' avoid obstructing the legend. Users shold choose between -1 and 1. 
#' 
#' fieldofview is set to a default of 0, which is an isometric projection.
#' Increasing it alters the degree of parallax in the perspective view, up to a
#' maximum of 179 degrees.
#' 
#' colors will support any vector of 8 colors, in any coloration scheme. Default
#' draws from the hsv color space to evenly space color information, however user
#' can supply a list of RGB values, character strings, or integers in place.
#'
#' @import
#' rgl
#'
#' @export
#' OPC3d


OPC3d <- function (OPC_Output_Object, 
                   binColors = hsv(h=(seq(10, 290, 40)/360), s=0.9, v=0.85),
                   patchOutline = FALSE, outlineColor = "black", maskDiscard = FALSE,
                   legend = TRUE, legendScale= 1, legendTextCol = "black",
                   legendLineCol = "black", leftOffset = 1, fieldofview = 0) 
{
  plyFile <- OPC_Output_Object$plyFile
  bins <- plyFile$Directional_Bins
  BinCount <- as.numeric(length(unique(plyFile$Directional_Bins)))
  BlackPatch <- NULL
  for (i in 1:BinCount) {
    Bin <- which(bins == i)
    bins[Bin] <- binColors[i]
    if (maskDiscard == TRUE) {
      if(OPC_Output_Object$Parameters$Minimum_Area==0){
        PatchList <- unlist(OPC_Output_Object$Patches[i], 
                            recursive = F)
        SmallPatch <- names(which(lapply(PatchList, length) < 
                                    OPC_Output_Object$Parameters$Minimum_Faces))
        Discarded <- as.numeric(unlist(PatchList[SmallPatch]))
        BlackPatch <- c(BlackPatch, Discarded)
      }
      if(OPC_Output_Object$Parameters$Minimum_Area>0){
        AreaList <- as.vector(OPC_Output_Object$Patch_Details[[i]][,2])
        MinAreaPercentage <- sum(OPC_Output_Object$plyFile$Face_Areas)*
          OPC_Output_Object$Parameters$Minimum_Area
        SmallPatchList <- which(AreaList < MinAreaPercentage)
        Discarded <- as.numeric(unlist(OPC_Output_Object$Patches[[i]][SmallPatchList]))
      }
      BlackPatch <- c(BlackPatch, Discarded)
    }
  }
  colormatrix <- bins
  if (maskDiscard == TRUE) {
    colormatrix[BlackPatch] <- "#000000"
  }
  colormatrix <- rep(colormatrix, 3)
  colormatrix <- matrix(colormatrix, nrow = 3, byrow = T)
  open3d()
  par3d(windowRect = c(100, 100, 900, 900))
  if (patchOutline == TRUE) {
    for (i in 1:BinCount) {
      Orientation <- OPC_Output_Object$Patches[i]
      PatchCount <- as.numeric(length(Orientation[[1]]))
      for (j in 1:PatchCount) {
        Patch <- Orientation[[1]][j]
        Patch <- as.numeric(Patch[[1]])
        Faces <- t(plyFile$it[, Patch])
        fnum <- length(Faces[, 1])
        vorder <- vector("list", fnum)
        for (i in 1:fnum) {vorder[[i]] <- unlist(sort(Faces[i, ]))}
        edges <- vector("list", fnum)
        for (i in 1:fnum) {
          Ordered <- vorder[[i]]
          G1 <- Ordered[1]
          G2 <- Ordered[2]
          G3 <- Ordered[3]
          ED1 <- paste(G1, G2, sep = "_")
          ED2 <- paste(G1, G3, sep = "_")
          ED3 <- paste(G2, G3, sep = "_")
          edges[[i]] <- paste(ED1, ED2, ED3, sep = ",")
        }
        for (i in 1:fnum) {edges[[i]] <- unlist(strsplit(edges[[i]], ","))}
        string <- unlist(edges)
        edgeframe <- data.frame(names = string)
        UniqueEdge <- aggregate(edgeframe, list(edgeframe$names), FUN = length)
        PatchEdge <- subset(UniqueEdge, UniqueEdge$names == 1)
        EdgeVerts <- as.numeric(unlist(strsplit(as.character(unlist(PatchEdge$Group.1)), "_")))
        EdgeCoords <- plyFile$vb[1:3, EdgeVerts]
        segments3d(t(EdgeCoords), color = outlineColor, 
                   lwd = 1.25, shininess = 120)
      }
    }
  }
  shade3d(plyFile, color = colormatrix, shininess = 110)
  if (legend == TRUE) {
    if(legendScale <= 0){stop("legendScale must be a positive number")}
    if(legendScale > 1.05){
      warning("legendScale greater than 1.05 will restrict legend visibility")
    }
    Fills <- rep("#FFFFFF", BinCount)
    for (i in 1:BinCount) {
      Fills[i] <- binColors[i]
    }
    molaR_bgplot(OPC_Legend(binColors=Fills, binNumber = BinCount, maskDiscard = maskDiscard,
                            size = legendScale, textCol=legendTextCol, lineCol=legendLineCol))
  }
  if (leftOffset > 1) {warning("Left offset greater than 1 may restrict mesh visibility")}
  if (leftOffset < -1) {warning("Left offset less than -1 may restrict mesh visibility")}
  rgl.viewpoint(fov = fieldofview)
  ZView <- par3d("observer")[3]
  XView <- leftOffset * ZView *0.055
  observer3d(XView, 0, ZView)
}