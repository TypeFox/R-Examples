
plot3Dfunrgl <- function(funcname, ...)  {
  form <- c(names(formals(funcname)), namespersp, "cex", "")
  dots <- list(...)
  plot <- dots$plot
  rep <- dots$rgltype == "rep"
  if (length(rep) == 0)
     rep <- FALSE
  if (! is.null(dots$rgltype)) {
    if (dots$rgltype == "new") {
      dots$new <- TRUE
      dots$add <- FALSE
    } else if (dots$rgltype == "add") {
      dots$add <- TRUE
      dots$new <- FALSE
    } else if (dots$rgltype == "rep") {
      dots$new <- FALSE
      dots$add <- TRUE
      ids <- rgl.ids()
      toreplace <- NULL
      if (funcname%in% c("persp3D", "slice3D","slicecont3D","isosurf3D",
        "surf3D","spheresurf3D","image3D")) toreplace <- "surface"
      else if (funcname%in% c("ribbon3D", "hist3D","box3D","rect3D"))
        toreplace <- "quads"
      else if (funcname%in% c("scatter3D", "points3D", "voxel3D"))
        toreplace <- "points"
      else if (funcname%in% c("lines3D", "segments3D", "border3D", "contour3D"))
        toreplace <- "lines"
      else if (funcname == "text3D")
        toreplace <- "text"
      plist <- getplist()
      plist$imgnr<-plist$img<-plist$segm<-plist$pt<-plist$CIpt<-plist$labels <- NULL
      setplist(plist)
    }
    dots$rgltype <- NULL
  }
  if (is.null(plot))
    plot <- TRUE
  dots$plot <- FALSE
  if (is.null(dots$add))
    dots$add <- FALSE
  nd <- names(dots)
  irgl <- unique(c(which(!nd %in% form ),
            which(nd %in% c("lighting", "smooth"))))
  rgldots <- c(dots[irgl], add = dots$add)

  dots[irgl] <- NULL
  if(dots$add)
    rgldots$new <- FALSE
  plot3D:::refresh(FALSE)
  do.call(funcname, dots)
  if (rep)
    rgl.pop(type = "shapes", id = ids[ids$type == toreplace, 1])

#  plot3D:::refresh(TRUE)
  if (plot) 
    do.call ("plotrgl", rgldots) 
}

persp3Drgl <- function(...)
  plot3Dfunrgl("persp3D", ...)
#  Persp3Drgl(...)         # faster
ribbon3Drgl <- function(...)
  plot3Dfunrgl("ribbon3D", ...)
hist3Drgl <- function(...)
  plot3Dfunrgl("hist3D", ...)
scatter3Drgl <- function(...)
  plot3Dfunrgl("scatter3D", ...)
points3Drgl <- function(...)
  plot3Dfunrgl("points3D", ...)
lines3Drgl <- function(...)
  plot3Dfunrgl("lines3D", ...)
slice3Drgl <- function(...)
  plot3Dfunrgl("slice3D", ...)
slicecont3Drgl <- function(...)
  plot3Dfunrgl("slicecont3D", ...)
isosurf3Drgl <- function(...)
  plot3Dfunrgl("isosurf3D", ...)
voxel3Drgl <- function(...)
  plot3Dfunrgl("voxel3D", ...)
triangle3Drgl <- function(...)
  plot3Dfunrgl("triangle3D", ...)
surf3Drgl <- function(...)
  plot3Dfunrgl("surf3D", ...)
spheresurf3Drgl <- function(...)
  plot3Dfunrgl("spheresurf3D", ...)
segments3Drgl <- function(...)
  plot3Dfunrgl("segments3D", ...)
box3Drgl <- function(...)
  plot3Dfunrgl("box3D", ...)
border3Drgl <- function(...)
  plot3Dfunrgl("border3D", ...)
rect3Drgl <- function(...)
  plot3Dfunrgl("rect3D", ...)
text3Drgl <- function(...)
  plot3Dfunrgl("text3D", ...)
image3Drgl <- function(...)
  plot3Dfunrgl("image3D", ...)
contour3Drgl <- function(...)
  plot3Dfunrgl("contour3D", ...)
  