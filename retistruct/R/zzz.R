.onLoad <- function(lib, pkg) {
  options(TF.col="darkcyan")
  options(TB.col="darkcyan")
  options(V.col="darkcyan")
  options(stitch.col="cyan")
  options(V.stitch.col="darkcyan")
  options(grid.maj.col="black")
  options(grid.min.col="gray")
  options(outline.col="black")
  options(contour.levels=c(5, 25, 50, 75, 95))
  options(max.proj.dim=400)
  options(show.sphere=TRUE)
  options(retistruct.print.pdf.width=6)
}
