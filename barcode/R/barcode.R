"barcode" <- function(x,
                      outer.margins=list(bottom=unit(2, "lines"),
                                         left=unit(2, "lines"),
                                         top=unit(2, "lines"),
                                         right=unit(2, "lines")),
                      horizontal=TRUE,
                      xlim=NULL,
                      nint=0,
                      main="",
                      xlab="",
                      labelloc=TRUE,
                      axisloc=TRUE,
                      labelouter=FALSE,
                      newpage=TRUE,
                      fontsize=9,
                      ptsize=unit(0.25, "char"),
                      ptpch=1,
                      bcspace=NULL,
                      use.points=FALSE,
                      buffer=0.02,
                      log=FALSE,
                      outerbox=TRUE
              ) {

  if (!require(grid)) stop("library(grid) is required and unavailable.\n\n")
  if (!require(lattice)) stop("library(lattice) is required and unavailable.\n\n")

  if (!is.null(labelloc)) {
    if (labelloc=="right" | labelloc=="top") labelloc <- FALSE
    if (labelloc=="left" | labelloc=="bottom") labelloc <- TRUE
  }
  if (!is.null(axisloc)) {
    if (axisloc=="right" | axisloc=="top") axisloc <- FALSE
    if (axisloc=="left" | axisloc=="bottom") axisloc <- TRUE
  }

  if (is.vector(x) && !is.list(x)) x <- list(x)
  if (is.null(names(x))) names(x) <- as.character(1:length(x))

  if (is.matrix(x)) x <- as.data.frame(x)

  if (newpage) grid.newpage()
  grid.text(main, 0.5, unit(1, "npc")-unit(1,"lines"), gp=gpar(fontface="bold"))

  # If there is an axis at the top as well as a title, and the labels
  # are inside the viewport, make extra space:
  if (!is.null(axisloc) && !axisloc && main!="" && !labelouter)
    outer.margins$top <- outer.margins$top + unit(2, "lines")

  xlaboffset <- unit(2.5, "lines")
  if (!is.null(axisloc) && xlab!="" && !labelouter) {
    if (axisloc) {
      if (horizontal) outer.margins$bottom <- outer.margins$bottom + unit(1.5, "lines")
      else outer.margins$top <- outer.margins$top + unit(1.5, "lines")
    } else {
      if (horizontal) outer.margins$top <- outer.margins$top + unit(1.5, "lines")
      else {
        outer.margins$left <- outer.margins$left + unit(1.5, "lines")
        outer.margins$top <- outer.margins$top
        outer.margins$right <- outer.margins$right - unit(1.5, "lines")
      }
    }
  }
    
  if (horizontal) {
    thisangle <- 0
    thisjust <- c("left", "bottom")
  } else {
    thisangle <- 90
    thisjust <- c("left", "top")
    pushViewport(viewport(x=0, y=0,
                          width=convertHeight(unit(1, "npc"), "inches"),
                          height=convertWidth(unit(1, "npc"), "inches"),
                          just=c("left", "bottom")))
  }
  if (labelouter) {
    outer.margins=list(bottom=unit(0, "lines"),
                                       left=unit(0, "lines"),
                                       top=unit(0, "lines"),
                                       right=unit(0, "lines"))
  }

  vp.main <- viewport(x=outer.margins$left,
                      y=outer.margins$bottom,
                      width=unit(1, "npc")-outer.margins$right-outer.margins$left,
                      height=unit(1, "npc")-outer.margins$top-outer.margins$bottom,
                      just=thisjust, angle=thisangle,
                      name="main", clip="off")
  pushViewport(vp.main)
  if (outerbox) grid.rect()

  barcode.panel(x, horizontal=horizontal, nint=nint, xlim=xlim,
                labelloc=labelloc, labelouter=labelouter,
                fontsize=fontsize, ptsize=ptsize, bcspace=bcspace,
                use.points=use.points, xlab=xlab, xlaboffset=xlaboffset,
                axisloc=axisloc, buffer=0.02, log=log)
  popViewport(1)
  if (!horizontal) popViewport(1)

}


