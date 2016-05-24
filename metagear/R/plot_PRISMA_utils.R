
# base constructor for generic phase box
phaseBox <- function(aLabel, 
                     x = 0.5, 
                     y = 0.5,
                     width = 25) {
  
  theLabels <- strwrap(aLabel, width)
  nLabel <- length(theLabels) 
  phaseVP <- viewport(x = x, 
                      y = y, 
                      width = max(stringWidth(theLabels)) + unit(12, "mm"), 
                      height = unit(nLabel + 2, "lines"))
  pushViewport(phaseVP)
    grid.roundrect()
    grid.text(theLabels, x = 0.5, y = unit(nLabel:1 + 0.5, "lines"))
  popViewport()
}

# generic constructor for phase box
phaseGrob <- function(aLabel, 
                      x = 0.5, 
                      y = 0.5, 
                      width = 25) {
  
  grob(aLabel = aLabel, x = x, y = y, width = width, cl = "box")
}

# used to plot the phase box on a grid
drawDetails.box <- function(x, 
                            ...) {
  
  phaseBox(x$aLabel, x$x, x$y, x$width)
}

# helpers to set phase characteristics
xDetails.box <- function(x, 
                         theta) {
  
  theLabels <- strwrap(x$aLabel, x$width)
  height <- unit(length(theLabels) + 2, "lines")
  width <- unit(12, "mm") + max(stringWidth(theLabels))
  grobX(roundrectGrob(x = x$x, y = x$y, width = width, height = height), theta)
}

yDetails.box <- function(x, 
                         theta) {
  
  theLabels <- strwrap(x$aLabel, x$width)
  height <- unit(length(theLabels) + 2, "lines")
  width <- unit(12, "mm") + max(stringWidth(theLabels))
  grobY(roundrectGrob(x = x$x, y = x$y, width = width, height = height), theta)
}

# generates all the phases with no X or Y specified
getPhaseGrobs <- function(aPhaseVector, colWidth = 40) {
  lapply(aPhaseVector, function(x) phaseGrob(x, width = colWidth))
}

# constructor for arrows that link phases
connectPhases <- function (parentPhase, 
                           childPhase) {
  
  grid.move.to(grobX(parentPhase, "north"), grobY(parentPhase, "south"))
  grid.line.to(grobX(childPhase, "south"), grobY(childPhase, "north") + unit(0.2, "mm"), 
               arrow = arrow(type = "closed", length = unit(4, "mm")), 
               gp = gpar(fill = "black"))
}

# constructor for double arrows that link start phase to a single daughter phase
marryPhases <- function (parentPhaseLeft, 
                         childPhase, 
                         parentPhaseRight) {
  
  grid.curve(grobX(parentPhaseLeft, "north"),
             grobY(parentPhaseLeft, "south"),
             grobX(childPhase, "south") - unit(5, "mm"),
             grobY(childPhase, "north"),
             inflect = TRUE,
             arrow = arrow(type = "closed",
                           angle = 30,
                           length = unit(4, "mm")),
             gp = gpar(fill = "black", lwd = 1))
  
  grid.curve(grobX(parentPhaseRight, "north"),
             grobY(parentPhaseRight, "south") ,
             grobX(childPhase, "south") + unit(5, "mm"),
             grobY(childPhase, "north") ,
             inflect = TRUE, curvature = -1,
             arrow = arrow(type = "closed",
                           angle = 30,
                           length = unit(4, "mm")),
             gp = gpar(fill = "black", lwd = 1)) 
}

# constructor for arrow that links exclusion (rightmost) phases
excludePhase <- function (parentPhase, 
                          exludedPhase) {
  
  grid.move.to(grobX(parentPhase, "east"), grobY(parentPhase, "east"))
  grid.line.to(grobX(exludedPhase, "west"), grobY(exludedPhase, "east"), 
               arrow = arrow(type = "closed", length = unit(4, "mm")), 
               gp = gpar(fill = "black"))
}

# generate simplified phase scheme
getPhaseScheme <- function(aPhaseVector) {
  aScheme <- lapply(aPhaseVector, function(x) {
    if(grepl("START_PHASE: ", x)) return("S") 
    if(grepl("EXCLUDE_PHASE: ", x)) return("E")
    return("P")
  })
  return(unlist(aScheme))
}

# generate main-line phase scheme (used to gauge length of flow chart)
getMainScheme <- function(aScheme) {
  mainScheme <- c(aScheme[1], aScheme[!(aScheme %in% "S" | aScheme %in% "E")])
  return (mainScheme)
}

# extracts phase labels without phase-definitions
getPhaseClean <- function(aPhaseVector) {
  aCleanVector <- lapply(aPhaseVector, function(x) {
    if(grepl("START_PHASE: ", x)) return(sub("START_PHASE: ", "", x)) 
    if(grepl("EXCLUDE_PHASE: ", x)) return(sub("EXCLUDE_PHASE: ", "", x))
    return(x)
  })
  return(unlist(aCleanVector))
}