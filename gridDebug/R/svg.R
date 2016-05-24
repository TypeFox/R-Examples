
##########################
# gridTreeTips()

# Add tooltip attributes to a graph node on the DL
garnishNodes <- function(elt) {
    if (inherits(elt, "gTree")) {
        garnishGrob(elt,
                    "pointer-events"="all",
                    onmouseover=paste("showTooltip(evt, '",
                      gsub("\n", "", elt$children[[2]]$label), "')",
                      sep=""),
                    onmouseout="hideTooltip()")
    } else {
        elt
    }
}

# User interface
gridTreeTips <- function(filename="Rplots.svg", ..., grid=TRUE) {
    if (!grid)
        stop("Can only add tooltips if scene tree is drawn using grid")
    gridTree(..., grid=grid)
    grid.DLapply(garnishNodes)
    grid.script(filename=system.file("js", "graphtips.js",
                  package="gridDebug"))
    gridToSVG(filename)
}


##########################
# grobBrowser()

# Add tooltip attributes to a grob on the DL
garnishAllGrobs <- function(elt) {
    if (inherits(elt, "grob")) {
        garnishGrob(elt,
                    onmousemove=paste("showTooltip(evt, '",
                      gsub("\n", " ", elt$name), "')",
                      sep=""),
                    onmouseout="hideTooltip()")
    } else {
        elt
    }
}

grobBrowser <- function(filename="Rplots.svg") {
    grid.DLapply(garnishAllGrobs)
    grid.script(filename=system.file("js", "grobtips.js",
                  package="gridDebug"))
    gridToSVG(filename)
}

