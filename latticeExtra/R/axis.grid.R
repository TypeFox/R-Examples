

axis.grid <-
    function(side = c("top", "bottom", "left", "right"),
             ..., ticks = c("default", "yes", "no"),
             scales, components, line.col)
{
    side <- match.arg(side)
    ticks <- match.arg(ticks)
    scales.tck <- switch(side,
                         left = , bottom = scales$tck[1], 
                         right = , top = scales$tck[2])
    ## only draw major ticks (those corresponding to labels)
    comps.major <- components
    mycomps <- components[[side]]
    if (is.list(mycomps)) {
        lab <- as.character(mycomps$labels$labels)
        if (any(lab != "")) {
            tck <- mycomps$ticks$tck
            if (any(tck * scales.tck != 0)) {
                tck <- rep(tck, length = length(lab))
                comps.major[[side]]$ticks$tck <- ifelse(lab == "", NA, tck)
            }
        }
    } else {
        ticks <- "no"
    }
    ## use axis.text for ticks because axis.line$col might be transparent
    axis.text <- trellis.par.get("axis.text")
    axis.default(side, scales = scales, ticks = ticks,
                 components = comps.major, ...,
                 line.col = axis.text$col)
    ## now draw grid lines corresponding to axis ticks.
    ## can only do this with the bottom and right sides;
    ## otherwise the strip viewports are current, not panel.
    if (side %in% c("top", "left"))
        return()
    if (scales$draw == FALSE)
        return()
    ref.line <- trellis.par.get("reference.line")
    if (side == "bottom") {
        tck <- abs(mycomps$ticks$tck)
        panel.refline(v = mycomps$ticks$at, 
                      lwd = ref.line$lwd * tck,
                      alpha = ref.line$alpha * tck / max(tck, na.rm = TRUE))
    }
    if (side == "right") {
        if (!is.list(mycomps))
            mycomps <- components[["left"]]
        tck <- abs(mycomps$ticks$tck)
        panel.refline(h = mycomps$ticks$at, 
                      lwd = ref.line$lwd * tck,
                      alpha = ref.line$alpha * tck / max(tck, na.rm = TRUE))
    }
}
