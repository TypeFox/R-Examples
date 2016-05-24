

panel.scaleArrow <-
    function(x = unit(0:1, "npc"), y = unit(0:1, "npc"),
             default.units = "npc",
             digits = 0, append = "", label = NULL,
             angle = 30, length = 0.5, unit = "char",
             type = "open", ends = "both",
             ...,
             col = add.line$col, fill = col, alpha = add.line$alpha,
             lty = add.line$lty, lwd = add.line$lwd,
             col.text = add.text$col, alpha.text = add.text$alpha)
{
    add.line <- trellis.par.get("add.line")
    add.text <- trellis.par.get("add.text")
    if (!is.unit(x)) x <- unit(x, default.units)
    if (!is.unit(y)) y <- unit(y, default.units)
    x <- rep(x, length.out = 2)
    y <- rep(y, length.out = 2)
    xnat <- convertX(x, "native", valueOnly = TRUE)
    ynat <- convertY(y, "native", valueOnly = TRUE)
    panel.arrows(xnat[1], ynat[1], xnat[2], ynat[2], ends = ends,
                 angle = angle, length = length, unit = unit, type = type,
                 col = col, fill = fill, alpha = alpha, lty = lty, lwd = lwd)
    #grid.lines(x = x, y = y, default.units = default.units,
    #           arrow = arrow, gp = gp)
    #xnat <- convertX(x, "native", valueOnly = TRUE)
    #ynat <- convertX(y, "native", valueOnly = TRUE)
    d <- sqrt(diff(xnat)^2 + diff(ynat)^2)
    d <- round(d, digits = digits)
    if (is.null(label)) label <- paste(d, append, sep = "")
    panel.text(x = 0.5 * (xnat[1] + xnat[2]), y = 0.5 * (ynat[1] + ynat[2]),
              labels = label, col = col.text, alpha = alpha.text, ...)
}
