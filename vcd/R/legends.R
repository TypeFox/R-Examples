legend_resbased <- function(fontsize = 12,
                            fontfamily = "",
                            x = unit(1, "lines"),
                            y = unit(0.1, "npc"),
                            height = unit(0.8, "npc"),
                            width = unit(0.7, "lines"),
			    digits = 2,
			    check_overlap = TRUE,
                            text = NULL,
                            steps = 200,
                            ticks = 10,
                            pvalue = TRUE,
                            range = NULL) {

  if(!is.unit(x)) x <- unit(x, "native")
  if(!is.unit(y)) y <- unit(y, "npc")
  if(!is.unit(width)) width <- unit(width, "lines")
  if(!is.unit(height)) height <- unit(height, "npc")

  function(residuals, shading, autotext) {
    res <- as.vector(residuals)

    if(is.null(text)) text <- autotext
    p.value <- attr(shading, "p.value")
    legend <- attr(shading, "legend")

    if (all(residuals == 0)) {
      pushViewport(viewport(x = x, y = y, just = c("left", "bottom"),
                            default.units = "native",
                            height = height, width = width))
      grid.lines(y = 0.5)
      grid.text(0, x = unit(1, "npc") + unit(0.8, "lines"),  y = 0.5,
                gp = gpar(fontsize = fontsize, fontfamily = fontfamily))
      warning("All residuals are zero.")

    } else {
        if (is.null(range))
            range <- range(res)
        if (length(range) != 2)
            stop("Range must have length two!")
        if (is.na(range[1]))
            range[1] <- min(res)
        if (is.na(range[2]))
            range[2] <- max(res)

        pushViewport(viewport(x = x, y = y, just = c("left", "bottom"),
                              yscale = range, default.units = "native",
                              height = height, width = width))


      if(is.null(legend$col.bins)) {
        col.bins <- seq(range[1], range[2], length = steps)
        at <- NULL
      } else {
        col.bins <- sort(unique(c(legend$col.bins, range)))
        col.bins <- col.bins[col.bins <= range[2] & col.bins >= range[1]]
        at <- col.bins
      }
      y.pos <- col.bins[-length(col.bins)]
      y.height <- diff(col.bins)

      grid.rect(x = unit(rep.int(0, length(y.pos)), "npc"),
                y = y.pos,
                height = y.height, default.units = "native",
                gp = gpar(fill = shading(y.pos + 0.5 * y.height)$fill, col = 0),
                just = c("left", "bottom"))

      grid.rect(gp = gpar(fill = "transparent"))

      if(is.null(at))
            at <- seq(from = head(col.bins, 1), to = tail(col.bins, 1),
                      length = ticks)
      lab <- format(round(at, digits = digits), nsmall = digits)
      tw <- lab[which.max(nchar(lab))]

      ## if(is.null(at))
      ##   at <- seq(from = head(col.bins, 1), to = tail(col.bins, 1), length = ticks)
      ## tw <- paste(rep("4", digits), collapse = "")
      ## if (any(trunc(at) != at))
      ##   tw <- paste(tw, ".", sep = "")
      ## if (any(at < 0))
      ##   tw <- paste(tw, "-", sep = "")

      grid.text(format(signif(at, digits = digits)),
                x = unit(1, "npc") + unit(0.8, "lines") + unit(1, "strwidth", tw),
                y = at,
                default.units = "native",
                just = c("right", "center"),
                gp = gpar(fontsize = fontsize, fontfamily = fontfamily),
                check.overlap = check_overlap)
      grid.segments(x0 = unit(1, "npc"), x1 = unit(1,"npc") + unit(0.5, "lines"),
                    y0 = at, y1 = at, default.units = "native")

    }

    popViewport(1)

    grid.text(text, x = x, y = unit(1, "npc") - y + unit(1, "lines"),
              gp = gpar(fontsize = fontsize, fontfamily = fontfamily,
              lineheight = 0.8),
              just = c("left", "bottom")
              )
    if(!is.null(p.value) && pvalue) {
      grid.text(paste("p-value =\n", format.pval(p.value), sep = ""),
                x = x,
                y = y - unit(1, "lines"),
                gp = gpar(fontsize = fontsize, fontfamily = fontfamily,
                lineheight = 0.8),
                just = c("left", "top"))
    }
  }
}
class(legend_resbased) <- "grapcon_generator"

legend_fixed <- function(fontsize = 12,
                         fontfamily = "",
                         x = unit(1, "lines"),
                         y = NULL,
                         height = NULL,
                         width = unit(1.5, "lines"),
                         steps = 200,
			 digits = 1,
                         space = 0.05,
                         text = NULL,
                         range = NULL) {

  if(!is.unit(x)) x <- unit(x, "native")
  if(!is.unit(y) && !is.null(y)) y <- unit(y, "npc")
  if(!is.unit(width)) width <- unit(width, "lines")
  if(!is.unit(height) && !is.null(height)) height <- unit(height, "npc")

  function(residuals, shading, autotext) {
    res <- as.vector(residuals)

    if(is.null(text)) text <- autotext

    if (is.null(y)) y <- unit(1, "strwidth", text) + unit(1, "lines")
    if (is.null(height)) height <- unit(1, "npc") - y

    pushViewport(viewport(x = x, y = y, just = c("left", "bottom"),
                          yscale = c(0,1), default.units = "npc",
                          height = height, width = width))

    p.value <- attr(shading, "p.value")
    legend <- attr(shading, "legend")

        if (is.null(range))
            range <- range(res)
        if (length(range) != 2)
            stop("Range must have length two!")
        if (is.na(range[1]))
            range[1] <- min(res)
        if (is.na(range[2]))
            range[2] <- max(res)

    if(is.null(legend$col.bins)) {
      col.bins <- seq(range[1], range[2], length = steps)
    } else {
      col.bins <- sort(unique(c(legend$col.bins, range)))
      col.bins <- col.bins[col.bins <= range[2] & col.bins >= range[1]]
    }
    l <- length(col.bins)
    y.height <- (1 - (l - 2) * space) / (l - 1)
    y.pos <- cumsum(c(0, rep(y.height + space, l - 2)))
    res <- col.bins[-l] + diff(col.bins) / 2

    grid.rect(x = unit(rep.int(0, length(y.pos)), "npc"),
              y = y.pos,
              height = y.height,
              default.units = "npc",
              gp = shading(res),
              just = c("left", "bottom"))
    numbers <- format(col.bins, nsmall = digits, digits = digits)
    wid <- unit(1, "strwidth", format(max(abs(col.bins)),
                                      nsmall = digits, digits = digits))
    grid.text(numbers[-l],
              x = unit(1, "npc") + unit(0.6, "lines") + wid,
              y = y.pos, gp = gpar(fontsize = fontsize, fontfamily = fontfamily),
              default.units = "npc", just = c("right", "bottom"))
    grid.text(numbers[-1],
              x = unit(1, "npc") + unit(0.6, "lines") + wid,
              y = y.pos + y.height, gp = gpar(fontsize = fontsize, fontfamily = fontfamily),
              default.units = "npc", just = c("right", "top"))
    wid2 <- unit(1, "strwidth", format(max(abs(trunc(col.bins))))) +
      unit(0.3, "strwidth", ".")
    grid.segments(x0 = unit(1, "npc") + wid2 + unit(0.6, "lines"),
                  x1 = unit(1, "npc") + wid2 + unit(0.6, "lines"),
                  y0 = unit(y.pos, "npc") + 1.5 * unit(1, "strheight", "-44.4"),
                  y1 = unit(y.pos + y.height, "npc") - 1.5 * unit(1, "strheight", "-44.4")
                  )

    popViewport(1)
    grid.text(text, x = x + 0.5 * width, y = 0,
              gp = gpar(fontsize = fontsize, fontfamily = fontfamily,
              lineheight = 0.8),
              just = c("left", "top"),
              rot = 90
              )
  }
}
class(legend_fixed) <- "grapcon_generator"
