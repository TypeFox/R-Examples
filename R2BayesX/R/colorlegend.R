colorlegend <- function(color = NULL, ncol = NULL, x = NULL, breaks = NULL, 
  pos = "center", shift = 0.02, side.legend = 1L, side.ticks = 1L, range = NULL, lrange = NULL, 
  width = 0.4, height = 0.06, scale = TRUE, xlim = NULL, ylim = NULL, plot = NULL, full = FALSE,
  add = FALSE, col.border = "black", lty.border = 1L, lwd.border = 1L, ticks = TRUE, 
  at = NULL, col.ticks = "black", lwd.ticks = 1L, lty.ticks = 1L, length.ticks = 0.3, 
  labels = NULL, distance.labels = 0.8, col.labels = "black", cex.labels = 1L, 
  digits = 2L, swap = FALSE, symmetric = TRUE, xpd = NULL,
  title = NULL, side.title = 2, shift.title = c(0, 0), ...)
{
  args <- list(...)
  if(is.null(xlim)) {
    if(add)
      xlim <- par("usr")[1L:2L]
    else
      xlim <- c(0L, 1L)
  }
  if(is.null(ylim)) {
    if(add)
      ylim <- par("usr")[3L:4L]
    else
      ylim <- c(0L, 1L)
  }
  if(!side.legend %in% c(1L, 2L)) {
    warning("argument side.legend is specified wrong, set to default!")
    side.legend <- 1L
  }
  if(full) {
    scale <- FALSE
    pos = c(0, 0)
    par(xaxs = "i")
    par(yaxs = "i")
    if(side.legend < 2L) {
      width <- diff(xlim)
      height <- diff(ylim)
    } else {
      height <- diff(xlim)
      width <- diff(ylim)
    }
  }
  if(is.null(plot) || plot == TRUE) {
    plot <- TRUE
    graphics::plot.default(xlim, ylim, type = "n", xlab = "", ylab = "", 
      axes = FALSE, xlim = xlim, ylim = ylim, asp = NA)
  } else plot <- FALSE
  if(is.null(xpd))
    xpd <- FALSE
  if(xpd)
    par(xpd = xpd)
  pos2 <- NULL
  postxt <- c("bottomleft", "topleft", "topright", "bottomright",
    "left", "right", "top", "bottom", "center")
  poscheck <- pmatch(pos, postxt)
  if(all(!is.na(poscheck)) && length(poscheck) > 0) {
    pos2 <- postxt[pmatch(pos, postxt)]
    pos <- c(0, 0)
  }
  if(is.null(pos)) {
    pos <- c(0.35, 0.15)
    if(side.legend < 2L)
      pos <- rev(pos)   
  }
  limits <- list(xlim, ylim)
  pos <- opos <- c(min(limits[[1L]], na.rm = TRUE) + pos[1L] * diff(limits[[1L]]), 
    min(limits[[2L]], na.rm = TRUE) + pos[2L] * diff(limits[[2L]])) 
  if(side.legend > 1L)
    limits <- rev(limits)
  if(scale) {
    width <- width * diff(limits[[1L]])
    height <- height * diff(limits[[2L]])
  }
  if(side.legend > 1L) {
    wi <- width
    width <- height
    height <- wi
  }
  if(full)
    shift <- 0
  shift <- rep(shift, length.out = 2)
  if(is.null(pos2)) {
    xlim <- range(c(pos[1L], pos[1L] + width, pos[1L] + width, pos[1L]))
    ylim <- range(c(pos[2L], pos[2L], pos[2L] + height, pos[2L] + height))
  } else {
    pos2 <- dopos(pos2, limits, width, height, side.legend, shift)
    xlim <- pos2$xlim
    ylim <- pos2$ylim
  }
  if(!is.null(x)) {
    if(is.null(lrange)) {
      if(is.null(range)) {      
        lrange <- range(x, na.rm = TRUE)
        if(symmetric) {
          mar <- max(abs(lrange))
          lrange <- c(0 - mar, mar)
        }
      } else lrange <- range
    }
    x <- unique(na.omit(sort(x)))
  } else { 
    if(is.null(range))
      range <- xlim
    if(is.null(lrange))
      lrange <- range
  }
  if(is.null(color))
    color <- grDevices::gray.colors
  args$col <- color
  args$ncol <- ncol
  args$data <- x
  args$range <- range
  args$breaks <- breaks
  args$swap <- swap
  args$symmetric <- symmetric
  pal <- do.call(make_pal, delete.args(make_pal, args))
  if(plot || add) {
    if(!is.null(lrange)) {
      if(min(lrange) > min(pal$breaks)) 
        pal$breaks[pal$breaks <= min(lrange)] <- min(lrange)
      if(max(lrange) < max(pal$breaks))
        pal$breaks[pal$breaks >= max(lrange)] <- max(lrange)
    }
    br <- c(min(pal$breaks, lrange), pal$breaks, max(pal$breaks, lrange))
    cl <- c(head(pal$colors, 1L), pal$colors, tail(pal$colors, 1L))
    obs2legend <- function(x, xr) ((x - lrange[1L]) / diff(lrange)) * diff(xr) + xr[1L]
    if(side.legend < 2L) {
      graphics::rect(obs2legend(head(br, -1L), xlim), ylim[1L], obs2legend(tail(br, -1L), xlim),
        ylim[2L], col = cl, border = cl, xpd = xpd, lwd = 0.01)
    } else {
      graphics::rect(xlim[1L], obs2legend(head(br, -1L), ylim), xlim[2L], 
        obs2legend(tail(br, -1L), ylim), col = cl, border = cl, xpd = xpd, lwd = 0.01)
    }
    graphics::rect(xlim[1L], ylim[1L], xlim[2L], ylim[2L], 
      border = col.border, lwd = lwd.border, lty = lty.border, xpd = xpd)
    dl <- TRUE
    if(!is.null(labels) && labels == FALSE)
      dl <- FALSE
    if(ticks || dl) {
      if(is.null(at)) {
        at <- pal$breaks
        if(abs(diff(lrange / max(lrange))) / length(at) < 0.2)
          at <- seq(min(lrange), max(lrange), length.out = 3L)
      }
      if(is.null(labels))
        labels <- round(at, digits = digits)
      if(side.legend < 2L) {
        at <- obs2legend(at, xlim)
        length.ticks <- length.ticks * height
        if(any(at > max(xlim))) 
          at[at > max(xlim)] <- max(xlim)
        if(any(at < min(xlim)))
          at[at < min(xlim)] <- min(xlim)
      } else {
        at <- obs2legend(at, ylim)
        length.ticks <- length.ticks * width
        if(any(at > max(ylim))) 
          at[at > max(ylim)] <- max(ylim)
        if(any(at < min(ylim)))
          at[at < min(ylim)] <- min(ylim)
      }
      at <- unique(at)
      if(side.ticks > 1L)
        length.ticks <- (-1) * length.ticks
      nat <- length(at)
      lwd.ticks <- rep(lwd.ticks, length.out = nat)
      lty.ticks <- rep(lty.ticks, length.out = nat)
      col.ticks <- rep(col.ticks, length.out = nat)
      col.labels <- rep(col.labels, length.out = nat)
      cex.labels <- rep(cex.labels, length.out = nat)
      labels <- rep(labels, length.out = nat)
      if(!full) {
        for(i in 1L:nat) {
          if(side.legend < 2L) {
            if(ticks) {
              graphics::lines(c(at[i], at[i]), c(ylim[side.ticks], ylim[side.ticks] - length.ticks),
                lwd = lwd.ticks[i], lty = lty.ticks[i], col = col.ticks[i])
            }
            if(dl) {
              graphics::text(at[i], ylim[side.ticks] - length.ticks - (distance.labels * length.ticks * 2),
                labels = labels[i], col = col.labels[i], cex = cex.labels[i], pos = 1, ...)
            }
          } else {
            if(ticks) {
              graphics::lines(c(xlim[side.ticks], xlim[side.ticks] - length.ticks), c(at[i], at[i]),
                lwd = lwd.ticks[i], lty = lty.ticks[i], col = col.ticks[i]) 
            }
            if(dl) {
              graphics::text(xlim[side.ticks] - length.ticks - (distance.labels * length.ticks * 2), 
                at[i], labels = labels[i], col = col.labels[i], cex = cex.labels[i],
                pos = if(side.ticks < 2L) 2 else 4, ...)
            }
          }
        }
      } else {
        if(side.legend < 2L) {
          if(side.ticks < 2L) where <- 1L else where <- 3L
        } else {
          if(side.ticks < 2L) where <- 2L else where <- 4L
        }
      axis(where, at = at, labels = labels, col = col.labels, 
        tick = ticks, lty = lty.ticks, col.ticks = col.ticks, 
        lwd.ticks = lwd.ticks, cex.axis = cex.labels[1])
      }
    }
    if(!is.null(title)) {
      if(length(shift.title) < 2)
        shift.title <- c(shift.title, 0)
      if(!full) {
        xp <- xlim[1L] + shift.title[1] * diff(range(xlim)) + diff(range(xlim)) / 2
        yp <- ylim[2L] + shift.title[2] * diff(range(ylim))
        text(if(side.legend < 2) xp else yp,
          if(side.legend < 2) yp else xp, title, pos = 3,
          srt = if(side.legend == 2) 270 else 0, cex = cex.labels, xpd = xpd)
      } else {
        mtext(title, side = side.title)
      }
    }
  }

  return(invisible(pal))
}

