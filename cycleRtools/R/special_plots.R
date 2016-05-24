## Data inputs (x) are assumed to be cycleRdata objects.

Wbal_plot <- function(x, xvar, CP) {
  if (all(is.na(x$Wexp.kJ)) && !is.null(CP))
    x[, "Wexp.kJ"] <- Wbal(x, "timer.s", "power.smooth.W", CP, noisy = FALSE)
  else if (all(is.na(x$Wexp.kJ)) && is.null(CP)) {
    plot.new()
    warning("No W' balance data to plot", call. = FALSE); return()
  }

  plot(
    x = x[, xvar], y = x[, "Wexp.kJ"],
    type = "l", col = "red", lwd = 1.7,
    ylab = "W' expended (kJ)", xlab = "",
    col.lab = "red", col.main = "red",
    main = paste(
      "Max W' Expended (kJ):",
      round(max(x[, "Wexp.kJ"], na.rm = TRUE), digits = 2)
    ),
    ylim = rev(extendrange(x[, "Wexp.kJ"]))  # Invert axis.
  )
}

pwr_plot <- function(x, xvar, laps, CP) {
  col <- ifelse(laps, "lap", "blue")
  xP  <- paste("xPower:",
               round(pwr_tf(x$power.smooth.W, 4), 2))
  smth_plot(data = x, xvar, "power.W", "power.smooth.W", colour = col,
            ylab = "Power (W)", xlab = "",
            col.lab = "blue", col.main = "blue",
            main = xP, character.only = TRUE)
  if (!is.null(CP))
    abline(h = CP, lty = 2, col = "red", lwd = 1.5) # Red line seems apt.
}

elev_plot <- function(x, xvar) {

  climbing <- sum(x[x[, "delta.elev"] > 0, "delta.elev"], na.rm = TRUE)
  xvals    <- x[, xvar]
  yvals    <- x[, "elevation.m"]
  ylim     <- extendrange(yvals)

  plot(xvals, yvals, type = "l", lwd = 3,
       ylab = "Elevation (m)", xlab = "",
       ylim = ylim,
       col.lab = "green", col.main = "green",
       main = paste("Climbing:", round(climbing, 2), "m"))

  # Draw polygons.
  yvals[is.na(yvals)] <- ylim[1]

  if (any(is.na(xvals))) {
    d <- cbind.data.frame(x = xvals, y = yvals)
    lapply(split(d, na_split(d[, "x"]))[-1], function(vals) {
      area <- list(
        x = c(vals[1, "x"], vals[, "x"], vals[nrow(vals), "x"]),
        y = c(ylim[1], vals[, "y"], ylim[1])
      )
      polygon(area$x, area$y, col = "green", border = NA)
    })
  } else {
    area <- list(
      x = c(xvals[1], xvals, xvals[length(xvals)]),
      y = c(ylim[1], yvals, ylim[1])
    )
    polygon(area$x, area$y, col = "green", border = NA)
  }
}
