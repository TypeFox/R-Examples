setMethod("plot", 
          signature(x = "RCLRMIX", y = "missing"),
function(x,
  y,
  nrow = 1,
  ncol = 1,
  cex = 0.8,
  fg = "black",
  lty = "solid",
  lwd = 1,
  pty = "m",
  tcl = 0.5,
  plot.cex = 0.8,
  plot.pch = 19, ...)
{
  if (missing(x)) {
    stop(sQuote("x"), " object of class RCLRMIX is requested!", call. = FALSE)
  }

  if (!is.wholenumber(nrow)) {
    stop(sQuote("nrow"), " integer is requested!", call. = FALSE)
  }

  if (nrow < 1) {
    stop(sQuote("nrow"), " must be greater than 0!", call. = FALSE)
  }

  if (!is.wholenumber(ncol)) {
    stop(sQuote("ncol"), " integer is requested!", call. = FALSE)
  }

  if (ncol < 1) {
    stop(sQuote("ncol"), " must be greater than 0!", call. = FALSE)
  }

  d <- ncol(x@x@Dataset[[x@pos]])

  nrow <- max(1, nrow)
  ncol <- max(1, ncol)

  N <- d * (d - 1) / 2
  
  opar <- par(mfrow = c(nrow, ncol),
    cex = cex,
    cex.axis = 1.0,
    fg = fg,
    lty = lty,
    lwd = lwd,
    mar = c(1.2, 1.2, 1.2, 1.2),
    oma = c(1.2, 0.2, 0.2, 0.2),
    pty = pty,
    tcl = tcl, ...)

  par(oma = c(1 + 0.2, 0.2, 0.2, 0.2))

  ey <- as.matrix(x@x@Dataset[[x@pos]]); ep <- as.numeric(x@Zp) - 1
  
  error <- is.error(x@Zt, x@Zp)
  
  ramp <- colorRamp(colors = c("magenta", "blue", "cyan", "green", "yellow", "red"),
    space = "rgb",
    interpolate = "linear")
    
  s <- length(levels(x@Zp))
  
  zlim <- c(0, s - 1); zmax <- zlim[2]
      
  plot.col <- rgb(ramp(ep / zmax), maxColorValue = 255);
  
  if (.Device == "tikz output") {
    if (sum(error) == 0) {
      legend <- paste("$", 1:s, "$", sep = "")
      legend.col <- rgb(ramp(zlim[1]:zlim[2] / zmax), maxColorValue = 255)
      legend.pch <- rep(plot.pch, s)
    }
    else {
      legend <- c(paste("$", 1:s, "$", sep = ""), "$\\mathrm{Error}$")
      legend.col <- c(rgb(ramp(zlim[1]:zlim[2] / zmax), maxColorValue = 255), "black")
      legend.pch <- c(rep(plot.pch, s), 1)
    }
  }
  else {
    if (sum(error) == 0) {
      legend <- paste(bquote(.(1:s)), sep = "")
      legend.col <- rgb(ramp(zlim[1]:zlim[2] / zmax), maxColorValue = 255)
      legend.pch <- rep(plot.pch, s)
    }
    else {
      legend <- c(paste(bquote(.(1:s)), sep = ""), "Error")
      legend.col <- c(rgb(ramp(zlim[1]:zlim[2] / zmax), maxColorValue = 255), "black")
      legend.pch <- c(rep(plot.pch, s), 1)
    }  
  }     
  
  if (N > 0) {
    figno <- 0

    for (i in 1:(d - 1)) {
      for (j in (i + 1):d) {
        plot(x = ey[which(error != 1), i],
          y = ey[which(error != 1), j],
          type = "p",
          main = "",
          sub = "",
          xlab = "",
          ylab = "",
          xlim = range(ey[, i]),
          ylim = range(ey[, j]),
          col = plot.col[which(error != 1)],
          axes = FALSE,
          lwd = 1,
          cex = plot.cex,
          pch = plot.pch)
          
        points(x = ey[which(error == 1), i],
          y = ey[which(error == 1), j],
          type = "p",
          main = "",
          sub = "",
          xlab = "",
          ylab = "",
          col = "black",
          lwd = 1,
          cex = plot.cex * 2,
          pch = 1)          

        box(col = fg, lty = "solid", lwd = 1)

        axis(side = 3,
          outer = FALSE,
          lty = "solid",
          lwd = 1,
          hadj = 0.5,
          padj = 1.0)

        axis(side = 2,
          outer = FALSE,
          lty = "solid",
          lwd = 1,
          hadj = 0.5,
          padj = 1.0)

        if (.Device == "tikz output") {
          text <- paste("$y_{", i, "}$", "$\\; - \\;$", "$y_{", j, "}$", sep = "")
        }
        else {
          text <- bquote(y[.(i)] - y[.(j)])
        }

        mtext(text = text,
          side = 1,
          line = 0,
          outer = FALSE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)
          
        figno <- figno + 1

        if ((figno == nrow * ncol) || ((i == d - 1) && (j == d))) {
          par(fig = c(0, 1, 0, 1), 
            oma = c(0, 0, 0, 0), 
            mar = c(0, 0, 0, 0), 
            new = TRUE)
          
          plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
          
          legend("bottom", 
            legend = legend,
            col = legend.col,
            lty = 0,
            pch = legend.pch,
            bty = "n",
            cex = 1.0,
            y.intersp = 0,
            horiz = TRUE,
            inset = c(0, 0),
            xpd = TRUE)
  
          par(mfrow = c(nrow, ncol),
            cex = cex,
            cex.axis = 1.0,
            fg = fg,
            lty = lty,
            lwd = lwd,
            mar = c(1.2, 1.2, 1.2, 1.2),
            oma = c(1.2, 0.2, 0.2, 0.2),
            pty = pty,
            tcl = tcl, ...)
    
          par(oma = c(1 + 0.2, 0.2, 0.2, 0.2))

          figno <- 0
        }
      }
    }
  }
  else {
    plot(x = ey[which(error != 1), 1],
      y = ep[which(error != 1)] + 1,
      type = "p",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      xlim = range(ey[, 1]),
      ylim = range(ep + 1),      
      col = plot.col[which(error != 1)],
      axes = FALSE,
      lwd = 1,
      cex = plot.cex,
      pch = plot.pch)
          
    points(x = ey[which(error == 1), 1],
      y = ep[which(error == 1)] + 1,
      type = "p",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      col = "black",
      lwd = 1,
      cex = plot.cex * 2,
      pch = 1)       

    box(col = fg, lty = "solid", lwd = 1)

    axis(side = 3,
      outer = FALSE,
      lty = "solid",
      lwd = 1,
      hadj = 0.5,
      padj = 1.0)

    axis(side = 2,
      outer = FALSE,
      lty = "solid",
      lwd = 1,
      hadj = 0.5,
      padj = 1.0)

    if (.Device == "tikz output") {
      text <- paste("$y_{1}$", "$\\; - \\;$", "$Z_{p}(y_{1})$", sep = "")
    }
    else {
      text <- bquote(y[1] - Z[p](y[1]))
    }

    mtext(text = text,
      side = 1,
      line = 0,
      outer = FALSE,
      adj = 0.5,
      padj = 0.2,
      cex = cex)

    par(fig = c(0, 1, 0, 1), 
      oma = c(0, 0, 0, 0), 
      mar = c(0, 0, 0, 0), 
      new = TRUE)
          
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
          
    legend("bottom", 
      legend = legend,
      col = legend.col,
      lty = 0,
      pch = legend.pch,
      bty = "n",
      cex = 1.0,
      y.intersp = 0,
      horiz = TRUE,
      inset = c(0, 0),
      xpd = TRUE)
      
    par(mfrow = c(nrow, ncol),
      cex = cex,
      cex.axis = 1.0,
      fg = fg,
      lty = lty,
      lwd = lwd,
      mar = c(1.2, 1.2, 1.2, 1.2),
      oma = c(1.2, 0.2, 0.2, 0.2),
      pty = pty,
      tcl = tcl, ...)
      
    par(oma = c(1 + 0.2, 0.2, 0.2, 0.2))
  }
  
  rm(list = ls()[!(ls() %in% c("opar"))])

  invisible(opar)
}) # plot         

setMethod("plot", 
          signature(x = "RCLRMVNORM", y = "missing"),
function(x,
  y,
  nrow = 1,
  ncol = 1,
  cex = 0.8,
  fg = "black",
  lty = "solid",
  lwd = 1,
  pty = "m",
  tcl = 0.5,
  plot.cex = 0.8,
  plot.pch = 19, ...)
{
  if (missing(x)) {
    stop(sQuote("x"), " object of class RCLRMVNORM is requested!", call. = FALSE)
  }

  if (!is.wholenumber(nrow)) {
    stop(sQuote("nrow"), " integer is requested!", call. = FALSE)
  }

  if (nrow < 1) {
    stop(sQuote("nrow"), " must be greater than 0!", call. = FALSE)
  }

  if (!is.wholenumber(ncol)) {
    stop(sQuote("ncol"), " integer is requested!", call. = FALSE)
  }

  if (ncol < 1) {
    stop(sQuote("ncol"), " must be greater than 0!", call. = FALSE)
  }

  d <- ncol(x@x@Dataset[[x@pos]])

  nrow <- max(1, nrow)
  ncol <- max(1, ncol)

  N <- d * (d - 1) / 2
  
  opar <- par(mfrow = c(nrow, ncol),
    cex = cex,
    cex.axis = 1.0,
    fg = fg,
    lty = lty,
    lwd = lwd,
    mar = c(1.2, 1.2, 1.2, 1.2),
    oma = c(1.2, 0.2, 0.2, 0.2),
    pty = pty,
    tcl = tcl, ...)

  par(oma = c(1 + 0.2, 0.2, 0.2, 0.2))

  ey <- as.matrix(x@x@Dataset[[x@pos]]); ep <- as.numeric(x@Zp) - 1
  
  error <- is.error(x@Zt, x@Zp)
  
  ramp <- colorRamp(colors = c("magenta", "blue", "cyan", "green", "yellow", "red"),
    space = "rgb",
    interpolate = "linear")
    
  s <- length(levels(x@Zp))
  
  zlim <- c(0, s - 1); zmax <- zlim[2]
      
  plot.col <- rgb(ramp(ep / zmax), maxColorValue = 255);
  
  if (.Device == "tikz output") {
    if (sum(error) == 0) {
      legend <- paste("$", 1:s, "$", sep = "")
      legend.col <- rgb(ramp(zlim[1]:zlim[2] / zmax), maxColorValue = 255)
      legend.pch <- rep(plot.pch, s)
    }
    else {
      legend <- c(paste("$", 1:s, "$", sep = ""), "$\\mathrm{Error}$")
      legend.col <- c(rgb(ramp(zlim[1]:zlim[2] / zmax), maxColorValue = 255), "black")
      legend.pch <- c(rep(plot.pch, s), 1)
    }
  }
  else {
    if (sum(error) == 0) {
      legend <- paste(bquote(.(1:s)), sep = "")
      legend.col <- rgb(ramp(zlim[1]:zlim[2] / zmax), maxColorValue = 255)
      legend.pch <- rep(plot.pch, s)
    }
    else {
      legend <- c(paste(bquote(.(1:s)), sep = ""), "Error")
      legend.col <- c(rgb(ramp(zlim[1]:zlim[2] / zmax), maxColorValue = 255), "black")
      legend.pch <- c(rep(plot.pch, s), 1)
    }  
  }    
  
  if (N > 0) {
    figno <- 0

    for (i in 1:(d - 1)) {
      for (j in (i + 1):d) {
        plot(x = ey[which(error != 1), i],
          y = ey[which(error != 1), j],
          type = "p",
          main = "",
          sub = "",
          xlab = "",
          ylab = "",
          xlim = range(ey[, i]),
          ylim = range(ey[, j]),
          col = plot.col[which(error != 1)],
          axes = FALSE,
          lwd = 1,
          cex = plot.cex,
          pch = plot.pch)
          
        points(x = ey[which(error == 1), i],
          y = ey[which(error == 1), j],
          type = "p",
          main = "",
          sub = "",
          xlab = "",
          ylab = "",
          col = "black",
          lwd = 1,
          cex = plot.cex * 2,
          pch = 1)          

        box(col = fg, lty = "solid", lwd = 1)

        axis(side = 3,
          outer = FALSE,
          lty = "solid",
          lwd = 1,
          hadj = 0.5,
          padj = 1.0)

        axis(side = 2,
          outer = FALSE,
          lty = "solid",
          lwd = 1,
          hadj = 0.5,
          padj = 1.0)

        if (.Device == "tikz output") {
          text <- paste("$y_{", i, "}$", "$\\; - \\;$", "$y_{", j, "}$", sep = "")
        }
        else {
          text <- bquote(y[.(i)] - y[.(j)])
        }

        mtext(text = text,
          side = 1,
          line = 0,
          outer = FALSE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)
          
        figno <- figno + 1

        if ((figno == nrow * ncol) || ((i == d - 1) && (j == d))) {
          par(fig = c(0, 1, 0, 1), 
            oma = c(0, 0, 0, 0), 
            mar = c(0, 0, 0, 0), 
            new = TRUE)
          
          plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
          
          legend("bottom", 
            legend = legend,
            col = legend.col,
            lty = 0,
            pch = legend.pch,
            bty = "n",
            cex = 1.0,
            y.intersp = 0,
            horiz = TRUE,
            inset = c(0, 0),
            xpd = TRUE)
  
          par(mfrow = c(nrow, ncol),
            cex = cex,
            cex.axis = 1.0,
            fg = fg,
            lty = lty,
            lwd = lwd,
            mar = c(1.2, 1.2, 1.2, 1.2),
            oma = c(1.2, 0.2, 0.2, 0.2),
            pty = pty,
            tcl = tcl, ...)
    
          par(oma = c(1 + 0.2, 0.2, 0.2, 0.2))

          figno <- 0
        }
      }
    }
  }
  else {
    plot(x = ey[which(error != 1), 1],
      y = ep[which(error != 1)] + 1,
      type = "p",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      xlim = range(ey[, 1]),
      ylim = range(ep + 1),      
      col = plot.col[which(error != 1)],
      axes = FALSE,
      lwd = 1,
      cex = plot.cex,
      pch = plot.pch)
          
    points(x = ey[which(error == 1), 1],
      y = ep[which(error == 1)] + 1,
      type = "p",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      col = "black",
      lwd = 1,
      cex = plot.cex * 2,
      pch = 1)       

    box(col = fg, lty = "solid", lwd = 1)

    axis(side = 3,
      outer = FALSE,
      lty = "solid",
      lwd = 1,
      hadj = 0.5,
      padj = 1.0)

    axis(side = 2,
      outer = FALSE,
      lty = "solid",
      lwd = 1,
      hadj = 0.5,
      padj = 1.0)

    if (.Device == "tikz output") {
      text <- paste("$y_{1}$", "$\\; - \\;$", "$Z_{p}(y_{1})$", sep = "")
    }
    else {
      text <- bquote(y[1] - Z[p](y[1]))
    }

    mtext(text = text,
      side = 1,
      line = 0,
      outer = FALSE,
      adj = 0.5,
      padj = 0.2,
      cex = cex)

    par(fig = c(0, 1, 0, 1), 
      oma = c(0, 0, 0, 0), 
      mar = c(0, 0, 0, 0), 
      new = TRUE)
          
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
          
    legend("bottom", 
      legend = legend,
      col = legend.col,
      lty = 0,
      pch = legend.pch,
      bty = "n",
      cex = 1.0,
      y.intersp = 0,
      horiz = TRUE,
      inset = c(0, 0),
      xpd = TRUE)
      
    par(mfrow = c(nrow, ncol),
      cex = cex,
      cex.axis = 1.0,
      fg = fg,
      lty = lty,
      lwd = lwd,
      mar = c(1.2, 1.2, 1.2, 1.2),
      oma = c(1.2, 0.2, 0.2, 0.2),
      pty = pty,
      tcl = tcl, ...)
      
    par(oma = c(1 + 0.2, 0.2, 0.2, 0.2))
  }
  
  rm(list = ls()[!(ls() %in% c("opar"))])

  invisible(opar)
}) # plot
