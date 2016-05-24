setMethod("plot", 
          signature(x = "REBMIX", y = "missing"),
function(x,
  y,
  pos = 1,
  what = c("density"),
  nrow = 1,
  ncol = 1,
  npts = 200,
  n = 200,
  cex = 0.8,
  fg = "black",
  lty = "solid",
  lwd = 1,
  pty = "m",
  tcl = 0.5,
  plot.cex = 0.8,
  plot.pch = 19,
  contour.drawlabels = FALSE,
  contour.labcex = 0.8,
  contour.method = "flattest",
  contour.nlevels = 12, ...)
{
  if (missing(x)) {
    stop(sQuote("x"), " object of class REBMIX is requested!", call. = FALSE)
  }

  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(pos) <- 1

  if ((pos < 1) || (pos > nrow(x@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(x@summary), "!", call. = FALSE)
  }
  
  if (!is.character(what)) {
    stop(sQuote("what"), " character vector is requested!", call. = FALSE)
  }

  what <- match.arg(what, .rebmix.plot$what, several.ok = TRUE)  

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

  if (!is.wholenumber(npts)) {
    stop(sQuote("npts"), " integer is requested!", call. = FALSE)
  }

  if (npts < 1) {
    stop(sQuote("npts"), " must be greater than 0!", call. = FALSE)
  }
  
  if (!is.wholenumber(n)) {
    stop(sQuote("n"), " integer is requested!", call. = FALSE)
  }

  if (n < 1) {
    stop(sQuote("n"), " must be greater than 0!", call. = FALSE)
  }  

  ni <- ncol(x@summary)

  Theta <- .extractThetaA(x@w[[pos]], x@Theta[[pos]])

  d <- nrow(Theta)

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

  C <- x@summary[pos, "Preprocessing"]

  if (.Device == "tikz output") {
    item <- list()

    item[[1]] <- "$\\textrm{Dataset}$"
    item[[2]] <- "$\\; = \\;$"
    item[[3]] <- paste("$\\textrm{", x@summary[pos, "Dataset"], "}$, ", sep = "")

    item[[4]] <- "$\\textrm{Preprocessing}$"
    item[[5]] <- "$\\; = \\;$"
    item[[6]] <- paste("$\\textrm{", x@summary[pos, "Preprocessing"], "}$, ", sep = "")

    item[[7]] <- "$\\textrm{Restraints}$"
    item[[8]] <- "$\\; = \\;$"
    item[[9]] <- paste("$\\textrm{", x@summary[pos, "Restraints"], "}$, ", sep = "")

    item[[10]] <- "$D$"
    item[[11]] <- ""
    item[[12]] <- ""

    item[[13]] <- "$c_{\\mathrm{max}}$"
    item[[14]] <- "$\\; = \\;$"
    item[[15]] <- paste("$", x@summary[pos, "cmax"], "$, ", sep = "")

    item[[16]] <- "$a_{\\mathrm{r}}$"
    item[[17]] <- "$\\; = \\;$"
    item[[18]] <- paste("$", as.number(x@summary[pos, "ar"]), "$, ", sep = "")

    item[[19]] <- "$c$"
    item[[20]] <- "$\\; = \\;$"
    item[[21]] <- paste("$", x@summary[pos, "c"], "$, ", sep = "")

    item[[22]] <- ""
    item[[23]] <- ""
    item[[24]] <- ""

    if (C == .rebmix$Preprocessing[1]) {
      item[[25]] <- "$v$"
      item[[26]] <- "$\\; = \\;$"
      item[[27]] <- paste("$", x@summary[pos, "v/k"], "$, ", sep = "")
    }
    else
    if (C == .rebmix$Preprocessing[2]) {
      item[[25]] <- "$v$"
      item[[26]] <- "$\\; = \\;$"
      item[[27]] <- paste("$", x@summary[pos, "v/k"], "$, ", sep = "")
    }
    else
    if (C == .rebmix$Preprocessing[3]) {
      item[[25]] <- "$k$"
      item[[26]] <- "$\\; = \\;$"
      item[[27]] <- paste("$", x@summary[pos, "v/k"], "$, ", sep = "")
    }

    item[[28]] <- paste("$\\mathrm{", x@summary[pos, "Criterion"], "}$", sep = "")
    item[[29]] <- "$\\; = \\;$"
    item[[30]] <- paste("$", as.number(x@summary[pos, "IC"]), "$, ", sep = "")

    item[[31]] <- "$\\mathrm{log}\\, L$"
    item[[32]] <- "$\\; = \\;$"
    item[[33]] <- paste("$", as.number(x@summary[pos, "logL"]), "$.", sep = "")

    i <- 1; legend <- list(); legend[[i]] <- item[[1]]

    for (j in c(2:9, 13:21, 25:33)) {
      legendwidth <- strwidth(paste(legend[[i]], item[[j]], sep = ""), units = "figure", cex = 1.0)

      if (legendwidth > ncol) {
        i <- i + 1; legend[[i]] <- item[[j]]
      }
      else {
        legend[[i]] <- paste(legend[[i]], item[[j]], sep = "")
      }
    }
  }
  else {
    item <- list()

    item[[1]] <- "Dataset"
    item[[2]] <- " = "
    item[[3]] <- paste(x@summary[pos, "Dataset"], ", ", sep = "")

    item[[4]] <- "Preprocessing"
    item[[5]] <- " = "
    item[[6]] <- paste(x@summary[pos, "Preprocessing"], ", ", sep = "")

    item[[7]] <- "Restraints"
    item[[8]] <- " = "
    item[[9]] <- paste(x@summary[pos, "Restraints"], ", ", sep = "")

    item[[10]] <- "D"
    item[[11]] <- ""
    item[[12]] <- ""

    item[[13]] <- bquote(c[max])
    item[[14]] <- " = "
    item[[15]] <- paste(x@summary[pos, "cmax"], ", ", sep = "")

    item[[16]] <- bquote(a[r])
    item[[17]] <- " = "
    item[[18]] <- paste(as.number(x@summary[pos, "ar"]), ", ", sep = "")

    item[[19]] <- "c"
    item[[20]] <- " = "
    item[[21]] <- paste(x@summary[pos, "c"], ", ", sep = "")

    item[[22]] <- ""
    item[[23]] <- ""
    item[[24]] <- ""

    if (C == .rebmix$Preprocessing[1]) {
      item[[25]] <- "v"
      item[[26]] <- " = "
      item[[27]] <- paste(x@summary[pos, "v/k"], ", ", sep = "")
    }
    else
    if (C == .rebmix$Preprocessing[2]) {
      item[[25]] <- "v"
      item[[26]] <- " = "
      item[[27]] <- paste(x@summary[pos, "v/k"], ", ", sep = "")
    }
    else
    if (C == .rebmix$Preprocessing[3]) {
      item[[25]] <- "k"
      item[[26]] <- " = "
      item[[27]] <- paste(x@summary[pos, "v/k"], ", ", sep = "")
    }

    item[[28]] <- as.character(x@summary[pos, "Criterion"])
    item[[29]] <- " = "
    item[[30]] <- paste(as.number(x@summary[pos, "IC"]), ", ", sep = "")

    item[[31]] <- "log L"
    item[[32]] <- " = "
    item[[33]] <- paste(as.number(x@summary[pos, "logL"]), ".", sep = "")

    i <- 1; legend <- list(); legend[[i]] <- bquote(.(item[[1]]))

    for (j in c(2:9, 13:21, 25:33)) {
      legendwidth <- strwidth(bquote(paste(.(legend[[i]]), .(item[[j]]), sep = "")), units = "figure", cex = 1.0)

      if (legendwidth > ncol) {
        i <- i + 1; legend[[i]] <- item[[j]]
      }
      else {
        legend[[i]] <- bquote(paste(.(legend[[i]]), .(item[[j]]), sep = ""))
      }
    }
  }

  par(oma = c(length(legend) + 0.2, 0.2, 0.2, 0.2))

  Dataset <- as.character(x@summary[pos, "Dataset"])

  ey <- as.matrix(x@Dataset[[which(names(x@Dataset) == x@summary[pos, "Dataset"])]])

  y0 <- array(data = 0.0, dim = d, dimnames = NULL)
  h <- array(data = 0.0, dim = d, dimnames = NULL)

  lim <- array(data = 0.0, dim = c(2, d), dimnames = NULL)

  py <- list(d)

  Variables <- match.arg(x@Variables, .rebmix$Variables, several.ok = TRUE)
  
  pdf <- match.arg(x@pdf, .rebmix$pdf, several.ok = TRUE)
  
  for (i in 1:d) {
    if (C == .rebmix$Preprocessing[1]) {
      k <- as.numeric(x@summary[pos, "v/k"])
      y0[i] <- as.numeric(x@summary[pos, paste("y0", if (d > 1) i, sep = "")])
      h[i] <- as.numeric(x@summary[pos, paste("h", if (d > 1) i, sep = "")])

      lim[, i] <- range(ey[, i], finite = TRUE)
    }
    else
    if (C == .rebmix$Preprocessing[2]) {
      h[i] <- as.numeric(x@summary[pos, paste("h", if (d > 1) i, sep = "")])

      lim[, i] <- range(ey[, i], finite = TRUE)
    }
    else
    if (C == .rebmix$Preprocessing[3]) {
      k <- as.numeric(x@summary[pos, "v/k"])

      h[i] <- as.numeric(x@summary[pos, paste("h", if (d > 1) i, sep = "")])

      lim[, i] <- range(ey[, i], finite = TRUE)
    }
    else {
      lim[, i] <- c(0.0, 1.0)
    }

    if (Variables[i] == .rebmix$Variables[2]) {
      py[[i]] <- seq(from = lim[1, i], to = lim[2, i], by = 1.0)
    }
    else {
      py[[i]] <- seq(from = lim[1, i], to = lim[2, i], length.out = npts)
    }
  }

  w <- as.numeric(x@w[[pos]])

  if (N > 0) {
    if (any(match(.rebmix.plot$what[1], what, nomatch = 0))) {  
      figno <- 0
    
      ramp <- colorRamp(colors = c("magenta", "blue", "cyan", "green", "yellow", "red"),
        space = "rgb",
        interpolate = "linear")

      for (i in 1:(d - 1)) {
        for (j in (i + 1):d) {
          if (C == .rebmix$Preprocessing[1]) {
            edens <- .densHistogram.xy(k, ey[, i], ey[, j], y0[i], y0[j], h[i], h[j], Variables[i], Variables[j], pdf[i], pdf[j])
          }
          else
          if (C == .rebmix$Preprocessing[2]) {
            edens <- .densParzenWindow.xy(ey[, i], ey[, j], h[i], h[j], n)
          }
          else
          if (C == .rebmix$Preprocessing[3]) {
            edens <- .densKNearestNeighbour.xy(ey[, i], ey[, j], k, h[i], h[j], n)
          }

          pdens <- outer(py[[i]], py[[j]], ".dfmix.xy", w, Theta[i,], Theta[j,])
 
          zlim <- range(edens$z, finite = TRUE); zmax <- max(zlim[2], pdens)

          zlim <- zlim / zmax

          plot(x = edens$x,
            y = edens$y,
            type = "p",
            main = "",
            sub = "",
            xlab = "",
            ylab = "",
            col = rgb(ramp(edens$z / zmax), maxColorValue = 255),
            axes = FALSE,
            lwd = 1,
            cex = plot.cex,
            pch = plot.pch)

          if ((Variables[i] == .rebmix$Variables[2]) && (Variables[j] == .rebmix$Variables[2])) {
            points(x = rep(py[[i]], length(py[[j]])),
              y = rep(py[[j]], each = length(py[[i]])),
              type = "p",
              xlab = "",
              ylab = "",
              col = rgb(ramp(pdens / zmax), maxColorValue = 255),
              lwd = 1,
              cex = plot.cex * 0.5,
              pch = plot.pch)
          }
          else
          if ((Variables[i] == .rebmix$Variables[2]) && (Variables[j] == .rebmix$Variables[1])) {
            for (l in 1:length(py[[i]])) {
              tx <- rep(py[[i]][l], length(py[[j]]))
              ty <- py[[j]]

              s <- 1:(length(tx)-1)

              segments(x0 = tx[s],
                y0 = ty[s],
                x1 = tx[s + 1],
                y1 = ty[s + 1],
                xlab = "",
                ylab = "",
                col = rgb(ramp((pdens[l, s] + pdens[l, s + 1]) / zmax / 2.0), maxColorValue = 255),
                cex = plot.cex)
            }
          }
          else
          if ((Variables[i] == .rebmix$Variables[1]) && (Variables[j] == .rebmix$Variables[2])) {
            for (l in 1:length(py[[j]])) {
              tx <- py[[i]]
              ty <- rep(py[[j]][l], length(py[[i]]))

              s <- 1:(length(tx)-1)

              segments(x0 = tx[s],
                y0 = ty[s],
                x1 = tx[s + 1],
                y1 = ty[s + 1],
                xlab = "",
                ylab = "",
                col = rgb(ramp((pdens[s, l] + pdens[s + 1, l]) / zmax / 2.0), maxColorValue = 255),
                cex = plot.cex)
            }
          }
          else {
            levels <- 10^seq(from = log(zlim[1]), to = log(zlim[2]), length.out = contour.nlevels)

            contour(x = py[[i]],
              y = py[[j]],
              z = pdens / zmax,
              levels = levels,
              xlim = lim[, i],
              ylim = lim[, j],
              zlim = zlim,
              labcex = contour.labcex, drawlabels = contour.drawlabels, method = contour.method,
              axes = FALSE, frame.plot = FALSE,
              col = rgb(ramp(levels), maxColorValue = 255),
              add = TRUE)
          }

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
            for (l in 1:length(legend)) {
              mtext(text = legend[[l]],
                side = 1,
                line = l - 1,
                outer = TRUE,
                adj = 0.5,
                padj = 0.2,
                cex = cex)
            }

            figno <- 0
          }
        }
      }
    }
    
    if (any(match(.rebmix.plot$what[6], what, nomatch = 0))) {  
      figno <- 0
    
      ramp <- colorRamp(colors = c("magenta", "blue", "cyan", "green", "yellow", "red"),
        space = "rgb",
        interpolate = "linear")

      for (i in 1:(d - 1)) {
        for (j in (i + 1):d) {
          edist <- .dist.xy(ey[, i], ey[, j], n)

          pdist <- outer(py[[i]], py[[j]], ".pfmix.xy", w, Theta[i,], Theta[j,])
 
          zlim <- range(edist$z, finite = TRUE); zmax <- max(zlim[2], pdist)

          zlim <- zlim / zmax

          plot(x = edist$x,
            y = edist$y,
            type = "p",
            main = "",
            sub = "",
            xlab = "",
            ylab = "",
            col = rgb(ramp(edist$z / zmax), maxColorValue = 255),
            axes = FALSE,
            lwd = 1,
            cex = plot.cex,
            pch = plot.pch)

          if ((Variables[i] == .rebmix$Variables[2]) && (Variables[j] == .rebmix$Variables[2])) {
            points(x = rep(py[[i]], length(py[[j]])),
              y = rep(py[[j]], each = length(py[[i]])),
              type = "p",
              xlab = "",
              ylab = "",
              col = rgb(ramp(pdist / zmax), maxColorValue = 255),
              lwd = 1,
              cex = plot.cex * 0.5,
              pch = plot.pch)
          }
          else
          if ((Variables[i] == .rebmix$Variables[2]) && (Variables[j] == .rebmix$Variables[1])) {
            for (l in 1:length(py[[i]])) {
              tx <- rep(py[[i]][l], length(py[[j]]))
              ty <- py[[j]]

              s <- 1:(length(tx)-1)

              segments(x0 = tx[s],
                y0 = ty[s],
                x1 = tx[s + 1],
                y1 = ty[s + 1],
                xlab = "",
                ylab = "",
                col = rgb(ramp((pdist[l, s] + pdist[l, s + 1]) / zmax / 2.0), maxColorValue = 255),
                cex = plot.cex)
            }
          }
          else
          if ((Variables[i] == .rebmix$Variables[1]) && (Variables[j] == .rebmix$Variables[2])) {
            for (l in 1:length(py[[j]])) {
              tx <- py[[i]]
              ty <- rep(py[[j]][l], length(py[[i]]))

              s <- 1:(length(tx)-1)

              segments(x0 = tx[s],
                y0 = ty[s],
                x1 = tx[s + 1],
                y1 = ty[s + 1],
                xlab = "",
                ylab = "",
                col = rgb(ramp((pdist[s, l] + pdist[s + 1, l]) / zmax / 2.0), maxColorValue = 255),
                cex = plot.cex)
            }
          }
          else {
            levels <- 10^seq(from = log(zlim[1]), to = log(zlim[2]), length.out = contour.nlevels)

            contour(x = py[[i]],
              y = py[[j]],
              z = pdist / zmax,
              levels = levels,
              xlim = lim[, i],
              ylim = lim[, j],
              zlim = zlim,
              labcex = contour.labcex, drawlabels = contour.drawlabels, method = contour.method,
              axes = FALSE, frame.plot = FALSE,
              col = rgb(ramp(levels), maxColorValue = 255),
              add = TRUE)
          }

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
            for (l in 1:length(legend)) {
              mtext(text = legend[[l]],
                side = 1,
                line = l - 1,
                outer = TRUE,
                adj = 0.5,
                padj = 0.2,
                cex = cex)
            }

            figno <- 0
          }
        }
      }
    }
  }
  else {
    if (any(match(.rebmix.plot$what[1], what, nomatch = 0))) {
      if (C == .rebmix$Preprocessing[1]) {
        edens <- .densHistogram.x(k, ey[, 1], y0[1], h[1], Variables[1], pdf[1])
      }
      else
      if (C == .rebmix$Preprocessing[2]) {
        edens <- .densParzenWindow.x(ey[, 1], h[1], n)
      }
      else
      if (C == .rebmix$Preprocessing[3]) {
        edens <- .densKNearestNeighbour.x(ey[, 1], k, h[1], n)
      }

      pdens <- .dfmix.x(py[[1]], w, Theta[1,])

      ylim <- c(0.0, max(edens$y, pdens))

      plot(x = edens$x,
        y = edens$y,
        type = "p",
        main = "",
        sub = "",
        xlab = "",
        ylab = "",
        ylim = ylim,
        col = "black",
        axes = FALSE,
        lwd = 1,
        cex = plot.cex,
        pch = plot.pch)

      points(x = py[[1]],
        y = pdens,
        type = "l",
        col = "black")

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
        text <- paste("$y_{1}$", "$\\; - \\;$", "$f(y_{1})$", sep = "")
      }
      else {
        text <- bquote(y[1] - f(y[1]))
      }

      mtext(text = text,
        side = 1,
        line = 0,
        outer = FALSE,
        adj = 0.5,
        padj = 0.2,
        cex = cex)
        
      for (l in 1:length(legend)) {
        mtext(text = legend[[l]],
          side = 1,
          line = l - 1,
          outer = TRUE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)
      }        
    }
    
    if (any(match(.rebmix.plot$what[6], what, nomatch = 0))) {
      edist <- .dist.x(ey[, 1], n)

      pdist <- .pfmix.x(py[[1]], w, Theta[1,])

      ylim <- c(0.0, max(edist$y, pdist))

      plot(x = edist$x,
        y = edist$y,
        type = "p",
        main = "",
        sub = "",
        xlab = "",
        ylab = "",
        ylim = ylim,
        col = "black",
        axes = FALSE,
        lwd = 1,
        cex = plot.cex,
        pch = plot.pch)

      points(x = py[[1]],
        y = pdist,
        type = "l",
        col = "black")

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
        text <- paste("$y_{1}$", "$\\; - \\;$", "$F(y_{1})$", sep = "")
      }
      else {
        text <- bquote(y[1] - F(y[1]))
      }

      mtext(text = text,
        side = 1,
        line = 0,
        outer = FALSE,
        adj = 0.5,
        padj = 0.2,
        cex = cex)
        
      for (l in 1:length(legend)) {
        mtext(text = legend[[l]],
          side = 1,
          line = l - 1,
          outer = TRUE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)
      }        
    }
  }
  
  m <- nrow * ncol * ceiling(N / nrow / ncol) - N  
  
  if (any(match(.rebmix.plot$what[2], what, nomatch = 0))) {
    for (i in 1:d) {
      if (C == .rebmix$Preprocessing[1]) {
        edens <- .densHistogram.x(k, ey[, i], y0[i], h[i], Variables[i], pdf[i])
      }
      else
      if (C == .rebmix$Preprocessing[2]) {
        edens <- .densParzenWindow.x(ey[, i], h[i], n)
      }
      else
      if (C == .rebmix$Preprocessing[3]) {
        edens <- .densKNearestNeighbour.x(ey[, i], k, h[i], n)
      }

      pdens <- .dfmix.x(py[[i]], w, Theta[i,])

      ylim <- c(0.0, max(edens$y, pdens))

      plot(x = edens$x,
        y = edens$y,
        type = "p",
        main = "",
        sub = "",
        xlab = "",
        ylab = "",
        ylim = ylim,
        col = "black",
        axes = FALSE,
        lwd = 1,
        cex = plot.cex,
        pch = plot.pch)

      points(x = py[[i]],
        y = pdens,
        type = "l",
        col = "black")

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
        text <- paste("$y_{", i, "}$", "$\\; - \\;$", "$f(y_{", i, "})$", sep = "")
      }
      else {
        text <- bquote(y[.(i)] - f(y[.(i)]))
      }

      mtext(text = text,
        side = 1,
        line = 0,
        outer = FALSE,
        adj = 0.5,
        padj = 0.2,
        cex = cex)
      
      if (m <= 0) {
        for (l in 1:length(legend)) {
          mtext(text = legend[[l]],
            side = 1,
            line = l - 1,
            outer = TRUE,
            adj = 0.5,
            padj = 0.2,
            cex = cex)
        } 
      
        m <- nrow * ncol - 1       
      }
      else {
        m <- m - 1
      }
    }
  }
  
  if (any(match(.rebmix.plot$what[3], what, nomatch = 0))) {
    ylim <- range(x@opt.IC[[pos]], finite = TRUE)
  
    plot(x = x@opt.c[[pos]],
      y = x@opt.IC[[pos]],
      type = "o",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      ylim = ylim,
      col = "black",
      axes = FALSE,
      lwd = 1,
      cex = plot.cex,
      pch = plot.pch)

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
      text <- paste("$c$", "$\\; - \\;$", item[[28]], sep = "")
    }
    else {
      text <- bquote(c - .(item[[28]]))
    }

    mtext(text = text,
      side = 1,
      line = 0,
      outer = FALSE,
      adj = 0.5,
      padj = 0.2,
      cex = cex)
      
    if (m <= 0) {
      for (l in 1:length(legend)) {
        mtext(text = legend[[l]],
          side = 1,
          line = l - 1,
          outer = TRUE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)
      } 
      
      m <- nrow * ncol - 1       
    }
    else {
      m <- m - 1
    }      
  }
  
  if (any(match(.rebmix.plot$what[4], what, nomatch = 0))) {
    ylim <- range(x@opt.logL[[pos]], finite = TRUE)
  
    plot(x = x@opt.c[[pos]],
      y = x@opt.logL[[pos]],
      type = "o",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      ylim = ylim,
      col = "black",
      axes = FALSE,
      lwd = 1,
      cex = plot.cex,
      pch = plot.pch)

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
      text <- paste("$c$", "$\\; - \\;$", item[[31]], sep = "")
    }
    else {
      text <- bquote(c - .(item[[31]]))
    }

    mtext(text = text,
      side = 1,
      line = 0,
      outer = FALSE,
      adj = 0.5,
      padj = 0.2,
      cex = cex)
      
    if (m <= 0) {
      for (l in 1:length(legend)) {
        mtext(text = legend[[l]],
          side = 1,
          line = l - 1,
          outer = TRUE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)
      } 
      
      m <- nrow * ncol - 1       
    }
    else {
      m <- m - 1
    }      
  }  
  
  if (any(match(.rebmix.plot$what[5], what, nomatch = 0))) {
    ylim <- range(x@opt.D[[pos]], finite = TRUE)
  
    plot(x = x@opt.c[[pos]],
      y = x@opt.D[[pos]],
      type = "o",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      ylim = ylim,
      col = "black",
      axes = FALSE,
      lwd = 1,
      cex = plot.cex,
      pch = plot.pch)

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
      text <- paste("$c$", "$\\; - \\;$", item[[10]], sep = "")
    }
    else {
      text <- bquote(c - .(item[[10]]))
    }

    mtext(text = text,
      side = 1,
      line = 0,
      outer = FALSE,
      adj = 0.5,
      padj = 0.2,
      cex = cex)
      
    if (m <= 0) {
      for (l in 1:length(legend)) {
        mtext(text = legend[[l]],
          side = 1,
          line = l - 1,
          outer = TRUE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)
      } 
      
      m <- nrow * ncol - 1       
    }
    else {
      m <- m - 1
    }      
  }
  
  if (any(match(.rebmix.plot$what[7], what, nomatch = 0))) {
    ylim <- range(x@all.IC[[pos]], finite = TRUE)
  
    plot(x = x@all.K[[pos]],
      y = x@all.IC[[pos]],
      type = "o",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      ylim = ylim,
      col = "black",
      axes = FALSE,
      lwd = 1,
      cex = plot.cex,
      pch = plot.pch)

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
      text <- paste("$K$", "$\\; - \\;$", item[[28]], sep = "")
    }
    else {
      text <- bquote(K - .(item[[28]]))
    }

    mtext(text = text,
      side = 1,
      line = 0,
      outer = FALSE,
      adj = 0.5,
      padj = 0.2,
      cex = cex)
      
    if (m <= 0) {
      for (l in 1:length(legend)) {
        mtext(text = legend[[l]],
          side = 1,
          line = l - 1,
          outer = TRUE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)
      } 
      
      m <- nrow * ncol - 1       
    }
    else {
      m <- m - 1
    }      
  }
  
  rm(list = ls()[!(ls() %in% c("opar"))])

  invisible(opar)
}) ## plot

setMethod("plot", 
          signature(x = "REBMVNORM", y = "missing"),
function(x,
  y,
  pos = 1,
  what = c("density"),
  nrow = 1,
  ncol = 1,
  npts = 200,
  n = 200,
  cex = 0.8,
  fg = "black",
  lty = "solid",
  lwd = 1,
  pty = "m",
  tcl = 0.5,
  plot.cex = 0.8,
  plot.pch = 19,
  contour.drawlabels = FALSE,
  contour.labcex = 0.8,
  contour.method = "flattest",
  contour.nlevels = 12, ...)
{
  if (missing(x)) {
    stop(sQuote("x"), " object of class REBMIX is requested!", call. = FALSE)
  }

  if (!is.wholenumber(pos)) {
    stop(sQuote("pos"), " integer is requested!", call. = FALSE)
  }
  
  length(pos) <- 1

  if ((pos < 1) || (pos > nrow(x@summary))) {
    stop(sQuote("pos"), " must be greater than 0 and less or equal than ", nrow(x@summary), "!", call. = FALSE)
  }
  
  if (!is.character(what)) {
    stop(sQuote("what"), " character vector is requested!", call. = FALSE)
  }

  what <- match.arg(what, .rebmix.plot$what, several.ok = TRUE)  

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

  if (!is.wholenumber(npts)) {
    stop(sQuote("npts"), " integer is requested!", call. = FALSE)
  }

  if (npts < 1) {
    stop(sQuote("npts"), " must be greater than 0!", call. = FALSE)
  }
  
  if (!is.wholenumber(n)) {
    stop(sQuote("n"), " integer is requested!", call. = FALSE)
  }

  if (n < 1) {
    stop(sQuote("n"), " must be greater than 0!", call. = FALSE)
  }  

  ni <- ncol(x@summary)

  Theta <- .extractThetaB(x@w[[pos]], x@Theta[[pos]])

  d <- length(Theta[[1]]$pdf)

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

  C <- x@summary[pos, "Preprocessing"]

  if (.Device == "tikz output") {
    item <- list()

    item[[1]] <- "$\\textrm{Dataset}$"
    item[[2]] <- "$\\; = \\;$"
    item[[3]] <- paste("$\\textrm{", x@summary[pos, "Dataset"], "}$, ", sep = "")

    item[[4]] <- "$\\textrm{Preprocessing}$"
    item[[5]] <- "$\\; = \\;$"
    item[[6]] <- paste("$\\textrm{", x@summary[pos, "Preprocessing"], "}$, ", sep = "")

    item[[7]] <- "$\\textrm{Restraints}$"
    item[[8]] <- "$\\; = \\;$"
    item[[9]] <- paste("$\\textrm{", x@summary[pos, "Restraints"], "}$, ", sep = "")

    item[[10]] <- "$D$"
    item[[11]] <- ""
    item[[12]] <- ""

    item[[13]] <- "$c_{\\mathrm{max}}$"
    item[[14]] <- "$\\; = \\;$"
    item[[15]] <- paste("$", x@summary[pos, "cmax"], "$, ", sep = "")

    item[[16]] <- "$a_{\\mathrm{r}}$"
    item[[17]] <- "$\\; = \\;$"
    item[[18]] <- paste("$", as.number(x@summary[pos, "ar"]), "$, ", sep = "")

    item[[19]] <- "$c$"
    item[[20]] <- "$\\; = \\;$"
    item[[21]] <- paste("$", x@summary[pos, "c"], "$, ", sep = "")

    item[[22]] <- ""
    item[[23]] <- ""
    item[[24]] <- ""

    if (C == .rebmix$Preprocessing[1]) {
      item[[25]] <- "$v$"
      item[[26]] <- "$\\; = \\;$"
      item[[27]] <- paste("$", x@summary[pos, "v/k"], "$, ", sep = "")
    }
    else
    if (C == .rebmix$Preprocessing[2]) {
      item[[25]] <- "$v$"
      item[[26]] <- "$\\; = \\;$"
      item[[27]] <- paste("$", x@summary[pos, "v/k"], "$, ", sep = "")
    }
    else
    if (C == .rebmix$Preprocessing[3]) {
      item[[25]] <- "$k$"
      item[[26]] <- "$\\; = \\;$"
      item[[27]] <- paste("$", x@summary[pos, "v/k"], "$, ", sep = "")
    }

    item[[28]] <- paste("$\\mathrm{", x@summary[pos, "Criterion"], "}$", sep = "")
    item[[29]] <- "$\\; = \\;$"
    item[[30]] <- paste("$", as.number(x@summary[pos, "IC"]), "$, ", sep = "")

    item[[31]] <- "$\\mathrm{log}\\, L$"
    item[[32]] <- "$\\; = \\;$"
    item[[33]] <- paste("$", as.number(x@summary[pos, "logL"]), "$.", sep = "")

    i <- 1; legend <- list(); legend[[i]] <- item[[1]]

    for (j in c(2:9, 13:21, 25:33)) {
      legendwidth <- strwidth(paste(legend[[i]], item[[j]], sep = ""), units = "figure", cex = 1.0)

      if (legendwidth > ncol) {
        i <- i + 1; legend[[i]] <- item[[j]]
      }
      else {
        legend[[i]] <- paste(legend[[i]], item[[j]], sep = "")
      }
    }
  }
  else {
    item <- list()

    item[[1]] <- "Dataset"
    item[[2]] <- " = "
    item[[3]] <- paste(x@summary[pos, "Dataset"], ", ", sep = "")

    item[[4]] <- "Preprocessing"
    item[[5]] <- " = "
    item[[6]] <- paste(x@summary[pos, "Preprocessing"], ", ", sep = "")

    item[[7]] <- "Restraints"
    item[[8]] <- " = "
    item[[9]] <- paste(x@summary[pos, "Restraints"], ", ", sep = "")

    item[[10]] <- "D"
    item[[11]] <- ""
    item[[12]] <- ""

    item[[13]] <- bquote(c[max])
    item[[14]] <- " = "
    item[[15]] <- paste(x@summary[pos, "cmax"], ", ", sep = "")

    item[[16]] <- bquote(a[r])
    item[[17]] <- " = "
    item[[18]] <- paste(as.number(x@summary[pos, "ar"]), ", ", sep = "")

    item[[19]] <- "c"
    item[[20]] <- " = "
    item[[21]] <- paste(x@summary[pos, "c"], ", ", sep = "")

    item[[22]] <- ""
    item[[23]] <- ""
    item[[24]] <- ""

    if (C == .rebmix$Preprocessing[1]) {
      item[[25]] <- "v"
      item[[26]] <- " = "
      item[[27]] <- paste(x@summary[pos, "v/k"], ", ", sep = "")
    }
    else
    if (C == .rebmix$Preprocessing[2]) {
      item[[25]] <- "v"
      item[[26]] <- " = "
      item[[27]] <- paste(x@summary[pos, "v/k"], ", ", sep = "")
    }
    else
    if (C == .rebmix$Preprocessing[3]) {
      item[[25]] <- "k"
      item[[26]] <- " = "
      item[[27]] <- paste(x@summary[pos, "v/k"], ", ", sep = "")
    }

    item[[28]] <- as.character(x@summary[pos, "Criterion"])
    item[[29]] <- " = "
    item[[30]] <- paste(as.number(x@summary[pos, "IC"]), ", ", sep = "")

    item[[31]] <- "log L"
    item[[32]] <- " = "
    item[[33]] <- paste(as.number(x@summary[pos, "logL"]), ".", sep = "")

    i <- 1; legend <- list(); legend[[i]] <- bquote(.(item[[1]]))

    for (j in c(2:9, 13:21, 25:33)) {
      legendwidth <- strwidth(bquote(paste(.(legend[[i]]), .(item[[j]]), sep = "")), units = "figure", cex = 1.0)

      if (legendwidth > ncol) {
        i <- i + 1; legend[[i]] <- item[[j]]
      }
      else {
        legend[[i]] <- bquote(paste(.(legend[[i]]), .(item[[j]]), sep = ""))
      }
    }
  }

  par(oma = c(length(legend) + 0.2, 0.2, 0.2, 0.2))

  Dataset <- as.character(x@summary[pos, "Dataset"])

  ey <- as.matrix(x@Dataset[[which(names(x@Dataset) == x@summary[pos, "Dataset"])]])

  y0 <- array(data = 0.0, dim = d, dimnames = NULL)
  h <- array(data = 0.0, dim = d, dimnames = NULL)

  lim <- array(data = 0.0, dim = c(2, d), dimnames = NULL)

  py <- list(d)

  Variables <- match.arg(x@Variables, .rebmix$Variables, several.ok = TRUE)
  
  pdf <- match.arg(x@pdf, .rebmix$pdf, several.ok = TRUE)
  
  for (i in 1:d) {
    if (C == .rebmix$Preprocessing[1]) {
      k <- as.numeric(x@summary[pos, "v/k"])
      y0[i] <- as.numeric(x@summary[pos, paste("y0", if (d > 1) i, sep = "")])
      h[i] <- as.numeric(x@summary[pos, paste("h", if (d > 1) i, sep = "")])

      lim[, i] <- range(ey[, i], finite = TRUE)
    }
    else
    if (C == .rebmix$Preprocessing[2]) {
      h[i] <- as.numeric(x@summary[pos, paste("h", if (d > 1) i, sep = "")])

      lim[, i] <- range(ey[, i], finite = TRUE)
    }
    else
    if (C == .rebmix$Preprocessing[3]) {
      k <- as.numeric(x@summary[pos, "v/k"])

      h[i] <- as.numeric(x@summary[pos, paste("h", if (d > 1) i, sep = "")])

      lim[, i] <- range(ey[, i], finite = TRUE)
    }
    else {
      lim[, i] <- c(0.0, 1.0)
    }

    if (Variables[i] == .rebmix$Variables[2]) {
      py[[i]] <- seq(from = lim[1, i], to = lim[2, i], by = 1.0)
    }
    else {
      py[[i]] <- seq(from = lim[1, i], to = lim[2, i], length.out = npts)
    }
  }

  w <- as.numeric(x@w[[pos]])

  if (N > 0) {
    if (any(match(.rebmix.plot$what[1], what, nomatch = 0))) {  
      figno <- 0
    
      ramp <- colorRamp(colors = c("magenta", "blue", "cyan", "green", "yellow", "red"),
        space = "rgb",
        interpolate = "linear")

      for (i in 1:(d - 1)) {
        for (j in (i + 1):d) {
          if (C == .rebmix$Preprocessing[1]) {
            edens <- .densHistogram.xy(k, ey[, i], ey[, j], y0[i], y0[j], h[i], h[j], Variables[i], Variables[j], pdf[i], pdf[j])
          }
          else
          if (C == .rebmix$Preprocessing[2]) {
            edens <- .densParzenWindow.xy(ey[, i], ey[, j], h[i], h[j], n)
          }
          else
          if (C == .rebmix$Preprocessing[3]) {
            edens <- .densKNearestNeighbour.xy(ey[, i], ey[, j], k, h[i], h[j], n)
          }

          pdens <- outer(py[[i]], py[[j]], ".dfmvnorm.xy", w, Theta, i, j)
          
          zlim <- range(edens$z, finite = TRUE); zmax <- max(zlim[2], pdens)

          zlim <- zlim / zmax

          plot(x = edens$x,
            y = edens$y,
            type = "p",
            main = "",
            sub = "",
            xlab = "",
            ylab = "",
            col = rgb(ramp(edens$z / zmax), maxColorValue = 255),
            axes = FALSE,
            lwd = 1,
            cex = plot.cex,
            pch = plot.pch)

          if ((Variables[i] == .rebmix$Variables[2]) && (Variables[j] == .rebmix$Variables[2])) {
            points(x = rep(py[[i]], length(py[[j]])),
              y = rep(py[[j]], each = length(py[[i]])),
              type = "p",
              xlab = "",
              ylab = "",
              col = rgb(ramp(pdens / zmax), maxColorValue = 255),
              lwd = 1,
              cex = plot.cex * 0.5,
              pch = plot.pch)
          }
          else
          if ((Variables[i] == .rebmix$Variables[2]) && (Variables[j] == .rebmix$Variables[1])) {
            for (l in 1:length(py[[i]])) {
              tx <- rep(py[[i]][l], length(py[[j]]))
              ty <- py[[j]]

              s <- 1:(length(tx)-1)

              segments(x0 = tx[s],
                y0 = ty[s],
                x1 = tx[s + 1],
                y1 = ty[s + 1],
                xlab = "",
                ylab = "",
                col = rgb(ramp((pdens[l, s] + pdens[l, s + 1]) / zmax / 2.0), maxColorValue = 255),
                cex = plot.cex)
            }
          }
          else
          if ((Variables[i] == .rebmix$Variables[1]) && (Variables[j] == .rebmix$Variables[2])) {
            for (l in 1:length(py[[j]])) {
              tx <- py[[i]]
              ty <- rep(py[[j]][l], length(py[[i]]))

              s <- 1:(length(tx)-1)

              segments(x0 = tx[s],
                y0 = ty[s],
                x1 = tx[s + 1],
                y1 = ty[s + 1],
                xlab = "",
                ylab = "",
                col = rgb(ramp((pdens[s, l] + pdens[s + 1, l]) / zmax / 2.0), maxColorValue = 255),
                cex = plot.cex)
            }
          }
          else {
            levels <- 10^seq(from = log(zlim[1]), to = log(zlim[2]), length.out = contour.nlevels)

            contour(x = py[[i]],
              y = py[[j]],
              z = pdens / zmax,
              levels = levels,
              xlim = lim[, i],
              ylim = lim[, j],
              zlim = zlim,
              labcex = contour.labcex, drawlabels = contour.drawlabels, method = contour.method,
              axes = FALSE, frame.plot = FALSE,
              col = rgb(ramp(levels), maxColorValue = 255),
              add = TRUE)
          }

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
            for (l in 1:length(legend)) {
              mtext(text = legend[[l]],
                side = 1,
                line = l - 1,
                outer = TRUE,
                adj = 0.5,
                padj = 0.2,
                cex = cex)
            }

            figno <- 0
          }
        }
      }
    }
    
    if (any(match(.rebmix.plot$what[6], what, nomatch = 0))) {  
      figno <- 0
    
      ramp <- colorRamp(colors = c("magenta", "blue", "cyan", "green", "yellow", "red"),
        space = "rgb",
        interpolate = "linear")

      for (i in 1:(d - 1)) {
        for (j in (i + 1):d) {
          edist <- .dist.xy(ey[, i], ey[, j], n)

          pdist <- outer(py[[i]], py[[j]], ".pfmvnorm.xy", w, Theta, i, j)          
 
          zlim <- range(edist$z, finite = TRUE); zmax <- max(zlim[2], pdist)

          zlim <- zlim / zmax

          plot(x = edist$x,
            y = edist$y,
            type = "p",
            main = "",
            sub = "",
            xlab = "",
            ylab = "",
            col = rgb(ramp(edist$z / zmax), maxColorValue = 255),
            axes = FALSE,
            lwd = 1,
            cex = plot.cex,
            pch = plot.pch)

          if ((Variables[i] == .rebmix$Variables[2]) && (Variables[j] == .rebmix$Variables[2])) {
            points(x = rep(py[[i]], length(py[[j]])),
              y = rep(py[[j]], each = length(py[[i]])),
              type = "p",
              xlab = "",
              ylab = "",
              col = rgb(ramp(pdist / zmax), maxColorValue = 255),
              lwd = 1,
              cex = plot.cex * 0.5,
              pch = plot.pch)
          }
          else
          if ((Variables[i] == .rebmix$Variables[2]) && (Variables[j] == .rebmix$Variables[1])) {
            for (l in 1:length(py[[i]])) {
              tx <- rep(py[[i]][l], length(py[[j]]))
              ty <- py[[j]]

              s <- 1:(length(tx)-1)

              segments(x0 = tx[s],
                y0 = ty[s],
                x1 = tx[s + 1],
                y1 = ty[s + 1],
                xlab = "",
                ylab = "",
                col = rgb(ramp((pdist[l, s] + pdist[l, s + 1]) / zmax / 2.0), maxColorValue = 255),
                cex = plot.cex)
            }
          }
          else
          if ((Variables[i] == .rebmix$Variables[1]) && (Variables[j] == .rebmix$Variables[2])) {
            for (l in 1:length(py[[j]])) {
              tx <- py[[i]]
              ty <- rep(py[[j]][l], length(py[[i]]))

              s <- 1:(length(tx)-1)

              segments(x0 = tx[s],
                y0 = ty[s],
                x1 = tx[s + 1],
                y1 = ty[s + 1],
                xlab = "",
                ylab = "",
                col = rgb(ramp((pdist[s, l] + pdist[s + 1, l]) / zmax / 2.0), maxColorValue = 255),
                cex = plot.cex)
            }
          }
          else {
            levels <- 10^seq(from = log(zlim[1]), to = log(zlim[2]), length.out = contour.nlevels)

            contour(x = py[[i]],
              y = py[[j]],
              z = pdist / zmax,
              levels = levels,
              xlim = lim[, i],
              ylim = lim[, j],
              zlim = zlim,
              labcex = contour.labcex, drawlabels = contour.drawlabels, method = contour.method,
              axes = FALSE, frame.plot = FALSE,
              col = rgb(ramp(levels), maxColorValue = 255),
              add = TRUE)
          }

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
            for (l in 1:length(legend)) {
              mtext(text = legend[[l]],
                side = 1,
                line = l - 1,
                outer = TRUE,
                adj = 0.5,
                padj = 0.2,
                cex = cex)
            }

            figno <- 0
          }
        }
      }
    }
  }
  else {
    if (any(match(.rebmix.plot$what[1], what, nomatch = 0))) {
      if (C == .rebmix$Preprocessing[1]) {
        edens <- .densHistogram.x(k, ey[, 1], y0[1], h[1], Variables[1], pdf[1])
      }
      else
      if (C == .rebmix$Preprocessing[2]) {
        edens <- .densParzenWindow.x(ey[, 1], h[1], n)
      }
      else
      if (C == .rebmix$Preprocessing[3]) {
        edens <- .densKNearestNeighbour.x(ey[, 1], k, h[1], n)
      }

      pdens <- .dfmvnorm.x(py[[1]], w, Theta, 1)      

      ylim <- c(0.0, max(edens$y, pdens))

      plot(x = edens$x,
        y = edens$y,
        type = "p",
        main = "",
        sub = "",
        xlab = "",
        ylab = "",
        ylim = ylim,
        col = "black",
        axes = FALSE,
        lwd = 1,
        cex = plot.cex,
        pch = plot.pch)

      points(x = py[[1]],
        y = pdens,
        type = "l",
        col = "black")

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
        text <- paste("$y_{1}$", "$\\; - \\;$", "$f(y_{1})$", sep = "")
      }
      else {
        text <- bquote(y[1] - f(y[1]))
      }

      mtext(text = text,
        side = 1,
        line = 0,
        outer = FALSE,
        adj = 0.5,
        padj = 0.2,
        cex = cex)
        
      for (l in 1:length(legend)) {
        mtext(text = legend[[l]],
          side = 1,
          line = l - 1,
          outer = TRUE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)
      }        
    }
    
    if (any(match(.rebmix.plot$what[6], what, nomatch = 0))) {
      edist <- .dist.x(ey[, 1], n)

      pdist <- .pfmvnorm.x(py[[1]], w, Theta, 1)

      ylim <- c(0.0, max(edist$y, pdist))

      plot(x = edist$x,
        y = edist$y,
        type = "p",
        main = "",
        sub = "",
        xlab = "",
        ylab = "",
        ylim = ylim,
        col = "black",
        axes = FALSE,
        lwd = 1,
        cex = plot.cex,
        pch = plot.pch)

      points(x = py[[1]],
        y = pdist,
        type = "l",
        col = "black")

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
        text <- paste("$y_{1}$", "$\\; - \\;$", "$F(y_{1})$", sep = "")
      }
      else {
        text <- bquote(y[1] - F(y[1]))
      }

      mtext(text = text,
        side = 1,
        line = 0,
        outer = FALSE,
        adj = 0.5,
        padj = 0.2,
        cex = cex)
        
      for (l in 1:length(legend)) {
        mtext(text = legend[[l]],
          side = 1,
          line = l - 1,
          outer = TRUE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)
      }        
    }
  }
  
  m <- nrow * ncol * ceiling(N / nrow / ncol) - N  
  
  if (any(match(.rebmix.plot$what[2], what, nomatch = 0))) {
    for (i in 1:d) {
      if (C == .rebmix$Preprocessing[1]) {
        edens <- .densHistogram.x(k, ey[, i], y0[i], h[i], Variables[i], pdf[i])
      }
      else
      if (C == .rebmix$Preprocessing[2]) {
        edens <- .densParzenWindow.x(ey[, i], h[i], n)
      }
      else
      if (C == .rebmix$Preprocessing[3]) {
        edens <- .densKNearestNeighbour.x(ey[, i], k, h[i], n)
      }

      pdens <- .dfmvnorm.x(py[[i]], w, Theta, i)

      ylim <- c(0.0, max(edens$y, pdens))

      plot(x = edens$x,
        y = edens$y,
        type = "p",
        main = "",
        sub = "",
        xlab = "",
        ylab = "",
        ylim = ylim,
        col = "black",
        axes = FALSE,
        lwd = 1,
        cex = plot.cex,
        pch = plot.pch)

      points(x = py[[i]],
        y = pdens,
        type = "l",
        col = "black")

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
        text <- paste("$y_{", i, "}$", "$\\; - \\;$", "$f(y_{", i, "})$", sep = "")
      }
      else {
        text <- bquote(y[.(i)] - f(y[.(i)]))
      }

      mtext(text = text,
        side = 1,
        line = 0,
        outer = FALSE,
        adj = 0.5,
        padj = 0.2,
        cex = cex)
      
      if (m <= 0) {
        for (l in 1:length(legend)) {
          mtext(text = legend[[l]],
            side = 1,
            line = l - 1,
            outer = TRUE,
            adj = 0.5,
            padj = 0.2,
            cex = cex)
        } 
      
        m <- nrow * ncol - 1       
      }
      else {
        m <- m - 1
      }
    }
  }
  
  if (any(match(.rebmix.plot$what[3], what, nomatch = 0))) {
    ylim <- range(x@opt.IC[[pos]], finite = TRUE)
  
    plot(x = x@opt.c[[pos]],
      y = x@opt.IC[[pos]],
      type = "o",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      ylim = ylim,
      col = "black",
      axes = FALSE,
      lwd = 1,
      cex = plot.cex,
      pch = plot.pch)

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
      text <- paste("$c$", "$\\; - \\;$", item[[28]], sep = "")
    }
    else {
      text <- bquote(c - .(item[[28]]))
    }

    mtext(text = text,
      side = 1,
      line = 0,
      outer = FALSE,
      adj = 0.5,
      padj = 0.2,
      cex = cex)
      
    if (m <= 0) {
      for (l in 1:length(legend)) {
        mtext(text = legend[[l]],
          side = 1,
          line = l - 1,
          outer = TRUE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)
      } 
      
      m <- nrow * ncol - 1       
    }
    else {
      m <- m - 1
    }      
  }
  
  if (any(match(.rebmix.plot$what[4], what, nomatch = 0))) {
    ylim <- range(x@opt.logL[[pos]], finite = TRUE)
  
    plot(x = x@opt.c[[pos]],
      y = x@opt.logL[[pos]],
      type = "o",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      ylim = ylim,
      col = "black",
      axes = FALSE,
      lwd = 1,
      cex = plot.cex,
      pch = plot.pch)

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
      text <- paste("$c$", "$\\; - \\;$", item[[31]], sep = "")
    }
    else {
      text <- bquote(c - .(item[[31]]))
    }

    mtext(text = text,
      side = 1,
      line = 0,
      outer = FALSE,
      adj = 0.5,
      padj = 0.2,
      cex = cex)
      
    if (m <= 0) {
      for (l in 1:length(legend)) {
        mtext(text = legend[[l]],
          side = 1,
          line = l - 1,
          outer = TRUE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)
      } 
      
      m <- nrow * ncol - 1       
    }
    else {
      m <- m - 1
    }      
  }  
  
  if (any(match(.rebmix.plot$what[5], what, nomatch = 0))) {
    ylim <- range(x@opt.D[[pos]], finite = TRUE)
  
    plot(x = x@opt.c[[pos]],
      y = x@opt.D[[pos]],
      type = "o",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      ylim = ylim,
      col = "black",
      axes = FALSE,
      lwd = 1,
      cex = plot.cex,
      pch = plot.pch)

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
      text <- paste("$c$", "$\\; - \\;$", item[[10]], sep = "")
    }
    else {
      text <- bquote(c - .(item[[10]]))
    }

    mtext(text = text,
      side = 1,
      line = 0,
      outer = FALSE,
      adj = 0.5,
      padj = 0.2,
      cex = cex)
      
    if (m <= 0) {
      for (l in 1:length(legend)) {
        mtext(text = legend[[l]],
          side = 1,
          line = l - 1,
          outer = TRUE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)
      } 
      
      m <- nrow * ncol - 1       
    }
    else {
      m <- m - 1
    }      
  }
  
  if (any(match(.rebmix.plot$what[7], what, nomatch = 0))) {
    ylim <- range(x@all.IC[[pos]], finite = TRUE)
  
    plot(x = x@all.K[[pos]],
      y = x@all.IC[[pos]],
      type = "o",
      main = "",
      sub = "",
      xlab = "",
      ylab = "",
      ylim = ylim,
      col = "black",
      axes = FALSE,
      lwd = 1,
      cex = plot.cex,
      pch = plot.pch)

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
      text <- paste("$K$", "$\\; - \\;$", item[[28]], sep = "")
    }
    else {
      text <- bquote(K - .(item[[28]]))
    }

    mtext(text = text,
      side = 1,
      line = 0,
      outer = FALSE,
      adj = 0.5,
      padj = 0.2,
      cex = cex)
      
    if (m <= 0) {
      for (l in 1:length(legend)) {
        mtext(text = legend[[l]],
          side = 1,
          line = l - 1,
          outer = TRUE,
          adj = 0.5,
          padj = 0.2,
          cex = cex)
      } 
      
      m <- nrow * ncol - 1       
    }
    else {
      m <- m - 1
    }      
  }
  
  rm(list = ls()[!(ls() %in% c("opar"))])

  invisible(opar)
}) # plot         
