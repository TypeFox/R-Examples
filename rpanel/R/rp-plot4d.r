#     Plots of two covariates coloured by a response variable and 
#     animated by a third covariate.

rp.plot4d <- function(x, z, y, model,
                  col.palette, col.breaks, col.labels, 
                  hscale = 1, vscale = hscale, panel = TRUE,
                  x1lab, x2lab, zlab, ylab,
                  background.plot = NULL, foreground.plot = NULL, 
                  z.window = "normal", z.window.pars = c(min(z), sd(z)/5),
                  coords = rep(NA, 2), radius = 0.05, location.plot = TRUE,
                  retain.location.plot = FALSE, eqscplot = FALSE,
                  location.plot.type = "histogram") {

   if (eqscplot & !require(MASS)) {
      cat("eqscplot requires the MASS package which is not available.\n")
      eqscplot <- FALSE
   }
   
   draw.plot <- function(panel) {
      with(panel, {
      	
      	 z0  <- z.window.pars["location"]
      	 zsd <- z.window.pars["width"]
      	 
         if (panel.plot) par(mar = c(3, 3, 1, 1) + 0.1, mgp = c(1.5, 0.2, 0), tcl = -0.2)
         if (eqscplot)
            eqscplot(x, type = "n", xlab = x1lab, ylab = x2lab)
         else
            plot(x, type = "n", xlab = x1lab, ylab = x2lab)
         
         if (is.list(model)) {
         	mt     <- dim(model$y)[3]
            m.ind  <- (z0 - min(model$z)) / (diff(range(model$z)) / (mt - 1))
            m.low  <- min(1 + floor(m.ind), mt)
            m.high <- min(1 + m.low, mt)
            m.p    <- 1 + m.ind - m.low
            if (m.low >= 1 & m.high <= mt) {
               my     <- (1 - m.p) * model$y[ , , m.low] + m.p * model$y[ , , m.high]
               my.ind <- !is.na(my)
              if (require(akima)) {
                  m.in  <- as.matrix(expand.grid(model$x[ , 1], model$x[ , 2]))
                  ngrd  <- 150
                  mx    <- cbind(seq(min(model$x[ , 1]), max(model$x[ , 1]), length = ngrd),
                                 seq(min(model$x[ , 2]), max(model$x[ , 2]), length = ngrd))
                  grd   <- interp(m.in[c(my.ind), 1], m.in[c(my.ind), 2], c(my)[c(my.ind)],
                                  mx[ , 1], mx[ , 2])
#                  d1.in <- model$x[2, 1] - model$x[1, 1]
#                  d2.in <- model$x[2, 2] - model$x[1, 2]
#                  d1.x  <- diff(range(model$x[ , 1])) / (ngrd - 1)
#                  d2.x  <- diff(range(model$x[ , 2])) / (ngrd - 1)
#                  nn1   <- min(1 + floor(mx[ , 1] / model$x[ , 1]), nrow(model$x))
#                  nn2   <- min(1 + floor(mx[ , 2] / model$x[ , 2]), ncol(model$x))
                  nn   <- floor((0:(ngrd - 1)) * (nrow(model$x) - 1) / (ngrd - 1))
                  nn   <- pmin(1 + nn, nrow(model$x) - 1)
                  nn   <- as.matrix(expand.grid(nn, nn))
                  ind  <- is.na(my[nn] + my[nn + 1] + my[cbind(nn[,1] + 1, nn[,2])] +
                                                      my[cbind(nn[,1], nn[,2] + 1)])
                  ind   <- matrix(ind, nrow = ngrd)
                  my    <- grd$z
                  my[ind] <- NA
               }
               else
                  mx <- model$x
               brks[is.infinite(brks) & (brks > 0)] <- max(y, model$y, na.rm = TRUE) + 1
               brks[is.infinite(brks) & (brks < 0)] <- min(y, model$y, na.rm = TRUE) - 1
               image(mx[ , 1], mx[ , 2], my, breaks = brks, col = col.palette, add = TRUE)
#               print(levelplot(my, 
#                     row.values = model$x[, 1], column.values = model$x[ , 2],
#                     cuts = length(col.palette), colorkey = FALSE, at = brks,
#                     region = TRUE, col.regions = col.palette))
            }
         }

         if (is.function(background.plot)) background.plot()

         if (is.list(model)) z.window <- "uniform"
         zsd1  <- if (zsd >= 1.49 * sdz) 4 * sdz else zsd
         alpha <- exp(-0.5 * (z - z0)^2 / zsd1^2)
         ord   <- order(alpha)
         if (z.window == "normal")
            clr <- hsv(clr[1, ], clr[2, ] * alpha, clr[3, ])
         else if (z.window == "uniform") {
            clr <- hsv(clr[1, ], clr[2, ], clr[3, ])
            ord <- ord[abs(z[ord] - z0) < 2 * zsd1]
            if (is.list(model)) points(x[ord, 1], x[ord, 2])
         }
         points(x[ord, 1], x[ord, 2], pch = 16, col = clr[ord])

         if (is.function(foreground.plot)) foreground.plot()

         if (all(!is.na(coords))) {
         	dr1          <- diff(range(panel$x[ , 1]))
   	        dr2          <- diff(range(panel$x[ , 2]))
   	        if (eqscplot) {
   	           dr1 <- max(dr1, dr2)
   	           dr2 <- max(dr1, dr2)
   	        }
            lines(coords[1] + circle[ , 1] * radius * dr1, 
                  coords[2] + circle[ , 2] * radius * dr2)
         }
      })
      panel
   }
   
   draw.key <- function(panel) {
   	  if (panel$missing.y) return(panel)
      if (is.factor(panel$y)) {
      	 par(mar = c(3, 0, 1, 0) + 0.1)
      	 plot(0:1, type = "n", axes = FALSE, xlab = "", ylab = "")
      	 for (i in 1:length(levels(panel$y)))
      	    text(1, 1 - i * 1.5 * strheight("A"), levels(panel$y)[i],
      	         col = panel$col.palette[i], pos = 4, offset = 0)
      	 # legend("topleft", levels(panel$y), # col = panel$col.palette,
         #    text.col = panel$col.palette)
      }
      else {
         rp.colour.key(panel$col.palette, panel$col.labels, par.mar = c(3, 1, 1, 1.5) + 0.1,
             natural = panel$natural)
         mtext(ylab, side = 2, line = 0.1, font = 1)
      }
      panel
   }
   
   draw.band <- function(panel) {
      with(panel, {
      	 z0      <- z.window.pars["location"]
      	 zsd     <- z.window.pars["width"]
      	 z0      <- z.window.pars[1]
      	 zsd     <- z.window.pars[2]
         mar.old <- par()$mar
         # par(mar = c(0, 3, 2, 1) + 0.1, mgp = c(1, 0.2, 0), tcl = -0.2)
         # plot(range(z), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "",
         #       xaxs = "i", yaxs = "i")
         par(mar = c(0, 3, 2, 1) + 0.1, mgp = c(1.5, 0.2, 0), tcl = -0.2)
         plot(range(z), c(0, 1), type = "n", axes = FALSE, yaxs = "i", xlab = " ", ylab = " ")
         zsd1  <- if (zsd >= 1.49 * sdz) 4 * sdz else zsd
         if (is.list(model)) z.window <- "uniform"
         if (z.window == "normal") {
            nrect <- 100
            zvec  <- seq(par()$usr[1], par()$usr[2], length = nrect + 1)
            zmid  <- (zvec[-nrect + 1] + zvec[-1]) / 2
            alpha <- exp(-0.5 * (zmid - z0)^2 / zsd1^2)
            clr   <- rgb2hsv(col2rgb("green"))
            clr   <- hsv(rep(clr[1, ], nrect), clr[2, ] * alpha, rep(clr[3, ], nrect))
            rect(zvec[-(nrect + 1)], 0, zvec[-1], 1, col = clr, border = NA)
         }
         else if (z.window == "uniform")
            rect(z0 - 2 * zsd1, 0, z0 + 2 * zsd1, 1, col = "green", border = NA)
         axis(3, font.main = 1)
         mtext(zlab, line = 1, font = 1)
         box()
         par(mar = mar.old)
      })
      panel
   }

   draw.location <- function(panel) {
      with(panel, {
         if (missing.y | is.factor(y)) {
            par(mar = c(3, 3, 1, 1) + 0.1, mgp = c(1.5, 0.2, 0), tcl = -0.2)
            if (missing.y) y <- factor(rep(1, length(z)))
            if (length(z[locind]) > 0) {
               if (nlevels(y) > 1 & require(lattice)) {
                  if (location.plot.type == "histogram")
                     print(histogram( ~ z[locind] | y[locind], xlab = zlab, xlim = range(z),
                                  layout = c(1, nlevels(y)), type = "count"))
                  else
                     print(densityplot( ~ z[locind] | y[locind], xlab = zlab, xlim = range(z),
                                  layout = c(1, nlevels(y)), type = "count"))
               }
               else {
               	  if (location.plot.type == "histogram") {
                     hist(z[locind], main = "", xlab = zlab, xlim = range(z))
                     box()
                  }
                  else
                     print(densityplot( ~ z[locind], xlab = zlab, xlim = range(z), 
                            type = "count"))                    
              }
            }
            else {
               plot(z, y, type = "n", axes = FALSE, xlab = zlab, ylab = "")
               axis(1)
               box()
            }
         }
         else {
            par(mar = c(3, 0.2, 1, 1) + 0.1, mgp = c(1.5, 0.2, 0), tcl = -0.2)
            plot(z, y, type = "n", axes = FALSE, xlab = zlab, ylab = "")
            axis(2, labels = FALSE)
            axis(1)
            box()
            points(z[locind], y[locind])
         }
      })
      # if (all(!is.na(panel$coords)))
      #    rp.text.change(panel, "textpane",
      #       paste("\n              location: (",
      #             signif(panel$coords[1], 5), ", ",
      #             signif(panel$coords[2], 5), ")", sep = ""))
      #       #       ")        radius: ", signif(panel$radius, 5), sep = ""))
      # else
      #    rp.text.change(panel, "textpane", "\n\n")
      panel
   }
   
   click <- function(panel, x, y) {
   	  dr1          <- diff(range(panel$x[ , 1]))
   	  dr2          <- diff(range(panel$x[ , 2]))
   	  if (eqscplot) {
   	     dr1 <- max(dr1, dr2)
   	     dr2 <- max(dr1, dr2)
   	  }
      d.pts        <- ((panel$x[ , 1] - x) / dr1)^2 + ((panel$x[ , 2] - y) / dr2)^2
      panel$locind <- which(d.pts <= panel$radius^2)
      panel$coords <- c(x, y)
      panel$zsdold <- panel$z.window.pars["width"]
      panel$z.window.pars["width"] <- panel$sdz * 4
      rp.control.put(panel$panelname, panel)
      rp.tkrreplot(panel, plot)
      rp.tkrreplot(panel, band)
      if (panel$location.plot.showing)
         rp.tkrreplot(panel, location)
      else {
         # rp.text(panel, "\n\n", name = "textpane", grid = "plots",
         #    row = 0, column = 1 + as.numeric(!panel$missing.y),
         #    sticky = "news", background = "white", fontsize = 12)
         rp.tkrplot(panel, location, draw.location, grid = "plots",
            row = 1, column = 1 + as.numeric(!panel$missing.y),
            hscale = panel$hscale, vscale = panel$vscale, background = "white")
      }
      panel$location.plot.showing <- TRUE
      panel
   }
   
   drag <- function(panel, x, y) {
   	  dr1          <- diff(range(panel$x[ , 1]))
   	  dr2          <- diff(range(panel$x[ , 2]))
   	  if (eqscplot) {
   	     dr1 <- max(dr1, dr2)
   	     dr2 <- max(dr1, dr2)
   	  }
      d.pts        <- ((panel$x[ , 1] - x) / dr1)^2 + ((panel$x[ , 2] - y) / dr2)^2
      panel$locind <- which(d.pts <= panel$radius^2)
      panel$coords <- c(x, y)
      rp.control.put(panel$panelname, panel)
      rp.tkrreplot(panel, plot)
      rp.tkrreplot(panel, location)
      panel
   }
   
   release <- function(panel, x, y) {
      if (!panel$retain.location.plot) {
         panel$z.window.pars["width"] <- panel$zsdold
         panel$locind <- integer(0)
         panel$coords <- rep(NA, 2)
         rp.control.put(panel$panelname, panel)
         rp.widget.dispose(panel, "location")
         rp.tkrreplot(panel, plot)
         # rp.tkrreplot(panel, location)
         rp.tkrreplot(panel, band)
         # rp.widget.dispose(panel, textpane)
         panel$location.plot.showing <- FALSE
      }
      panel
   }

   plot.3d <- function(panel) {
      with(panel, {
         rp.plot3d(x[ , 1], x[ , 2], z, col = colour, xlab = x1lab, ylab = x2lab, zlab = zlab)
      })
      panel
   }

   redraw4d <- function(panel) {
   	  if (panel$location.plot.showing) {
   	     if (panel$retain.location.plot)
            rp.do(panel, click, panel$coords[1], panel$coords[2])
         else {
            rp.do(panel, release, panel$coords[1], panel$coords[2])
            panel$location.plot.showing <- FALSE
            panel$coords <- rep(NA, 2)
         }
   	  }
   	  else {
         rp.tkrreplot(panel, plot)
         rp.tkrreplot(panel, band)
      }
      panel
   }
   
   xlab <- deparse(substitute(x))
   if (missing(x1lab)) {
      if (!is.null(colnames(x)[1]))
         x1lab <- colnames(x)[1]
      else
         x1lab <- paste(xlab, "1", sep = "-")
   }
   if (missing(x2lab)) {
      if (!is.null(colnames(x)[2]))
         x2lab <- colnames(x)[2]
      else
         x2lab <- paste(xlab, "2", sep = "-")
   }
   if (missing(zlab)) zlab <- deparse(substitute(z))
   missing.y <- missing(y)
   if (!missing.y) {
      if (missing(ylab)) ylab <- deparse(substitute(y))
   }
   else {
      y <- factor(rep(1, length(z)))
      ylab <- ""
   }
   panel.flag <- panel
   
   if (is.data.frame(x)) x <- as.matrix(x)
   if (!is.matrix(x) && ncol(x) == 2) stop("x should be a two-column matrix")
   w   <- cbind(x, z)
   if (!missing.y) w <- cbind(w, y)
   ind <- apply(w, 1, function(x) any(is.na(x)))
   x   <- x[!ind, ]
   z   <- z[!ind]
   if (!missing.y) y <- y[!ind]

   if (missing(model)) model <- NULL
   if (!is.null(model)) {
      if (!is.list(model)) {
         cat("model$y is not a list and will not be used.\n")
         model <- NULL
      }      
      else if (length(dim(model$y)) != 3) {
         cat("model$y is not a three-dimensional array and will not be used.\n")
         model$y <- NULL
      }
   }

   brks    <- NA
   natural <- NA
   missing.col.labels <- missing(col.labels)
   if (missing.col.labels) col.labels <- NA
   key <- 0.25
   if (missing.y) {
      ind         <- rep(1, length(z))
      col.palette <- "blue"
   }
   else if (is.factor(y)) {
      if (missing(col.palette) || all(is.na(col.palette)))
         col.palette <- topo.colors(nlevels(y))
      ind <- as.numeric(y)
   }
   else {
      if (missing(col.palette) || all(is.na(col.palette)))
         col.palette <- topo.colors(20)
      if (!missing(col.breaks)) {
         if (length(col.breaks) != length(col.palette) + 1)
            stop("the length of col.breaks should be length(col.palette) + 1.")
      	 brks <- col.breaks
      }
      else {
         rng  <- if (!is.list(model)) range(y) else range(y, model$y, na.rm = TRUE)
         del  <- 0.04 * diff(rng)
         brks <- seq(rng[1] - del, rng[2] + del, length = length(col.palette) + 1)
         # brks <- seq(rng[1], rng[2], length = length(col.palette) + 1)
      }
      natural <- missing.col.labels
      if (natural) col.labels <- brks
      ind     <- cut(y, brks, labels = FALSE)
      key     <- 0.15
   }
   colour <- col.palette[ind]
   clr    <- col2rgb(colour)
   clr    <- rgb2hsv(clr)
   # rad    <- mean(apply(x, 2, function(w) diff(range(w)))) / 20
   theta  <- seq(0, 2 * pi, length = 50)
   circle <- matrix(c(cos(theta), sin(theta)), ncol = 2)
   n      <- length(z)
      
   names(z.window.pars) <- c("location", "width")
   
   if (panel.flag) {
      panel <- rp.control(x = x, y = y, z = z, missing.y = missing.y,
                  x1lab = x1lab, x2lab = x2lab, ylab = ylab, zlab = zlab,
                  model = model, brks = brks,
                  col.palette = col.palette, col.labels = col.labels, natural = natural,
                  coords = rep(NA, 2), radius = radius, circle = circle, n = n, sdz = sd(z),
                  z.window = z.window, z.window.pars = z.window.pars,
                  colour = colour, clr = clr, hscale = hscale, vscale = vscale, 
                  location.plot.showing = FALSE,
                  retain.location.plot = retain.location.plot, eqscplot = eqscplot,
                  locind = integer(0),
                  background.plot = background.plot, foreground.plot = foreground.plot,
                  panel.plot = TRUE)
      rp.grid(panel, "controls", row = 0, column = 0, sticky = "n", background = "grey")
      rp.grid(panel, "plots",    row = 0, column = 1, background = "white")
      rp.tkrplot(panel, band, draw.band,
                hscale = hscale, vscale = 0.12 * vscale,
                grid = "plots", row = 0, column = 0, background = "white")
      if (location.plot)
         rp.tkrplot(panel, plot, draw.plot, click, drag, release,
                hscale = hscale, vscale = vscale,
                grid = "plots", row = 1, column = 0, background = "white")
      else
         rp.tkrplot(panel, plot, draw.plot,
                hscale = hscale, vscale = vscale,
                grid = "plots", row = 1, column = 0, background = "white")
      if (!missing.y)
         rp.tkrplot(panel, key,  draw.key, hscale = key * hscale, vscale = vscale,
                grid = "plots", row = 1, column = 1, background = "white")
      rp.slider(panel, z.window.pars, c(min(z), sd(z) / 20), c(max(z), sd(z) * 1.5), redraw4d,
                labels = c("centre", "width"),
                title = paste(zlab, "window"),
                grid = "controls", row = 0, column = 0)
      if (is.null(model))
         rp.radiogroup(panel, z.window, c("normal", "uniform"), action = redraw4d,
                grid = "controls", row = 1, column = 0, title = paste(zlab, "window shape"))
      if (location.plot) {
         rp.slider(panel, radius, 0.05 / 5,  0.05 * 5, redraw4d, 
                labels = paste(xlab, "window radius"),
                grid = "controls", row = 2, column = 0)
         rp.checkbox(panel, retain.location.plot, redraw4d, "Retain location plot",
                grid = "controls", row = 3, column = 0)
         if (is.factor(y) & require(lattice))
            rp.radiogroup(panel, location.plot.type, c("histogram", "density"), action = redraw4d,
                grid = "controls", row = 4, column = 0, title = "location plot type")
      }
      if (require(rgl))
         rp.button(panel, plot.3d, "3D plot", grid = "controls",
                row = 4 + as.numeric(require(lattice)), column = 0)
   }
   else {
      panel <- list(x = x, y = y, z = z, missing.y = missing.y,
                  x1lab = x1lab, x2lab = x2lab, ylab = ylab, zlab = zlab,
                  model = model, brks = brks, natural = natural,
                  col.palette = col.palette, brks = brks, col.labels = col.labels,
                  coords = coords, radius = radius, circle = circle, n = n, sdz = sd(z), 
                  colour = colour, clr = clr, eqscplot = eqscplot,
                  z.window = z.window, z.window.pars = z.window.pars,
                  background.plot = background.plot, foreground.plot = foreground.plot,
                  panel.plot = FALSE)
      draw.plot(panel)
   }
   
   invisible()
}

rp.spacetime <- function(space, time, y, model, 
                  col.palette, col.breaks, col.labels,
                  hscale = 1, vscale = hscale, panel = TRUE,
                  x1lab, x2lab, zlab, ylab, background.plot = NULL, foreground.plot = NULL,
                  time.window = "normal",
                  time.window.pars = c(min(time), sd(time)/5),
                  coords = rep(NA, 2), radius = 0.05, location.plot = TRUE,
                  retain.location.plot = FALSE, eqscplot = TRUE,
                  location.plot.type = "histogram") {

   xlab <- deparse(substitute(space))
   if (!is.null(colnames(space)[1]))
      x1lab <- colnames(space)[1]
   else
      x1lab <- paste(xlab, "1", sep = "-")
   if (!is.null(colnames(space)[2]))
      x2lab <- colnames(space)[2]
   else
      x2lab <- paste(xlab, "2", sep = "-")
   if (missing(zlab)) zlab <- deparse(substitute(time))
   missing.y <- missing(y)
   if (!missing.y) {
      if (missing(ylab)) ylab <- deparse(substitute(y))
   }
   else {
   	  y        <- jitter(rep(1, length(time)))
   	  ylab     <- ""
   }
   if (missing(model)) model <- NULL

   rp.plot4d(space, time, y, model,
                  col.palette = col.palette, col.breaks = col.breaks, col.labels = col.labels,
                  hscale = hscale, vscale = vscale, panel = panel,
                  x1lab = x1lab, x2lab = x2lab, zlab = zlab, ylab = ylab,
                  background.plot = background.plot, foreground.plot = foreground.plot,
                  z.window = time.window, z.window.pars = time.window.pars,
                  coords = coords, radius = radius, location.plot = location.plot,
                  retain.location.plot = retain.location.plot, eqscplot = eqscplot,
                  location.plot.type = location.plot.type)
}
