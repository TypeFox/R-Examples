#     Animated bubble plot

rp.bubbleplot <- function(x, y, year, size, col, col.palette = topo.colors(20),
                     interpolate = FALSE, fill.in = FALSE, labels = rownames(x),
                     hscale = 1, vscale = hscale) {
   
   size.label <- deparse(substitute(size))
   col.label  <- deparse(substitute(col))
   
   if (!require(rpanel)) stop("the rpanel package is required.")

   size.missing <- missing(size)
   if (size.missing)    size <- matrix(1, ncol = ncol(x), nrow = nrow(x))
   if (missing(col))    col  <- matrix("lightblue", ncol = ncol(x), nrow = nrow(x))
   if (is.vector(size)) size <- matrix(size, nrow = nrow(x), ncol = ncol(x))
   if (is.vector(col))  col  <- matrix(col,  nrow = nrow(x), ncol = ncol(x))
   scl  <- as.matrix(size)
   clr  <- as.matrix(col)
   if (is.numeric(clr)) {
      clr.brks <- seq(min(clr, na.rm = TRUE), max(clr, na.rm = TRUE), length = 21)
      del <- 0.001 * diff(range(clr, na.rm = TRUE))
      clr.brks[1] <- clr.brks[1] - del
      clr.brks[length(clr.brks)] <- clr.brks[length(clr.brks)] + del
   }
   xlab <- deparse(substitute(x))
   ylab <- deparse(substitute(y))
   
   # Remove those cases where no data are available
   ind  <- apply(x, 1, function(y) all(is.na(y))) | apply(x, 1, function(y) all(is.na(y)))
   x    <-    x[-ind, ]
   y    <-    y[-ind, ]
   size <- size[-ind, ]
   scl  <-  scl[-ind, ]
   clr  <-  clr[-ind, ]

   # Compute the scaling factors (cex)
   if (length(which(c(scl) < 0))) stop("the size information must be positive.")
   scl            <- 15 * sqrt(scl) / max(sqrt(scl), na.rm = TRUE)
   # scl[scl < 0.5] <- 0.5
   
   # This function fills in gaps of missing data with the largest previous value
   if (fill.in) {
      fn <- function(x) {
            for (i in 2:length(x)) {
            j <- which(!is.na(x[1:(i-1)]))
            if (is.na(x[i]) & length(j) > 0) x[i] <- x[max(j)]
         }
         x
      }
      x <- t(apply(x, 1, fn))
      y <- t(apply(y, 1, fn))
   }
                 
   bubble.draw <- function(panel) {
      with(panel, {
         # currently assumed to be years and integers
         if (interpolate) {
         	i    <- which(year == floor(year.ind))
         	p    <- year.ind - floor(year.ind)
            xi   <- (1 - p) *    x[ , i] + p *    x[ , min(i + 1, ncol(x))]
            yi   <- (1 - p) *    y[ , i] + p *    y[ , min(i + 1, ncol(x))]
            sizi <- (1 - p) * size[ , i] + p * size[ , min(i + 1, ncol(x))]
            scli <- (1 - p) *  scl[ , i] + p *  scl[ , min(i + 1, ncol(x))]
            if (is.numeric(clr)) {
               clri <- (1 - p) * clr[ , i] + p * clr[ , min(i + 1, ncol(x))]
               clri <- col.palette[findInterval(clri, clr.brks)]
            }
            else
               clri <- clr[ , i]
         }
         else {
            i    <- which(year == round(year.ind))
         	xi   <- x[ , i]
         	yi   <- y[ , i]
         	sizi <- scl[ , i]
         	scli <- scl[ , i]
         	clri <- clr[ , i]
         }
         plot(xi, yi, type = "n", xlab = xlab, ylab = ylab,
            xlim = range(x, na.rm = TRUE), ylim = range(y, na.rm = TRUE))
         xticks <- par()$xaxp[1] + (0:round(par()$xaxp[3])) * 
                                       (par()$xaxp[2] - par()$xaxp[1]) / par()$xaxp[3]
         yticks <- par()$yaxp[1] + (0:round(par()$yaxp[3])) * 
                                       (par()$yaxp[2] - par()$yaxp[1]) / par()$yaxp[3]
         text(mean(par()$usr[1:2]), mean(par()$usr[3:4]), year[which(year == round(year.ind))],
              col = "lightgrey", cex = 5)
         # text(x[,i], y[,i], rownames(x))
         segments(xticks, par()$usr[3], xticks, par()$usr[4], col = "grey")
         segments(par()$usr[1], yticks, par()$usr[2], yticks, col = "grey")
         ind <- rev(order(scli))
         points(xi[ind], yi[ind], cex = scli[ind], col = "black", bg = clri[ind], pch = 21)
         if (country != "none") {
         	id <- match(country, labels)
            points(xi[id], yi[id], cex = scli[id], col = "red", bg = "red", pch = 21)
         }
         if (all(!is.na(coords))) {
         	dst <- sqrt(((xi - coords[1])/par()$cxy[1])^2 + ((yi - coords[2])/par()$cxy[1])^2)
         	id  <- which((dst / scli) < 0.3)
         	if (length(id) > 0) {
         	   # xsgn <- sign(coords[1] - mean(par()$usr[1:2]))
         	   # if (xsgn == 0) xsgn <- 1
         	   # xpos <- coords[1] - xsgn * diff(par()$usr[1:2]) / 8
         	   xpos <- mean(par()$usr[1:2])
         	   ysgn <- sign(coords[2] - mean(par()$usr[3:4]))
         	   if (ysgn == 0) ysgn <- 1
         	   ypos <- coords[2] - ysgn * diff(par()$usr[3:4]) / 8
         	   legend(xpos, ypos, paste(labels[id], ": population ", sizi[id], sep = ""),
         	          fill = clri[id], xjust = 0.5, yjust = 1 - ysgn)
         	   # text(xi[id], yi[id], labels[id])
         	}
         }
         # points(xi[ind], yi[ind], cex = scli[ind], col = clr[ind])
         # for (j in ind) {
         #    points(xi[j], yi[j], cex = scl[j], col = clr[j], pch = 16)
         #    points(xi[j], yi[j], cex = scl[j])
         # }
         mtext(paste("The areas of the points are proportional to", size.label), line = 2)
      })
      panel
   } 

   bubble.redraw <- function(panel) {
      rp.tkrreplot(panel, plot)
      panel
   }

   key.draw <- function(panel) {
   	  par.mar <- par()$mar
   	  p2      <- par.mar[2]
   	  par(mar = par.mar * c(1, 0, 1, 1) + c(0, p2 %% floor(p2), 0, 0))
      rp.colour.key(panel$col.palette, panel$clr.brks)
      mtext(panel$col.label, side = 4, line = 1.5, font = 1)
      par(mar = par.mar)
      panel
   }

   bubble.movie.start <- function(panel) {
      panel$movie.stop <- FALSE
   	  # rp.control.put(panel$panelname, panel)      
   	  rp.timer(panel, 1, bubble.movie.call, function(panel) !panel$movie.stop)
      panel
   }
   
   bubble.movie.stop <- function(panel) {
      panel$movie.stop <- TRUE
      panel
   }
   
   bubble.movie.call <- function(panel) {
      if (panel$year.ind < max(panel$year)) {
         panel$year.ind <- panel$year.ind + (max(panel$year) - min(panel$year)) / 30
         panel$year.ind <- min(panel$year.ind, max(panel$year))
   	     # rp.control.put(panel$panelname, panel)
         rp.tkrreplot(panel, plot)
         rp.slider.change(panel, "slider", panel$year.ind)
      }
      else
         panel$movie.stop <- TRUE
      panel
   }
   
   click <- function(panel, x, y) {
   	  panel$coords <- c(x, y)
   	  # rp.control.put(panel$panelname, panel)
      rp.tkrreplot(panel, plot)
      panel
   }
   
   release <- function(panel, x, y) {
   	  panel$coords <- rep(NA, 2)
   	  # rp.control.put(panel$panelname, panel)
      rp.tkrreplot(panel, plot)
      panel
   }
   
   panel <- rp.control(x = x, y = y, year = year, year.ind = min(year), coords = rep(NA, 2),
                       scl = scl, clr = clr, clr.brks = clr.brks, col.palette = col.palette,
                       size.label = size.label, col.label = col.label,
                       interpolate = interpolate, country = "none", movie.stop = FALSE)
#   rp.grid(panel, "plot",      row = 0, column = 0)
#   rp.grid(panel, "key",       row = 0, column = 1)
#   rp.grid(panel, "listbox",   row = 0, column = 2)
#   rp.grid(panel, "controls1", row = 1, column = 0)
#   rp.grid(panel, "controls2", row = 1, column = 1)
   rp.tkrplot(panel, plot, bubble.draw, hscale = hscale, vscale = vscale,
      row = 0, column = 0, action = click, mousedrag = click, mouseup = release,
      background = "white")
   rp.tkrplot(panel, key,  key.draw, hscale = 0.15 * hscale, vscale = vscale,
      row = 0, column = 1, background = "white")
   rp.slider(panel, year.ind, min(year), max(year), bubble.redraw, labels = "Year", name = "slider",
      row = 1, column = 0)
   rp.button(panel, bubble.movie.start, "movie",
      row = 1, column = 1, sticky = "news")
   rp.button(panel, bubble.movie.stop, "movie stop",
      row = 1, column = 2, sticky = "news")
   rp.listbox(panel, country, labels = c("none", labels), rows = round(31 * vscale),
      action = bubble.redraw,
      row = 0, column = 2, title = "Country")

   invisible()
}
