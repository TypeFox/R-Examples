rp.surface <- function(surface, covariance, x1grid, x2grid, x, y, Display = "persp",
                       hscale = 1, vscale = hscale, Speed = 5, ntime = 10, ninterp = 50,
                       zlim = NULL, col.palette = topo.colors(100)) {

   if (!require(tkrplot)) stop("the tkrplot package is not available.")
   if (!require(akima))   stop("the akima package is not available.")
   # if (!require(lattice)) stop("the lattice package is not available.")
   # if (!require(rgl))     stop("the rgl package is not available.")
   
   x1lab <- deparse(substitute(x1grid))
   x2lab <- deparse(substitute(x2grid))
   ylab  <- deparse(substitute(surface))

   #   This needs to allow different grid lengths in different directions.
   n1grid  <- length(x1grid)
   n2grid  <- length(x2grid)
   xgrid   <- as.matrix(expand.grid(x1grid, x2grid))
   eig     <- eigen(covariance)
   e.vals  <- pmax(eig$values, 0)
   e.vecs  <- eig$vectors
   se.fit  <- sqrt(diag(covariance))
   surface <- c(surface)
   
   if (is.null(zlim))
      zlim <- range(surface - 3 * se.fit, surface + 3 * se.fit, na.rm = TRUE)
   brks <- seq(zlim[1], zlim[2], length = length(col.palette) + 1)

   col.fn <- function(panel) {
      rp.colour.key(panel$col.palette, brks, par.mar = c(3, 0, 1, 1.5) + 0.1, margin = TRUE)
      if (all(!is.na(panel$coords))) lines(rep(-0.5, 2), panel$cipos, lwd = 3)
      panel
   }

   #   Set up a finer grid and mask for the image plots
   xg1    <- seq(min(x1grid), max(x1grid), length = ninterp)
   xg2    <- seq(min(x2grid), max(x2grid), length = ninterp)
   xg     <- as.matrix(expand.grid(xg1, xg2))
   ind    <- apply(xg, 1, function(x) which.min((xgrid[ , 1] - x[1])^2 + (xgrid[ , 2] - x[2])^2))
   mask   <- as.numeric(!is.na(surface[ind]))
   mask[mask == 0] <- NA
   mask   <- matrix(mask, ncol = ninterp)
   ind.na <- !is.na(surface)

   rp.surface.draw <- function(panel) {
   	
      with(panel, {
      	   	
      surf <- surface
      if (animation) surf <- surf + e.isurf
      
      if (Display == "image") {
         surf <- interp(xgrid[ind.na, 1], xgrid[ind.na, 2], surf[ind.na], xg1, xg2)$z * mask
         par(mar = c(3, 3, 1, 0) + 0.1, mgp = c(1.5, 0.2, 0), tcl = -0.2)
         image(xg1, xg2, surf, zlim = zlim, col = col.palette, xlab = x1lab, ylab = x2lab)
      }
         # print(contourplot(c(surf) ~ xgrid[,1]*xgrid[,2], region = TRUE))
         # filled.contour(x1grid, x2grid, surf, zlim = zlim, color.palette = heat.colors)
      else if (Display == "persp") {
      	 surf <- matrix(surf, ncol = n2grid)
      	 ng1  <- 1:(n1grid - 1)
      	 ng2  <- 1:(n2grid - 1)
      	 clr  <- array(c(surf[ng1, ng2], surf[ng1 + 1, ng2], surf[ng1, ng2 + 1],
      	          surf[ng2 + 1, ng2 + 1]), dim = c(n1grid - 1, n2grid - 1, 4))
      	 ind  <- apply(clr, 1:2, function(x) if (length(which(is.na(x))) > 1) NA else 1)
      	 clr  <- apply(clr, 1:2, function(x) mean(x, na.rm = TRUE)) * ind
         clr  <- cut(c(clr), brks, labels = FALSE, include.lowest = TRUE)
         clr  <- col.palette[clr]
         par(mar = c(3, 3, 1, 0.5) + 0.1, mgp = c(1.5, 0.2, 0), tcl = -0.2)
         persp(x1grid, x2grid, surf,
                zlim = zlim, theta = theta, phi = phi, ticktype = "detailed", col = clr,
                d = 4, xlab = x1lab, ylab = x2lab, zlab = ylab)
      }
      
      # Draw the confidence interval against the colour key
      if (all(!is.na(coords)) & Display == "image" & !animation) {
         ind   <- !is.na(surface)
         pred1 <- interp(xgrid[ind, 1], xgrid[ind, 2], surface[ind], coords[1], coords[2])$z
         se1   <- interp(xgrid[ind, 1], xgrid[ind, 2], se.fit[ind],  coords[1], coords[2])$z
         # Is pretty ok here?
         col.range <- range(pretty(range(pred1), 20))
         ci <- c(pred1 - 2 * se1, pred1 + 2 * se1)
         cipos <- c((ci[1] - col.range[1]) / diff(col.range),
                    (ci[2] - col.range[1]) / diff(col.range))
         # cipos <- min(x2grid) + cipos * diff(range(x2grid))
         # Will the weights below always work?  Change to image and separate colour key.
         # xpos <- 0.28 * min(x1grid) + 0.72 * max(x1grid)
         panel$cipos <- ci
         rp.control.put(panel$panelname, panel)
         rp.tkrreplot(panel, key)
      }

#      if (animation == "none") {
#         scaling <- rp.plot3d(x[,1], y, x[,2], ylim = range(y, e.sim),
#            xlab = "Longitude", ylab = "Score1", zlab = "Latitude")
#         save(scaling, file = "scaling.dmp")
#         xgrid <- cbind(x1grid, x2grid)
#         sm.surface3d(xgrid, e.sim[ , , 1],
#                    scaling = scaling, zlim = range(e.sim))
#      }
      
      if (animation == "rgl") {
#         scaling <- rp.plot3d(x[,1], y, x[,2], ylim = range(y, e.sim))
         load("scaling.dmp")
         xgrid <- cbind(x1grid, x2grid)
#         sm.surface3d(xgrid, e.sim[ , , 1],
#                    scaling = scaling, zlim = range(e.sim))
         for (i in 2:nsim) {
            for (wt in seq(0, 1, length = 15)) {
               par3d(skipRedraw = TRUE)
               pop3d()
               pop3d()
               sm.surface3d(xgrid,
                    (1 - wt) * e.sim[ , , i - 1] + wt * e.sim[ , , i],
                    scaling = scaling, zlim = range(e.sim))
               par3d(skipRedraw = FALSE)
               # Sys.sleep(0.1)
            }  
         }
      }
      
      })
      
      panel
   }
   
   rp.surface.redraw <- function(panel) {
      rp.tkrreplot(panel, plot)
      panel
   }
   
   rp.key.redraw <- function(panel) {
      rp.tkrreplot(panel, key)
      panel
   }
   
   animate <- function(panel) {
      panel$animation  <- !panel$animation
      if (panel$animation) {
         panel$e.sim      <- panel$e.vecs %*% diag(sqrt(panel$e.vals)) %*%
                                     rnorm(length(panel$surface))
         panel$e.sim.old  <- rep(0, length(panel$surface))
         panel$isurf      <- 1
   	     rp.control.put(panel$panelname, panel)      
   	     rp.timer(panel, 1, animation.call, function(panel) panel$animation)
         panel$e.isurf    <- rep(0, length(panel$surface))
  	     rp.control.put(panel$panelname, panel)      
         rp.tkrreplot(panel, plot)
      }
      panel
   }
   
   animation.call <- function(panel) {
   	  Sys.sleep(0.01 + panel$Speed / 100)
   	  if (panel$isurf == ntime + 1) {
         panel$e.sim.old <- panel$e.sim
         panel$e.sim     <- panel$e.vecs %*% diag(sqrt(panel$e.vals)) %*%
                                  rnorm(length(panel$surface))
   	  	 panel$isurf     <- 1
   	  }
   	  wt            <- panel$isurf / ntime
   	  wt1           <- sqrt(wt^2 + (1 - wt)^2)
      panel$e.isurf <- panel$e.sim.old * (1 - wt) / wt1 + panel$e.sim * wt / wt1
      panel$isurf   <- panel$isurf + 1
   	  rp.control.put(panel$panelname, panel)
      rp.tkrreplot(panel, plot)
      panel
   }
   
   mouse <- function(panel, x, y) {
   	  panel$coords <- c(x, y)
   	  rp.control.put(panel$panelname, panel)
   	  rp.tkrreplot(panel, plot)
      panel
   }

   release <- function(panel, x, y) {
   	  panel$coords <- c(NA, NA)
   	  rp.control.put(panel$panelname, panel)
   	  rp.tkrreplot(panel, plot)
      panel
   }

   panel <- rp.control(x = x, y = y, surface = c(surface), se.fit = se.fit,
                       e.vecs = e.vecs, e.vals = e.vals,
                       x1grid = x1grid, x2grid = x2grid, xgrid = xgrid,
                       n1grid = n1grid, n2grid = n2grid,
                       zlim = zlim, coords = rep(NA, 2), theta = -30, phi = 40,
                       Display = Display, animation = FALSE,
                       Speed = Speed, ntime = ntime, ninterp = ninterp,
                       xg1 = xg1, xg2 = xg2, mask = mask, ind.na = ind.na,
                       x1lab = x1lab, x2lab = x2lab, ylab = ylab,
                       brks = brks, col.palette = col.palette, col.fn = col.fn)
   rp.grid(panel, "controls", row = 0, column = 0, sticky = "n")
   rp.grid(panel, "plot",     row = 0, column = 1, background = "white")
   rp.grid(panel, "key",      row = 0, column = 2, background = "white")
   rp.tkrplot(panel, plot, rp.surface.draw, mouse, mouse, release,
              hscale = hscale, vscale = vscale,
              grid = "plot", row = 0, column = 0, background = "white")
   rp.tkrplot(panel, key, col.fn, hscale = 0.15 * hscale, vscale = vscale,
              grid = "key", row = 0, column = 0, background = "white")
   rp.radiogroup(panel, Display, c("persp", "image"), action = rp.surface.redraw,
              grid = "controls", row = 0, column = 0)
   rp.slider(panel, theta, -180, 180, rp.surface.redraw, "persp angle 1",
              grid = "controls", row = 1, column = 0)
   rp.slider(panel, phi,      0,  90, rp.surface.redraw, "persp angle 2",
              grid = "controls", row = 2, column = 0)
   rp.button(panel, animate, "Animate: on/off",
              grid = "controls", row = 3, column = 0)
   rp.doublebutton(panel, Speed, 0.95, log = TRUE, action = rp.surface.redraw,
              grid = "controls", row = 5, column = 0)
   
}
