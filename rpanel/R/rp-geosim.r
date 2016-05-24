# To be done:
#    - the variogram plot no longer disappears when the checkboxes are not checked.
# Display the values of the anisotropy parameters in the plot title.

rp.geosim <- function(max.Range = 0.5, max.pSill = 1, max.Nugget = 1, max.Kappa = 10,
                      max.aniso.ratio = 5,
                      min.ngrid = 10, max.ngrid = 25, hscale = NA, vscale = hscale,
                      col.palette = terrain.colors(40)) {

field.new <- function(panel) {
	
      r <- panel$Range / (2 * sqrt(panel$kappa))
      panel$family <- "matern"
      # Old code which matched the surface and data grids - no longer required.
      # cgrid <- ceiling((panel$smgrid - 1) / (panel$ngrid - 1))
      # panel$smgrid <- 1 + cgrid * (panel$ngrid - 1)
      pars.changed <- (panel$ngrid != panel$ngrid.old) | (panel$Range != panel$Range.old) |
                      (panel$pSill != panel$pSill.old) | (panel$kappa != panel$kappa.old)
      warn <- options()$warn
      options(warn = -1)
      if (!panel$points.only | pars.changed)
         panel$fieldsm <- grf(panel$smgrid^2, grid = "reg", cov.model = panel$family,
            cov.pars = c(panel$pSill, r), nugget = 0, messages = FALSE, kappa = panel$kappa,
            aniso.pars = c(panel$aniso.angle, panel$aniso.ratio))
      panel$fieldnug <- grf(panel$ngrid^2, grid = "reg", cov.model = "pure.nugget",
         cov.pars = c(panel$Nugget, 0), messages = FALSE)
      options(warn = warn)
      # igrid <- seq(1, panel$smgrid, by = cgrid)
      igrid <- 1 + round(((1:panel$ngrid) - 1) * (panel$smgrid - 1) / (panel$ngrid - 1))
      igrid <- as.matrix(expand.grid(igrid, igrid))
      panel$fieldsm$data <- matrix(panel$fieldsm$data, ncol = panel$smgrid)
      panel$data <- panel$fieldsm$data[igrid] + panel$fieldnug$data
      
      panel$ngrid.old  <- panel$ngrid
      panel$Range.old  <- panel$Range
      panel$pSill.old  <- panel$pSill
      panel$Nugget.old <- panel$Nugget
      panel$kappa.old  <- panel$kappa
   
   if (!panel$first)
      rp.do(panel, graphics.update)
   else
      panel$first <- FALSE
   
   panel
   }
   
graphics.update <- function(panel) {
	
   rp.tkrreplot(panel, plot1)
   
   if (any(panel$vgm.checks)) {
      if (!panel$vgm.present) {
      	  rp.tkrplot(panel, plot2, vario.update, hscale = panel$hscale, vscale = panel$vscale, 
             grid = "rightplot", row = 0, column = 0)
         panel$vgm.present <- TRUE
         }
      else
         rp.tkrreplot(panel, plot2)
      }
      
   if (("rgl plot" %in% names(panel$display.checks)) && (panel$display.checks["rgl plot"])) {
      w  <- panel$data
      x  <- panel$fieldnug$coords[, 1]
      y  <- panel$fieldnug$coords[, 2]
      if ("rgl.id" %in% names(panel))
         try.out <- try(rgl.set(panel$rgl.id), silent = TRUE)
      else
         try.out <- "try-error"
      if (is.null(try.out)) {
      	 ind <- (rgl.ids()$type == "points")
         if (any(ind))  pop3d(id = rgl.ids()$id[ind])
      	 ind <- (rgl.ids()$type == "surface")
         if (any(ind))  pop3d(id = rgl.ids()$id[ind])
         }
      else
         panel$scaling <- rp.plot3d(y, w, x, col = "red",
             ylim = panel$z.range, ylab = "z", xlab = "x", zlab = "y", type = "n")
      if (panel$display.checks["points"]) {
         a <- panel$scaling(y, w, x)
         panel$points.id <- points3d(a$z, a$y, a$x, col = "red", size = 3)
         }
      if (panel$display.checks["surface"]) {
         z  <- panel$fieldsm$data
         brks <- seq(panel$z.range[1], panel$z.range[2], length = length(panel$col.palette) + 1)
   	  	 brks[1] <- min(brks[1], min(z) - 1)
   	  	 brks[length(brks)] <- max(brks[length(brks)], max(z) + 1)
         clr  <- cut(z, brks, labels = FALSE)
         cols <- panel$col.palette[clr]
         x    <- seq(0, 1, length = panel$smgrid)
         y    <- seq(0, 1, length = panel$smgrid)
         a    <- panel$scaling(x, z, y)
         panel$surface.id <- rgl.surface(a$x, a$z, a$y, col = cols)
         }
      bg3d("white")
      panel$rgl.id  <- rgl.cur()
      }
   else if (panel$rgl.old) {
      try.out <- try(rgl.set(panel$rgl.id), silent = TRUE)
      if (is.null(try.out)) rgl.close()
      }
   
   if (("rgl plot" %in% names(panel$display.checks)) && (panel$display.checks["rgl plot"]))
      panel$rgl.old <- panel$display.checks["rgl plot"]
   
   panel
   }
   
cont.update <- function(panel) {
   mar.old <- par()$mar
   with(panel, {
   	  layout(matrix(1:2, ncol = 2), widths = c(0.84, 0.16))
   	  g <- seq(0, 1, length = smgrid)
   	  par(mar = c(5, 4, 4, 0) + 0.1)
   	  plot(g, g, type = "n", xlab = "x", ylab = "y", xaxs = "i", yaxs = "i")
   	  if (display.checks["surface"]) {
         image(x = g, y = g, z = matrix(fieldsm$data, nrow = smgrid), add = TRUE,
                     zlim = z.range, col = col.palette)
         box()
   	     }
   	  if (display.checks["points"]) {
   	  	 brks <- seq(z.range[1], z.range[2], length = length(col.palette) + 1)
   	  	 brks[1] <- min(brks[1], min(data) - 1)
   	  	 brks[length(brks)] <- max(brks[length(brks)], max(data) + 1)
         clr  <- cut(data, brks, labels = FALSE)
   	     g    <- seq(0, 1, length = ngrid)
         points(rep(g, ngrid), rep(g, each = ngrid))
         points(rep(g, ngrid), rep(g, each = ngrid), col = col.palette[clr], pch = 16)
         points(rep(g, ngrid), rep(g, each = ngrid))
   	     }
      title(paste("Range =", Range, "  Partial sill =", pSill, "  Nugget=", Nugget), 
                    line = 2, cex = 1)
      title(paste("Kappa =", kappa, "  Data grid =", ngrid, "x", ngrid),   line = 1, cex = 1)
      if (aniso.ratio > 1)
         title(paste("Anisotropy: ratio =", round(aniso.ratio, 2), 
                                " angle =", round(aniso.angle, 2)), line = 0, cex = 1)

   	  par(mar = c(5, 2, 4, 2) + 0.1)
      rp.colour.chart(col.palette, z.range)
      layout(1)
      })
   par(mar = mar.old)
   panel
   }

vario.update <- function(panel) {
   if (panel$aniso.ratio > 1) {
      plot(1:10, type = "n", axes = FALSE, xlab = "", ylab = "")
      return(panel)
   }
   with(panel, {
      vg  <- pSill * (1 - matern(xa, Range / (2 * sqrt(kappa)), kappa))
      fld <- data
      g   <- seq(0, 1, length = ngrid)
      vga <- variog(as.geodata(cbind(fieldnug$coords, fld)), messages=FALSE,
                    max.dist = 0.7)
      plot(vga$v ~ vga$u, xlim = c(0, 0.7), ylim = c(0, 2), type = "n",
         xlab = "Distance", ylab = "Semivariogram")
      if (vgm.checks["sample"])
         points(vga$v ~ vga$u, col = "red")
      if (vgm.checks["true"]) {
         if (Range > 0) {
            lines(xa, vg + Nugget, col = "blue")
            abline(h = pSill + Nugget, lty = 2)
            ind <- min(which(vg >= 0.95 * pSill))
            if (length(ind) > 0) abline(v = xa[ind], lty = 2)
            }
         else
            lines(xa, rep(pSill + Nugget, length(xa)), col = "blue")
      }
      title("Semivariogram", line = 2, cex = 1)
      title("True (blue)    Sample (red)", line = 1, cex = 1)
      })
   panel
   }
   
   if (!require(tkrplot))
      stop("the tkrplot package is not available.")
   if (!require(geoR))
      stop("the geoR package is not available.")
   if (!require(RandomFields))
      stop("the RandomFields package is not available.")
      
   if (is.na(hscale)) {
      if (.Platform$OS.type == "unix") hscale <- 1.2
      else                             hscale <- 1.4
      }
   if (is.na(vscale)) 
      vscale <- hscale
      
   if (min.ngrid < 8) {
      cat("min.ngrid reset to 8.\n")
      min.ngrid <- 8
      }
   if (max.ngrid < min.ngrid) {
      cat("max.ngrid reset to", min.ngrid, "\n")
      max.ngrid <- min.ngrid
      }
   if (max.Range <= 0) {
      cat("max.Range reset to 0.5\n")
      max.Range <- 0.5
      }
   if (max.pSill <= 0) {
      cat("max.pSill reset to 1\n")
      max.pSill <- 1
      }
   if (max.Nugget <= 0) {
      cat("max.Nugget reset to 1\n")
      max.Nugget <- 1
      }
   if (max.Kappa <= 0.5) {
      cat("max.Kappa reset to 10\n")
      max.Kappa <- 10
      }
      
   display.checks.init <- c(TRUE, FALSE)
   checks.lbls <- c("surface", "points")
   if (require(rgl)) {
      display.checks.init <- c(display.checks.init, FALSE)
      checks.lbls <- c(checks.lbls, "rgl plot")
      }
   names(display.checks.init) <- checks.lbls
      
   panel <- rp.control("Spatial correlation",
      first = TRUE,
      xa = seq(0, 0.7, by = 0.01), smgrid = 25, z.range = c(-5, 5),
      ngrid = 15, Range = 0.1, pSill = 0.5, Nugget = 0, kappa = 4,
      ngrid.old = 15, Range.old = 0.05, pSill.old = 0.5, kappa.old = 4,
      aniso.angle = 0, aniso.ratio = 1,
      hscale = hscale, vscale = vscale, rgl.old = FALSE, vgm.present = FALSE,
      display.checks = display.checks.init,
      vgm.checks = c(true = FALSE, sample = FALSE),
      points.only = FALSE, points.id = NA, surface.id = NA, col.palette = col.palette)
   rp.do(panel, field.new)
   rp.grid(panel, "controls",  row = 0, column = 0)
   rp.grid(panel, "leftplot",  row = 0, column = 1, background = "white")
   rp.grid(panel, "rightplot", row = 0, column = 2, background = "white")
   rp.tkrplot(panel, plot1, cont.update, hscale = hscale, vscale = vscale, 
            grid = "leftplot", row = 0, column = 0, sticky = "ew")
   rp.button(panel, field.new, "New simulation",     
            grid = "controls", row = 1, column = 0, sticky = "ew")
   rp.slider(panel, Range, 0, max.Range,  field.new, "Range",     
            grid = "controls", row = 2, column = 0, sticky = "ew")
   rp.slider(panel, pSill, 0, max.pSill,  field.new, "Partial sill",     
            grid = "controls", row = 3, column = 0, sticky = "ew")
   rp.slider(panel, Nugget, 0, max.Nugget, field.new, "Nugget",    
            grid = "controls", row = 4, column = 0, sticky = "ew")
   rp.checkbox(panel, display.checks, graphics.update, checks.lbls, title = "Display", 
       grid = "controls", row = 5, column = 0, sticky = "ew")
   rp.checkbox(panel, vgm.checks, graphics.update, c("true", "sample"), title = "Variogram", 
       grid = "controls", row = 6, column = 0, sticky = "ew")
   rp.checkbox(panel, points.only, title = "Sample points only", 
       grid = "controls", row = 7, column = 0)
   rp.slider(panel, ngrid, min.ngrid, max.ngrid, field.new, "Data grid", 
       resolution = 1, grid = "controls", row = 8, column = 0, sticky = "ew")
   rp.slider(panel, kappa, 0.5, max.Kappa, field.new, "Kappa",     
       grid = "controls", row = 9, column = 0, sticky = "ew")
   rp.slider(panel, aniso.angle, 0, pi, field.new, "Anisotropy angle",
       grid = "controls", row = 10, column = 0, sticky = "ew")
   rp.slider(panel, aniso.ratio, 1, max.aniso.ratio, field.new, "Anisotropy ratio",     
       grid = "controls", row = 11, column = 0, sticky = "ew")
   rp.do(panel, graphics.update)
   }

rp.colour.chart <- function(cols, zlim)  {
   ngrid <- length(cols)
   plot(0:1, zlim, type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n",
        xlab = "", ylab = "")
   axis(4)
   xvec <- rep(0, ngrid)
   yvec <- seq(zlim[1], zlim[2], length = ngrid + 1)
   rect(xvec, yvec[-length(yvec)], xvec + 1, yvec[-1], col = cols, border = NA)
   box()
   }
