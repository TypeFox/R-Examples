# What scale of measurement is appropriate?
# Consider trend, cov.pars and nugget settings.
# Remove redundant items from mururoa.list. (trendm, setRF)
# A yes/no dialogue box for the "Take sample" button would be good.

rp.mururoa <- function(hscale = NA, col.palette = rev(heat.colors(40)), col.se = "blue",
                 file = NA, parameters = NA) {

mururoa.points <- function(panel) {
	
   if (panel$sample.taken) return(panel)
	
   if (panel$random.alignment) {
      panel$gx <- runif(1)
      panel$gy <- runif(1)
      }
   panel$random.alignment.old <- panel$random.alignment
	
   if (panel$stype == "Random") {
      wh          <- which(panel$mur.mat > 0)     # wh is a vector length 20301
      sind        <- sample(wh, panel$npts) - 1
      panel$sampx <- (sind %%  201) / 2           # convert pos. in vector to x, y
      panel$sampy <- (sind %/% 201) / 2           # coords using quotient-remainder
      }
      
   else if (panel$stype == "Grid") {
   	  gsp         <- sqrt(8522 / (4 * panel$npts))
      txc         <- round(2 * seq(panel$gx * gsp, 100, by = gsp)) / 2
      tyc         <- round(2 * seq(panel$gy * gsp,  50, by = gsp)) / 2
      gr          <- expand.grid(txc, tyc)
      a           <- which(panel$mur.mat[1 + 2 * txc, 1 + 2 * tyc] > 0)
      gr          <- gr[a, ]
      panel$sampx <- gr[ , 1]
      panel$sampy <- gr[ , 2]
      panel$txc   <- txc
      panel$tyc   <- tyc
      }
      
   else if (panel$stype == "Transect") {
   	  a  <- atan(panel$transect.angle)        # Work with axes rotated by angle a.
      sp <- sqrt(8522 / (4 * panel$npts))     # Get (average) point spacing required.
      tt <- round(44 * cos(a) / sp)           # Set no. of transects - cos(a) correction due to rotation
      gy <- panel$gy
      gx <- panel$gx
      if (panel$tsb == "Random")
        ty <- sample(1:89, tt)
      else {
        ty <- round(seq(1 + gy * 2 * sp / cos(a), 89, by = 2 * sp / cos(a)))
        ty <- ty[ty < 90]
        }
                                              # choose which transects from the 89 possibilities
      tt <- length(ty)                        # reset in case one has been dropped
      L  <- rep(0, tt)
      for (j in 1:tt) {                       # find relative lengths of trans.
        L[j] <- sum(panel$mur.tran[ , ty[j]]) # ty(j)-th column of mur.tran
        }
      L1   <- sum(L)
      grx  <- c()
      gry  <- c()
      grxp <- c()
      gryp <- c()
      if (panel$tsw == "Systematic") {
        sp <- L1 / (2 * cos(a) * panel$npts)                  # Adjust spacing to reflect
        for (j in 1:tt) {                                     # total transect length.
          x1   <- seq(gx * sp, 103, by = sp)                       # coords in rotated system
          x    <- cos(a) * x1 - panel$transect.angle * (-10 + (ty[j] - 1) / 2) # coords in usual system
          # x    <- cos(a) * x1 - 0.22 * (-10 + (ty[j] - 1) / 2) # coords in usual system
          xp   <- x
          x    <- round(2 * x) / 2                             # align to grid
          x    <-  x[x  >= 0 & x  <= 100]                         # remove out-of-range x
          xp   <- xp[xp >= 0 & xp <= 100]
          x    <-  x[panel$mur.tran[1 + 2 * x,  ty[j]] == 1]     # remove pts outside Mururoa
          xp   <- xp[panel$mur.tran[1 + 2 * xp, ty[j]] == 1]
          grx  <- c(grx,  x)
          grxp <- c(grxp, xp)
          gry  <- c(gry,  -10 + (ty[j] - 1) / 2 + round(panel$transect.angle * x  * 2) / 2)
          gryp <- c(gryp, -10 + (ty[j] - 1) / 2 + round(panel$transect.angle * xp * 2) / 2)
          # gry  <- c(gry,  -10 + (ty[j] - 1) / 2 + round(0.22 * x  * 2) / 2)
          # gryp <- c(gryp, -10 + (ty[j] - 1) / 2 + round(0.22 * xp * 2) / 2)
          }
        }
      else {
        for (j in 1:tt) {
          x <- 0:200/2                                    # candidate x-coords
          x <- x[panel$mur.tran[ , ty[j]] == 1]           # pts inside Mururoa
          samp.size <- min(length(x), round(panel$npts * L[j] / L1))
          if (samp.size >= 1) {
            xpick <- sample(x, samp.size)                 # take suitable sample
            grx   <- c(grx, xpick)
            gry   <- c(gry, -10 + (ty[j] - 1) / 2 + round(panel$transect.angle * xpick * 2) / 2)
            # gry   <- c(gry, -10 + (ty[j] - 1) / 2 + round(0.22 * xpick * 2) / 2)
            }
          }
        }
      panel$sampx <- grx
      panel$sampy <- gry
      panel$plotx <- grxp
      panel$ploty <- gryp
      panel$tyc   <- -10 + (ty - 1) / 2               # store y-intercepts for later
      }
      
   if (panel$sampling.started) {
   	  rp.control.put(panel$panelname, panel)
      rp.tkrreplot(panel, plot1)
   }
   
   panel 
   }

mururoa.draw <- function(panel) {
	
  ind <- which(panel$trend.setting == c("cte", "1st", "2nd"))

  with(panel, {
  	
  plot(murx, mury, type = "n", asp = 1,
          xlab = "easting", ylab = "northing", main = paste("Number of points = ", length(sampx)))
  polygon(murx, mury, col = col.inside, density = -1, border = NA)
  
  if (sample.taken) {
  	
    if (display.options["predicted surface"])
       image(seq(0, 100,len = 201), seq(0, 50, len = 101), kpred[,,ind], 
              add = TRUE, col = col.palette, zlim = zlim)

    if (display.options["prediction s.e."])
       contour(seq(0, 100,len = 201), seq(0, 50, len = 101), kse[,,ind], add = TRUE,
              levels = pretty(range((kse[,,ind] * mask)[kse[,,ind] > 0], na.rm = TRUE), 10), col = col.se)
    
    if (display.options["points"]) {
   	   brks <- seq(zlim[1], zlim[2], length = length(col.palette) + 1)
       clr  <- cut(sz[ , 3], brks, labels = FALSE, include.lowest = TRUE, right = FALSE)
       points(sz[ , 1], sz[ , 2], col = col.palette[clr], pch = 16)
       points(sz[ , 1], sz[ , 2])
       }
       
    if (mururoa.true) 
       image(seq(0, 100,len = 201), seq(0, 50, len = 101), true.surface,
              add = TRUE, col = col.palette, zlim = zlim)
       
    }  
    
  else if (panel$sampling.started) {
  
     sx <- sampx
     sy <- sampy
     if (stype == "Grid") {
        x0 <- txc
        y0 <- tyc
        segments(x0,  0,  x0, 50, lty = 2, col = col.points)
        segments( 0, y0, 100, y0, lty = 2, col = col.points)
        }
     if (stype == "Transect") {
        y0 <- tyc
        segments(0, y0, 100, y0 + round(100 * panel$transect.angle), lty = 2, col = col.points)
        # segments(0, y0, 100, y0 + 22, lty = 2, col = col.points)
        if (tsb == "Systematic" & tsw == "Systematic") {
     	   # The following code makes sure that the plotted points lie on the
     	   # transect, rather than using the grid positions.
           # points(plotx, ploty, pch = 16, col = col.points)
           ang <- transect.angle
           for (i in 1:length(sx)) {
              dd  <- (cos(ang) * sx[i] + sin(ang) * (sy[i] - y0))
              pt  <- cbind(dd * cos(ang), dd * sin(ang) + y0)
              dst <- (pt[ , 1] - sx[i])^2 + (pt[ , 2] - sy[i])^2
              ind <- min(which(dst == min(dst)))
              points(pt[ind, 1], pt[ind, 2], pch = 16, col = col.points)
              }
           }
        }
     if (!(stype == "Transect" & tsb == "Systematic" & tsw == "Systematic"))
        points(round(2 * sx) / 2, round(2 * sy) / 2, pch = 16, col = col.points)
     axis(1, at = seq(0, 100, by = 10))
     }
     
     usr <- par()$usr
     polygon(c(murx, murx[1], usr[1],  usr[1], usr[2], usr[2], usr[1], usr[1]), 
             c(mury, mury[1], mury[1], usr[4], usr[4], usr[3], usr[3], mury[1]),
             col = col.outside, density = -1, border = NA)
     box()
     lines(c(murx, murx[1]), c(mury, mury[1]), col = col.border, lwd = 3)

  })
    
  panel
  }
  
mururoa.predict <- function(panel) {

  ind <- which(panel$trend.setting == c("cte", "1st", "2nd"))
      
  if (!panel$prediction.computed[ind]) {
     sz   <- panel$sz
     z0 <- sz
     fg <- as.geodata(z0)                            # convert to format usable by geoR
     vg2 <- variog(fg, trend = panel$trend.setting, max.dist = 60, messages = FALSE)  # quadratic trend
     vfit <- variofit(vg2, ini.cov.pars = c(max(vg2$v), 20), nugget = min(vg2$v), cov.model = "matern", 
                      kappa = 4, messages = FALSE)
     xx <- 0:60
     yy <- 100 + 625 * (xx / 20 - 0.5 * (xx / 30)^3)
     yy[32:61] <- 725
     pred.grid <- expand.grid(0.25 + 0:199/2, 0.25 + 0:99/2)      # offset from samples
     KC  <- krige.control(obj.model = vfit, trend.d = panel$trend.setting, trend.l = panel$trend.setting)
     shh <- output.control(messages = FALSE)
     kc  <- try(krige.conv(fg, locations = pred.grid, krige = KC, output = shh), silent = TRUE)
     if (class(kc) == "try-error" & panel$trend.setting == "cte") {
        rp.messagebox("There are numerical problems in producing predictions with this model.  Try fitting a linear or quadratic trend function.")
        return(panel)
     }
     kcm <- matrix(kc$predict, nrow = 200)
     pg2 <- expand.grid(0:200 / 2, 0:100 / 2)               # use 'real' grid to find RSS
     KC2 <- krige.control(obj.model = vfit, trend.d = panel$trend.setting, trend.l = panel$trend.setting)
     kc2 <- krige.conv(fg, locations = pg2, krige = KC2, output = shh)
     panel$kpred[ , , ind] <- matrix(kc2$predict, ncol = 101) * panel$mur.mat
     panel$kse[ , , ind]   <- sqrt(matrix(kc2$krige.var, ncol = 101))
     panel$prediction.computed[ind] <- TRUE
     }
     
  panel$zlim  <- range(panel$zlim,  panel$kpred[ , ,ind] * panel$mask,
                       panel$true.surface, na.rm = TRUE)
  # panel$zlim  <- range(panel$sz[ , 3],  panel$kpred[,,1] * panel$mask,
  #                      panel$kpred[,,2] * panel$mask,
  #                      panel$kpred[,,3] * panel$mask, panel$true.surface, na.rm = TRUE)
  rp.control.put(panel$panelname, panel)
  rp.do(panel, mururoa.colour.chart.redraw)
  rp.do(panel, mururoa.samp.redraw)

  panel
  }
  
mururoa.colour.reset <- function(panel) {
  ind <- which(panel$trend.setting == c("cte", "1st", "2nd"))
  panel$zlim  <- range(panel$sz[ , 3],  panel$kpred[ , ,ind] * panel$mask,
                       panel$true.surface, na.rm = TRUE)
  rp.control.put(panel$panelname, panel)
  rp.do(panel, mururoa.colour.chart.redraw)
  rp.do(panel, mururoa.samp.redraw)
  panel
  }
  
mururoa.samp <- function(panel) {

  if (!panel$sampling.started) return(panel)
	
  panel$sample.taken <- TRUE
  x <- panel$sampx
  y <- panel$sampy
  x <- round(2 * x) / 2               # round to half-integer
  y <- round(2 * y) / 2
  x1 <- x[which(panel$mur.mat[1 + 2 * x + 201 * 2 * y] > 0)]
  y1 <- y[which(panel$mur.mat[1 + 2 * x + 201 * 2 * y] > 0)]
  dup <- c()                          # remove duplicate points
  xn <- length(x1)
  for (i in 1:(xn - 1)) {
    for (j in seq(i + 1, xn)) {
      if (x1[i] == x1[j] & y1[i] == y1[j]) dup <- c(dup, j)
    }
  }
  if (length(dup) > 0) x1 <- x1[-dup]
  if (length(dup) > 0) y1 <- y1[-dup]
  #  if (length(x1)>50) {
  #    x1 <- x1[1:50]                               # [restrict to 50 sampling points]
  #    y1 <- y1[1:50] }
  a <- 1 + 2 * x1 + 201 * 2 * y1                  # find right place in trend etc. vectors
  ptsm <- as.matrix(expand.grid(seq(0, 100, by = 0.5), seq(0, 50, by = 0.5)))
  warn <- options()$warn
  options(warn = -1)
  nug  <- grf(nrow(ptsm), ptsm, cov.model = "pure.nugget", cov.pars = c(panel$nugget, 0), 
                  nugget = 0, messages = FALSE)
  options(warn = warn)
  z <- panel$field[a] + c(panel$trendmurmat)[a] + nug$data[a]
  sz <- cbind(x1, y1, z)                          # return 3 columns with x,y,Z values
  panel$sz <- sz
  mask <- panel$mur.mat
  mask[mask == 0] <- NA
  panel$mask <- mask
  panel$true.surface <- panel$trendmurmat + panel$fieldm
  panel$zlim <- range(panel$sz[ , 3], c(panel$true.surface))
  rp.control.put(panel$panelname, panel)

  rp.tkrplot(panel, plot1a, mururoa.colour.chart, 
        hscale = 0.2, vscale = panel$hscale * 0.7, 
        grid = "plot1a", row = 0, column = 0, background = "white")
  rp.checkbox(panel, display.options, labels = c("points", "predicted surface", "prediction s.e."),
        action = mururoa.predict, title = "Display", 
        grid = "controls2", row = 0, column = 0)
  rp.radiogroup(panel, trend.setting, c("cte", "1st", "2nd"), 
        labels = c("constant", "linear", "quadratic"),
        action = mururoa.predict, title = "Trend",
        grid = "controls2", row = 1, column = 0, sticky = "ew")
  rp.checkbox(panel, mururoa.true, mururoa.samp.redraw, 
        labels = "true surface", grid = "controls2", row = 2, column = 0)
  rp.button(panel, mururoa.colour.reset, "Reset colour scale", 
        grid = "controls2", row = 3, column = 0, sticky = "ew")
  rp.do(panel, mururoa.samp.redraw)
       
  if (!is.na(panel$file)) {
     mururoa.data <- data.frame(x = x1, y = y1, z = z)
     save(mururoa.data, file = panel$file)
     }

  panel
  }
  
mururoa.colour.chart <- function(panel) {
  par(mar = c(5, 1, 4, 2) + 0.1)
  rp.colour.chart(panel$col.palette, panel$zlim)
  panel
  }
  
mururoa.samp.redraw <- function(panel) {
   rp.control.put(panel$panelname, panel)
   rp.tkrreplot(panel, plot1)
   panel
   }
  
mururoa.colour.chart.redraw <- function(panel) {
   rp.control.put(panel$panelname, panel)
   rp.tkrreplot(panel, plot1a)
   panel
   }
  
mururoa.help <- function(panel) {
   eval(parse(text = "help(rp.mururoa)"), envir = .GlobalEnv)
   panel
   }
   
mururoa.blank <- function(panel) {
   panel$sampling.started <- FALSE
   panel
   }

mururoa.start <- function(panel) {
   panel$sampling.started <- TRUE
   panel
   }
  
   if (!require(tkrplot))      stop("The tkrplot package is not available.")
   if (!require(geoR))         stop("The geoR package is not available.")
   if (!require(RandomFields)) stop("the RandomFields package is not available.")

   if (is.list(parameters)) {
   	  nms <- names(parameters)
      for (i in 1:length(nms))
         mururoa.list[nms[i]] <- parameters[nms[i]]
      }
      
   if (is.na(hscale)) {
      if (.Platform$OS.type == "unix") hscale <- 1
      else                             hscale <- 1.4
      }
	
   ptsm        <- as.matrix(expand.grid(seq(0, 100, by = 0.5), seq(0, 50, by = 0.5)))
   trendmurmat <- apply(ptsm, 1, mururoa.list$trend.fn)
   trendmurmat <- matrix(trendmurmat, ncol = 101)
#   warn <- options()$warn
#   options(warn = -1)
#   field <- grf(nrow(ptsm), ptsm, cov.model = "matern", cov.pars = mururoa.list$cov.pars, 
#                  kappa = 4, nugget = 0, messages = FALSE)
#   options(warn = warn)

   ptsm  <- as.matrix(expand.grid(seq(0, 2, length = 201), seq(0, 1, length = 101)))
   covp  <- mururoa.list$cov.pars
   field <- GaussRF(ptsm, grid = FALSE, model = "matern", 
                     param = c(NA, covp[1], 0, covp[2] / 100, 4))
   
   panel <- rp.control("Sampling at Mururoa",
     stype = "Random", npts = 25, gsp = 10, 
     cov.pars = mururoa.list$cov.pars, nugget = mururoa.list$nugget,
     murx = mururoa.list$murx, mury = mururoa.list$mury, 
     msx.big  = mururoa.list$msx.big, msy.big  = mururoa.list$msy.big, 
     mnx.big  = mururoa.list$mnx.big, mny.big  = mururoa.list$mny.big,
     mur.mat  = mururoa.list$mur.mat, mur.tran = mururoa.list$mur.tran,
     trendmurmat = trendmurmat, trend.setting = "cte",
     kpred = array(dim = c(201, 101, 3)), kse = array(dim = c(201, 101, 3)),
     prediction.computed = rep(FALSE, 3),
     col.inside = mururoa.list$col.inside, col.outside = mururoa.list$col.outside,
     col.points = mururoa.list$col.points, col.border = mururoa.list$col.border,
     sample.taken = FALSE, hscale = hscale, col.palette = col.palette, col.se = col.se,
     transect.angle = mururoa.list$transect.angle,
     sampx = c(), sampy = c(), gx = 0, gy = 0, numt = 5, tsb = "Random", tsw = "Random",
     random.alignment = FALSE, random.alignment.old = FALSE, file = file,
     txc = c(), tyc = c(), diffm = 0, sampling.started = FALSE,
     # field = field$data, fieldm = matrix(field$data, nrow = 201),
     field = field, fieldm = matrix(field, nrow = 201),
     display.options = c("points" = TRUE, "predicted surface" = FALSE, "prediction s.e." = FALSE),
     mururoa.true = FALSE, trend.setting = "constant")
     
   rp.grid(panel, "controls1", row = 0, column = 0)
   rp.grid(panel, "plot1",     row = 0, column = 1, background = "white")
   rp.grid(panel, "plot1a",    row = 0, column = 2, background = "white")
   rp.grid(panel, "controls2", row = 0, column = 3)
   
   rp.tkrplot(panel, plot1, mururoa.draw, 
     hscale = hscale, vscale = hscale * 0.7,
     grid = "plot1", row = 0, column = 0, background = "white")
   
   rp.radiogroup(panel, stype, c("Random", "Grid", "Transect"),
     title = "Sample type", action = mururoa.points, 
     grid = "controls1", row = 0, column = 0, sticky = "ew")
   rp.radiogroup(panel, tsb, c("Random", "Systematic"),
     title = "Between transect spacing", action = mururoa.points, 
     grid = "controls1", row = 1, column = 0, sticky = "ew")
   rp.radiogroup(panel, tsw, c("Random", "Systematic"),
     title = "Within transect spacing", action = mururoa.points, 
     grid = "controls1", row = 2, column = 0, sticky = "ew")
   rp.doublebutton(panel, gx, range = c(0, 1), step = 0.05, action = mururoa.points,
     title = "Grid/Transect x-align", grid = "controls1", row = 3, column = 0)
   rp.doublebutton(panel, gy, range = c(0, 1), step = 0.05, action = mururoa.points,
     title = "Grid/Transect y-align", grid = "controls1", row = 4, column = 0)
   rp.checkbox(panel, random.alignment, mururoa.points, 
     title = "Random x,y-alignment", grid = "controls1", row = 5, column = 0)
   rp.slider(panel, npts, 10, 100, log = TRUE, action = mururoa.points,
     title = "Number of points", grid = "controls1", row = 6, column = 0, sticky = "ew")
   rp.button(panel, mururoa.points, "New sampling positions",
     grid = "controls1", row = 7, column = 0)
   rp.do(panel, mururoa.blank)
   rp.do(panel, mururoa.points)
   rp.do(panel, mururoa.start)
   
   rp.text(panel, " \n \n \n \n \n ", grid = "controls2", row = 0, column = 0)
   rp.button(panel, mururoa.samp, "Take sample", grid = "controls2", row = 2, column = 0)
   }

