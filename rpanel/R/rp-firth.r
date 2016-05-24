# What scale of measurement is appropriate?
# Consider trend, cov.pars and nugget settings.
# Remove redundant items from firth.list (setRF).

rp.firth <- function(hscale = NA, col.palette = rev(heat.colors(40)), col.se = "blue",
                 file = NA, parameters = NA) {

firth.points <- function(panel) {
	
   if (panel$random.alignment) {
      panel$gx <- runif(1)
      panel$gy <- runif(1)
      }
   panel$random.alignment.old <- panel$random.alignment

  if (panel$stype=="Random") {
     wh          <- which(panel$strat > 0)
     sind        <- sample(wh, panel$npts) - 1
     panel$sampx <- (sind %% 201) 
     panel$sampy <- (sind %/% 201)
     }

  if (panel$stype == "Systematic") {
     gsp <- sqrt(6041 / (panel$npts))	# find (average) grid spacing
     txc <- round(seq(panel$gx * gsp, 200, by = gsp))
     tyc <- round(seq(panel$gy * gsp,  60, by = gsp))
     if (panel$sp.type == "Random x") {
        grx <- c()
        gry <- c()
        ycoords <- tyc
        nt <- length(ycoords)
        for (i in 1:nt) {
          xcoords <- which(panel$strat[ , ycoords[i] + 1] > 0) - 1
          Ly <- length(xcoords)   
     	  # find which pts on transect are in Firth using strat matrix as look-up
          if (Ly > 0) grx <- c(grx, sample(xcoords, size = round(Ly / gsp)))
	      # take suitably sized sample
          if (Ly > 0) gry <- c(gry, rep(ycoords[i], round(Ly / gsp)))
          }
       panel$sampx <- grx
       panel$sampy <- gry
       }
    else {
       gr <- expand.grid(txc, tyc)
       a <- which(panel$strat[txc + 1, tyc + 1] > 0)
       gr <- gr[a, ] 
       panel$sampx <- gr[ , 1]
       panel$sampy <- gr[ , 2]
       }
     panel$txc <- txc
     panel$tyc <- tyc
     }

  if (panel$stype == "Stratified") {
     alloc <- rep(0, 4)
     if (panel$str.type == "Equal")
        alloc <- rep(round(panel$npts / 4), 4)
     else
       alloc <- round(panel$npts * panel$str.size / 100)
     sind1  <- sample(which(panel$strat == 1), alloc[1])
     sind2  <- sample(which(panel$strat == 2), alloc[2])
     sind3  <- sample(which(panel$strat == 3), alloc[3])
     sind4  <- sample(which(panel$strat == 4), alloc[4])
     sind   <- c(sind1, sind2, sind3, sind4)
     panel$sampx <- (sind %% 201) - 1
     panel$sampy <- (sind %/% 201)
     }

  if (panel$sampling.started) {
     rp.control.put(panel$panelname, panel)
     rp.tkrreplot(panel, plot1)
  }
  
  panel
  }

firth.predict <- function(panel) {
	
      ind <- which(panel$trend.setting == c("cte", "1st", "2nd", "stratum"))
      
      if (!panel$prediction.computed[ind]) {
         sz   <- panel$sz
         fg   <- as.geodata(sz)
         vg2  <- variog(fg, trend = "cte", max.dist = 120, messages = FALSE)
         shh  <- output.control(messages = FALSE)
      	 if (ind == 4) {
            a    <- 1 + sz[ , 1] + 201 * sz[ , 2]
            sgd  <- factor(panel$strat[a])
            fg   <- as.geodata(data.frame(sz[,1], sz[,2], sz[,3], sgd), 
                      covar.col = 4, covar.names = "stratum")
      	 	fit  <- likfit(fg, ini.cov.pars = c(max(vg2$v), 20), trend = ~stratum)
            pg2  <- expand.grid(0:200, 0:60)
            sgp  <- c(panel$strat)
            # indg <- (sgp %in% as.numeric(levels(sgd)))
            # pg2  <- pg2[indg, ]
            # sgp  <- sgp[indg]
            sgp[sgp == 0] <- 1
            sgp  <- factor(sgp)
            KC   <- krige.control(obj.model = fit, trend.d = ~stratum, trend.l = ~sgp)
            kc   <- krige.conv(fg, locations = pg2, krige = KC, output = shh)
      	    }
      	 else {
            vfit <- variofit(vg2, ini.cov.pars = c(max(vg2$v), 20), nugget = 10, cov.model = "matern",
                             kappa = 4, messages = FALSE)
            xx   <- 0:120
            yy   <- 0.04 + 0.25 * (xx / 20 - 0.5 * (xx / 30)^3)
            yy[32:121] <- 0.29
            pred.grid <- expand.grid(0.5 + 0:199, 0.5 + 0:59) 	# offset from samples
            KC   <- krige.control(obj.model = vfit, trend.d = panel$trend.setting, 
                           trend.l = panel$trend.setting)
            kc   <- krige.conv(fg, locations = pred.grid, krige = KC, output = shh)
            kcm  <- matrix(kc$predict, nrow = 200)
            pg2  <- expand.grid(0:200, 0:60) # use 'real' grid to find RSS
            KC   <- krige.control(obj.model = vfit, trend.d = panel$trend.setting, 
                         trend.l = panel$trend.setting)
            kc   <- try(krige.conv(fg, locations = pg2, krige = KC , output = shh), silent = TRUE)
            if (class(kc) == "try-error" & panel$trend.setting == "cte") {
               rp.messagebox("There are numerical problems in producing predictions with this model.  Try fitting a linear or quadratic trend function, or a stratum effect.")
               return(panel)
            }
      	 }
         panel$kpred[ , , ind] <- matrix(kc$predict, ncol = 61)
         panel$kse[ , , ind]   <- sqrt(matrix(kc$krige.var, ncol = 61))
         panel$prediction.computed[ind] <- TRUE
      }
      panel$zlim  <- range(panel$zlim,  panel$kpred[ , ,ind] * panel$mask,
                       panel$true.surface, na.rm = TRUE)
#      panel$zlim  <- range(panel$sz[ , 3],
#                           panel$kpred[,,1] * panel$mask, panel$kpred[,,2] * panel$mask,
#                           panel$kpred[,,3] * panel$mask, panel$kpred[,,4] * panel$mask,
#                           panel$true.surface, na.rm = TRUE)
      rp.control.put(panel$panelname, panel)
      rp.do(panel, firth.colour.chart.redraw)
      rp.do(panel, firth.samp.redraw)
      
      panel
      }

firth.draw <- function(panel) {
	
   ind <- which(panel$trend.setting == c("cte", "1st", "2nd", "stratum"))
   
   with(panel, {
   	
   plot(c(10, 10), type = "n", asp = 1, xaxs = "i", xlim = c(0, 200), ylim = c(0, 60),
        xlab = "easting", ylab = "northing", xaxp = c(0, 200, 10))
   polygon(rosx, rosy, col = hsv(0.08, 0.3, 0.9))
   polygon(p1xa, p1ya, col = hsv(0.16, 0.3, 1))
   polygon(p3xa, p3ya, col = hsv(0.12, 0.3, 1))
   polygon(p2xa, p2ya, col = hsv(0.12, 0.3, 1))
   polygon(p4x,  p4y,  col = hsv(0.08, 0.3, 0.8))
   polygon(p5xa, p5ya, col = hsv(0.16, 0.3, 1))
   polygon(p6xa, p6ya, col = hsv(0.12, 0.3, 1))
   polygon(p7xa, p7ya, col = hsv(0.16, 0.3, 1))
   polygon(p8xa, p8ya, col = hsv(0.12, 0.3, 1))
   
   if (sample.taken) {
   	
     if (display.options["predicted surface"])
        image(seq(0, 200, len = 201), seq(0, 60, len = 61), kpred[,,ind], 
               add = TRUE, col = col.palette, zlim = zlim)

     if (display.options["prediction s.e."])
        contour(seq(0, 200, len = 201), seq(0, 60, len = 61), kse[,,ind], add = TRUE,
               levels = pretty(range((kse[,,ind] * mask)[kse[,,ind] > 0], na.rm = TRUE), 10), col = col.se)
     
     if (firth.true)
        image(seq(0, 200, len = 201), seq(0, 60, len = 61), true.surface,
               add = TRUE, col = col.palette, zlim = zlim)
        
     if (display.options["predicted surface"] | firth.true) {
        polygon(rosx, rosy)
        polygon(p1xa, p1ya)
        polygon(p2xa, p2ya)
        polygon(p3xa, p3ya)
        polygon(p4x,  p4y)
        polygon(p5xa, p5ya)
        polygon(p6xa, p6ya)
        polygon(p7xa, p7ya)
        polygon(p8xa, p8ya)
        }

     if (display.options["points"] & !firth.true) {
        brks <- seq(zlim[1], zlim[2], length = length(col.palette) + 1)
        clr  <- cut(sz[ , 3], brks, labels = FALSE, include.lowest = TRUE, right = FALSE)
        points(sz[ , 1], sz[ , 2], col = col.palette[clr], pch = 16)
        points(sz[ , 1], sz[ , 2])
        }
       
      }
      
   else if (panel$sampling.started) {
   
      sx <- round(as.numeric(sampx))
      sy <- round(as.numeric(sampy))
      if (stype == "Systematic") {
         x0 <- panel$txc
         y0 <- panel$tyc
         segments(x0,  0,  x0, 60, lty = 2)
         segments( 0, y0, 200, y0, lty = 2) 
         }
      points(sx, sy, pch = 19)
      mtext(paste("Number of points = ",length(sx)), line = 3)
      a <- 1 + sx + 201 * sy
      mtext(bquote(paste("Allocation: ", n[1]==.(sum(strat[a]==1)), ", ",
        n[2]==.(sum(strat[a]==2)), ", ", n[3]==.(sum(strat[a]==3)), ", ",
        n[4]==.(sum(strat[a]==4)))), line = 2)
      mtext(bquote(paste("[Proportional:  ", n[1]==.(round(length(sx)*str.size[1]/100)),
        ", ", n[2]==.(round(length(sx)*str.size[2]/100)),
        ", ", n[3]==.(round(length(sx)*str.size[3]/100)),
        ", ", n[4]==.(round(length(sx)*str.size[4]/100)), "]")), line = 1)
      }
      
   usr <- par()$usr
   polygon(c(rsx[1:123], rsx[123], rsx[1]), c(rsy[1:123], rep(usr[3], 2)), 
              col = col.outside, density = -1, border = NA)
   polygon(c(rnx[1:138], rnx[138], rnx[1]), c(rny[1:138], rep(usr[4], 2)), 
              col = col.outside, density = -1, border = NA)
   box()
   polygon(rosx, rosy)

   })
   
   panel
   }

firth.colour.reset <- function(panel) {
  ind <- which(panel$trend.setting == c("cte", "1st", "2nd", "stratum"))
  panel$zlim  <- range(panel$sz[ , 3],  panel$kpred[ , ,ind] * panel$mask, 
                       panel$true.surface, na.rm = TRUE)
  rp.control.put(panel$panelname, panel)
  rp.do(panel, firth.colour.chart.redraw)
  rp.do(panel, firth.samp.redraw)
  panel
  }
  
firth.samp <- function(panel) {
	
   if (!panel$sampling.started) return(panel)
	
   panel$sample.taken <- TRUE
   x   <- round(panel$sampx)
   y   <- round(panel$sampy)
   x1  <- x[which(panel$strat[1 + x + 201 * y] > 0)]
   y1  <- y[which(panel$strat[1 + x + 201 * y] > 0)]
   dup <- c() 				# remove duplicate points
   xn  <- length(x1)
   for (i in 1:(xn - 1)) {
     for (j in seq(i + 1, xn)) {
       if (x1[i] == x1[j] & y1[i] == y1[j]) dup <- c(dup, j)
       }
     }
   if (length(dup) > 0) x1 <- x1[-dup]
   if (length(dup) > 0) y1 <- y1[-dup]
   #  if (length(x1)>50) {
   #    x1 <- x1[1:50] 			 # [restrict to 50 sampling points]
   #    y1 <- y1[1:50] }
   a <- 1 + x1 + 201 * y1		       # find right place in trend etc. vectors
   warn <- options()$warn
   options(warn = -1)
   nug <- grf(nrow(panel$pts), panel$pts, cov.model = "pure.nugget",	# add nugget effect
              nugget = 0, cov.pars = c(panel$nugget, 0), messages = FALSE)	# 'sampling error'
   options(warn = warn)
   # z <- panel$field[a] + panel$trend[a] + nug$data[a] + panel$strat.sm[a] / 2
   z <- panel$field[a] + panel$trend[a] + nug$data[a] + panel$strat.effect[c(panel$strat)[a]]
   panel$sz        <- cbind(x1, y1, z) 		# return 3 columns with x,y,z values
   mask            <- panel$strat
   mask[mask >  0] <- 1
   mask[mask == 0] <- NA
   panel$mask      <- mask
   panel$true.surface <- matrix(panel$trend, ncol = 61) + matrix(panel$field, nrow = 201) +
                             apply(panel$strat, 1:2, 
                                     function(x) if (x > 0) panel$strat.effect[x] else 0)
   panel$zlim <- range(panel$sz[ , 3], c(panel$true.surface))
   rp.control.put(panel$panelname, panel)
   
   rp.tkrplot(panel, plot1a, firth.colour.chart, 
         hscale = 0.2, vscale = panel$hscale * 0.7, 
         grid = "plot1a", row = 0, column = 0, background = "white")
   rp.checkbox(panel, display.options,
         labels = c("points", "predicted surface", "prediction s.e."),
         action = firth.predict, title = "Display", 
         grid = "controls2", row = 0, column = 0, sticky = "ew")
   rp.radiogroup(panel, trend.setting, c("cte", "1st", "2nd", "stratum"), 
         labels = c("constant", "linear", "quadratic", "stratum"),
         action = firth.predict, title = "Trend", 
         grid = "controls2", row = 1, column = 0, sticky = "ew")
   rp.checkbox(panel, firth.true, firth.samp.redraw,
         labels = "true surface", grid = "controls2", row = 2, column = 0)
   rp.button(panel, firth.colour.reset, "Reset colour scale", 
         grid = "controls2", row = 3, column = 0, sticky = "ew")
   rp.control.put(panel$panelname, panel)
   rp.do(panel, firth.samp.redraw)
       
   if (!is.na(panel$file)) {
      firth.data <- data.frame(x = x1, y = y1, z = z)
      save(firth.data, file = panel$file)
      }

   panel
   }
   
   firth.colour.chart <- function(panel) {
      par(mar = c(5, 1, 4, 2) + 0.1)
      rp.colour.chart(panel$col.palette, panel$zlim)
      panel
      }
  
   firth.samp.redraw <- function(panel) {
	  rp.tkrreplot(panel, plot1)
      panel
      }
  
   firth.colour.chart.redraw <- function(panel) {
      rp.tkrreplot(panel, plot1a)
      panel
      }
  
   firth.blank <- function(panel) {
      panel$sampling.started <- FALSE
      panel
      }

   firth.start <- function(panel) {
      panel$sampling.started <- TRUE
      panel
      }
  
   if (!require(tkrplot))      stop("The tkrplot package is not available.")
   if (!require(geoR))         stop("The geoR package is not available.")
   if (!require(RandomFields)) stop("the RandomFields package is not available.")
   if (is.list(parameters)) {
   	  nms <- names(parameters)
      for (i in 1:length(nms))
         firth.list[nms[i]] <- parameters[nms[i]]
      }
      
   if (is.na(hscale)) {
      if (.Platform$OS.type == "unix") hscale <- 1
      else                             hscale <- 1.4
      }

#   warn <- options()$warn
#   options(warn = -1)
#   field <- grf(nrow(pts), grid = pts, cov.model = "matern", cov.pars = firth.list$cov.pars, 
#                  kappa = 4, nugget = 0, messages = FALSE)
#   options(warn = warn)
   pts   <- as.matrix(expand.grid(seq(0, 2, length = 201), seq(0, 0.6, length = 61)))
   covp  <- firth.list$cov.pars
   field <- GaussRF(pts, grid = FALSE, model = "matern", 
                     param = c(NA, covp[1], 0, covp[2] / 100, 4))
   pts   <- firth.list$pts
   trend <- apply(pts, 1, firth.list$trend.fn)

   panel <- rp.control("Sampling in a firth",
                stype = "Random", npts = 25, gsp = 10,
                sampx = c(), sampy = c(), gx = 0, gy = 0, tbw = 0, file = file,
                cov.pars = firth.list$cov.pars, nugget = firth.list$nugget,
                strat.effect = firth.list$strat.effect,
                display.options = c("points" = TRUE, "predicted surface" = FALSE, 
                                     "prediction s.e." = FALSE),
                trend.setting = "cte", firth.true = FALSE,
                sample.taken = FALSE, hscale = hscale, col.palette = col.palette, col.se = col.se,
                str.type = "Proportional", sp.type = "Equal", pts = pts,
                txc = round(seq(0, 200, by = sqrt(6041 / 25))),
                tyc = round(seq(0,  60, by = sqrt(6041 / 25))),
                diffm = 0, kvf = 0, 
                kpred = array(dim = c(201, 61, 4)), kse = array(dim = c(201, 61, 4)),
                prediction.computed = rep(FALSE, 4), sampling.started = FALSE,
                random.alignment = FALSE, random.alignment.old = FALSE, 
                field = field, trend = trend, col.outside = "palegreen4",
                rnx.big = firth.list$rnx.big, rny.big = firth.list$rny.big,
                rsx.big = firth.list$rsx.big, rsy.big = firth.list$rsy.big,
                rosx = firth.list$rosx, rosy = firth.list$rosy,
                rnx  = firth.list$rnx,  rny  = firth.list$rny,
                rsx  = firth.list$rsx,  rsy  = firth.list$rsy,
                p1xa = firth.list$p1xa, p1ya = firth.list$p1ya,
                p2xa = firth.list$p2xa, p2ya = firth.list$p2ya,
                p3xa = firth.list$p3xa, p3ya = firth.list$p3ya,
                p4x  = firth.list$p4x,  p4y  = firth.list$p4y,
                p5xa = firth.list$p5xa, p5ya = firth.list$p5ya,
                p6xa = firth.list$p6xa, p6ya = firth.list$p6ya,
                p7xa = firth.list$p7xa, p7ya = firth.list$p7ya,
                p8xa = firth.list$p8xa, p8ya = firth.list$p8ya,
                str.size = firth.list$str.size, strat   = firth.list$strat,
                strat.sm = firth.list$strat.sm, strat2  = firth.list$strat2,
                stratk   = firth.list$stratk,   stratk2 = firth.list$stratk2)
   rp.grid(panel, "controls1", row = 0, column = 0)
   rp.grid(panel, "controls2", row = 0, column = 3)
   rp.grid(panel, "plot1",     row = 0, column = 1, background = "white")
   rp.grid(panel, "plot1a",    row = 0, column = 2, background = "white")

   rp.tkrplot(panel, plot1, firth.draw, 
      hscale = hscale, vscale = hscale * 0.7,
      grid = "plot1", row = 0, column = 0, background = "white")

   rp.radiogroup(panel, stype, c("Random","Systematic","Stratified"), 
      title = "Sample type", action = firth.points,
      grid = "controls1", row = 0, column = 0, sticky = "ew")
   rp.radiogroup(panel, sp.type, c("Equal","Random x"), 
      title = "Systematic spacing type", action = firth.points, 
      grid = "controls1", row = 1, column = 0, sticky = "ew")
   rp.radiogroup(panel, str.type, c("Proportional","Equal"), 
      title = "Allocation type (stratified)", action = firth.points, 
      grid = "controls1", row = 2, column = 0, sticky = "ew")
   rp.button(panel, firth.points, title = "New Sampling Points", 
      grid = "controls1", row = 3, column = 0, sticky = "ew")
   rp.slider(panel, npts, 10, 100, action = firth.points,
      title = "Number of points", log = FALSE, 
      grid = "controls1", row = 4, column = 0, sticky = "ew")
   rp.doublebutton(panel, gx, 0.05, range = c(0, 1), action = firth.points, title = "Grid x-align", 
      grid = "controls1", row = 5, column = 0)
   rp.doublebutton(panel, gy, 0.05, range = c(0, 1), action = firth.points, title = "Grid y-align", 
      grid = "controls1", row = 6, column = 0)
   rp.checkbox(panel, random.alignment, firth.points, 
     labels = "Random x,y-alignment", grid = "controls1", row = 7, column = 0)
   rp.do(panel, firth.blank)
   rp.do(panel, firth.points)
   rp.do(panel, firth.start)

   rp.text(panel, " \n \n \n \n \n ", grid = "controls2", row = 0, column = 0)
   rp.button(panel, firth.samp, "Take sample", grid = "controls2", row = 2, column = 0)
   }
