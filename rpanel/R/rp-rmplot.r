#   An rpanel function for rmplot

rp.rmplot <- function(y, id = NA, timept = NA, fac = NA, type = "all",
                          xlab = NA, ylab = NA, xlabels = NA, add = FALSE,
                          lwd = 1, col = NA, lty = NA, panel = TRUE, 
                          panel.plot = TRUE, hscale = NA, vscale = hscale, ...) {

   if (is.na(hscale)) {
      if (.Platform$OS.type == "unix") hscale <- 1
      else                             hscale <- 1.4
      }
   if (is.na(vscale)) 
      vscale <- hscale

rmplot <- function(y, id = NA, timept = NA, fac = NA, type = "all", 
                   add = FALSE, xlab = NA, ylab = NA, xlabels = NA, 
                   lwd = 1, col = NA, lty = NA, ...) {

   if (!is.na(fac) && !is.factor(fac)) stop("fac must be a factor.")

   if (is.matrix(y) | is.data.frame(y)) 
      result <- rmplot.table(y, fac, type, timept, add = add, xlab = xlab, ylab = ylab,
                          lwd = lwd, col = col, lty = lty, xlabels = xlabels, ...)
   
   else {

      if (all(is.na(fac)))    fac <- factor(rep(1, length(y)))
      if (!is.factor(id))     stop("id  must be a factor.")
      if (all(is.na(timept))) stop("timept must be specified.")

      n     <- length(id)
      ids   <- unique(id)
      times <- sort(unique(timept))
      nid   <- length(ids)
      nt    <- length(times)

      ymat   <- matrix(NA, nrow = nid, ncol = nt)
      id.fac <- factor(rep(NA, nid), levels = levels(fac))
      for (i in (1:n)) {
        rowind               <- (1:nid)[id[i] == ids]
        colind               <- (1:nt)[timept[i] == times]
        ymat[rowind, colind] <- y[i]
        if (!is.na(id.fac[rowind]) & (fac[i] != id.fac[rowind]))
           stop("Mismatch between id and fac.")
        id.fac[rowind] <- fac[i]
        }

      result <- rmplot.table(ymat, times = times, fac = id.fac, type = type, 
                      xlab = xlab, ylab = ylab, xlabels = xlabels, ...)
      # result <- list(id = ids, ymat = ymat, times = times, fac = id.fac)
      }

   invisible(result)

   }

rmplot.table <- function(y, fac = NA, type = "all", times = (1:ncol(y)),
          xlab = xlab, ylab = ylab, add = FALSE, 
          lwd = 1, col = NA, lty = NA, xlabels = NA, ...) {

   nt      <- ncol(y)
   if (all(is.na(fac))) fac <- factor(rep(1, nrow(y)))
   nlevels <- length(levels(fac))   
   nfac    <- nlevels
   if (all(is.na(times))) times <- 1:ncol(y)
   if (all(is.na(fac))) {
      fac <- factor(rep(1, nrow(y)))
      if (is.na(col)) col <- 1
      if (is.na(lty)) lty <- 1
      }
   else {
      if (!is.factor(fac)) fac <- factor(fac)
      if (all(is.na(col))) col <- 1:nlevels
      if (all(is.na(lty))) lty <- 1:nlevels
      }
   if (all(is.na(xlabels))) xlabels <- as.character(times)

#--------------------------------------------------------------     
#       Plot individual profiles
#--------------------------------------------------------------     

   if(type == "all") {
      if (!add) plot(range(times), range(y, na.rm = TRUE),
                     type = "n", xlab = xlab, ylab = ylab, axes = FALSE, ...)
      for(i in 1:nrow(y)) {
         if (!is.na(fac[i])) {
            # ifac <- (1:nlevels)[fac[i] == levels(fac)]
            ifac <- as.numeric(fac[i])
            if (!all(is.na(y[i,]))) {
            lines(times, y[i,  ], col = col[ifac], lty = lty[ifac], lwd = lwd)
            #   Plot isolated points
            yi <- c(NA, y[i, ], NA)
            ind <- (diff(diff(is.na(yi))) == 2)
            if (any(ind)) points(times[ind], y[i, ind], col = col[ifac])
            }
         }
      }
   }

#--------------------------------------------------------------
#       Calculate means and standard errors, if necessary
#--------------------------------------------------------------     

   if(!(type == "all")) {
      m  <- matrix(NA, nfac, nt)
      se <- matrix(NA, nfac, nt)
      for(i in 1:nfac) {
         # yi      <- y[fac == facs[i], ]
         yi      <- y[fac == levels(fac)[i], ]
         fun     <- function(y) length(y[!is.na(y)])
         ni      <- apply(yi, 2, fun)
         ni      <- max(ni, 1)
         m[i,  ] <- apply(yi, 2, mean, na.rm = TRUE)
         se[i, ] <- apply(yi, 2,   sd, na.rm = TRUE) / sqrt(ni)
         }
      }
   shift <- ((1:nfac) - 0.5 - nfac / 2) * diff(range(times)) / 200

#--------------------------------------------------------------
#       Plot mean profiles
#--------------------------------------------------------------     
        
   if(type == "mean") {
      if (!add) plot(range(times),
                     # range(m),
                     range(m, m - 2 * se, m + 2 * se, na.rm = TRUE),
                     type = "n", xlab = xlab, ylab = ylab, axes = FALSE, ...)
      for(i in 1:nfac) {
         # ifac <- (1:nlevels)[facs[i] == levels(fac)]
         ifac <- i
         lines(times + shift[i], m[i,  ], col = col[ifac], lty = lty[ifac])
         if (nfac > 1) text(times + shift[i], m[i,  ], ifac, col = col[ifac])
         }
      }
     
#--------------------------------------------------------------     
#       Plot mean profiles + 2 s.e. bars, possibly also with band
#--------------------------------------------------------------  
   
   if(type == "band") {
      if (nlevels < 2)
         cat("Band cannot be drawn with only one group.\n")
      else if (nlevels > 2)
         cat("Band cannot be drawn with more than two groups.\n")
      else if (any(is.na(m + se)))
         cat("Band cannot be drawn with missing means or se's.\n")
      if (nlevels < 2 | nlevels > 2 | any(is.na(m + se)))
         type <- "mean+bar"
   }

   if(type == "band") {
      if (!add) plot(range(times),
                     range(m, m - 2 * se, m + 2 * se, na.rm = TRUE),
                     type = "n", xlab = xlab, ylab = ylab, axes = FALSE, ...)
      hi <- rep(0, nt)
      lo <- rep(0, nt)
      for(j in 1:nt) {
         av    <- mean(c(m[1, j], m[2, j]))
         width <- sqrt(se[1, j]^2 + se[2, j]^2)
         hi[j] <- av + width
         lo[j] <- av - width
         }
      polygon(c(times, rev(times)), c(hi, rev(lo)), 
               density = -1, border = NA, col = "lightgreen")
      }

   if ((type == "band") | (type == "mean+bar")) {
      if (type == "mean+bar" & !add)
         plot(range(times), 
              range(m, m - 2 * se, m + 2 * se, na.rm = TRUE), 
              type = "n", xlab = xlab, ylab =  ylab, axes = FALSE, ...)
      for (i in 1:nfac) {
         # ifac <- (1:nlevels)[facs[i] == levels(fac)]
         ifac <- i
         lines(times + shift[i], m[i,  ], col = col[ifac], lty = lty[ifac])
         if (nfac > 1) text( times + shift[i], m[i,  ], i, 
                            col = col[ifac])
         if (type == "mean+bar")
            segments(times + shift[i], m[i, ] - 2 * se[i, ],
                     times + shift[i], m[i, ] + 2 * se[i, ],
                        col = col[ifac], lty = lty[ifac])
         }
      }
#--------------------------------------------------------------     
#       Scan individual cases
#--------------------------------------------------------------        

   if (type == "scan") {
      choices <- c("Next", "Previous", "Select", "Complete", "Stop")
      par(col = 1, lty = 1)
      plot(range(times), range(y, na.rm = TRUE), type = "n",
                xlab = xlab, ylab = ylab, ...)
      rmplot.table(y = y, fac = fac, type = "mean", times = times,
                xlab = xlab, ylab = ylab, add = TRUE, axes = FALSE, ...)

      plot.profile <- function(i, times, y, fac, erase = FALSE) {
         if (!is.na(fac[i])) {
            ifac <- (1:nlevels)[fac[i] == levels(fac)]
            if (!all(is.na(y[i,]))) {
           if (erase) lines(times, y[i,  ], col = "white", lty = ifac)
           else       lines(times, y[i,  ], col = ifac,    lty = ifac)
           }
            }
         if (erase) title(paste("Case", i), cex = 0.5, col.main = "white")
         else       title(paste("Case", i), cex = 0.5)
      }

      i <- 1
      direction <- 1
      plot.profile(i, times, y, fac)
      choice <- menu(choices)
      while (choice < 5) {
         plot.profile(i, times, y, fac, erase = TRUE)
         if (choice == 2) direction <- -1 else direction <- 1
         i <- i + direction
         if(choice == 3) i <- scan(what = integer())
         if(i < 1)       i <- nrow(y)
         if(i > nrow(y)) i <- 1
         plot.profile(i, times, y, fac)
         if(choice < 4 | i == nrow(y)) choice <- menu(choices)
         }
      }

   if (type %in% c("all", "band", "mean", "mean+bar")) {
      axis(2)
      if (length(xlabels) < 100)
         axis(1, at = times, labels = xlabels)
      else
         axis(1)
      box()
      }

   results <- list(ymat = y, fac = fac, times = times)
   if (type %in% c("band", "mean", "mean+bar"))
      results$mean <- m
   if (type %in% c("band", "mean+bar"))
      results$se <- se
   
   invisible(results)
   }

#-------------------------------------------
#     Define the action function
#-------------------------------------------

rmplot.draw <- function(panel) {
   with(panel, {
      if (fac.showing) fac1 <- fac else fac1 <- NA
      if (data.range)
         rmplot(y, timept = timept, fac = fac1, type = type,
                    xlab = xlab, ylab = ylab, 
                    lwd = 1, col = col, lty = lty,
                    ylim = range(y, na.rm = TRUE), xlabels = xlabels)
      else
         rmplot(y, timept = timept, fac = fac1, type = type, 
                lwd = 1, col = col, lty = lty,
                xlab = xlab, ylab = ylab, xlabels = xlabels)
      if (case.showing) {
         rmplot(matrix(y[case,], nrow = 1), timept = timept, fac = fac1[case], type = "all",
                    lwd = 3, col = col, lty = lty,
                    xlab = xlab, ylab = ylab, xlabels = xlabels, add = TRUE)
         title(paste("case =", case))
      }
   })
   panel
}

   rmplot.redraw <- function(panel) {
      rp.tkrreplot(panel, plot)
      panel
   }

   if (all(is.na(fac))) {
      if (is.matrix(y) | is.data.frame(y))
         fac <- factor(rep(1, nrow(y)))
      else
         fac <- factor(rep(1, length(y)))
   }
   else {
      if (!is.factor(fac)) fac <- factor(fac)
      fac.temp <- fac[!is.na(fac)]
      if (length(levels(fac.temp)) != length(unique(fac.temp)))
         stop("data must be available for all factor levels.")
      }
   if (is.na(xlab))  
      xlab <- if (is.vector(y)) deparse(substitute(timept)) else "Time"
   if (is.na(ylab))
      ylab <- deparse(substitute(y))
   timept <- as.numeric(timept)

   if (panel) {
      result <- rmplot(y, id, timept, fac, type = "none")
      y      <- result$ymat
      fac    <- result$fac
      times  <- result$times
      rmplot.panel <- rp.control("Repeated measurements", y = y, 
                      timept = times, fac = fac, col = col, lty = lty, 
                      xlab = xlab, ylab = ylab, xlabels = xlabels,
                      fac.showing = length(levels(fac)) > 1, type = "all",
                      case.showing = FALSE, case = 1, data.range = FALSE)
      if (panel.plot) {
         rp.tkrplot(rmplot.panel, plot, rmplot.draw, pos = "right",
                hscale = hscale, vscale = vscale, background = "white")
         action.fn <- rmplot.redraw
         }
      else
         action.fn <- rmplot.draw
      plot.options <- c("all", "mean", "mean+bar")
      if (length(levels(fac)) == 2) plot.options <- c("all", "mean", "mean+bar", "band")
      rp.radiogroup(rmplot.panel, type, plot.options, title = "Display", action = action.fn)
      if (length(levels(fac)) > 1)
         rp.checkbox(rmplot.panel, fac.showing, action.fn, title = "Show groups")
      rp.checkbox(rmplot.panel, case.showing, action.fn, title = "Show cases")
      rp.doublebutton(rmplot.panel, case, 1, range = c(1, nrow(y)),
                      title = "Case", initval = 1, action = action.fn)
      rp.checkbox(rmplot.panel, data.range, action.fn,
                      title = "Use data range")
      rp.do(rmplot.panel, action.fn)
      }
   else {
      # fac.showing  <- (!all(is.na(fac)))
      # data.range   <- (type == "all")
      # rmplot.panel <- list(y = y, timept = times, fac = fac, lwd = lwd,
      #               col = col, lty = lty, xlab = xlab, ylab = ylab,
      #               fac.showing = fac.showing, data.range = data.range,
      #               case.showing = FALSE, xlabels = xlabels)   
      # rmplot.draw(rmplot.panel)
      rmplot(y, id, timept, fac, type, add,
                xlab, ylab, xlabels, lwd = 1, col = col, lty = lty, ...)
      }

   invisible()
}
