rp.likelihood <- function(loglik.fn, data, theta.low, theta.high,
                            form = "log-likelihood", hscale = NA, vscale = hscale) {

   if (is.na(hscale)) {
      if (.Platform$OS.type == "unix") hscale <- 1
      else                             hscale <- 1.4
      }
   if (is.na(vscale)) 
      vscale <- hscale

#   One-parameter likelihood plots

rp.loglik1 <- function(loglik.text, data, theta.low, theta.high,
                            form = "log-likelihood") {

   if (!require(tkrplot)) stop("the tkrplot library is required.")

   theta.range <- c(theta.low, theta.high)

   loglik1.plot <- function(panel) {
   	
      with(panel, {
         if (form == "likelihood") {
   	        fun.text <- paste("function(theta, data) exp(", panel$loglik.text, ")")
            if (display["relative scale"]) ylab <- "Relative likelihood"
               else                        ylab <- "Likelihood"
            }
         else {
   	        fun.text <- paste("function(theta, data)", panel$loglik.text)
            if (display["relative scale"]) ylab <- "Relative log-likelihood"
               else                        ylab <- "Log-likelihood"
            }

         warn <- options()$warn
         options(warn = -1)
         
         fnscale <- apply(matrix(mean(theta.range)), 1, eval(parse(text = fun.text)), data = data)
         mloglik <- optim(mean(theta.range), eval(parse(text = fun.text)), data = data, 
                         method = "BFGS", hessian = TRUE,
                         control = list(fnscale = -abs(fnscale)))
         convergence <- (mloglik$convergence == 0)
         if (!convergence) rp.messagebox("The mle could not be located.")
            else mle <- mloglik$par
         options(warn = warn)

         theta.grid <- seq(theta.range[1], theta.range[2], length = ngrid)
         ll         <- apply(as.matrix(theta.grid), 1, eval(parse(text = fun.text)), data = data)
         ll[is.nan(ll)] <- NA
         if (display["relative scale"]) {
            if (form == "likelihood") ll <- ll / mloglik$value
               else                   ll <- ll - mloglik$value
            }
         rx        <- range(theta.grid)
         ry        <- range(ll, na.rm = TRUE)
         if (!is.na(axis.active))
            rx <- range(rx, new.range)
         plot(ry ~ rx, type = "n", xaxs = "i", xlab = "theta", ylab = ylab)
         lines(theta.grid, ll, col = "blue", lwd = 2)
         if (!is.na(axis.active))
            abline(v = new.range, lty = 2)
      
         if (display["mle"] & convergence) {
            if (display["relative scale"]) {
               if (form == "likelihood") ymax <- 1
                  else                   ymax <- 0
               }
            else 
               ymax <- mloglik$value
            segments(mle, par()$usr[3], mle, ymax, lty = 2)
            points(mle, ymax, pch = 16, col = "red")
            title(paste("mle =", signif(mle, 5)), line = 2, cex.main = 1)
            }

         if (display["threshold proportion"] & convergence) {
         	 if (form == "likelihood") {
         	 	if (display["relative scale"]) thresh <- prop
         	 	   else                        thresh <- prop * mloglik$value
         	    }
         	 else {
         	 	if (display["relative scale"]) thresh <- log(prop)
         	 	   else                        thresh <- log(prop) + mloglik$value 	
         	    }
         	 abline(h = thresh, col = "brown")
            }
      
         if (display["quadratic approximation"] & 
                  (form == "log-likelihood") & convergence) {
            quad.fn <- function(x, mloglik) mloglik$value + 0.5 * c(mloglik$hessian) * (x - mle)^2
            quad    <- apply(as.matrix(theta.grid), 1, quad.fn, mloglik = mloglik)
            if (display["relative scale"]) quad <- quad - mloglik$value
            lines(theta.grid, quad, lty = 2, col = "red", lwd = 2)
            }
      
         if (display["ci"] & 
                  (form == "log-likelihood") & convergence) {
            if (display["relative scale"]) ref <- -1.92
               else                        ref <- mloglik$value - 1.92
            abline(h = ref, lty = 2, col = "green")
         
            ind <- min(which(ll > ref))
            if (ind > 1) {
               ci.left  <- theta.grid[ind - 1] + (theta.grid[ind] - theta.grid[ind - 1]) * 
                                        (ref - ll[ind - 1]) / (ll[ind] - ll[ind - 1])
               segments(ci.left, par()$usr[3], ci.left, ref, lty = 2, col = "grey")
               }
            else
               ci.left <- theta.grid[1]
            ind <- max(which(ll > ref))
            if (ind < ngrid) {
               ci.right  <- theta.grid[ind] + (theta.grid[ind + 1] - theta.grid[ind]) * 
                                        (ll[ind] - ref) / (ll[ind] - ll[ind + 1])
               segments(ci.right, par()$usr[3], ci.right, ref, lty = 2, col = "grey")
               }
            else
               ci.right <- theta.grid[ngrid]
            ypos <- par()$usr[3]
            segments(ci.left, ypos, ci.right, ypos, lwd = 2, col = "blue")
            if (ci.left > theta.grid[1] & ci.right < theta.grid[ngrid]) {
               ypos <- ypos + strheight("(")
               text(mean(c(ci.left, ci.right)), ypos,
                    paste("(", signif(ci.left, 5), ", ", signif(ci.right, 5), ")", sep = ""),
                    adj = rep(0.5, 2), col = "blue")
               }

            if (display["quadratic approximation"]) {
               ci.left  <- mle - 1.96 / sqrt(-mloglik$hessian)
               ci.right <- mle + 1.96 / sqrt(-mloglik$hessian)
               if (ci.left  >= theta.grid[1])
                  segments(ci.left,  par()$usr[3], ci.left,  ref, lty = 2, col = "grey")
               if (ci.right <= theta.grid[ngrid])
                  segments(ci.right, par()$usr[3], ci.right, ref, lty = 2, col = "grey")
               ypos <- ypos + 1.5 * strheight("(")
               segments(ci.left, ypos, ci.right, ypos, lwd = 2, col = "red")
               if (ci.left > theta.grid[1] & ci.right < theta.grid[ngrid]) {
                  ypos <- ypos + strheight("(")
                  text(mean(c(ci.left, ci.right)), ypos,
                       paste("(", signif(ci.left, 5), ", ", signif(ci.right, 5), ")", sep = ""),
                       adj = rep(0.5, 2), col = "red")
                  }
               }
            }

         title(paste(ylab, "function"), line = 3, col.main = "blue", cex.main = 1)
         title("Click and drag the vertical axes to adjust the plotting range.",
            line = 1, cex.main = 0.85)
         })
      panel
      }
   
   find.pt <- function(panel, x, y) {
      tol   <- diff(range(panel$theta.range)) / 25
      d.pts <- abs(panel$theta.range - x)
      axis.active   <- min(which(d.pts == min(d.pts)))
      if ((abs(panel$theta.range[axis.active] - x) < tol)) {
         panel$axis.active <- axis.active
         panel$new.range <- x
	     }
      rp.control.put(panel$panelname, panel)
      rp.tkrreplot(panel, plot1)
      panel
   }

   drag <- function(panel, x, y) {
      if (!is.na(panel$axis.active))
         panel$new.range <- x
      rp.control.put(panel$panelname, panel)
      rp.tkrreplot(panel, plot1)
      panel
   }
   
   release <- function(panel, x, y) {
	  if (!is.na(panel$axis.active))
         panel$theta.range[panel$axis.active] <- x
      panel$axis.active <- NA
      rp.control.put(panel$panelname, panel)
      rp.tkrreplot(panel, plot1)
      panel
   }

   loglik1.replot <- function(panel) {
      if (!is.na(panel$axis.active))   
         panel$theta.range[panel$axis.active] <- panel$new.range
      panel$axis.active <- NA
      rp.control.put(panel$panelname, panel)
      rp.tkrreplot(panel, plot1)
      panel
   }

   panel <- rp.control(data = data, ngrid = 50, loglik.text = loglik.text,
                       form = form, theta.range = theta.range, axis.active = NA,
                       display = c(mle = FALSE, "relative scale" = FALSE, 
                       "threshold proportion" = FALSE, ci = FALSE, 
                       "quadratic approximation" = FALSE), prop = 0.9)
   rp.tkrplot(panel, plot1, loglik1.plot, find.pt, drag, release, 
      hscale = hscale, vscale = vscale, pos = "right", background = "white")
   rp.textentry(panel, loglik.text, loglik1.replot,
      title = "Log-likelihood")
   rp.radiogroup(panel, form, c("likelihood", "log-likelihood"), action = loglik1.replot,
      title = "Form")
   rp.checkbox(panel, display, loglik1.replot, 
      c("mle", "relative scale", "threshold proportion", "ci", "quadratic approximation"))
   rp.slider(panel, prop, 0.01, 1, loglik1.replot, "Threshold proportion", 
      resolution = 0.01, showvalue = TRUE)
   rp.do(panel, loglik1.replot)
}

#   Two-parameter likelihood plots

rp.loglik2 <- function(loglik.text, data, theta.low, theta.high) {
	
   if (!require(rgl)) stop("the rgl library is required.")
      
   theta1.range <- c(theta.low[1], theta.high[1])
   theta2.range <- c(theta.low[2], theta.high[2])

   loglik2.plot <- function(panel) {

      fun.text <- paste("function(theta, data)", panel$loglik.text)

      warn <- options()$warn
      options(warn = -1)
      panel$mloglik <- optim(c(mean(panel$theta1.range), mean(panel$theta2.range)), 
                             eval(parse(text = fun.text)), data = panel$data,
                             control = list(fnscale = -1), hessian = TRUE)
      if (panel$mloglik$convergence > 0) {
         rp.messagebox("The mle cannot be identified.")
         convergence <- FALSE
         }
      else
         convergence <- TRUE
      options(warn = warn)
      
      panel$theta1.low  <- as.numeric(panel$ranges["theta1.low"])
      panel$theta1.high <- as.numeric(panel$ranges["theta1.high"])
      panel$theta2.low  <- as.numeric(panel$ranges["theta2.low"])
      panel$theta2.high <- as.numeric(panel$ranges["theta2.high"])
      theta1       <- seq(panel$theta1.low, panel$theta1.high, length = panel$ngrid)
      theta2       <- seq(panel$theta2.low, panel$theta2.high, length = panel$ngrid)
      thetamat     <- as.matrix(expand.grid(theta1, theta2))
      loglik.mat <- apply(thetamat, 1, eval(parse(text = fun.text)), data = panel$data)
      loglik.mat <- matrix(loglik.mat, ncol = panel$ngrid)
      threshold  <- min(apply(loglik.mat, 1, max))
      ind             <- (loglik.mat <= threshold)
      loglik.mat[ind] <- threshold
      ints <- seq(min(loglik.mat), max(loglik.mat), length = 51)
      ind  <- findInterval(c(loglik.mat), ints)
      clr  <- terrain.colors(50)[ind]
      scaling <- rp.plot3d(range(theta1), range(loglik.mat), range(theta2), 
                     xlab = "theta1", ylab = "log-likelihood", zlab = "theta2",
                     type = "n", new.window = FALSE)
      panel$loglik.mat    <- loglik.mat
      panel$threshold     <- threshold
      panel$theta1        <- theta1
      panel$theta2        <- theta2
      panel$scaling       <- scaling
      panel$clr.surface   <- clr
      panel$n.add         <- 0
      panel               <- loglik2.plot.add(panel)
      panel
   }

   loglik2.plot.add <- function(panel) {

      if (panel$n.add > 0) {
         for (i in 1:panel$n.add) pop3d()
         panel$n.add <- 0
      }

      with(panel, {
      	
      	 if (display["transparent"]) alpha.surface <- rep(0.5, ngrid^2)
            else                     alpha.surface <- rep(  1, ngrid^2)
         ind <- (loglik.mat <= threshold)
         alpha.surface[c(ind)] <- 0
         ax <- scaling(theta1,  theta1,  theta1)$x
         ay <- scaling(loglik.mat, loglik.mat, loglik.mat)$y
         az <- scaling(theta2, theta2, theta2)$z
         pop3d()
         rgl.surface(ax, az, ay, alpha = alpha.surface, col = clr.surface)
         material3d(alpha = 1)
         
         if (display["mle"]) {
           if ((mloglik$par[1] >= theta1.low) & (mloglik$par[1] <= theta1.high) & 
                (mloglik$par[2] >= theta2.low) & (mloglik$par[2] <= theta2.high)) {
               a <- scaling(mloglik$par[1], mloglik$value, mloglik$par[2])
               points3d(a$x, a$y, a$z, col = "red", size = 5)
               }
            }
            
         if (display["ci"]) {
            a <- scaling(mloglik$par[1], mloglik$value, mloglik$par[2])
            if ((mloglik$value - 3 > min(loglik.mat)) & 
                (mloglik$value - 3 < max(loglik.mat))) {
               a <- scaling(c(theta1.low, theta1.low, theta1.high, theta1.high),
                      rep(mloglik$value - 3, 4), 
                      c(theta2.high, theta2.low, theta2.low, theta2.high))
               quads3d(a$x, a$y, a$z, col = "red", size = 5)
               }
            }
            
         if (display["quadratic"] & ("hessian" %in% names(mloglik))) {
            theta1    <- seq(theta1.low, theta1.high, length = ngrid)
            theta2    <- seq(theta2.low, theta2.high, length = ngrid)
            thetamat  <- as.matrix(expand.grid(theta1, theta2))
            quad.fn   <- function(theta, mloglik) mloglik$value + 
                 0.5 * as.vector((theta - mloglik$par) %*% mloglik$hessian %*% (theta - mloglik$par))
            quad.mat  <- apply(thetamat, 1, quad.fn, mloglik = mloglik)
            quad.mat  <- matrix(quad.mat, ncol = ngrid)
            ind       <- (quad.mat <= threshold)
            quad.mat[ind] <- threshold
            ints <- seq(min(quad.mat), max(quad.mat), length = 51)
            ind  <- findInterval(c(quad.mat), ints)
            clr  <- terrain.colors(50)[ind]
            
      	    if (display["transparent"]) alpha.surface <- rep(0.5, ngrid^2)
               else                     alpha.surface <- rep(  1, ngrid^2)
            ind <- (quad.mat <= threshold)
            alpha.surface[c(ind)] <- 0
            ax   <- scaling(theta1,  theta1,  theta1)$x
            ay   <- scaling(quad.mat, quad.mat, quad.mat)$y
            az   <- scaling(theta2, theta2, theta2)$z
            rgl.surface(ax, az, ay, alpha = alpha.surface, col = clr, front = "line", back = "line")
            material3d(alpha = 1)
            }
         })
      panel$n.add <- 0
      if (panel$display["mle"] & (panel$mloglik$par[1] >= panel$theta1.low)  & 
                                 (panel$mloglik$par[1] <= panel$theta1.high) & 
                                 (panel$mloglik$par[2] >= panel$theta2.low)  & 
                                 (panel$mloglik$par[2] <= panel$theta2.high))
         panel$n.add <- panel$n.add + 1
      if (panel$display["ci"] & (panel$mloglik$value - 3 > min(panel$loglik.mat)) & 
                                (panel$mloglik$value - 3 < max(panel$loglik.mat)))   
         panel$n.add <- panel$n.add + 1
      if (panel$display["quadratic"] & ("hessian" %in% names(panel$mloglik))) 
         panel$n.add <- panel$n.add + 1
      panel
      }

   open3d()
   bg3d(color = c("white", "black"))
      
   panel <- rp.control(data = data, ngrid = 50, loglik.text = loglik.text,
                       theta1.range = theta1.range, theta2.range = theta2.range,
                       theta1.low  = theta1.range[1], theta1.high = theta1.range[2],
                       theta2.low  = theta2.range[1], theta2.high = theta2.range[2])
   rp.textentry(panel, loglik.text, loglik2.plot, title = "Log-likelihood")
   rp.textentry(panel, ranges, loglik2.plot,
                 labels = c("theta1.low", "theta1.high", "theta2.low", "theta2.high"),
                 initval = c(theta1.range[1], theta1.range[2], theta2.range[1], theta2.range[2]))
   rp.checkbox(panel,  display,  loglik2.plot.add, c("mle", "ci", "quadratic", "transparent"))
   rp.do(panel, loglik2.plot)
   }

   if (is.function(loglik.fn))
      loglik.text <- paste(deparse(substitute(loglik.fn)), "(theta, data)", sep = "")
   else
      loglik.text <- loglik.fn

   if (length(theta.low) == 1)
      rp.loglik1(loglik.text, data, theta.low, theta.high, form)
   else if (length(theta.low) == 2)
      rp.loglik2(loglik.text, data, theta.low, theta.high)
   else
      stop("theta.low must be a vector of length 1 or 2.")
      
   invisible()
   }
   
