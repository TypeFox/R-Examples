"sm.variogram" <- function(x, y, h,
                       df.se = "automatic", max.dist = NA, original.scale = TRUE, 
                       varmat = FALSE, ...) {

   type     <- "binned"
   bin.type <- "log"
   type.se  <- "smooth-monotonic-original"
   
   if ("geodata" %in% class(x)) {
      y <- x$data
      x <- x$coords
      }
   
   opt     <- sm.options(list(...))
   data    <- sm.check.data(x = x, y = y, ...)
   x       <- data$x
   y       <- data$y
   n       <- data$nobs
   ndim    <- data$ndim
   opt     <- data$options
   # rawdata <- list(x = x, y = y, nbins = opt$nbins, nobs = n, ndim = ndim)
   
   model <- opt$model
   replace.na(opt, band, (model != "none"))
   if (model == "isotropic") {
      replace.na(opt, ngrid, 20)
      replace.na(opt, display, "image")
      }
   if (model == "stationary") {
      replace.na(opt, ngrid, 12)
      }
   else {
      replace.na(opt, ngrid,   100)
      replace.na(opt, display, "means")
      }
   if      (model == "stationary") replace.na(opt, df, 20)
   else if (model == "isotropic")  replace.na(opt, df, 12)
   else                            replace.na(opt, df, 6)
   replace.na(opt, band,    (model != "none"))
   replace.na(opt, test,    (model != "none"))
   if (is.null(opt$weights.penalty)) opt$weights.penalty <- NA
   replace.na(opt, se,      TRUE)
   replace.na(opt, col,     "black")
   # replace.na(opt, weights.penalty,      FALSE)
   opt$weight.penalty <- FALSE
   if (model == "isotropic")  replace.na(opt, nbins, 10)
   if (model == "stationary") replace.na(opt, nbins,  6)
   if (ndim == 1) {
      x <- matrix(x, ncol = 1)
      replace.na(opt, nbins, 100)
      if (opt$nbins == 0) opt$nbins <- 100
      }
   else {
      replace.na(opt, nbins, ceiling((n * (n - 1) / 2)^(1/3)))
      if (opt$nbins == 0) opt$nbins <- ceiling((n * (n - 1) / 2)^(1/3))
      }

   x1mat <- matrix(rep(x[, 1], n), ncol = n, byrow = TRUE)
   if (ndim == 2) {
      x2mat <- matrix(rep(x[, 2], n), ncol = n)
      hmat  <- sqrt((t(x1mat) - x1mat)^2 + (t(x2mat) - x2mat)^2)
      }
   else
      hmat <- abs(t(x1mat) - x1mat)
   dmat  <- matrix(rep(y, n), ncol = n)
   dmat  <- t(dmat) - dmat
   dmat0 <- dmat^2
   dmat  <- sqrt(abs(dmat))
   imat  <- matrix(rep(1:n, n), ncol = n, byrow = TRUE)
   ind   <- as.vector(t(imat) - imat)
   hall  <- (as.vector(hmat))[ind > 0]
   dall0 <- (as.vector(dmat0))[ind > 0]
   dall  <- (as.vector(dmat))[ind > 0]
   i1    <- (as.vector(imat))[ind > 0]
   i2    <- (as.vector(t(imat)))[ind > 0]
   ipair <- cbind(i1, i2)
   
   # if (!is.na(max.dist)) {
   #    ind   <- (hall <= max.dist)
   #    hall  <- hall[ind]
   #    dall  <- dall[ind]
   #    ipair <- ipair[ind, ]
   #    }
   
   results <- list(distance = hall, sqrtdiff = dall, ipair = ipair)
   
#------------------------------------------------------------
#                   Bin the differences
#------------------------------------------------------------

   if (!(model %in% c("isotropic", "stationary"))) {
      
      if (bin.type == "regular") {
         bins   <- binning(hall, dall, nbins = opt$nbins)
         hh     <- bins$x
         dd     <- bins$means
         wts    <- bins$x.freq
         breaks <- bins$breaks
         empse  <- sqrt(bins$devs / (wts - 1)) / sqrt(wts)
         empse[wts == 1] <- 0
         igp    <- as.vector(cut(hall, breaks, labels = FALSE))
         ibin   <- match(igp, sort(unique(igp)))
      }
      else if (bin.type == "balanced") {
         test.pts <- sort(unique(signif(results$distance, 7)))
         test.pts <- test.pts[-length(test.pts)] + diff(test.pts) / 2
         test.pts <- c(test.pts, max(results$distance) + 1)
         breaks   <- -1
         for (i in 1:length(test.pts)) {
	        ind <- ((results$distance > max(breaks)) & (results$distance <= test.pts[i]))
	        if (sum(ind) >= length(results$distance) / opt$nbins) breaks <- c(breaks, test.pts[i])
         	}
         nbrks         <- length(breaks)
         breaks[nbrks] <- max(breaks[nbrks], max(results$distance))
         igp           <- as.numeric(cut(results$distance, breaks, labels = FALSE))
         ibin          <- match(igp, sort(unique(igp)))
         breaks[1]     <- 0
         dd            <- tapply(dall, ibin, mean)
         # hh            <- breaks[-1] - diff(breaks) / 2
         hh            <- tapply(results$distance, ibin, mean)
         wts           <- table(ibin) 
      }
      else if (bin.type == "unique") {
         hh   <- sort(unique(signif(hall, 6)))
         ibin <- match(signif(hall, 6), hh)
         wts  <- table(ibin)
         dd   <- tapply(dall, ibin, mean)
      }
      else if (bin.type == "log") {
         #  This doesn't handle 0 distances.  They should be added as a separate bin.
         breaks <- -1
         ind    <- (hall < 2 * .Machine$double.eps)
         nzero  <- length(which(ind))
         if (nzero > 4) breaks <- c(breaks, 2 * .Machine$double.eps)
         else           ind <- rep(FALSE, length(hall))
         breaks <- c(breaks, exp(min(log(hall[!ind])) + 
                            (1:opt$nbins) * diff(range(log(hall[!ind]))) / opt$nbins))
         nbrks  <- length(breaks)
         breaks[nbrks] <- breaks[nbrks] + 1
         igp    <- as.numeric(cut(results$distance, breaks), labels = FALSE)
         ibin   <- match(igp, sort(unique(igp)))
         dd0    <- as.vector(tapply(dall0, ibin, mean))
         dd     <- as.vector(tapply(dall,  ibin, mean))
         hh     <- as.vector(tapply(hall,  ibin, mean))
         wts    <- as.vector(table(ibin))
         breaks[nbrks] <- breaks[nbrks] - 1
         breaks[1] <- 0
      }

      results$distance.mean <- hh
      results$sqrtdiff.mean <- dd
      results$sqdiff.mean   <- dd0
      results$weights       <- wts
      results$ibin          <- ibin
      results$breaks        <- breaks
      results$nbins         <- length(hh)
      
      nbins <- length(results$sqrtdiff.mean)
      if (!is.numeric(df.se)) df.se <- round(0.8 * nbins)


#------------------------------------------------------------
#       Construct the estimate
#------------------------------------------------------------

      if (type == "binned") {
         ev        <- hh
         gamma.hat <- dd
      }
      else {
         if (missing(h)) h <- h.select(hh, dd, weights = wts, 
                            nbins = 0, df = df.se, method = opt$method)
         replace.na(opt, eval.points, seq(min(hall), max(hall), length = opt$ngrid))
         ev        <- opt$eval.points
         W         <- sm.weight(hh, ev, h, weights = wts, options = opt)
         gamma.hat <- as.vector(W %*% dd)
      }

#------------------------------------------------------------
#       Find the standard errors of the estimated values
#------------------------------------------------------------

      if (opt$se & (model == "none")) {
         if (is.na(max.dist)) max.dist <- max(hh) + 1
         # if (type.se == "true")
            # gamma.hat.V <- 1 - cov.spatial(results$distance.mean, cov.pars = cov.pars, kappa = kappa)
            # #    True gamma, evaluated at the observations
            # # gamma.hat <- 1 - cov.spatial(results$distance, cov.pars = c(1, phi))
         if (type.se == "binned")
            gamma.hat.V <- (dd / 0.977741)^4
         if (type.se == "cressie")
            gamma.hat.V <- 0.5 * dd^4 / (0.457 + 0.494 / wts)
         if (type.se == "smooth-monotonic-original") {
            sm.model            <- ps.normal(hh, 0.5 * dd0, df = df.se, weights = wts, 
                                       eval.points = hh, increasing = TRUE,
                                       display = "none")
            gamma.hat.V         <- sm.model$estimate
            # results$gamma.hat.V <- gamma.hat.V
         }
         if (type.se == "smooth-original") {
            sm.model            <- ps.normal(hh, 0.5 * dd0, df = df.se, weights = wts, 
                                       eval.points = hh, increasing = FALSE,
                                       display = "none")
            gamma.hat.V         <- sm.model$estimate
            # results$gamma.hat.V <- gamma.hat.V
         }
         if (type.se %in% c("smooth", "smooth-monotonic", "smooth-w", "smooth-monotonic-w")) {
            # ind                 <- (hh <= max.dist)
            # gamma.hat.V         <- rep(0, length(hh))
            # gamma.hat.V[ind]    <- sm.regression(hh, dd, weights = wts, eval.points = hh[ind], 
            #                          display = "none")$estimate
            # gamma.hat.V[!ind]   <- gamma.hat.V[hh == max(hh[ind])]
            inc <- (type.se %in% c("smooth-monotonic", "smooth-monotonic-w"))
            wp  <- (type.se %in% c("smooth-w", "smooth-monotonic-w"))
            sm.model            <- ps.normal(hh, dd, df = df.se, weights = wts, 
                                       eval.points = hh, increasing = inc, kappa = 1e8,
                                       weights.penalty = wp, display = "none")
            gamma.hat.V         <- sm.model$estimate
            # results$gamma.hat.V <- gamma.hat.V
            gamma.hat.V         <- (gamma.hat.V / 0.977741)^4
         }
         if (type.se %in% c("monotonic", "monotonic-original")) {
         	B     <- diag(nbins)
            D1    <- diff(diag(nbins), diff = 1)
            ddx   <- if (type.se == "monotonic") dd else dd0
            beta  <- ddx
            delta <- 1
            while (delta > 1e-5) {
               v <- as.numeric(diff(beta) <= 0)
               B1    <- solve(diag(wts) + 10000 * t(D1) %*% diag(v) %*% D1)
               beta.old <- beta
               beta  <- as.vector(B1 %*% diag(wts) %*% ddx)
               delta <- sum((beta - beta.old)^2) / sum(beta.old^2)
            }
            gamma.hat.V         <- beta
            # results$gamma.hat.V <- gamma.hat.V
            if (type.se == "monotonic") gamma.hat.V <- (gamma.hat.V / 0.977741)^4
         }
         # if (type.se == "matern") {
            # vgm.emp   <- variog(coords = x, data = y, estimator.type = "modulus", breaks = breaks,
                                   # messages = FALSE)
            # gamma.hat.V <- variofit(vgm.emp, ini = c(var(y), 0.2), fix.nugget = TRUE, fix.kappa = TRUE,
                                   # max.dist = max.dist, messages = FALSE)
            # results$cov.pars  <- gamma.hat.V$cov.pars 
            # results$kappa     <- gamma.hat.V$kappa
            # # plot(vgm.emp)
            # # lines(gamma.hat)
            # gamma.hat.V <- gamma.hat.V$cov.pars[1] - 
                # cov.spatial(results$distance.mean, cov.pars = gamma.hat.V$cov.pars, 
                # kappa = gamma.hat.V$kappa)
         # }
         #    Smoothing - small df and tilted to small distances
         # h  <- h.select(hh, dd, weights = wts, df = 4)
         # hw <- exp(hh^2) / mean(exp(hh^2))
         # gamma.hat <- sm.regression(hh, dd, df = 4, weights = wts, h.weights = hw,
         #                    eval.points = hh, display = "none")$estimate
         # gamma.hat <- (gamma.hat / 0.977741)^4
      
         gamma.hat.V <- pmax(gamma.hat.V, 0)
         
         if (type == "binned") {
            se <- rep(0, nbins)
# LAST MODIFICATION
# DIAG.COV.BIN.FUN
           # for (i in 1:nbins) {
           #    # ib    <- sort(unique(results$ibin))[i]
	     #      se[i] <- sqrt(cov.bin.fun(i, i, results, gamma.hat.V))
           # }
            result <- as.vector(matrix(1.0, nrow=length(gamma.hat.V), ncol = 1))
            rho.n  <- 50

            output <- .Fortran("diag_cov_bin_fun",
                        as.integer(length(gamma.hat.V)),
                        as.integer(nrow(results$ipair)),
                        as.integer(rho.n),
                        as.integer(results$ibin),
                        matrix(as.integer(results$ipair), ncol = 2),
                        as.double(gamma.hat.V),
                        res = as.double(result), PACKAGE = "sm")
            se     <- sqrt(output$res)
         }
         if (type != "binned" | varmat) {
            V <- matrix(0, nrow = nbins, ncol = nbins)

# LAST MODIFICATION
# FULL.COV.BIN.FUN
           # for (i in 1:nbins) {
	     #      # if (opt$verbose > 0) cat(i, "")
           #    for (j in i:nbins){
           #      # ib      <- sort(unique(results$ibin))[i]	        
           #       # jb      <- sort(unique(results$ibin))[j]	        
	     #         V[i, j] <- cov.bin.fun(i, j, results, gamma.hat.V)
	     #         if (j > i) V[j, i] <- V[i, j]
           #    }
           # }
	     # V <- full.cov.bin.fun(results,gamma.hat.V)
         result <- matrix(1.0, nrow=length(gamma.hat.V), ncol=length(gamma.hat.V))
         rho.n  <- 50
         output <- .Fortran("full_cov_bin_fun",
            as.integer(length(gamma.hat.V)),
            as.integer(nrow(results$ipair)),
            as.integer(rho.n),
            as.integer(results$ibin),
            matrix(as.integer(results$ipair), ncol = 2),
            as.double(gamma.hat.V),
            res=as.double(result), PACKAGE = "sm")
         V <- matrix(data=output$res, nrow=length(gamma.hat.V), ncol=length(gamma.hat.V), byrow= TRUE)

            # if (opt$verbose > 0) cat("\n")
            # save(V, file = "V.dmp")
            # load("V.dmp")
            if (type != "binned") se <- sqrt(diag(W %*% V %*% t(W)))
         }
      }
   }

#------------------------------------------------------------
#                     Test of independence
#------------------------------------------------------------

   if (model == "independent") {
   	
      vv    <- 0.1724
      cv    <- 0.03144
      Sigma <- table(c(igp, igp), c(i1, i2))
      Sigma <- Sigma %*% t(Sigma)
      Sigma <- cv * (Sigma - diag(2 * wts))
      Sigma <- Sigma / outer(wts, wts)
      Sigma <- diag(vv / wts) + Sigma

      if (opt$test | opt$band) {
      	 h    <- h.select(hh, hh, weights = wts, df = opt$df)
         W    <- sm.weight(hh, hh, h, weights = wts, options = opt)
         est  <- W %*% dd
         r0   <- sum(wts * (dd - mean(dall))^2)
         r1   <- sum(wts * (dd - est)^2)
         tobs <- (r0 - r1) / r1
         nb   <- length(hh)
         A    <- matrix(rep(wts / sum(wts), nb), ncol = nb, byrow = TRUE)
         A    <- t(diag(nb) - A) %*% diag(wts) %*% (diag(nb) - A)
         A    <- A - (1 + tobs) * t(diag(nb) - W) %*% diag(wts) %*% (diag(nb) - W)
         pval <- p.quad.moment(A, Sigma, 0, 0)
         if (opt$verbose > 0) 
            cat("Test of spatial independence: p = ",round(pval, 3), "\n")
         results$h <- h
         results$p <- pval
      }
         
   }

#------------------------------------------------------------
#                           Plots
#------------------------------------------------------------

   if ((opt$display != "none") & !(model %in% c("isotropic", "stationary"))) {
   	
      fn <- if (original.scale) function(x) (x / 0.977741)^4 else I
      
      if (opt$display %in% c("bins", "means")) {
         xx <- hh
         yy <- fn(dd)
      }
      else {
         xx <- hall
         yy <- if (original.scale) dall^4 else dall
      }
      
   	  if (!opt$add) {
   	  	
      	 replace.na(opt, xlab, "Distance")
      	 replace.na(opt, ylab, if (original.scale) " Squared difference" else "Square-root abs. difference")
         replace.na(opt, xlim, range(xx))

         r <- yy
         if (model == "independent") {
            replace.na(opt, eval.points, seq(min(hall), max(hall), length = opt$ngrid))
            ev  <- opt$eval.points
            W   <- sm.weight(hh, ev, h, weights = wts, options = opt)
            est <- c(W %*% dd)
            r <- c(r, fn(est))
            if (opt$band | opt$se) {
               sigmahat <- sqrt(var(y))
               nmeans   <- length(wts)
               V        <- matrix(rep(wts / sum(wts), length(ev)), ncol = nmeans, byrow = TRUE)
               se.band  <- sigmahat * sqrt(diag((W - V) %*% Sigma %*% t(W - V)))
               r        <- c(r, fn(mean(dall) + 2 * se.band), fn(mean(dall) - 2 * se.band))
            }
         }
         if (model == "none" & opt$se)
            r <- c(r, fn(gamma.hat - 2 * se), fn(gamma.hat + 2 * se))
         replace.na(opt, ylim, range(r))
         xlm <- opt$xlim
         if (!is.na(max.dist)) xlm[2] <- max.dist
         plot(xx, yy, xlab = opt$xlab, ylab = opt$ylab, xlim = xlm, ylim = opt$ylim, type = "n")
   	  }

      if (model == "independent" & opt$band)
         polygon(c(ev, rev(ev)), fn(c(mean(dall) + 2 * se.band, rev(mean(dall) - 2 * se.band))),
                 border = FALSE, col = opt$col.band)

      points(xx, yy, col = opt$col.points, pch = opt$pch)

      if (model == "none" & opt$se)
         segments(hh, fn(dd - 2 * se), hh, fn(dd + 2 * se), col = opt$col.points)

      if (model == "independent")
         lines(ev, fn(est), col = opt$col, lty = opt$lty)

#      if (type.se == "smooth") {
#         lines(ev, gamma.hat, col = opt$col, lty = opt$lty)
#         if (opt$se) {
#            lines(ev, gamma.hat + 2 * se, lty = 2, col = opt$col)
#            lines(ev, gamma.hat - 2 * se, lty = 2, col = opt$col)
#         }
#      }
      
   }

#------------------------------------------------------------
#                   Test of isotropy
#------------------------------------------------------------

   if (model == "isotropic") {

      imat   <- matrix(rep(1:n, n), ncol = n, byrow = TRUE)
      ind    <- as.vector(t(imat) - imat)
      amat   <- atan2(t(x2mat) - x2mat, t(x1mat) - x1mat)
      amat[amat < 0] <- amat[amat < 0] + pi
      angles <- (as.vector(amat))[ind > 0]
      
      # Remember to handle the case of 0 distances.
      centres1 <- seq(0, pi, length = opt$nbins + 1)
      centres1 <- centres1[-1] - diff(centres1) / 2
      centres2 <- c(0, exp(min(log(hall)) + (1:opt$nbins) * diff(range(log(hall))) / opt$nbins))
      centres2 <- centres2[-1] - diff(centres2) / 2
      centres  <- as.matrix(expand.grid(centres1, centres2))

      # bins  <- binning(cbind(angles, hall), dall, nbins = opt$nbins)
      identify.grid <- function(x, centres) {
         d      <- (x[1] - centres[, 1])^2 + (x[2] - centres[, 2])^2
         ind.pt <- which(d == min(d))
         gpt    <- centres[ind.pt, ]
         if (length(ind.pt) > 1) ind.pt <- ind.pt[gpt[, 1] == min(gpt[, 1])]
         if (length(ind.pt) > 1) ind.pt <- ind.pt[gpt[, 2] == min(gpt[, 2])]
         ind.pt
         }
      ibin  <- apply(cbind(angles, hall), 1, identify.grid, centres)   
      ibin  <- match(ibin, unique(ibin))
      dd    <- as.vector(tapply(dall,   ibin, mean))
      dd0   <- as.vector(tapply(dall0,  ibin, mean))
      hh    <- as.vector(tapply(hall,   ibin, mean))
      ang   <- as.vector(tapply(angles, ibin, mean))
      wts   <- as.vector(table(ibin))
      nbins <- length(unique(ibin))

      # sm.regression(cbind(hh, ang), dd, df = 12, weights = wts, nbins = 0, display = "rgl",
      #         col.points = "red", alpha = 0.7, size = 2, period = c(NA, pi))
         
      h <- h.select(cbind(hh, ang), dd, weights = wts, df = opt$df, nbins = 0, period = c(NA, pi))
      if (!is.numeric(df.se)) df.se <- round(0.8 * opt$nbins)
      
      results$distance.mean <- hh
      results$sqrtdiff.mean <- dd
      results$angles        <- angles
      results$angles.mean   <- ang
      results$weights       <- wts
      results$ibin          <- ibin
      results$h             <- h
      
      if (is.na(max.dist)) max.dist <- max(hh) + 1
      ind <- (hh <= max.dist)
      gamma.hat.V <- rep(0, length(hh))
      
      if (type.se == "binned")
         gg <- dd[ind]
      else if (type.se == "smooth")
         gg <- sm.regression(hh, dd, weights = wts, eval.points = hh[ind], 
                                     display = "none", nbins = 0)$estimate
      else if (type.se == "smooth-monotonic")
         gg <- ps.normal(hh, dd, df = df.se, weights = wts, eval.points = hh[ind],
                            increasing = TRUE, weights.penalty = FALSE, display = "none")$estimate
      else if (type.se == "smooth-monotonic-original")
         gg <- ps.normal(hh, 0.5 * dd0, df = df.se, weights = wts, 
                            negative = FALSE, increasing = TRUE,
                            eval.points = hh[ind], display = "none")$estimate
      
      gamma.hat.V[ind]    <- pmax(gg, 0)
      gamma.hat.V[!ind]   <- gamma.hat.V[hh == max(hh[ind])]
      # results$gamma.hat.V <- gamma.hat.V
      if (type.se != "smooth-monotonic-original")
         gamma.hat.V <- (gamma.hat.V / 0.977741)^4
               
#      plot(hh, dd)
#      plot(hh, 0.5 * dd0)
#      ps.normal(hh, 0.5 * dd0, df = df.se, weights = wts, negative = FALSE, increasing = TRUE)
#      print(cbind(hh, gamma.hat.V))
#      stop()
      
      V <- matrix(0, nrow = nbins, ncol = nbins)


# LAST MODIFICATION
# FULL.COV.BIN.FUN

#      if (opt$verbose > 0) cat(nbins, ": ")
#      for (i in 1:nbins) {
#	     if (opt$verbose > 0) cat(i, "")
#         for (j in i:nbins){
#            # ib      <- sort(unique(results$ibin))[i]	        
#            # jb      <- sort(unique(results$ibin))[j]	 
#	        V[i, j] <- cov.bin.fun(i, j, results, gamma.hat.V)
#	        if (j > i) V[j, i] <- V[i, j]
#            }
#         }

      result <- matrix(1.0, nrow=length(gamma.hat.V), ncol=length(gamma.hat.V))
      rho.n  <- 50
      output <- .Fortran("full_cov_bin_fun",
                         as.integer(length(gamma.hat.V)),
                         as.integer(nrow(results$ipair)),
                         as.integer(rho.n),
                         as.integer(results$ibin),
                         matrix(as.integer(results$ipair), ncol = 2),
                         as.double(gamma.hat.V),
                         res=as.double(result), PACKAGE = "sm")
      V <- matrix(data=output$res, nrow=length(gamma.hat.V), ncol=length(gamma.hat.V), byrow= TRUE)

      # if (opt$verbose > 0) cat("\n")
      
      opt$period <- c(NA, pi)
      
      if (opt$test | (opt$display != "none")) {
         mdl1  <- sm(dd ~ s(cbind(hh, ang), df = opt$df, period = opt$period), 
                     weights = wts, display = "none")
         df0   <- ceiling(sqrt(opt$df))
         mdl0  <- sm(dd ~ s(hh, df = df0), weights = wts, display = "none")
         # mdl0  <- sm(dd ~ s(hh, df = opt$df), weights = wts, display = "none")
         # mdl0  <- sm(dd ~ s(hh, df = opt$df / 2), weights = wts, display = "none")
         # mdl0  <- sm(dd ~ s(hh, df = opt$df / 3), weights = wts, display = "none")
         # mdl0  <- sm(dd ~ s(hh, df = 4), weights = wts, display = "none")
         # mdl0  <- sm(dd ~ s(hh, lambda = mdl$lambda[[1]][1] * (mdl$nseg[[1]][1] + 3)), 
         #                 weights = wts, display = "none")
      }
      
      if (opt$test) {
         
         # save(model, V, wts, dd, hh, ang, h, opt, file = "temp.dmp")
         
         S1  <- mdl1$B %*% mdl1$B1 %*% t(mdl1$B * wts)
         S0  <- mdl0$B %*% mdl0$B1 %*% t(mdl0$B * wts)  
         # est1  <- S0 %*% dd

#         S1   <- S0 - S1
#         V1   <- solve(V)
#         # V1   <- solve(S1 %*% V %*% t(S1))
#         ds   <- S1 %*% dd
#         tobs <- c(t(ds) %*% V1 %*% ds)
#         pval1 <- p.quad.moment.adjusted(t(S1) %*% V1 %*% S1, V, tobs)
#         
      	 # S0   <- sm.weight(hh, hh, h[1], weights = wts) 
         # S1   <- sm.weight2(cbind(hh, ang), cbind(hh, ang), h, weights = wts, options = opt)
         # cat("traces:", round(sum(diag(S0))), round(sum(diag(S1))), "\n")
#         est0 <- S0 %*% dd
#         plot(range(hh), range(est1, est0), type = "n")
#         points(hh, est0)
#         points(hh, est1, col = "green")
         # sm.regression(cbind(hh, ang), dd, h, weights = wts, options = opt, display = "image")

         # S0   <- t(diag(nbins) - S0) %*% V1 %*% (diag(nbins) - S0)
         # S1   <- t(diag(nbins) - S1) %*% V1 %*% (diag(nbins) - S1)
         # S0   <- t(diag(nbins) - S0) %*% (diag(nbins) - S0)
         # S1   <- t(diag(nbins) - S1) %*% (diag(nbins) - S1)
         # r0   <- c(dd %*% S0 %*% dd)
         # r1   <- c(dd %*% S1 %*% dd)
         # tobs <- (r0 - r1) / r1
         # S0   <- S0 - (1 + tobs) * S1
         # pval <- p.quad.moment(S0, V, 0, 0)
         S1   <- S0 - S1
         V1   <- solve(V)
         # V1   <- solve(S1 %*% V %*% t(S1))
         ds   <- S1 %*% dd
         tobs <- c(t(ds) %*% V1 %*% ds)
         pval <- p.quad.moment.adjusted(t(S1) %*% V1 %*% S1, V, tobs)
         
         # cat(round(pval, 3), round(pval1, 3), "\n")
         
#         mdl <- sm(dd ~ s(hh, df = 4) * s(ang, df = 4, period = pi), 
#                     weights = wts, display = "none")
#         # save(model, V, wts, file = "temp.dmp")
#         ind  <- c(mdl$b.ind[[2]], mdl$b.ind[[3]])
#         tobs <- sum(mdl$alpha[ind]^2)
#         I.i  <- rep(0, ncol(mdl$B))
#         I.i[ind] <- 1
#         A    <- mdl$B1 %*% t(mdl$B * wts)
#         pval <- p.quad.moment.adjusted(t(A) %*% diag(I.i) %*% A, V, tobs)
         
         if (opt$verbose > 0) cat("Test of isotropy: p = ", round(pval, 3), "\n")
         results$h <- h
         results$p <- pval
         # results$df0 <- df0
         }

      opt$covmat <- V
      opt$nbins  <- 0
      opt$alpha.mesh <- 0.7
      opt$test <- FALSE
      op <- opt
      replace.na(op, xlab, "Distance")
      replace.na(op, zlab, "Square-root difference")
      replace.na(op, ylab, "Angle")
      op$display <- "none"
      
      u     <- list(length = 2)
      for (j in 1:2)
         u[[j]] <- seq(mdl1$xrange[[1]][j, 1], mdl1$xrange[[1]][j, 2], length = opt$ngrid)
      U     <- as.matrix(expand.grid(u))
      mask  <- sm.mask(cbind(hh, ang), cbind(u[[1]], u[[2]]), mask.method = opt$mask.method)
      
      B1    <- ps.matrices(U, mdl1$xrange[[1]], 2, nseg = mdl1$nseg[[1]], period = mdl1$period[[1]])$B
      B1    <- cbind(1, B1)
      S1    <- B1 %*% mdl1$B1 %*% t(mdl1$B * wts)
      est1  <- matrix(S1 %*% dd, nrow = opt$ngrid) * mask
      
      B0    <- ps.matrices(as.matrix(U[ , 1]), mdl0$xrange[[1]], 1, nseg = mdl0$nseg[[1]])$B
      B0    <- cbind(1, B0)
      S0    <- B0 %*% mdl0$B1 %*% t(mdl0$B * wts)
      est0  <- matrix(S0 %*% dd, nrow = opt$ngrid) * mask
      
      stde  <- sqrt(diag((S1 - S0) %*% V %*% t(S1 - S0)))
      sdiff <- (est1 - est0) / matrix(stde, ncol = opt$ngrid)
      ev        <- u
      gamma.hat <- est1
      
      results$eval.points <- u
      results$estimate    <- est1
      results$sdiff       <- sdiff

      # surf <- sm.regression.2d(cbind(hh, ang), dd, h, model = "isotropic", weights = wts, 
      #                      rawdata = list(x = cbind(hh, ang), y = dd), options = op)
      # ngrid   <- nrow(surf$eval.points)
      # ev      <- rep(surf$eval.points[, 1], ngrid)
      # gamma.hat <- surf$estimate
      # X       <- sm.weight(hh, ev, h[1], weights = wts)
      # model.y <- matrix(as.vector(X %*% dd), ncol = opt$ngrid) * surf$estimate / surf$estimate

      if (opt$display == "rgl") {
      	 # This code needs to be updated
         # sm.surface3d(surf$eval.points, model.y, surf$scaling, col.mesh = "grey",
         #          alpha = 0.5, alpha.mesh = 0)
      }
      else if (opt$display == "image") {
      	 if (!require(akima)) stop("this option requires the akima package.")
         a     <- U
         a     <- rbind(a, cbind(a[ , 1], a[ , 2] + pi))
         a1    <- a[ , 1] * cos(a[ , 2])
         a2    <- a[ , 1] * sin(a[ , 2])
         b     <- rep(c(est1), 2)
         sdiff <- rep(c(sdiff), 2)
         ind   <- !is.na(b) & !duplicated(cbind(a1, a2))
         inte  <- interp(a1[ind], a2[ind], b[ind])
         ints  <- interp(a1[ind], a2[ind], sdiff[ind])
         cts   <- contourLines(ints)
         lvls  <- rep(0, length(cts))
         for (i in 1:length(cts)) lvls[i] <- cts[[i]]$level
         lvls  <- lvls[lvls <= -2]
         results$levels <- lvls
         filled.contour(inte, color.palette = topo.colors,
            plot.axes = {
               axis(1)
               axis(2)
               # op <- opt
               # op$display <- "none"
               # surf  <- sm.regression.2d(cbind(hh, ang), dd, h, model = "isotropic", weights = wts, 
               #             rawdata = list(x = cbind(hh, ang), y = dd), options = op)
               if (length(lvls) > 0) contour(ints, levels = lvls, add = TRUE, col = "red")
#         a    <- as.matrix(expand.grid(surf$eval.points[ , 1], surf$eval.points[ , 2]))
#         a    <- rbind(a, cbind(a[ , 1], a[ , 2] + pi))
#         a1   <- a[ , 1] * cos(a[ , 2])
#         a2   <- a[ , 1] * sin(a[ , 2])
#         b    <- rep(c(surf$estimate), 2)
#         ind  <- !is.na(b)
#         #image(interp(a1[ind], a2[ind], b[ind]))
#         #print(contourplot(interp(a1[ind], a2[ind], b[ind])$z))
#         filled.contour(interp(a1[ind], a2[ind], b[ind]), color.palette = topo.colors,
#            plot.axes = {
#               axis(1)
#               axis(2)
#               op <- opt
#               op$display <- "none"
#               surf  <- sm.regression.2d(cbind(hh, ang), dd, h, model = "isotropic", weights = wts, 
#                           rawdata = list(x = cbind(hh, ang), y = dd), options = op)
#               sdiff <- rep(c(surf$sdiff), 2)
#               cts   <- contourLines(interp(a1[ind], a2[ind], sdiff[ind]))
#               lvls  <- rep(0, length(cts))
#               for (i in 1:length(cts)) lvls[i] <- cts[[i]]$level
#               lvls  <- lvls[lvls < 0]
#               contour(interp(a1[ind], a2[ind], sdiff[ind]), levels = lvls, add = TRUE)
            }
         )
      }
      
      }

#------------------------------------------------------------
#                   Test of stationarity
#------------------------------------------------------------

   if (model == "stationary") {

      imat      <- matrix(rep(1:n, n), ncol = n, byrow = TRUE)
      ind       <- as.vector(t(imat) - imat)
      av1       <- (t(x1mat) + x1mat) / 2
      av2       <- (t(x2mat) + x2mat) / 2
      av1       <- (as.vector(av1))[ind > 0]
      av2       <- (as.vector(av2))[ind > 0]

      centres1  <- seq(min(av1), max(av1), length = opt$nbins + 1)
      centres2  <- seq(min(av2), max(av2), length = opt$nbins + 1)
      centres3  <- c(0, exp(min(log(hall)) + (1:opt$nbins) * diff(range(log(hall))) / opt$nbins))
      centres1  <- centres1[-1] - diff(centres1) / 2
      centres2  <- centres2[-1] - diff(centres2) / 2
      centres3  <- centres3[-1] - diff(centres3) / 2
      centres   <- as.matrix(expand.grid(centres1, centres2, centres3))

      identify.grid <- function(x, centres) {
         d      <- (x[1] - centres[, 1])^2 + (x[2] - centres[, 2])^2 + (x[3] - centres[, 3])^2
         ind.pt <- which(d == min(d))
         gpt    <- centres[ind.pt, ]
         if (length(ind.pt) > 1) ind.pt <- ind.pt[gpt[, 1] == min(gpt[, 1])]
         if (length(ind.pt) > 1) ind.pt <- ind.pt[gpt[, 2] == min(gpt[, 2])]
         if (length(ind.pt) > 1) ind.pt <- ind.pt[gpt[, 3] == min(gpt[, 3])]
         ind.pt
         }
      ibin  <- apply(cbind(av1, av2, hall), 1, identify.grid, centres)   
      ibin  <- match(ibin, unique(ibin))
      dd    <- as.vector(tapply(dall,  ibin, mean))
      dd0   <- as.vector(tapply(dall0, ibin, mean))
      hh    <- as.vector(tapply(hall,  ibin, mean))
      a1    <- as.vector(tapply(av1,   ibin, mean))
      a2    <- as.vector(tapply(av2,   ibin, mean))
      wts   <- as.vector(table(ibin))
      nbins <- length(unique(ibin))
      
      # h <- h.select(cbind(hh, ang), dd, weights = wts, df = opt$df, nbins = 0, period = c(NA, pi))
      if (!is.numeric(df.se)) df.se <- round(0.8 * opt$nbins)

      results$distance.mean <- hh
      results$sqrtdiff.mean <- dd
      # results$av1           <- av1
      # results$av2           <- av2
      # results$a1            <- a1
      # results$a2            <- a2
      results$weights       <- wts
      results$ibin          <- ibin
      
      xx <- cbind(x = a1, y = a2, Distance = hh)
      
      if (is.na(max.dist)) max.dist <- max(hh) + 1
      ind <- (hh <= max.dist)
      gamma.hat.V <- rep(0, length(hh))
      
      if (type.se == "binned")
         gg <- dd[ind]
      else if (type.se == "smooth")
         gg <- sm.regression(hh, dd, weights = wts, eval.points = hh[ind], 
                                     display = "none", nbins = 0)$estimate
      else if (type.se == "smooth-monotonic")
         gg <- ps.normal(hh, dd, df = df.se, weights = wts, eval.points = hh,
                            increasing = TRUE, weights.penalty = FALSE, display = "none")$estimate
      else if (type.se == "smooth-monotonic-original")
         gg <- ps.normal(hh, 0.5 * dd0, df = df.se, weights = wts, 
                            negative = FALSE, increasing = TRUE,
                            eval.points = hh[ind], display = "none")$estimate
      
      gamma.hat.V[ind]    <- pmax(gg, 0)
      gamma.hat.V[!ind]   <- gamma.hat.V[hh == max(hh[ind])]
      # results$gamma.hat.V <- gamma.hat.V
      if (type.se != "smooth-monotonic-original")
         gamma.hat.V <- (gamma.hat.V / 0.977741)^4
               

      V <- matrix(0, nrow = nbins, ncol = nbins)
      # if (opt$verbose > 0) cat(nbins, ": ")


# LAST MODIFICATION
# FULL.COV.BIN.FUN

#      for (i in 1:nbins) {
#	     if (opt$verbose > 0) cat(i, "")
#         for (j in i:nbins){
#	        V[i, j] <- cov.bin.fun(i, j, results, gamma.hat.V)
#	        if (j > i) V[j, i] <- V[i, j]
#          }
#       }

      result <- matrix(1.0, nrow=length(gamma.hat.V), ncol=length(gamma.hat.V))
      rho.n  <- 50
      output <- .Fortran("full_cov_bin_fun",
                         as.integer(length(gamma.hat.V)),
                         as.integer(nrow(results$ipair)),
                         as.integer(rho.n),
                         as.integer(results$ibin),
                         matrix(as.integer(results$ipair), ncol = 2),
                         as.double(gamma.hat.V),
                         res=as.double(result), PACKAGE = "sm")
      V <- matrix(data=output$res, nrow=length(gamma.hat.V), ncol=length(gamma.hat.V), byrow= TRUE)

      # if (opt$verbose > 0) cat("\n")

      model1 <- sm(dd ~ s(xx, df = opt$df), weights = wts, display = "none")
            
      if (opt$test) {
      	 df0 <- ceiling(opt$df^(1/3))
      	 # df0 <- ceiling(opt$df / 3)
      	 # df0 <- ceiling(opt$df / 2)
      	 model0 <- sm(dd ~ s(hh, df = df0), weights = wts, display = "none")
         # mdl0  <- sm(dd ~ s(hh, lambda = mdl$lambda[[1]][1] * (mdl$nseg[[1]][1] + 3)^2), 
         #                 weights = wts, display = "none")
      	 # model0 <- sm(dd ~ s(hh, df = 3), weights = wts, display = "none")
         S0     <- model0$B %*% model0$B1 %*% t(model0$B * wts)
      	 S1     <- model1$B %*% model1$B1 %*% t(model1$B * wts)
         V1     <- solve(V)
         S0     <- t(S0 - S1) %*% V1 %*% (S0 - S1)
         tobs   <- c(dd %*% S0 %*% dd)
         pval   <- p.quad.moment.adjusted(S0, V, tobs)
         if (opt$verbose > 0) 
            cat("Test of stationarity: p = ", round(pval, 3), "\n")
         results$p   <- pval
         # results$df0 <- df0
      }

      u     <- list(length = 3)
      for (j in 1:3)
         u[[j]] <- seq(model1$xrange[[1]][j, 1], model1$xrange[[1]][j, 2], length = opt$ngrid)
      U     <- as.matrix(expand.grid(u))      

      B1    <- ps.matrices(U, model1$xrange[[1]], 3, nseg = model1$nseg[[1]], period = model1$period[[1]])$B
      B1    <- cbind(1, B1)
      S1    <- B1 %*% model1$B1 %*% t(model1$B * wts)
      est1  <- S1 %*% dd
      
      B0    <- ps.matrices(as.matrix(U[ , 3]), model0$xrange[[1]], 1, nseg = model0$nseg[[1]])$B
      B0    <- cbind(1, B0)
      S0    <- B0 %*% model0$B1 %*% t(model0$B * wts)
      est0  <- S0 %*% dd
      
      stde  <- diag((S1 - S0) %*% V %*% t(S1 - S0))
      stde[stde < 0] <- NA
      stde  <- sqrt(stde)
      model1$sdiff <- array((est1 - est0) / stde, dim = rep(opt$ngrid, 3))

      # ev    <- u
      # gamma.hat <- est1

      opt$covmat <- V
      #  Add V into the estimated surface.
      #  Return the estimate from the fitted surface.
      # plot(model1)
      
      # replace.na(opt, xlab, "Distance")
      # replace.na(opt, zlab, "Square-root difference")
      # replace.na(opt, ylab, "Angle")
      gamma.hat <- model1$fitted
      ev <- xx
      # results$model <- model1
      
      names(u) <- c("x1", "x2", "distance")
      results$eval.points <- u
      results$estimate    <- array(est1, dim = rep(opt$ngrid, 3))
      results$sdiff       <- model1$sdiff

   }


#------------------------------------------------------------
#                   Return values
#------------------------------------------------------------

   # results$eval.points <- ev
   # results$estimate    <- gamma.hat
   # results$gamma.hat.V <- gamma.hat.V
   results$df    <- opt$df
      
   if (opt$se & model == "none")            results$se <- se
   if ((model == "independent") & opt$band) results$se.band <- se.band
   if (varmat | model %in% c("isotropic", "stationary")) results$V <- V
   invisible(results)
   
   }



p.quad.moment.adjusted <- function (A, Sigma, tobs) {
    B <- A %*% Sigma
    k1 <- sum(diag(B)) - tobs
    C <- B %*% B
    k2 <- 2 * sum(diag(C))
    k3 <- 8 * sum(diag(C %*% B))
    aa <- abs(k3/(4 * k2))
    bb <- (8 * k2^3)/k3^2
    cc <- k1 - aa * bb
    1 - pchisq(-cc/aa, bb)
    }

# hg <- function(rho) {
   # a    <- 3/4
   # b    <- 3/4
   # cc   <- 1/2
   # fn   <- gamma(a) * gamma(b) / gamma(cc)
   # facn <- 1
   # n    <- 1
   # fnold <- 0.1
   # while (abs(fn - fnold) / fnold > 0.0001 ) {
      # facn  <- facn * n
      # fnold <- fn
      # fn    <- fn + gamma(a + n) * gamma(b + n) * rho^n / (gamma(cc + n) * facn)
      # n     <- n + 1
      # }
   # fn * gamma(cc) / (gamma(a) * gamma(b))
   # }

# cor.sqrtabs <- function(rho) {
   # # library(Davies)
   # # gamma(0.75)^2 * ((1 - rho^2) * hypergeo(0.75, 0.75, 0.5, rho^2) - 1) / (sqrt(pi) - gamma(0.75)^2)
   # gamma(0.75)^2 * ((1 - rho^2) * hg(rho^2) - 1) / (sqrt(pi) - gamma(0.75)^2)
   # }
