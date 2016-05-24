
"sm.regression" <- function(x, y, h, design.mat = NA, model = "none",
         weights = NA, group = NA, ... ) {
         	
   if(!all(is.na(group))) 
      return(sm.ancova(x, y, group, h, model, weights=weights,...))

   x.name <- deparse(substitute(x))
   if (isMatrix(x)) x.names <- dimnames(x)[[2]]
   y.name <- deparse(substitute(y))

   opt     <- sm.options(list(...))
   data    <- sm.check.data(x = x, y = y, weights = weights, group = group, ...)
   x       <- data$x
   y       <- data$y
   weights <- data$weights
   group   <- data$group
   nobs    <- data$nobs
   ndim    <- data$ndim
   opt     <- data$options

   replace.na(opt, nbins, round((nobs > 500) * 8 * log(nobs) / ndim))
   rawdata <- list(x = x, y = y, nbins = opt$nbins, nobs = nobs, ndim = ndim)
   
   if (!((model %in% "none")      | (model %in% "no effect") | 
         (model %in% "no.effect") | (model %in% "linear")))
     stop("invalid setting for model argument.", call. = FALSE)
   if (model != "none") replace.na(opt, test, TRUE)

   if(missing(h))
     h <- h.select(x = x, y = y, weights = weights, ...)
   else
     {if(length(h) != ndim) stop("length(h) does not match size of x")}

   if(opt$nbins > 0) {
     if (!all(weights == 1) & opt$verbose > 0)
        cat("Warning: weights overwritten by binning\n")
     if (!all(is.na(opt$h.weights)))
        stop("use of h.weights is incompatible with binning - set nbins=0")
     bins         <- binning(x, y, nbins = opt$nbins)
     x            <- bins$x
     y            <- bins$means
     weights      <- bins$x.freq
     rawdata$devs <- bins$devs
     nx           <- length(y)
     }
   else 
     nx <- nobs
  
   replace.na(opt, h.weights, rep(1,nx))
   if (opt$panel && !require(rpanel)) {
      opt$panel <- FALSE
      cat("The rpanel package is not available.\n")
      }

   if (ndim == 1) {
      replace.na(opt, xlab,  x.name)
      replace.na(opt, ylab,  y.name)
      replace.na(opt, ngrid, 50)
      opt$period <- opt$period[1]
      if (opt$pch == ".") replace.na(opt, cex, 1)
         else replace.na(opt, cex,  2/log(rawdata$nobs))
      if (!opt$panel) 
         est <- sm.regression.1d(x, y, h, design.mat, 
                  model, weights, rawdata, options = opt)
      else {
         rp.smooth1(x, y, h, design.mat, model, weights, rawdata, opt)
         }
      }
   else {
      replace.na(opt, ngrid, 20)
      dimn <- x.names # dimnames(x)[[2]]
      name.comp<-if(!is.null(dimn) & !all(dimn=="")) dimn
             else {if(!is.null(attributes(x)$names)) attributes(x)$names
             else outer(x.name,c("[1]","[2]"),paste,sep="")}
      replace.na(opt, xlab, name.comp[1])
      replace.na(opt, ylab, name.comp[2])
      replace.na(opt, zlab, y.name)
      if (all(is.na(opt$period))) opt$period <- rep(NA, 2)
      if (!(length(opt$period) == 2))
         stop("the length of period should match the number of covariates.")
      if (opt$panel) 
         rp.smooth2(x, y, h, model, weights, rawdata, opt)
      else
         est <- sm.regression.2d(x, y, h, model, weights, rawdata, options = opt)            
      }
      
   if (opt$panel)
      invisible()
   else {
      est$data <- list(x = x, y = y, opt$nbins, freq = weights)
      est$call <- match.call()
      invisible(est)
      }
      
   }


"sm.regression.1d" <- function (x, y, h, design.mat = NA, model = "none", 
                    weights = rep(1, length(x)), rawdata, options = list()) {
    	
    opt <- sm.options(options)
    replace.na(opt, ngrid,      50)
    replace.na(opt, xlim,       range(rawdata$x))
    replace.na(opt, ylim,       range(rawdata$y))
    replace.na(opt, display,    "line")
    replace.na(opt, col,        "black")
    replace.na(opt, col.band,   "cyan")
    replace.na(opt, col.points, "black")
    replace.na(opt, se,         FALSE)
    hmult <- opt$hmult
    if (model == "none") {
        opt$band <- FALSE
        opt$test <- FALSE
        }
    else
        replace.na(opt, band, TRUE)
    band <- opt$band
    if (opt$add | opt$display %in% "none")
        opt$panel <- FALSE
    r <- list(x = NA, y = NA, model.y = NA, se = NA, sigma = NA,
           h = h * hmult, hweights = opt$h.weights, weights = weights)
    if (!opt$add & !(opt$display %in% "none"))
        plot(rawdata$x, rawdata$y, xlab = opt$xlab, ylab = opt$ylab,
           xlim = opt$xlim, ylim = opt$ylim, type = "n")
    if (!(opt$display %in% "none")) {
        opt1 <- opt
        opt1$test <- FALSE
        r <- smplot.regression(x, y, design.mat, h, r, model, weights,
            rawdata, options = opt1)
        }
    if (opt$test)
        rtest <- sm.regression.test(x, y, design.mat, h, model,
            weights, rawdata, options = opt)
    if (!(any(is.na(opt$eval.points))))
        r <- sm.regression.eval.1d(x, y, design.mat, h, model,
            weights, rawdata, options = opt)
    else if ((opt$display %in% "none") & (model == "none")) {
        opt$eval.points <- seq(min(x), max(x), length = opt$ngrid)
        r <- sm.regression.eval.1d(x, y, design.mat, h, model,
            weights, rawdata, options = opt)
        }
    if (opt$test)
        r <- list(eval.points = r$eval.points, estimate = r$estimate,
            model.y = r$model.y, se = r$se, sigma = r$sigma,
            h = r$h, hweights = r$hweights, weights = weights,
            model = rtest$model, p = rtest$p, q.squared=rtest$q.squared)
    r
    }


"smplot.regression" <- function (x, y, design.mat, h, r, model, weights, rawdata = list(),
    options = list(), ...) {
    	
    opt <- sm.options(options)
    rnew <- sm.regression.eval.1d(x, y, design.mat, h, model,
        weights = weights, rawdata = rawdata, options = opt)

    if (!any(is.na(r$x))) {
        if (opt$band) {
            upper <- r$model.y + 2 * r$se
            upper <- pmin(pmax(upper, par()$usr[3]), par()$usr[4])
            lower <- r$model.y - 2 * r$se
            lower <- pmin(pmax(lower, par()$usr[3]), par()$usr[4])
            polygon(c(r$eval.points, rev(r$eval.points)), c(lower, rev(upper)),
                    col = 0, border = 0)
            }
        if (opt$se | (opt$display %in% "se")) {
            upper <- r$estimate + 2 * r$se
            upper <- pmin(pmax(upper, par()$usr[3]), par()$usr[4])
            lower <- r$estimate - 2 * r$se
            lower <- pmin(pmax(lower, par()$usr[3]), par()$usr[4])
            lines(r$eval.points, upper, lty = 3, col = 0)
            lines(r$eval.points, lower, lty = 3, col = 0)
            }
        lines(r$eval.points, r$estimate, col = 0)
    }
    if (opt$band) {
        upper <- rnew$model.y + 2 * rnew$se
        upper <- pmin(pmax(upper, par()$usr[3]), par()$usr[4])
        lower <- rnew$model.y - 2 * rnew$se
        lower <- pmin(pmax(lower, par()$usr[3]), par()$usr[4])
        polygon(c(rnew$eval.points, rev(rnew$eval.points)), c(lower,
            rev(upper)), col = opt$col.band, border = 0)
    }
    lines(rnew$eval.points, rnew$estimate, lty = opt$lty, col = opt$col, lwd = opt$lwd)
    if ((model == "none") & (opt$se | (opt$display %in% "se"))) {
        upper <- rnew$estimate + 2 * rnew$se
        upper <- pmin(pmax(upper, par()$usr[3]), par()$usr[4])
        lower <- rnew$estimate - 2 * rnew$se
        lower <- pmin(pmax(lower, par()$usr[3]), par()$usr[4])
        lines(rnew$eval.points, upper, lty = 3, col = opt$col)
        lines(rnew$eval.points, lower, lty = 3, col = opt$col)
    }
    if (!opt$add)
        points(rawdata$x, rawdata$y, col = opt$col.points, pch = opt$pch,
               cex = opt$cex)
    box(col = 1, lty = 1)
    rnew
}


"sm.regression.2d" <- function (x, y, h, model = "none", weights = rep(1, length(y)), rawdata,
    options = list()) {

    opt <- sm.options(options)
    # Identify the setting of opt$col and other things follow. 
    if (!(model %in% "none")) {
       if (is.na(opt$col)) {
          if (!is.na(opt$se)) {
       	     if (opt$se) opt$col <- "se"
             }
          else if (!is.na(opt$band)) {
       	     if (opt$band) opt$col <- "se"
             }
          else
             opt$col <- "se"
          }
       }
    else {
       if (is.na(opt$col)) {
          if (!is.na(opt$se)) {
       	     if (opt$se) opt$col <- "se"
             }
          }
       }
    if (is.na(opt$col))  opt$col <- "green"
    if (opt$col == "se") opt$se <- TRUE
    surf.ids <- rep(NA, 2)
       	  
    if (!is.na(opt$band) && opt$band)          opt$col <- "se"
    if (!is.na(opt$col)  && opt$col == "se" && is.na(opt$se))   opt$se  <- TRUE
    replace.na(opt, h.weights, rep(1, length(y)))
    replace.na(opt, display,   "persp")
    replace.na(opt, band,      FALSE)
    replace.na(opt, col,       "green")
    replace.na(opt, se,        FALSE)
    
    replace.na(opt, ngrid, 20)
    if (any(is.na(opt$eval.points))) {
        replace.na(opt, xlim, range(x[, 1]))
        if (opt$display %in% "rgl") {
        	replace.na(opt, zlim, range(x[, 2]))
            replace.na(opt, eval.points,
                cbind(seq(opt$xlim[1], opt$xlim[2], length = opt$ngrid),
                      seq(opt$zlim[1], opt$zlim[2], length = opt$ngrid)))
            }
        else {
            replace.na(opt, ylim, range(x[, 2]))
            replace.na(opt, eval.points,
                cbind(seq(opt$xlim[1], opt$xlim[2], length = opt$ngrid),
                      seq(opt$ylim[1], opt$ylim[2], length = opt$ngrid)))
            }
        }
    else {
        replace.na(opt, xlim, range(opt$eval.points[, 1]))
        if (opt$display %in% "rgl")
            replace.na(opt, zlim, range(opt$eval.points[, 2]))
        else
            replace.na(opt, ylim, range(opt$eval.points[, 2]))
        if (opt$eval.grid &
                (any(diff(opt$eval.points[,1]) < 0) |
                 any(diff(opt$eval.points[,2]) < 0)))
            stop(paste("eval.points are not suitable for grid evaluation",
                       "(eval.grid = TRUE)."))
        }
        
    # sigma <- sm.sigma(rawdata$x, rawdata$y, nbins = opt$nbins)$estimate
    sigma <- sm.sigma(x, y, rawdata, weights = weights, nbins = 0)$estimate
    
    if (!opt$eval.grid) {
        w <- sm.weight2(x, opt$eval.points, h, weights = weights,
                        options = opt)
        est <- as.vector(w %*% y)
        model.y <- est
        if (opt$se) se <- diag(w %*% t(w)) * sigma
	    }
    else {
        est <- sm.regression.eval.2d(x, y, h, model, opt$eval.points,
                          opt$hull, weights, options = opt)
        x1grid  <- opt$eval.points[,1]
        x2grid  <- opt$eval.points[,2]
        ngrid   <- length(x1grid)
        evpts   <- cbind(rep(x1grid, ngrid), rep(x2grid, each = ngrid))
        model.y <- est
        est.col <- est
        n       <- length(y)
        
        if (opt$se) {
        	
           S    <- sm.weight2(x, evpts, h = h, weights = weights, options = opt)
           mask <- est / est
           
           if (model == "none") {
              se      <- matrix(diag(S %*% t(S)) * sigma, ncol = ngrid) * mask
              }
           else if (model == "no effect") {
              X <- matrix(rep(1, n), ncol = 1)
              Z <- matrix(rep(1, ngrid^2), ncol = 1)
              }
           else if (model == "linear") {
              X <- cbind(rep(1, n), x[,1] - mean(x[,1]), x[,2] - mean(x[,2]))
              Z <- cbind(rep(1, ngrid * ngrid), evpts[,1] - mean(x[,1]),
                                                evpts[,2] - mean(x[,2]))
              }
           if (model %in% c("no effect", "linear")) {
              X       <- Z %*% solve(t(X) %*% diag(weights) %*% X) %*% t(X) %*% diag(weights)
              model.y <- matrix(as.vector(X %*% y), ncol = ngrid)
              S       <- S - X
              }
              
           if (model == "isotropic") {
              # est    <- S %*% dd
              # S      <- sm.weight2(cbind(hh, ang), cbind(hh, ang), sp, weights = wts)
              X       <- sm.weight(x[, 1], evpts[, 1], h[1], weights = weights)
              S       <- S - X
              se      <- matrix(sqrt(diag(S %*% opt$covmat %*% t(S))), ncol = ngrid) * mask
              model.y <- matrix(as.vector(X %*% y), ncol = ngrid) * mask
              }
           else {
              se    <- matrix(sqrt(diag(S %*% t(S))), ncol = ngrid) * sigma * mask
              }
           sdiff <- (est - model.y) / se
              
           if (opt$display %in% "persp") {
              se    <- array(c(se[-ngrid, -ngrid], se[    -1, -ngrid], 
                               se[-ngrid,     -1], se[    -1,     -1]), 
                               dim = c(ngrid - 1, ngrid - 1, 4))
              se    <- apply(se, 1:2, function(x) 
                          if (length(which(is.na(x))) > 1) NA else mean(x, na.rm = TRUE))
              se    <- matrix(c(se), nrow = ngrid - 1, ncol = ngrid - 1)
              sdiff <- array(c(sdiff[-ngrid, -ngrid], sdiff[    -1, -ngrid],
                                sdiff[-ngrid,     -1], sdiff[    -1,     -1]),
                                dim = c(ngrid - 1, ngrid - 1, 4))
              sdiff <- apply(sdiff, 1:2, function(x) 
                          if (length(which(is.na(x))) > 1) NA else mean(x, na.rm = TRUE))
              sdiff <- matrix(c(sdiff), nrow = ngrid - 1, ncol = ngrid - 1)
              }
           }
        else
           se <- NA
           
        if (opt$col == "height") {
           if (opt$display %in% "persp") {
              est.col <- array(c(est.col[-ngrid, -ngrid], est.col[    -1, -ngrid],
                                 est.col[-ngrid,     -1], est.col[    -1,     -1]),
                                 dim = c(ngrid - 1, ngrid - 1, 4))
              est.col <- apply(est.col, 1:2, function(x) 
                          if (length(which(is.na(x))) > 1) NA else mean(x, na.rm = TRUE))
              }
           opt$col <- opt$col.palette[cut(c(est.col), length(opt$col.palette), labels = FALSE)]
           }
        else if (opt$col == "se") {
           if (model == "none") {
              opt$col <- opt$col.palette[cut(c(se), length(opt$col.palette), labels = FALSE)]
              }
           else {
           	  if (length(opt$col.palette) != length(opt$se.breaks) + 1)
           	     opt$col.palette <- 
           	                rev(rainbow(length(opt$se.breaks) + 1, start = 0/6, end = 4/6))
           	  opt$se.breaks <- c(min(sdiff, na.rm = TRUE) - 1, sort(opt$se.breaks), 
           	                     max(sdiff, na.rm = TRUE) + 1)
              opt$col <- opt$col.palette[cut(c(sdiff), opt$se.breaks, labels = FALSE)]
              }
           }
        if (length(opt$col) > 1) {
           if      (opt$display %in% "rgl") opt$col <- matrix(opt$col, ncol = ngrid)
           else if (opt$display %in% "rgl") opt$col <- matrix(opt$col, ncol = ngrid - 1)
           }
           
        if (opt$display %in% "rgl") {
           if (opt$col.mesh == "height")
              opt$col.mesh <- opt$col.palette[cut(c(est), length(opt$col.palette), labels = FALSE)]
           else if ((opt$col.mesh == "se") & (model == "none"))
              opt$col.mesh <- opt$col.palette[cut(c(se), length(opt$col.palette), labels = FALSE)]
           if (length(opt$col.mesh) > 1) opt$col.mesh <- matrix(opt$col.mesh, ncol = ngrid)
           }

        if (opt$display %in% "image") {
            replace.na(opt, zlim, range(est, na.rm = TRUE))
            image(x1grid, x2grid, est, col = opt$col.palette,
                  xlab = opt$xlab, ylab = opt$ylab,
                  xlim = opt$xlim, ylim = opt$ylim,  zlim = opt$zlim, add = opt$add)
            }
        else if (opt$display %in% "slice")
            contour(x1grid, x2grid, est,
                  xlab = opt$xlab, ylab = opt$ylab,
                  xlim = opt$xlim, ylim = opt$ylim,
                  lty = opt$lty, col = opt$col, add = opt$add)
        else if (opt$display %in% "persp") {
           replace.na(opt, zlim, range(est, na.rm = TRUE))
           if (length(opt$col) == 1) {
              if (opt$col == 1) opt$col <- "green"
              opt$col <- matrix(opt$col, nrow = ngrid, ncol = ngrid)
              }
            persp(x1grid, x2grid, est,
                  xlab = opt$xlab, ylab = opt$ylab, zlab = opt$zlab,
                  xlim = opt$xlim, ylim = opt$ylim, zlim = opt$zlim,
                  theta = opt$theta, phi = opt$phi,
                  ticktype = "detailed", col = c(opt$col), d = 4)
            }
        else if ((opt$display %in% "rgl") && (require(rgl) & require(rpanel))) {
            replace.na(opt, ylim, range(y, est, na.rm = TRUE))
            if (!opt$add)
               opt$scaling <- rp.plot3d(rawdata$x[, 1], rawdata$y, rawdata$x[, 2],
                         xlab = opt$xlab, ylab = opt$zlab, zlab = opt$ylab,
                         xlim = opt$xlim, ylim = opt$ylim, zlim = opt$zlim,
                         size = opt$size, col = opt$col.points)
            surf.ids <- sm.surface3d(cbind(x1grid, x2grid), est, opt$scaling, 
                   col = opt$col, col.mesh = opt$col.mesh, 
                   alpha = opt$alpha, alpha.mesh = opt$alpha.mesh,
                   lit = opt$lit)
            }
        }

    r <- list(eval.points = opt$eval.points, estimate = est, model.y = model.y,
        sigma = sigma, h = h * opt$hmult, hweights = opt$h.weights,
        weights = weights, scaling = opt$scaling, surf.ids = surf.ids)
    if (opt$se) {
       r$se <- se
       r$sdiff <- sdiff
    }
    if (model != "none" & opt$test) {
        rtest <- sm.regression.test(x, y, design.mat = NA, h,
                    model, weights, rawdata, opt)
        r$model     <- rtest$model
        r$p         <- rtest$p
        r$q.squared <- rtest$q.squared
        }
    r
}


"sm.surface3d" <- function(eval.points, surf, scaling, 
                              col = "green", col.mesh = "black", alpha = 0.7, alpha.mesh = 1, 
                              lit = TRUE, ...) {

      #     This function adds a surface to the current rgl plot.
      
      if (!is.function(scaling)) stop("a scaling must be specified.")
      
      if (all(is.na(col)))                  col   <- "green"
      if ((length(col) == 1) && (col == 1)) col   <- "green"

      ep <- eval.points
      if (is.matrix(ep) && ncol(ep) == 2) {
         ep1 <- ep[ , 1]
         ep2 <- ep[ , 2]
      }
      else if (is.list(ep) && length(ep) == 2) {
         ep1 <- ep[[1]]
         ep2 <- ep[[2]]
      }
      else
         stop("the form of eval.points in sm.surface3d is invalid.")

      ngrid1 <- length(ep1)
      ngrid2 <- length(ep2)         
      col      <- matrix(c(col),      nrow = ngrid1, ncol = ngrid2)
      col.mesh <- matrix(c(col.mesh), nrow = ngrid1, ncol = ngrid2)
      
      xg1    <- rep(ep1[-ngrid1], ngrid2 - 1)
      xg2    <- rep(ep1[     -1], ngrid2 - 1)
      xg3    <- rep(ep1[     -1], ngrid2 - 1)
      xg4    <- rep(ep1[-ngrid1], ngrid2 - 1)
      zg1    <- rep(ep2[-ngrid2], each = ngrid1 - 1)
      zg2    <- rep(ep2[-ngrid2], each = ngrid1 - 1)
      zg3    <- rep(ep2[     -1], each = ngrid1 - 1)
      zg4    <- rep(ep2[     -1], each = ngrid1 - 1)
      yg1    <- c(surf[-ngrid1, -ngrid2])
      yg2    <- c(surf[     -1, -ngrid2])
      yg3    <- c(surf[     -1,      -1])
      yg4    <- c(surf[-ngrid1,      -1])
      col1   <- c(col[-ngrid1, -ngrid2])
      col2   <- c(col[     -1, -ngrid2])
      col3   <- c(col[     -1,      -1])
      col4   <- c(col[-ngrid1,      -1])
      ind1   <- !is.na(yg1 + yg2 + yg3)
      ind2   <- !is.na(yg1 + yg3 + yg4)
      xg     <- c(c(rbind( xg1,  xg2,  xg3)[, ind1]), c(rbind( xg1,  xg3,  xg4)[, ind2]))
      yg     <- c(c(rbind( yg1,  yg2,  yg3)[, ind1]), c(rbind( yg1,  yg3,  yg4)[, ind2]))
      zg     <- c(c(rbind( zg1,  zg2,  zg3)[, ind1]), c(rbind( zg1,  zg3,  zg4)[, ind2]))
      colg   <- c(c(rbind(col1, col2, col3)[, ind1]), c(rbind(col1, col3, col4)[, ind2]))
      ind3   <- is.na(yg3) & !is.na(yg1 + yg2 + yg4)
      xg     <- c(  xg, c(rbind( xg1,  xg2,  xg4)[, ind3]))
      yg     <- c(  yg, c(rbind( yg1,  yg2,  yg4)[, ind3]))
      zg     <- c(  zg, c(rbind( zg1,  zg2,  zg4)[, ind3]))
      colg   <- c(colg, c(rbind(col1, col2, col4)[, ind3]))
      ind4   <- is.na(yg1) & !is.na(yg2 + yg3 + yg4)
      xg     <- c(  xg, c(rbind( xg2,  xg3,  xg4)[, ind4]))
      yg     <- c(  yg, c(rbind( yg2,  yg3,  yg4)[, ind4]))
      zg     <- c(  zg, c(rbind( zg2,  zg3,  zg4)[, ind4]))
      colg   <- c(colg, c(rbind(col2, col3, col4)[, ind4]))
      a      <- scaling(xg, yg, zg)
      id1    <- triangles3d(a$x, a$y, a$z, col = colg, alpha = alpha, lit = lit,...)
      
      xg1    <- rep(ep1[-ngrid1], ngrid2)
      xg2    <- rep(ep1[     -1], ngrid2)
      xg3    <- rep(ep1         , each = ngrid2 - 1)
      xg4    <- rep(ep1         , each = ngrid2 - 1)
      zg1    <- rep(ep2         , each = ngrid1 - 1)
      zg2    <- rep(ep2         , each = ngrid1 - 1)
      zg3    <- rep(ep2[-ngrid2], ngrid1)
      zg4    <- rep(ep2[     -1], ngrid1)
      yg1    <- c(surf[-ngrid1, ])
      yg2    <- c(surf[     -1, ])
      yg3    <- c(t(surf[     , -ngrid2]))
      yg4    <- c(t(surf[     ,      -1]))
      col1   <- c(col.mesh[-ngrid1,        ])
      col2   <- c(col.mesh[     -1,        ])
      col3   <- c(t(col.mesh[     , -ngrid2]))
      col4   <- c(t(col.mesh[     ,      -1]))
      ind1   <- !is.na(yg1 + yg2)
      ind2   <- !is.na(yg3 + yg4)
      xg     <- c(c(rbind( xg1,  xg2)[, ind1]), c(rbind( xg3,  xg4)[, ind2]))
      yg     <- c(c(rbind( yg1,  yg2)[, ind1]), c(rbind( yg3,  yg4)[, ind2]))
      zg     <- c(c(rbind( zg1,  zg2)[, ind1]), c(rbind( zg3,  zg4)[, ind2]))
      colg   <- c(c(rbind(col1, col2)[, ind1]), c(rbind(col3, col4)[, ind2]))
      a      <- scaling(xg, yg, zg)
      id2    <- segments3d(a$x, a$y, a$z, col = colg, alpha = alpha.mesh, lit = lit, ...)

      invisible(c(id1, id2))
      }
      
      
"sm.regression.eval.1d" <- function (x, y, design.mat, h, model = "none", weights = rep(1,
    length(x)), rawdata, options = list()) {
    	
    opt <- sm.options(options)

    replace.na(opt, band, FALSE)
    replace.na(opt, test, FALSE)
    replace.na(opt, ngrid, 50)
    replace.na(opt, eval.points, seq(min(x), max(x), length = opt$ngrid))
    if (missing(rawdata))
        rawdata <- list(x = x, y = y, nbins = 0)
    band <- opt$band
    test <- opt$test
    ngrid <- opt$ngrid
    h.weights <- opt$h.weights
    eval.points <- opt$eval.points
    w <- sm.weight(x, eval.points, h, weights = weights, options = opt)
    est <- as.vector(w %*% y)
    sig <- sm.sigma(x, y, rawdata = rawdata, weights = weights)$estimate
    n <- length(x)
    ne <- length(eval.points)
    if (model == "none") {
        model.y <- est
        se <- as.vector(sig * sqrt(((w^2) %*% (1/weights))))
    }
    else if ((model == "no.effect") | (model == "no effect")) {
        if (is.na(as.vector(design.mat)[1])) {
            X <- matrix(rep(1, n), ncol = 1)
            model.y <- rep(wmean(y, weights), ne)
        }
        else {
            X <- design.mat
            model.y <- rep(0, ne)
        }
        X <- diag(n) - X %*% solve(t(X) %*% diag(weights) %*%
            X) %*% t(X) %*% diag(weights)
        se <- sig * sqrt(diag(w %*% X %*% diag(1/weights) %*%
            t(w)))
    }
    else if (model == "linear") {
        e <- cbind(rep(1, ne), eval.points - mean(x))
        l <- cbind(rep(1, n), x - mean(x))
        l <- e %*% solve(t(l) %*% diag(weights) %*% l) %*% t(l) %*%
            diag(weights)
        model.y <- as.vector(l %*% y)
        se <- as.vector(sig * sqrt(((w - l)^2) %*% (1/weights)))
    }
    list(eval.points = eval.points, estimate = est, model.y = model.y,
        se = se, sigma = sig, h = h * opt$hmult, hweights = h.weights,
        weights = weights)
}


"sm.regression.eval.2d" <- function (x, y, h, model, eval.points, 
               hull = TRUE, weights, options = list()) {

    opt <- sm.options(options)
    hmult <- opt$hmult
    h.weights <- opt$h.weights
    n <- nrow(x)
    ngrid <- nrow(eval.points)
    wd1 <- matrix(rep(eval.points[, 1], n), ncol = n)
    wd1 <- wd1 - matrix(rep(x[, 1], ngrid), ncol = n, byrow = TRUE)
    wd2 <- matrix(rep(eval.points[, 2], n), ncol = n)
    wd2 <- wd2 - matrix(rep(x[, 2], ngrid), ncol = n, byrow = TRUE)
    wy <- matrix(rep(h.weights, ngrid), ncol = n, byrow = TRUE)
    if (!is.na(opt$period[1])) 
       w1 <- exp(cos(2 * pi * wd1 / opt$period[1]) / (h[1] * hmult * wy))
    else
       w1 <- exp(-0.5 * (wd1 / (h[1] * hmult * wy))^2)
    w1 <- w1 * matrix(rep(weights, ngrid), ncol = n, byrow = TRUE)
    if (!is.na(opt$period[2])) 
       w2 <- exp(cos(2 * pi * wd2 / opt$period[2]) / (h[2] * hmult * wy))
    else
       w2 <- exp(-0.5 * (wd2 / (h[2] * hmult * wy))^2)
    wy <- matrix(rep(y, ngrid), ncol = n, byrow = TRUE)
    if ((opt$poly.index == 0) | (sum(is.na(opt$period)) == 0))
        est <- w1 %*% t(w2 * wy)/(w1 %*% t(w2))
    else if ((opt$poly.index == 1) & (length(opt$period) == 2) & (sum(is.na(opt$period)) == 1)) {
        x1grid  <- opt$eval.points[,1]
        x2grid  <- opt$eval.points[,2]
        ngrid   <- length(x1grid)
        evpts   <- cbind(rep(x1grid, ngrid), rep(x2grid, each = ngrid))
        w       <- sm.weight2(x, evpts, h, weights = weights, options = opt)
        est     <- matrix(c(w %*% y), ncol = ngrid)
        }
    else {
        a11 <- w1 %*% t(w2)
        a12 <- (w1 * wd1) %*% t(w2)
        a13 <- w1 %*% t(w2 * wd2)
        a22 <- (w1 * wd1^2) %*% t(w2)
        a23 <- (w1 * wd1) %*% t(w2 * wd2)
        a33 <- w1 %*% t(w2 * wd2^2)
        d   <- a22 * a33 - a23^2
        b1  <- 1/(a11 - ((a12 * a33 - a13 * a23) * a12 + (a13 *
                a22 - a12 * a23) * a13)/d)
        b2  <- (a13 * a23 - a12 * a33) * b1/d
        b3  <- (a12 * a23 - a13 * a22) * b1/d
        c1  <- w1 %*% t(w2 * wy)
        c2  <- (w1 * wd1) %*% t(w2 * wy)
        c3  <- w1 %*% t(w2 * wy * wd2)
        est <- b1 * c1 + b2 * c2 + b3 * c3
    }
    if (hull) {
        hull.points <- x[order(x[, 1], x[, 2]), ]
        dh <- diff(hull.points)
        hull.points <- hull.points[c(TRUE, !((dh[, 1] == 0) & (dh[, 2] == 0))), ]
        hull.points <- hull.points[chull(hull.points), ]
        nh <- nrow(hull.points)
        gstep <- matrix(rep(eval.points[2, ] - eval.points[1, ], nh),
                        ncol = 2, byrow = TRUE)
        hp.start <- matrix(rep(eval.points[1, ], nh), ncol = 2,
            byrow = TRUE)
        hull.points <- hp.start + gstep * round((hull.points -
            hp.start)/gstep)
        hull.points <- hull.points[chull(hull.points), ]
        grid.points <- cbind(rep(eval.points[, 1], ngrid), rep(eval.points[,
            2], rep(ngrid, ngrid)))
        D <- diff(rbind(hull.points, hull.points[1, ]))
        temp <- D[, 1]
        D[, 1] <- D[, 2]
        D[, 2] <- (-temp)
        C <- as.vector((hull.points * D) %*% rep(1, 2))
        C <- matrix(rep(C, ngrid^2), nrow = ngrid^2, byrow = TRUE)
        D <- t(D)
        wy <- ((grid.points %*% D) >= C)
        wy <- apply(wy, 1, all)
        wy[wy] <- 1
        wy[!wy] <- NA
        wy <- matrix(wy, ncol = ngrid)
    }
    else {
        w1 <- (w1 > exp(-2))
        w2 <- (w2 > exp(-2))
        wy <- w1 %*% t(w2)
        wy[wy > 0] <- 1
        wy[wy == 0] <- NA
    }
    est <- est * wy
    invisible(est)
}


"sm.regression.test" <- function (x, y, design.mat = NA, h, model = "no.effect", 
          weights = rep(1, length(y)), rawdata, options = list())
{
    opt <- sm.options(options)
    if (length(dim(x)) > 0) {
        ndim <- 2
        n <- dim(x)[1]
        W <- sm.weight2(x, x, h, weights = weights, options = opt)
        S <- cbind(rep(1, n), x[, 1] - mean(x[, 1]), x[, 2] - mean(x[, 2]))
        }
    else {
        ndim <- 1
        n <- length(x)
        W <- sm.weight(x, x, h, weights = weights, options = opt)
        S <- cbind(rep(1, n), x - mean(x))
        }
    if ((model == "no.effect") | (model == "no effect")) {
        if (is.na(as.vector(design.mat)[1]))
            S <- matrix(rep(1, n), ncol = 1)
        else 
            S <- design.mat
        }
    if ((model == "linear") | (model == "no.effect") | (model == "no effect")) {
        S <- diag(n) - S %*% solve(t(S) %*% diag(weights) %*% S) %*% t(S) %*% diag(weights)
        W <- diag(n) - W
        W <- t(W) %*% diag(weights) %*% W
        e <- as.vector(S %*% y)
        r0 <- sum(weights * e^2) + sum(rawdata$devs)
        r1 <- as.numeric(t(e) %*% W %*% e) + sum(rawdata$devs)
        ts <- (r0 - r1)/r1
        p <- p.quad.moment(diag(weights) - (1 + ts) * W, S %*%
            diag(1/weights), ts, sum(weights) - length(weights))
        q.squared <- (r0-r1)/(r0 + sum(rawdata$devs))
        }
    if(opt$verbose > 0) 
      cat(paste("Test of", model,"model:  significance = ",round(p,3),"\n"))
    list(model = model, p = p, h = h * opt$hmult, hweights =
         opt$h.weights,  q.squared=q.squared)
         
}


"p.quad.moment" <- function (A, Sigma, tobs, ndevs) {
    B <- A %*% Sigma
    k1 <- sum(diag(B)) - tobs * ndevs
    C <- B %*% B
    k2 <- 2 * sum(diag(C)) + 2 * tobs^2 * ndevs
    k3 <- 8 * sum(diag(C %*% B)) - 8 * tobs^3 * ndevs
    aa <- abs(k3/(4 * k2))
    bb <- (8 * k2^3)/k3^2
    cc <- k1 - aa * bb
    1 - pchisq(-cc/aa, bb)
    }


"p.quad.moment.old" <- 
function(A, Sigma, cnst)
{
        B <- A %*% Sigma
        k1 <- sum(diag(B))
        C <- B %*% B
        k2 <- 2 * sum(diag(C))
        k3 <- 8 * sum(diag(C %*% B))
        aa <- abs(k3/(4 * k2))
        bb <- (8 * k2^3)/k3^2
        cc <- k1 - aa * bb
        # print(paste("Degrees of freedom = ",bb))
        # print(paste("Chisq = ",((cnst-cc)/aa-bb)/sqrt(2*bb)))
        1 - pchisq((cnst - cc)/aa, bb)
}


"sm.sigma" <- function (x, y, 
               rawdata = NA, weights = rep(1, length(y)), 
               diff.ord = 2, ci = FALSE, model = "none", h = NA, ...) {

    opt <- sm.options(list(...))
    if (!is.list(rawdata)) rawdata <- list(devs = 0)
    
	data    <- sm.check.data(x = x, y = y, weights = weights)
    x       <- data$x
    y       <- data$y
    weights <- data$weights
    n       <- data$nobs
    ndim    <- data$ndim
    
    replace.na(opt, nbins, round((n > 500) * 8 * log(n) / ndim))
    if(opt$nbins > 0) {
      if (!all(weights == 1) & opt$verbose > 0)
         cat("Warning: weights overwritten by binning\n")
      if (!all(is.na(opt$h.weights)))
         stop("use of h.weights is incompatible with binning - set nbins=0")
      bins         <- binning(x, y, nbins = opt$nbins)
      x            <- bins$x
      y            <- bins$means
      weights      <- bins$x.freq
      rawdata$devs <- bins$devs
      n            <- length(y)
      }

    if (ndim == 2) return(sm.sigma2(x, y, 
           rawdata = rawdata, weights = weights, ci = ci, model = model, h = h,
           options = opt))

    if (diff.ord == 1) {
        yd <- diff(y[order(x)])
        ww <- 1/weights[order(x)]
        wd <- ww[2:n] + ww[1:(n - 1)]
        ssq <- sum(yd^2/wd) + sum(rawdata$devs)
        sig <- sqrt(ssq / (sum(weights) - 1))
        }
    else {
        yy <- y[order(x)]
        xx <- sort(x)
        xx1 <- diff(xx)
        xx2 <- diff(xx, lag = 2)
        a <- xx1[-1]/xx2
        b <- xx1[-(n - 1)]/xx2
        a[xx2 == 0] <- 0.5
        b[xx2 == 0] <- 0.5
        ww <- weights[order(x)]
        cc <- a^2/ww[1:(n - 2)] + b^2/ww[3:n] + 1/ww[2:(n - 1)]
        eps <- yy[1:(n - 2)] * a + yy[3:n] * b - yy[2:(n - 1)]
        ssq <- sum(eps^2/cc) + sum(rawdata$devs)
        sig <- sqrt(ssq / (sum(ww) - 2))
    }
    invisible(list(estimate = sig))
}


"sm.sigma2" <- function(x, y, 
            rawdata = list(devs = 0), weights = rep(1, length(y)),
            stand = "local", cross = FALSE, ci = FALSE,
			simple = FALSE, model = "none", h = NA, strip = FALSE,
			display = "none", options = list()) {
    	
  opt <- sm.options(options)

  n  <- length(y)
  x1 <- x[, 1]
  x2 <- x[, 2]
  if (any(is.na(x1 + x2 + y))) stop("Missing data not allowed in sm.sigma2")
  
  if (strip) {
    ch <- chull(x1, x2)
    x1 <- x1[-ch]
    x2 <- x2[-ch]
    y  <-  y[-ch]
    }
  X <- cbind(x1, x2)

  #  Repeated values can cause difficulties with zero nearest neighbour
  #  distances, so use the unique values.
  
  if (all(weights == rep(1, length(y))) & (nrow(unique(X)) < nrow(X))) {
    X <- paste(as.character(X[,1]), as.character(X[,2]), sep = ",")
    weights <- table(X)
    rawdata$devs <- tapply(y, factor(X), function(x) sum((x - mean(x))^2))
    y <- tapply(y, factor(X), mean)
    X <- sort(unique(X))
    X <- paste(X, collapse = ",")
    X <- paste("c(", X, ")", sep = "")
    X <- eval(parse(text = X))
    X <- matrix(X, ncol = 2, byrow = TRUE)
    }
  
  # if (simple) {
  #   S <- sm.weight2.nn(cbind(x1, x2), cross = cross)
  #   }
  # else {
  	nx <- nrow(X)
    d  <- matrix(rep(X[,1], nx), ncol = nx)
    D  <- (d - t(d))^2
    d  <- matrix(rep(X[,2], nx), ncol = nx)
    D  <- sqrt(D / var(X[,1]) + ((d - t(d))^2) / var(X[,2]))
    nn <- function(x) {sort(x)[5]}
    hw <- apply(D, 1, nn) / 2
    hh <- c(sqrt(var(X[,1])), sqrt(var(X[,2])))
    S  <- sm.weight2(X, X, hh, cross = cross, weights = weights, 
               options = list(h.weights = hw, hw.eval = TRUE))
  #   }

  A <- (diag(nx) - S)

  if (stand=="local") {
     A   <- t(A) %*% diag(1/diag(A %*% t(A))) %*% A
     adf <- n
     }
  else {
     A   <- t(A) %*% A
     adf <- sum(diag(A))
     }
  
  sigmahat <- as.numeric(((t(y) %*% A %*% y) + sum(rawdata$devs)) / sum(weights))
  if (sigmahat <= 0) {
     if (opt$verbose > 0) warning(" estimated standard deviation is 0")
     sigmahat <- 0
     }
  sigmahat <- sqrt(sigmahat)
  # sigmahat <- as.numeric(sqrt((t(y) %*% A %*% y) / adf))
  result   <- list(estimate = sigmahat, qmat = A / adf)

  if (model == "constant") {

    if (any(is.na(h))) stop("Smoothing parameter missing")
    
    P          <- sqrt(diag(1/diag(A %*% t(A)))) %*% A
    pseudo.res <- P %*% y
    svcomp     <- sqrt(abs(pseudo.res))

    S  <- sm.weight2(X, X, h, options = list())
    S  <- diag(n) - S
    S  <- t(S) %*% S
    L  <- matrix(1/n, ncol = n, nrow=n)
    r0 <- as.numeric(t(svcomp) %*% (diag(n) - L) %*% svcomp)
    r1 <- as.numeric(t(svcomp) %*% S %*% svcomp)
    ts <- (r0 - r1)/r1
    P  <- (P %*% t(P))^2
    diag(P) <- 0
    P  <- 5.545063 * ((1-P) * hyperg(P) - 1)
    diag(P) <- 1
    p  <- p.quad.moment(diag(n) - L - (1 + ts) * S, P, ts, 0)
    # cat(paste("Test of constant variance:  significance = ", round(p, 3),
    # 		"\n"))

    if (!(display == "none")) {
      surface <- sm.regression(X, svcomp, h, display = "none")
      contour(surface$eval.points[,1], surface$eval.points[,2],
    			surface$estimate)
      }

    result$pseudo.residuals <- pseudo.res
    result$p <- p
    }

  #  Confidence interval

  ie <- NA
  if (ci) {
    B  <- A / adf
    k1 <- sum(diag(B))
    C  <- B %*% B
    k2 <- 2 * sum(diag(C))
    k3 <- 8 * sum(diag(C %*% B))
    aa <- abs(k3 / (4 * k2))
    bb <- (8 * k2^3) / k3^2
    cc <- k1 - aa * bb
    q  <- qchisq(c(0.025, 0.975), bb)
    ie <- rev(sigmahat/sqrt(cc + q * aa))
    result$ci <- ie
    }

  invisible(result)

  }


"hyperg" <- function(z) {
  a  <- 0.75
  b  <- 0.75
  cc <- 0.5
  hg    <- gamma(a) * gamma(b) / gamma(cc)
  hgold <- 0.5
  n     <- 0
  lfac  <- 0
  while (max((hg - hgold)/hgold) > 0.001) {
    n     <- n + 1
    lfac  <- lfac + log(n)
    hgold <- hg
    hg    <- hg + exp(n * log(z) - lfac +
    			lgamma(a + n) + lgamma(b + n) - lgamma(cc + n))
    }
  hg <- hg * gamma(cc) / (gamma(a) * gamma(b))
  hg
  }
  
  
"sm.sigma2.compare" <- function(x1, y1, x2, y2) {

  data   <- sm.check.data(x1, y1)
  x1     <- data$x
  y1     <- data$y
  n1     <- data$nobs
  ndim1  <- data$ndim
  data   <- sm.check.data(x2, y2)
  x2     <- data$x
  y2     <- data$y
  n2     <- data$nobs
  ndim2  <- data$ndim
  if (!(ndim1 == 2 & ndim2 == 2)) 
     stop("x1 and x2 should be two-column matrices.")

  sig  <- sm.sigma2(x1, y1, options = list(nbins = 0))
  est1 <- sig$estimate
  A1   <- sig$qmat

  sig  <- sm.sigma2(x2, y2, options = list(nbins = 0))
  est2 <- sig$estimate
  A2   <- sig$qmat

  Fobs <- est1^2 / est2^2
  Fobs <- max(Fobs, 1/Fobs)
  Al <- rbind(A1, matrix(0, nrow=n2, ncol=n1))
  Ar <- rbind(matrix(0, nrow=n1, ncol=n2), -Fobs * A2)
  A  <- cbind(Al, Ar)
  p  <- 2 * p.quad.moment(A, diag(n1 + n2), 0, 0)
  # cat(paste("Test of equality of variances: p =", round(p, 3), "\n"))

  invisible(p)
  }


"sm.weight" <- function (x, eval.points, h, cross = FALSE, weights = rep(1, length(x)),
          options = list()) {

    if (!exists(".sm.Options")) stop("cannot find .sm.Options")
    opt <- sm.options(options)
    replace.na(opt, hmult, 1)
    replace.na(opt, h.weights, rep(1, length(x)))
    replace.na(opt, poly.index, 1)
    poly.index <- opt$poly.index
    h.weights <- opt$h.weights
    hmult     <- opt$hmult
    period    <- opt$period[1]
    n  <- length(x)
    ne <- length(eval.points)
    wd <- matrix(rep(eval.points, rep(n, ne)), ncol = n, byrow = TRUE)
    wd <- wd - matrix(rep(x, ne), ncol = n, byrow = TRUE)
    w <- matrix(rep(h.weights, ne), ncol = n, byrow = TRUE)
    if (!is.na(opt$period)) 
       w <- exp(cos(2 * pi * wd / period) / (h * hmult * w))
    else
       w <- exp(-0.5 * (wd / (h * hmult * w))^2)
    w <- w * matrix(rep(weights, ne), ncol = n, byrow = TRUE)
    if (cross)
        diag(w) <- 0
    if ((poly.index == 0) | (!is.na(period))) {
        den <- w %*% rep(1, n)
        w   <- w / matrix(rep(den, n), ncol = n)
        }
    else {
        s0 <- w %*% rep(1, n)
        s1 <- (w * wd) %*% rep(1, n)
        s2 <- (w * wd^2) %*% rep(1, n)
        w  <- w * (matrix(rep(s2, n), ncol = n) - wd * matrix(rep(s1, n), ncol = n))
        w  <- w / (matrix(rep(s2, n), ncol = n) * matrix(rep(s0, n), ncol = n) - 
                    matrix(rep(s1, n), ncol = n)^2)
        }
    }


"sm.weight2" <- function (x, eval.points, h, cross = FALSE, weights = rep(1, nrow(x)),
    options = list()) {

    opt <- sm.options(options)
    if (all(is.na(opt$period))) opt$period <- rep(NA, 2)
    replace.na(opt, hmult, 1)
    replace.na(opt, h.weights, rep(1, nrow(x)))
    replace.na(opt, poly.index, 1)
    poly.index <- opt$poly.index
    h.weights  <- opt$h.weights
    hmult      <- opt$hmult
    n          <- nrow(x)
    ne         <- nrow(eval.points)

    wd1 <- matrix(rep(eval.points[, 1], rep(n, ne)), ncol = n, byrow = TRUE)
    wd1 <- wd1 - matrix(rep(x[, 1], ne), ncol = n, byrow = TRUE)
    if (("hw.eval" %in% names(opt)) & (opt$hw.eval = TRUE))
       hw  <- matrix(rep(h.weights, n),  ncol = n)
    else
       hw  <- matrix(rep(h.weights, ne), ncol = n, byrow = TRUE)
    if (!is.na(opt$period[1])) 
       w <- exp(cos(2 * pi * wd1 / opt$period[1]) / (h[1] * hmult * hw))
    else
       w <- exp(-0.5 * (wd1 / (h[1] * hmult * hw))^2)
    wd2 <- matrix(rep(eval.points[, 2], rep(n, ne)), ncol = n, byrow = TRUE)
    wd2 <- wd2 - matrix(rep(x[, 2], ne), ncol = n, byrow = TRUE)
    if (!is.na(opt$period[2])) 
       w <- w * exp(cos(2 * pi * wd2 / opt$period[2]) / (h[2] * hmult * hw))
    else
       w <- w * exp(-0.5 * (wd2 / (h[2] * hmult * hw))^2)
    w <- w * matrix(rep(weights, ne), ncol = n, byrow = TRUE)
    
    if (cross)
        diag(w) <- 0

    if ((opt$poly.index == 0) | (sum(is.na(opt$period)) == 0)) {
        den <- w %*% rep(1, n)
        w   <- w/matrix(rep(den, n), ncol = n)
        }
    else if ((opt$poly.index == 1) & (sum(is.na(opt$period)) == 1)) {
    	if (is.na(opt$period[2])) wd1 <- wd2
    	s0 <- w %*% rep(1, n)
        s1 <- (w * wd1) %*% rep(1, n)
        s2 <- (w * wd1^2) %*% rep(1, n)
        w  <- w * (matrix(rep(s2, n), ncol = n) - 
                 wd1 * matrix(rep(s1, n), ncol = n))
        w  <- w / (matrix(rep(s2, n), ncol = n) * 
                   matrix(rep(s0, n), ncol = n) - matrix(rep(s1, n), ncol = n)^2)
        }
    else {
        a11 <- w %*% rep(1, n)
        a12 <- (w * wd1) %*% rep(1, n)
        a13 <- (w * wd2) %*% rep(1, n)
        a22 <- (w * wd1^2) %*% rep(1, n)
        a23 <- (w * wd1 * wd2) %*% rep(1, n)
        a33 <- (w * wd2^2) %*% rep(1, n)
        d   <- a22 * a33 - a23^2
        b1  <- 1/(a11 - ((a12 * a33 - a13 * a23) * a12 + (a13 *
                  a22 - a12 * a23) * a13)/d)
        b2  <- (a13 * a23 - a12 * a33) * b1/d
        b3  <- (a12 * a23 - a13 * a22) * b1/d
        wt  <- matrix(rep(b1, n), ncol = n)
        wt  <- wt + matrix(rep(b2, n), ncol = n) * wd1
        wt  <- wt + matrix(rep(b3, n), ncol = n) * wd2
        w   <- wt * w
    }
    w
}


"sm.sigweight" <- function(x, weights) {

  if (is.vector(x)) {
  
    n   <- length(x)
    xx  <- sort(x)
    xx1 <- diff(xx)
    xx2 <- diff(xx, lag = 2)

    a <- xx1[-1]/xx2
    b <- xx1[-(n-1)]/xx2
    a[xx2==0] <- 0.5
    b[xx2==0] <- 0.5
    c <- sqrt(a^2/weights[1:(n-2)] + b^2/weights[3:n] +
                1/weights[2:(n-1)])

    A <- cbind(rep(0,n-2), diag(-1/c), rep(0,n-2)) +
            cbind(diag(a/c), rep(0,n-2), rep(0,n-2)) +
            cbind(rep(0,n-2), rep(0,n-2), diag(b/c))
    A <- rbind(rep(0,n), A, rep(0,n))
    
    A <- t(A) %*% A
    }
    
  if (is.matrix(x)) {
    
    x1 <- x[,1]
    x2 <- x[,2]
    n  <- length(x1)
    X  <- cbind(x1, x2)

    d  <- matrix(rep(x1,n),ncol=n)
    D  <- (d-t(d))^2
    d  <- matrix(rep(x2,n),ncol=n)
    D  <- sqrt(D / var(x1) + ((d-t(d))^2) / var(x2))
    nn <- function(x) {sort(x)[5]}
    hw <- apply(D, 1, nn)/2
    h  <- c(sqrt(var(x1)), sqrt(var(x2)))
    S  <- sm.weight2(X, X, h, cross=T, weights=weights,
            options=list(h.weights=hw))

    A  <- (diag(n)-S)
    A  <- t(A) %*% diag(1/diag(A %*% t(A))) %*% A
    # A  <- A / sum(weights)
    
    }
    
  A

  }


"sig.trace" <- function (expn, hvec, ...) {

    opt <- sm.options(list(...))
    replace.na(opt, display, "line")
    expn.char <- paste(deparse(substitute(expn)), collapse = "")
    lead.char <- substring(expn.char, 1, nchar(expn.char) - 1)
    nvec <- length(hvec)
    pvec <- vector("numeric", length = nvec)
    for (i in 1:nvec) {
        extn.char <- paste(lead.char, ", h = ", as.character(hvec[i]), ")")
        result <- eval(parse(text = extn.char))
        pvec[i] <- result$p
        }
    if (!(opt$display == "none")) {
        plot(hvec, pvec, type = "l", ylim = c(0, max(pvec)),
            xlab = "Smoothing parameter, h", ylab = "p-value")
        if (max(pvec) >= 0.05)
            lines(range(hvec), c(0.05, 0.05), lty = 2)
        }
    invisible(list(h = hvec, p = pvec))
    }
