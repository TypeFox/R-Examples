
sm.discontinuity <- function(x, y, h, hd, ...) {

#		A test for the presence of one or more discontinuities

   x.name <- deparse(substitute(x))
   if (isMatrix(x)) x.names <- dimnames(x)[[2]]
   y.name <- deparse(substitute(y))

   opt <- sm.options(list(...))
   data    <- sm.check.data(x = x, y = y, ...)
   x       <- data$x
   y       <- data$y
   n       <- data$nobs
   ndim    <- data$ndim
   opt     <- data$options
   if (ndim > 2) x <- x[, 1:2]

   replace.na(opt, display, "lines")
   replace.na(opt, se,      TRUE)
   replace.na(opt, band,    TRUE)
   replace.na(opt, test,    TRUE)
   replace.na(opt, col,     "black")
   replace.na(opt, df,      5)
   if (ndim == 1) {
      replace.na(opt, ngrid, 100)
      replace.na(opt, xlab,  x.name)
      replace.na(opt, ylab,  y.name)
      replace.na(opt, xlim,  range(x))
      replace.na(opt, ylim,  range(y))
      if (length(opt$lty) == 1) 
         opt$lty <- c(opt$lty, opt$lty + 1)
      }
   else {
      replace.na(opt, ngrid, 21)
      dimn <- x.names
      name.comp<-if(!is.null(dimn) & !all(dimn == "")) dimn
             else {if (!is.null(attributes(x)$names)) attributes(x)$names
             else outer(x.name, c("[1]", "[2]"), paste, sep = "")}
      replace.na(opt, xlab, name.comp[1])
      replace.na(opt, ylab, name.comp[2])
      replace.na(opt, xlim, range(x[, 1]))
      replace.na(opt, ylim, range(x[, 2]))
      }

   if(missing(h)) 
      h <- h.select(x, y, ...)
   doublesmooth <- TRUE
   if (missing(hd)) {
      if (ndim == 1) {
         hd <- h * sqrt(0.25)
         h  <- h * sqrt(0.75)
         }
      else {
         hd <- h * sqrt(0.5)
         h  <- h * sqrt(0.5)
         }
      }
   else if (all(hd == rep(0, ndim)))
      doublesmooth <- FALSE
         
   if (ndim == 1) 
      result <- sm.discontinuity.1d(x, y, h, hd, doublesmooth, opt)
   else 
      result <- sm.discontinuity.2d(x, y, h, hd, doublesmooth, opt)
   
   result$h  <- h
   if (doublesmooth) result$hd <- hd
   invisible(result)
   }

sm.discontinuity.1d <- function(x, y, h, hd, doublesmooth, opt) {

   y <- y[order(x)]
   x <- sort(x)
   n <- length(x)

#		Define z to be the mid-points of distinct x values.
#		Restrict z so that there are at least 5 observations,
#		over more than one design point, on each size of
#		every value.

	z    <- x[c(1, diff(x)) > 0]
	nz   <- length(z)
	z    <- (z[1:(nz-1)] + z[2:nz]) / 2
	nz   <- length(z)
	flag <- rep(T, nz)
	for (i in 1:nz) {
	   left    <- x[x < z[i]]
	   right   <- x[x > z[i]]
	   flag[i] <- (length(left) > 5 
		    & length(right) > 5
		    & length(diff(left)[diff(left)   > 0]) > 1
		    & length(diff(right)[diff(right) > 0]) > 1)
	}
	z  <- z[flag]
	nz <- length(z)

	ghat.left  <- vector("numeric", length = nz)
	ghat.right <- vector("numeric", length = nz)

    wd <- matrix(rep(z, rep(n, nz)), ncol = n, byrow = T)
    wd <- wd - matrix(rep(x, nz), ncol = n, byrow = T)
    w  <- exp(-.5 * (wd / h)^2)

	wl <-  w * (sign(wd) + 1) / 2
	s0 <-  wl          %*% rep(1, n)
	s1 <- (wl * wd)    %*% rep(1, n)
	s2 <- (wl * wd^2 ) %*% rep(1, n)
	wl <-  wl * (matrix(rep(s2, n), ncol = n) - wd * matrix(rep(s1, n), ncol = n))
	wl <-  wl / (matrix(rep(s2, n), ncol = n) * matrix(rep(s0, n), ncol = n)
			- matrix(rep(s1, n), ncol = n)^2)
	ghat.left <- wl %*% y

	wr <-  w * (1 - sign(wd)) / 2
	s0 <-  wr          %*% rep(1, n)
	s1 <- (wr * wd)    %*% rep(1, n)
	s2 <- (wr * wd^2 ) %*% rep(1, n)
	wr <-  wr * (matrix(rep(s2, n), ncol = n) - wd * matrix(rep(s1, n), ncol = n))
	wr <-  wr / (matrix(rep(s2, n), ncol = n) * matrix(rep(s0, n), ncol = n)
			- matrix(rep(s1, n), ncol = n)^2)
	ghat.right <- wr %*% y

	A    <- sm.sigweight(x, rep(1, length(x))) / (n - 2)
	w    <- wl - wr
    if (doublesmooth) {
       ws         <- sm.weight(z, z, hd)
       w          <- ws %*% w
       ghat.left  <- as.vector(ws %*% as.vector(ghat.left))
       ghat.right <- as.vector(ws %*% as.vector(ghat.right))
       }

	shat <- sqrt(as.vector(t(as.matrix(y)) %*% A %*% as.matrix(y)))
	s.e. <- as.vector(shat * sqrt((w^2) %*% rep(1, n)))
	ts   <- sum(((ghat.left - ghat.right) / s.e.) ^2 )
	A    <- t(w) %*% diag((shat / s.e.)^2) %*% w - A * ts
	p    <- p.quad.moment(A, diag(n), 0, 0)

	if (opt$display != "none") {
	   if (!opt$add)
	      plot(x, y, xlab = opt$xlab, ylab = opt$ylab, xlim = opt$xlim, ylim = opt$ylim, type = "n")
	   av   <- (ghat.left + ghat.right) / 2
	   if (opt$band & opt$se)
	      polygon(c(z, rev(z)), c(av + s.e., rev(av - s.e.)), col = opt$col.band, border = FALSE)
	   lines(z, ghat.left,  lty = opt$lty[1])
	   lines(z, ghat.right, lty = opt$lty[2])
	   points(x, y, col = opt$col.points, pch = opt$pch)
	   }

    if (opt$verbose > 0)
        cat("Test of continuity:  significance = ", round(p, 3), "\n")
    st.diff <- (ghat.left - ghat.right)/ s.e.
    diffmat <- cbind(z, round(st.diff, 2))[abs(st.diff) > 2.5,]
    #  The following line forces a matrix when there is only one row in diffmat.
    if (!is.matrix(diffmat)) diffmat <- matrix(diffmat, ncol = 2)
    if ((opt$verbose > 0) & (nrow(diffmat) > 0)) {
       cat("location  st.diff\n")
       for (i in 1:nrow(diffmat)) cat(diffmat[i, ], "\n")
       }
	
	invisible(list(p = p, sigma = shat, eval.points = z, st.diff = st.diff, diffmat = diffmat))

   }

sm.discontinuity.2d <- function(x, y, h, hd, doublesmooth, opt, 
              nangles = 4, trim = 1/6, hull = FALSE) {

   #	Discontinuity detection with two covariates

   n           <- nrow(x)
   del1        <- diff(range(x[,1])) * trim
   del2        <- diff(range(x[,2])) * trim
   x1grid      <- seq(min(x[,1]) + del1, max(x[,1]) - del1, length = opt$ngrid)
   x2grid      <- seq(min(x[,2]) + del2, max(x[,2]) - del2, length = opt$ngrid)
   ev.points   <- cbind(x1grid, x2grid)
   replace.na(opt, eval.points, ev.points)
   eval.points <- opt$eval.points
   ngrid       <- nrow(eval.points)
   weights     <- rep(1, n)

   wd1   <- matrix(rep(eval.points[, 1], n), ncol = n)
   wd1   <- wd1 - matrix(rep(x[, 1], ngrid), ncol = n, byrow = TRUE)
   wd2   <- matrix(rep(eval.points[, 2], n), ncol = n)
   wd2   <- wd2 - matrix(rep(x[, 2], ngrid), ncol = n, byrow = TRUE)
   w1    <- exp(-0.5 * (wd1 / h[1])^2)
   w1    <- w1 * matrix(rep(weights, ngrid), ncol = n, byrow = TRUE)
   w2    <- exp(-0.5 * (wd2 / h[2])^2)
   wy    <- matrix(rep(y, ngrid), ncol = n, byrow=TRUE)

   a11   <-  w1          %*% t(w2)
   a12   <- (w1 * wd1)   %*% t(w2)
   a13   <-  w1          %*% t(w2 * wd2)
   a22   <- (w1 * wd1^2) %*% t(w2)
   a23   <- (w1 * wd1)   %*% t(w2 * wd2)
   a33   <-  w1          %*% t(w2 * wd2^2)

   c1    <-  w1        %*% t(w2 * wy)
   c2    <- (w1 * wd1) %*% t(w2 * wy)
   c3    <-  w1        %*% t(w2 * wy * wd2)

   beta1 <- sm.regression.invert(a22,a12,a23,a11,a13,a33,c2,c1,c3)
   beta2 <- sm.regression.invert(a33,a23,a13,a22,a12,a11,c3,c2,c1)

   wmask <- matrix(1, nrow = ngrid, ncol = ngrid)

   if (hull) { 
     hull.points <- x[order(x[,1], x[,2]),]
     dh          <- diff(hull.points)
     hull.points <- hull.points[c(TRUE, !((dh[,1] == 0) & (dh[,2] == 0))),]
     hull.points <- hull.points[chull(hull.points), ]
     nh          <- nrow(hull.points)
     gstep       <- matrix(rep(eval.points[2, ] - eval.points[1,], nh),   
                       ncol = 2, byrow = TRUE)
     hp.start    <- matrix(rep(eval.points[1, ], nh), ncol = 2, byrow = TRUE)
     hull.points <- hp.start + gstep * round((hull.points - hp.start) / gstep)
     hull.points <- hull.points[chull(hull.points), ]
     grid.points <- cbind(rep(eval.points[, 1], ngrid),
                          rep(eval.points[, 2], rep(ngrid, ngrid)))
     D      <- diff(rbind(hull.points, hull.points[1, ]))
     temp   <- D[, 1]
     D[,1]  <- D[, 2]
     D[,2]  <- -temp
     C      <- as.vector((hull.points * D) %*% rep(1, 2))
     C      <- matrix(rep(C, ngrid^2), nrow = ngrid^2, byrow = TRUE)
     D      <- t(D)
     wmask         <- ((grid.points %*% D) >= C)
     wmask         <- apply(wmask, 1, all)
     wmask[wmask]  <- 1
     wmask[!wmask] <- NA
     wmask         <- matrix(wmask, ncol = ngrid)
     }

   w1 <- w1[rep(1:ngrid, ngrid), ]
   w2 <- w2[rep(1:ngrid, each = ngrid), ]
   ind.select <- function(i, w1, w2, selection) {
   	  iset1 <- selection[i, ]
   	  iset1 <- iset1[!is.na(iset1)]
   	  iset2 <- !selection[i, ]
   	  iset2 <- iset2[!is.na(iset2)]
   	  (length(iset1) > 4) && (length(iset2) > 4)
         sum((w1[i, iset1] > exp(-2)) & (w2[i, iset1] > exp(-2)), na.rm = TRUE) > 4 &&
         sum((w1[i, iset2] > exp(-2)) & (w2[i, iset2] > exp(-2)), na.rm = TRUE) > 4
      }

   beta1 <- beta1 * wmask
   beta2 <- beta2 * wmask
   var.angle <- atan2(beta2, beta1) + pi / 2
   var.angle <- as.vector(var.angle)

   sig2 <- sm.sigma(x, y)
   A    <- sig2$qmat
   shat <- sig2$estimate
   ts   <- 0
   B    <- matrix(0, nrow = n, ncol = n)

   for (ang in ((1:nangles) * pi / nangles)) {

      angle <- var.angle
      angle <- rep(ang, ngrid^2) * angle / angle
      ev.points      <- matrix(0, nrow=ngrid^2, ncol = 2)
      ev.points[, 1] <- rep(eval.points[, 1], ngrid)
      ev.points[, 2] <- rep(eval.points[, 2], rep(ngrid, ngrid))
      selection <- matrix(rep(cos(angle), n), ncol = n) *
               (matrix(rep(x[, 2], ngrid^2), ncol = n, byrow = TRUE)
               - matrix(rep(ev.points[, 2], n), ncol = n))
      selection <- selection - matrix(rep(sin(angle), n), ncol = n) *
               (matrix(rep(x[, 1], ngrid^2), ncol = n, byrow = TRUE)
               - matrix(rep(ev.points[, 1], n), ncol = n))

      selection <- (selection > 0)
      ind <- sapply(1:nrow(selection), ind.select, 
                w1 = w1, w2 = w2, selection = selection)
      ev.points <- ev.points[ind, ]
      selection <- selection[ind, ]

      selection[selection >  0] <-  1
      selection[selection <= 0] <- -1
      wl   <- sm.discon.weight2(x, ev.points, h, (1 + selection) / 2)
      wr   <- sm.discon.weight2(x, ev.points, h, (1 - selection) / 2)
      w    <- wl - wr
      if (doublesmooth) w <- sm.weight2(ev.points, ev.points, hd) %*% w
      dhat <- as.vector(w %*% y)
      s.e. <- as.vector(shat * sqrt((w^2) %*% rep(1,n)))
      ts   <- ts + sum((dhat / s.e.) ^2 )
      B    <- B + t(w) %*% diag((shat/s.e.)^2) %*% w
      }

   C <- B - A * ts
   p <- p.quad.moment(C, diag(n), 0, 0)

   #	Calculations for the reference band

   angle     <- var.angle
   ev.points <- as.matrix(expand.grid(eval.points[, 1], eval.points[, 2]))
   selection <- matrix(rep(cos(angle), n), ncol = n) *
               (matrix(rep(x[, 2], ngrid^2), ncol = n, byrow = TRUE)
               -matrix(rep(ev.points[, 2], n), ncol = n))
   selection <- selection - matrix(rep(sin(angle), n), ncol = n) *
               (matrix(rep(x[, 1], ngrid^2), ncol = n, byrow = TRUE)
               -matrix(rep(ev.points[, 1], n), ncol = n))
   selection <- (selection > 0)
   ind <- sapply(1:nrow(selection), ind.select,
             w1 = w1, w2 = w2, selection = selection)
   ev.points <- ev.points[ind, ]
   selection <- selection[ind, ]
   
   selection[selection >  0]  <-  1
   selection[selection <= 0]  <- -1
   wl   <- sm.discon.weight2(x, ev.points, h, (1 + selection) / 2)
   wr   <- sm.discon.weight2(x, ev.points, h, (1 - selection) / 2)
   w    <- wl - wr
   if (doublesmooth) w <- sm.weight2(ev.points, ev.points, hd) %*% w
   dhat <- as.vector(w %*% y)
   s.e. <- as.vector(shat * sqrt((w^2) %*% rep(1, n)))
   std  <- rep(NA, ngrid * ngrid)
   std[ind] <- dhat / s.e.
   std <- matrix(abs(std), ncol = ngrid)
   
   results <- list(p = p, sigma = shat, eval.points = eval.points,
	         	st.diff = std, angle = matrix(angle, ncol = ngrid))

   if (opt$display != "none") {
      if (!opt$add)
         plot(x[, 1], x[, 2], xlab = opt$xlab, ylab = opt$ylab, 
                xlim = opt$xlim, ylim = opt$ylim, pch = opt$pch, col = opt$col.points)
      mx <- max(std, na.rm = TRUE)
      if (mx > 2.5)
         contour(eval.points[, 1], eval.points[, 2], std, 
              levels = seq(2.5, mx, by = 0.5),
              lty = opt$lty, col = opt$col, add = TRUE)
      }

   if (opt$verbose)
      cat(paste("Test of continuity:  significance = ", round(p, 3)), "\n")

   invisible(results)

   }

#-------------------------------------------------------------

sm.regression.invert <- function(a11,a12,a13,a22,a23,a33,c1,c2,c3) {
   #	Creates local linear intercept or slopes with two covariates
   d     <- a22 * a33 - a23^2
   b1    <- 1 / (a11 - ((a12*a33 - a13*a23)*a12 + 
     			(a13*a22 - a12*a23)*a13)/d)
   b2    <- (a13*a23 - a12*a33) * b1 / d
   b3    <- (a12*a23 - a13*a22) * b1 / d
   est   <- b1 * c1 + b2 * c2 + b3 * c3
   invisible(est)
   }

#-------------------------------------------------------

sm.discon.weight2 <- function(x, eval.points, h, selection,
                                 weights = rep(1, nrow(x)))  {

#	Amended version of sm.weight2 which uses different points 
#	at each grid position

   n  <- nrow(x)
   ne <- nrow(eval.points)

   wd1 <- matrix(rep(eval.points[,1], rep(n, ne)), ncol = n, byrow = TRUE)
   wd1 <- wd1 - matrix(rep(x[,1],           ne),   ncol = n, byrow = TRUE)
   w   <- exp(-.5 * (wd1 / h[1])^2)
   wd2 <- matrix(rep(eval.points[,2], rep(n, ne)), ncol = n, byrow = TRUE)
   wd2 <- wd2 - matrix(rep(x[,2],            ne),  ncol = n, byrow = TRUE)
   w   <- w * exp(-.5 * (wd2 / h[2])^2)
   w   <- w * matrix(rep(weights, ne), ncol = n, byrow = TRUE)
   w   <- w * selection
     
   a11 <-  w              %*% rep(1, n)
   a12 <- (w * wd1      ) %*% rep(1, n)
   a13 <- (w * wd2      ) %*% rep(1, n)
   a22 <- (w * wd1^2    ) %*% rep(1, n)
   a23 <- (w * wd1 * wd2) %*% rep(1, n)
   a33 <- (w * wd2^2    ) %*% rep(1, n)

   d   <- a22 * a33 - a23^2
   
   b1  <- 1 / (a11 - ((a12*a33 - a13*a23)*a12 + (a13*a22 - a12*a23)*a13)/d)
   b2  <- (a13*a23 - a12*a33) * b1 / d
   b3  <- (a12*a23 - a13*a22) * b1 / d

   wt  <-      matrix(rep(b1, n), ncol = n)
   wt  <- wt + matrix(rep(b2, n), ncol = n) * wd1
   wt  <- wt + matrix(rep(b3, n), ncol = n) * wd2
   w   <- wt * w

   invisible(w)

   }
