sm.monotonicity <- function(x, y, N = rep(1, length(y)), h, type = "continuous", ...) {

   #   A test of monotonicity with nonparametric regression

   if (!(type %in% c("continuous", "binomial")))
      stop("only continuous and binomial data can be handed.")

   x.name <- deparse(substitute(x))
   y.name <- deparse(substitute(y))

   if (isMatrix(x)) {
      cat("Warning: only the first covariate has been used.\n")
      x <- x[, 1]
      }
      
   opt <- sm.options(list(...))
   data    <- sm.check.data(x = cbind(x, N), y = y, ...)
   x       <- data$x[, 1]
   y       <- data$y
   N       <- data$x[, 2]
   n       <- data$nobs
   ndim    <- data$ndim
   opt     <- data$options

   replace.na(opt, display, "lines")
   replace.na(opt, nboot,   200)
   replace.na(opt, col,     "black")
   replace.na(opt, df,      5)
   replace.na(opt, ngrid,   100)
   replace.na(opt, xlab,    x.name)
   replace.na(opt, ylab,    y.name)
   replace.na(opt, xlim,    range(x))
   if (type == "continuous")
      replace.na(opt, ylim, range(y))
   else
      replace.na(opt, ylim, range(y / N))
   ngrid   <- opt$ngrid
   nboot   <- opt$nboot
   display <- opt$display

   if ((type == "continuous")) {
      if (missing(h))
         h <- h.select(x, y, ...)
      h0 <- h
      }

#....................Find boundary h value...............................

   n      <- length(y)
   r      <- (max(x) - min(x))
   hstart <- r / 50
   hend   <- r / 2
   if (type == "binomial") {
     hstart <- hstart * 2
     hend   <- hend * 2
     }
   hstep  <- (hend - hstart) / (ngrid - 1)
   h      <- hstart

   shape <- shapesmooth(x, y, N, h, ngrid, type)
   if (shape == "increasing" | shape == "decreasing") {
     cat("The test cannot be performed as the smooth curve is already",
         "monotonic at the smallest value of h.\n")
     return()
     }

   while (shape == "non-monotonic" & h < hend) {
     h     <- h + hstep
     shape <- shapesmooth(x, y, N, h, ngrid, type)
     }

   if (shape == "flat") {
     cat("The test cannot be performed as the only monotonic shape identified is flat.\n")
     return()
     }
   else if (shape == "non-monotonic") {
     stop("The test cannot be performed as the smooth curves are non-monotonic at all values of h.\n")
     return()
     }

   if (shape == "increasing") article <- "an"
      else                    article <- "a"
   if (opt$verbose > 0)
      cat("The smallest h which produces", article, 
    		shape,"curve is",signif(h, digits = 5), "\n")


   #............Find residuals from which to bootstrap (continuous case)........

   if (type=="continuous") {
     smres <- sm.regression(x, y, h0, eval.points = x, display = "none")
     e       <- y - smres$estimate
     e       <- e - mean(e)
     if (opt$verbose > 0) {
        cat("Standard deviation of estimated errors is:",  
   			signif(sqrt(var(e)), digits = 5), "\n")
        cat("Smoothing parameter used for estimation of residuals is:",
			signif(h0, digits = 5), "\n")
        }
     }


   #.............Bootstrap data from monotonic estimator....................


   if (type=="continuous") 
      sm <- sm.regression(x, y, h, eval.points = x, display = "none")
   else if (type=="binomial")
      sm <- sm.binomial(x, y, N, h, eval.points = x, display = "none")
   if (display != "none") {
      if (!opt$add) {
         if (type == "continuous") yy <- y
            else               yy <- jitter(y / N, amount = 0)
         plot(x, yy, xlab = opt$xlab, ylab = opt$ylab, xlim = opt$xlim, ylim = opt$ylim,
                  col = opt$col.points, pch = opt$pch)
         }
      }
   p <- 0
   for (i in 1:nboot) {
     if (type=="continuous") ystar <- sm$estimate + sample(e, size = n, replace = TRUE)
        else             ystar <- rbinom(n, N, sm$estimate)
     shsm <- shapesmooth(x, ystar, N, h, ngrid, type)
     if (shsm == "non-monotonic") p <- p + 1
     if (shsm == "missing") cat("Warning: shape not identifiable.\n")
     if (display != "none") {
       if (shsm == "non-monotonic") clr <- "red"
          else                      clr <- "grey"
       if (type == "continuous")
          sm.regression(x, ystar, h, ngrid=ngrid, col = clr, add = TRUE)
       else {
          a <- sm.binomial(x, ystar, N, h, ngrid=ngrid, display = "none")
          lines(a$eval.points, a$estimate, col = clr, lty = opt$lty)
          }
       }
     }
   p <- p / nboot
   if (opt$verbose > 0) cat("Test of monotonicity: p =", round(p, 3), "\n")

   results <- list(p = p, hcrit = h)
   if (type == "continuous") results$h <- h0
   invisible(results)

   }


#------------------------------------------------------------------------

shapesmooth <- function(x, y, N = rep(1, length(y)), h, ngrid, type) {

   #      Identifies the shape of a nonparametric regression function

   if (type == "continuous")
      sm <- sm.regression(x, y, h, ngrid = ngrid, display = "none")
   else
      sm <- sm.binomial(x, y, N, h, ngrid = ngrid, display = "none")
   d  <- diff(sm$estimate)
   nplus  <- length(d[d >  0])
   nminus <- length(d[d <  0])
   nzero  <- length(d[d == 0])
   nsm    <- length(d)
   shapeind <- "missing"
   if (nplus  == nsm) shapeind <- "increasing"
   if (nminus == nsm) shapeind <- "decreasing"
   if (nzero  == nsm) shapeind <- "flat"
   if (nplus > 0 & nminus > 0) shapeind <- "non-monotonic"
   shapeind

   }
