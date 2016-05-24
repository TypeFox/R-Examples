`metaCDW` <-
function (dat, method="weibull6", xstart = 55,
           xmin = 0, xmax = 365,
           minpeak = 0.1, mincut = 0.382,
           quantile = 0.05, symmetric = FALSE,
           p0 = NULL, linint = -1, findpeak = TRUE, maxit=2000) {
   if (method == "weibull4" & is.null(p0)){p0 <- c(0.1, 50, 5, 100)}
   if (is.null(dat$flag)) dat$flag <- TRUE
   dat$sample <- as.factor(dat$sample)
   if (!(length(xstart) %in% c(1, nlevels(dat$sample)))) stop("length of xstart does not match number of samples")
   datnew  <- dat

   xstartframe   <- data.frame(sample = levels(datnew$sample), xstart = xstart)
   datnew  <- subset(dat, dat$flag == TRUE)
   rows    <- nrow(datnew)
   setflag <- NULL
   # find the relevant peaks
   if (findpeak == TRUE) {
     for (i in levels(datnew$sample)) {
       dattmp  <- subset(datnew, sample == i)
       if (nrow(dattmp) < 3) stop(paste("less than 3 valid points in sample`", i,"'"))
       thexstart <- xstartframe$xstart[xstartframe$sample == i]
       peaks   <- peakwindow(dattmp$x, dattmp$y,
                             xstart = thexstart, xmax = xmax,
                             minpeak = minpeak, mincut = mincut)
       smd     <- peaks$smd.indices
       flagtmp <- rep(FALSE, nrow(dattmp))
       flagtmp[smd] <- TRUE
       setflag <- c(setflag, flagtmp)
     }
     datnew <- subset(datnew, setflag == TRUE)
   }
   # calculate CWD
   if (method == "weibull6") {fitweibull <- fitweibull6}
   else if (method == "weibull4") {fitweibull <- fitweibull4}
   else {stop("method unknown")}
   metares <- NULL
   res     <- NULL
   j       <- 0
   for (i in levels(datnew$sample)) {
     cat("processing sample", i, "\t")
     j       <- j + 1
     dattmp  <- subset(datnew, sample == i)
     x       <- dattmp$x
     y       <- dattmp$y

     res[j]  <- list(fitweibull(x, y, p0, linint, maxit=maxit))
     cat(ifelse(res[[j]]$convergence ==0, "converged, ", " --      , "))
     cat("r2 =", round(res[[j]]$r2, 4),"\n")
     ## ---> Testing code
     card <- tryCatch({
         CDW(res[[j]], xmin, xmax, quantile, symmetric)
       },
         error = function(ex) {
           #print(ex$message)
           cat("Warning: no useful Weibull function found (you may try manual fit).\n")
           ## testing defaults!!! change this !!!
           list(x = rep(0, 3),
                y = rep(0 ,3),
                p = rep(0, 6)
           )
         }
       )
     ## <--- end testing code

     metarestmp <- data.frame(
       sample = i,
       tMid   = card$x[[1]],
       tBegin = card$x[[2]],
       tEnd   = card$x[[3]],
       yMid   = card$y[1] * res[[j]]$ymax,
       yBegin = card$y[2] * res[[j]]$ymax,
       yEnd   = card$y[3] * res[[j]]$ymax,
       p1     = card$p[1],
       p2     = card$p[2],
       p3     = card$p[3],
       p4     = card$p[4]
     )
     if (method == "weibull6") {
       metarestmp$p5 <- card$p[5]
       metarestmp$p6 <- card$p[6]
     }
     metarestmp$r2 <- res[[j]]$r2
     metares <-  rbind(metares, metarestmp)
   }
   ret <- list(metares = metares, weibullfits = res)
   class(ret) <- c("list", "cardiMetacdw")
   ret
}

