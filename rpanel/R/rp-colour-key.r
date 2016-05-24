rp.colour.key <- function(cols, brks, par.mar = c(5, 0, 4, 3) + 0.1, natural = TRUE,
                          margin = FALSE)  {
   ngrid <- length(cols)
   xvec  <- rep(0, ngrid)
    if (length(brks) == 2)
      brks <- seq(brks[1], brks[2], length = ngrid + 1)
   else if (length(brks) != ngrid + 1)
      stop("inappropriate length of brks in rp.colour.key.")
   if (!natural) {
      zlim      <- c(0, ngrid)
      brks.orig <- brks
      brks      <- 0:ngrid
      yaxs      <- "i"
   }
   else {
      zlim <- range(brks)
      yaxs <- "r"
   }
   par(mar = par.mar, mgp = c(1.5, 0.2, 0), tcl = -0.2)
   xrange <- if (margin) c(-1, 1) else c(0, 1)
   plot(xrange, zlim, type = "n", axes = FALSE, xaxs = "i", yaxs = yaxs, xlab = " ", ylab = " ")
   if (natural)
      axis(4)
   else {
      ticks <- pretty(0:ngrid)
      lbls  <- as.character(signif(brks.orig))[match(ticks, 0:ngrid)]
      lbls[lbls == "Inf"] <- NA
      axis(4, at = ticks, labels = lbls)
   }
   nbrks <- length(brks)
   brks[c(1, nbrks)] <- par()$usr[3:4]
   rect(xvec, brks[-nbrks], xvec + 1, brks[-1], col = cols, border = NA)
   lines(c(0, 0, 1, 1, 0), brks[c(1, nbrks, nbrks, 1, 1)])
   invisible()  
   }
