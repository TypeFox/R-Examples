colour.key <- function(cols, brks, mar.adj = TRUE, natural = TRUE)  {
   ngrid <- length(cols)
   xvec  <- rep(0, ngrid)
   if (length(brks) == 2)
      brks <- seq(brks[1], brks[2], length = ngrid + 1)
   else if (length(brks) != ngrid + 1)
      stop("inappropriate length of brks in colour.key.")
   if (!natural) {
      zlim <- c(0, ngrid)
      brks.orig <- brks
      brks <- 0:ngrid
   }
   if (mar.adj)
      plot(c(0,0.1), range(brks), type = "n", xaxs = "i", yaxs = "i", 
                xaxt = "n", yaxt = "n", xlab = "", ylab = "", mar = c(5, 0, 4, 2))
   else
      plot(c(0,0.1), zlim, type = "n", xaxs = "i", yaxs = "i", 
                xaxt = "n", yaxt = "n", xlab = "", ylab = "")
   if (natural)
      axis(4)
   else {
      ticks <- pretty(0:ngrid)
      print(ticks)
      print(as.character(signif(brks.orig))[match(ticks, 0:ngrid)])
      axis(4, at = ticks, labels = as.character(signif(brks.orig))[match(ticks, 0:ngrid)])
   }
   rect(xvec, brks[-length(brks)], xvec + 1, brks[-1], col = cols, border = NA)    
   invisible()  
   }
