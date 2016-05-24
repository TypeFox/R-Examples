`specplot` <-
function (cddews, scaling = "by.level", arrangement = c(3, 
    3), page = TRUE, dataname = "Image", verbose = FALSE, display = "persp", 
    reset = TRUE, wtitle = TRUE) 
{
    now <- proc.time()[1:2]
    if (class(cddews) != "cddews") 
        stop("Sorry! The quantity supplied to this package must be of class cddews.\n")
    if (display != "image" && display != "persp") 
        stop("Can only display these spectra using image and persp!")
    if (verbose == TRUE) 
        cat("Computing the wavelet spectrum \ncoefficients\n")
    enwdata <- cddews$S
    switch <- cddews$structure
    J <- dim(enwdata)[1]/3
    nr <- dim(enwdata)[2]
    nc <- dim(enwdata)[3]
    vert <- array(0, dim = c(J, nr, nc))
    horiz <- array(0, dim = c(J, nr, nc))
    diag <- array(0, dim = c(J, nr, nc))
    if (switch == "level") {
            oldpar <- par(mfrow = arrangement, pty = "s", oma = c(0.1, 
                0.1, 0.1, 0.1), ask = page)
            for (i in 1:J) {
                msub1 <- enwdata[(3 * i - 2), , ]
                msub2 <- enwdata[(3 * i - 1), , ]
                msub3 <- enwdata[(3 * i), , ]
                if (scaling == "by.level") {
                  limits <- range(range(msub1), range(msub2), 
                    range(msub3))
                }
                if (scaling == "global") {
                  limits <- range(range(enwdata))
                }
                if (display == "image") {
                  xlabstr <- paste("Level", -i, "(vertical)")
                  image(msub1, xlab = xlabstr)
		  if(wtitle==TRUE){
                  title("cddews of", dataname)
                  }
                  xlabstr <- paste("Level", -i, "(horizontal)")
                  if(wtitle==TRUE){
                  title("cddews of", dataname)
                  }
                  image(msub2, xlab = xlabstr)
                  xlabstr <- paste("Level", -i, "(diagonal)")
                  image(msub3, xlab = xlabstr)
                  if(wtitle==TRUE){
                  title("cddews of", dataname)
                  }
                }
                if (display == "persp") {
                  zlabstr <- paste("Power")
                  persp(msub1, theta=30, phi=30,  xlab = "R", ylab = "S", zlab = zlabstr, 
                    zlim = limits)
                  titlestr <- paste("(Corrected, unsmoothed) DDEWS of ", 
                    dataname, "\n [Level", -i, "(vertical)]")
		  if (wtitle==TRUE){
                  title(titlestr)
		  }
                  persp(msub2, theta = 30, phi = 30, xlab = "R", ylab = "S", zlab = zlabstr, 
                    zlim = limits)
                  titlestr <- paste("(Corrected, unsmoothed) DDEWS of ", 
                    dataname, "\n [Level", -i, "(horizontal)]")
                  if (wtitle==TRUE){
                  title(titlestr)
                  }
                  persp(msub3, theta = 30, phi = 30, xlab = "R", ylab = "S", zlab = zlabstr, 
                    zlim = limits)
                  titlestr <- paste("(Corrected, unsmoothed) DDEWS of ", 
                    dataname, "\n[Level", -i, "(diagonal)]")
                  if (wtitle==TRUE){
                  title(titlestr)
                  }
                }
            }
            par(oldpar)
    }
    if (switch == "direction") {
            oldpar <- par(mfrow = arrangement, pty = "s", oma = c(0.1, 
                0.1, 0.1, 0.1), ask = page)
            for (i in 1:J) {
                msub1 <- enwdata[i, , ]
                msub2 <- enwdata[(J + i), , ]
                msub3 <- enwdata[(2 * J + i), , ]
                if (scaling == "by.level") {
                  limits <- range(range(msub1), range(msub2), 
                    range(msub3))
                }
                if (scaling == "global") {
                  limits <- range(range(enwdata))
                }
                if (display == "image") {
                  xlabstr <- paste("Level", -i, "(vertical)")
                  image(msub1, xlab = xlabstr)
                  if(wtitle==TRUE){
                  title("cddews of", dataname)
                  }
                  xlabstr <- paste("Level", -i, "(horizontal)")
                  if(wtitle==TRUE){
                  title("cddews of", dataname)
                  }
                  image(msub2, xlab = xlabstr)
                  xlabstr <- paste("Level", -i, "(diagonal)")
                  image(msub3, xlab = xlabstr)
                  if(wtitle==TRUE){
                  title("cddews of", dataname)
                  }
                }
                if (display == "persp") {
                  zlabstr <- paste("Power")
                  persp(msub1, theta = 30, phi = 30, xlab = "R", ylab = "S", zlab = zlabstr, 
                    zlim = limits)
                  titlestr <- paste(dataname, "\n [Level", -i, "(vertical)]")
                  if (wtitle==TRUE){
                  title(titlestr)
                  }
                  persp(msub2, theta = 30, phi = 30, xlab = "R", ylab = "S", zlab = zlabstr, 
                    zlim = limits)
                  titlestr <- paste(dataname, "\n [Level", -i, "(horizontal)]")
                  if (wtitle==TRUE){
                  title(titlestr)
                  }
                  persp(msub3, theta = 30, phi = 30, xlab = "R", ylab = "S", zlab = zlabstr, 
                    zlim = limits)
                  titlestr <- paste(dataname, "\n[Level", -i, "(diagonal)]")
                  if (wtitle==TRUE){
                  title(titlestr)
                  }
                }
            }
            par(oldpar)
    }
}

