"statgen" <-
function (valuelist, xgrid = seq(0, 1, length = 21), ygrid = seq(0, 
    1, length = 21), binsize = 32, plot.m = FALSE, plot.v = FALSE, 
    plot.ks = FALSE, ptype = "persp") 
{
    afmean <- matrix(0, length(xgrid), length(ygrid))
    anmean <- matrix(0, length(xgrid), length(ygrid))
	frmean<-anmean
    afvar <- matrix(0, length(xgrid), length(ygrid))
    anvar <- matrix(0, length(xgrid), length(ygrid))
	frvar<-anvar
    afks <- matrix(0, length(xgrid), length(ygrid))
    anks <- matrix(0, length(xgrid), length(ygrid))
frks<-anks
    af <- valuelist[[1]]
    an <- valuelist[[2]]
    fr <- valuelist[[3]]
    afmean <- apply(af, c(1, 2), mean)
    anmean <- apply(an, c(1, 2), mean)
frmean <- apply(fr, c(1, 2), mean)
    afvar <- apply(af, c(1, 2), var) 
    anvar <- apply(an, c(1, 2), var) 
    frvar <- apply(fr, c(1, 2), var) * 1 * (binsize + 0.5)
    for (i in 1:length(xgrid)) {
        for (j in 1:length(ygrid)) {
            afks[i, j] <- ks.test(af[i, j, ], y = "pnorm", mean =afmean[i, j], sd = sqrt(afvar[i, j]))$statistic
            anks[i, j] <- ks.test(an[i, j, ], "pnorm", mean =anmean[i, j], sd = sqrt(anvar[i,j]))$statistic
		frks[i, j] <- ks.test(fr[i, j, ], "pnorm", mean = frmean[i, j], sd = sqrt(afvar[i,j]))$statistic
        }
    }
    if (plot.m == TRUE) {
        getOption("device")()
        if (ptype == "persp") {
            persp(xgrid, ygrid, abs(afmean - asymean(xgrid, ygrid, 
                binsize)), xlab = "p_1", ylab = "p_2", zlab = "", 
                ticktype = "detailed", theta = -45, phi = 15, 
                scale = FALSE, zlim = c(0, 0.5))
        }
        if (ptype == "contour") {
            contour(xgrid, ygrid, abs(afmean - asymean(xgrid, 
                ygrid, binsize)), xlab = "p_1", ylab = "p_2")
        }
    }
    if (plot.v == TRUE) {
        getOption("device")()
        if (ptype == "persp") {
            persp(xgrid, ygrid, afvar, xlab = "p_1", ylab = "p_2", 
                zlab = "Var(zeta(X_1,X_2))", ticktype = "detailed", 
                theta = -30, phi = 15, scale = FALSE, zlim = c(0, 
                  1.2))
        }
        if (ptype == "contour") {
            contour(xgrid, ygrid, afvar, xlab = "p_1", ylab = "p_2")
        }
        getOption("device")()
        plot(xgrid[2:(length(xgrid) - 1)], diag(afvar)[2:(length(xgrid) - 
            1)], type = "l", ylim = c(0, 1.2), xlab = "binomial probability", 
            ylab = "variance")
        lines(xgrid[2:(length(xgrid) - 1)], diag(anvar)[2:(length(xgrid) - 
            1)], lty = 2)
	  lines(xgrid[2:(length(xgrid) - 1)], diag(frvar)[2:(length(xgrid) - 
            1)], lty = 3)
        getOption("device")()
        plot(xgrid[2:(length(xgrid) - 1)], (1 - diag(afvar)[2:(length(xgrid) - 
            1)])^2, type = "l", ylim = c(0, 1), xlab = "binomial probability", 
            ylab = "residual")
        lines(xgrid[2:(length(xgrid) - 1)], (1 - diag(anvar)[2:(length(xgrid) - 
            1)])^2, lty = 2)
	  lines(xgrid[2:(length(xgrid) - 1)], (1 - diag(frvar)[2:(length(xgrid) - 
            1)])^2, lty = 3)

    }
    if (plot.ks == TRUE) {
        getOption("device")()
       if (ptype == "persp") {
            persp(xgrid, ygrid, anks - afks, xlab = "p_1", ylab = "p_2", 
                zlab = "", ticktype = "detailed", theta = -30, 
                phi = 15, scale = FALSE)
	  getOption("device")()
		persp(xgrid, ygrid, frks - afks, xlab = "p_1", ylab = "p_2", 
                zlab = "", ticktype = "detailed", theta = -30, 
                phi = 15, scale = FALSE)
	getOption("device")()
		 persp(xgrid, ygrid, anks - frks, xlab = "p_1", ylab = "p_2",
                zlab = "", ticktype = "detailed", theta = -30,
                phi = 15, scale = FALSE)

        }
        if (ptype == "contour") {
	getOption("device")()
            contour(xgrid, ygrid, anks - afks, xlab = "p_1", 
                ylab = "p_2")
	getOption("device")()
		contour(xgrid, ygrid, frks - afks, xlab = "p_1", 
                ylab = "p_2")
	getOption("device")()
	contour(xgrid, ygrid, anks - frks, xlab = "p_1",
                ylab = "p_2")
        }
    }
    return(list(afm = afmean, anm = anmean, frm=frmean, afv = afvar, anv = anvar, frv=frvar, 
        afk = afks, ank = anks, frk=frks))
}

