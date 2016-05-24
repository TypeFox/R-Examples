"onesamp" <-
function(dset, x = "unsprayed", y = "sprayed", xlab = NULL, ylab = NULL, 
	dubious = NULL, conv = NULL, dig = 2)
{
	if(!is.null(conv))
		dset <- round(dset * conv, 1)
	xlabel <- xlab
	ylabel <- ylab
	if(is.null(xlabel))
		xlabel <- x
	if(is.null(ylabel))
		ylabel <- y
	xname <- x
	yname <- y
	xv <- dset[, xname]
	yv <- dset[, yname]
	omit <- is.na(xv) | is.na(yv)
	if(!is.null(dubious))
		omit[dubious] <- TRUE
	ylim <- range(c(xv[!omit], yv[!omit]))
	xlim <- ylim
	plot(dset[!omit, xname], dset[!omit, yname], pch = 1, lwd = 1, xlab = 
		xlabel, ylab = ylabel, xlim = xlim, ylim = ylim)
	if(sum(omit) != 0) {
		points(dset[omit, xname], dset[omit, yname], pch = 4)
	}
	abline(0, 1)
	xmid <- mean(par()$usr[1:2])
	ymid <- mean(par()$usr[3:4])
	chw <- par()$cxy[1]
	chh <- par()$cxy[2]
	d <- dset[!omit, yname] - dset[!omit, xname]
	dbar <- mean(d)
	se <- sqrt(var(d)/length(d))
	xpos <- xmid - 3.0 * chw
	ypos <- ymid - 3.0*chh
	lines(c(xpos, xpos), c(ypos - se/2, ypos + se/2), lwd = 2)
	lines(c(xpos - chw/4, xpos + chw/4), rep(ypos + se/2, 2), lwd = 2)
	lines(c(xpos - chw/4, xpos + chw/4), rep(ypos - se/2, 2), lwd = 2)
	text(xpos + chw/2, ypos, paste("SED =", format(round(se, dig))), adj = 
		0)
	abline(dbar, 1, lty = 2)
	n <- dim(dset)[1]
	if(sum(omit) > 0)
		sex <- sqrt(var(dset[ - omit, xname])/n)
	else sex <- sqrt(var(dset[, xname])/n)
	if(sum(omit) > 0)
		sey <- sqrt(var(dset[ - omit, yname])/n)
	else sey <- sqrt(var(dset[, yname])/n)
	axis(3, at = c(xmid - sex/2, xmid + sex/2), labels = FALSE)
	axis(4, at = c(ymid - sey/2, ymid + sey/2), labels = FALSE)
	mtext(side = 3, line = 0.75, text = paste("SE =", format(round(sex, dig
		))), at = xmid, adj = 0.5)
	mtext(side = 4, line = 0.75, text = paste("SE =", format(round(sey, dig
		))), at = ymid, adj = 0.5, srt = 90)
	if(sum(omit) > 0)
		cat("\n", yname, sqrt(var(dset[ - omit, yname])), sqrt(var(dset[
			 - omit, xname])), "\n")
	else cat("\n", yname, sqrt(var(dset[, yname])), sqrt(var(dset[, xname])
			), sqrt(var(d)), "\n")
	if(sum(omit) > 0)
		r <- cor(dset[ - omit, yname], dset[ - omit, xname])
	else r <- cor(dset[, yname], dset[, xname])
	topleft <- par()$usr[c(1, 4)] + c(chw/4,  - chh/3)
	mtext(side=3,line=0.15, paste("r =", format(round(r, 4))), 
adj = 1.1, cex=0.75
		)
	print(t.test(d))
	invisible()
}
