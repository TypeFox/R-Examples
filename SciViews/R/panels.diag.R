## Panel function to use on the diagonal (univariate graphs)
## TODO: define fill colors differently

## Boxplot
panel.boxplot <- function (x, col = par("col"), box.col = "cornsilk", ...)
{
    ## Note: col is defined here, but unused, because otherwise redefining
	## col would cause an error about duplicated 'col' arguments to boxplot()!
	## further arguments to boxplot are allowed (try notch = TRUE ... not very
	## useful here, but just for test). Note that warnings are generates in
	## pairs() in case of a call with non-graphic arguments, or even, col.box =
	par(new = TRUE)
    boxplot(x, axes = FALSE, col = box.col, horizontal = TRUE,
		xlim = c(0.5, 2), ...)
}

## Density plot
panel.density <- function (x, adjust = 1, rug = TRUE, col = par("col"),
lwd = par("lwd"), line.col = col, line.lwd = lwd,...)
{
	## Further arguments to density() are allowed (see examples) but it generates
	## warnings in pairs()
	dens.x <- density(x, adjust = adjust, ...)
    lines(dens.x$x, min(x) + dens.x$y * diff(range(x)) / diff(range(dens.x$y)),
		col = line.col, lwd = line.lwd)
    if (isTRUE(rug))
		points(x, rep(min(x), length(x)), pch = "|", col = line.col)
}

## Histogram
panel.hist <- function (x, breaks = "Sturges", hist.col = "cornsilk",
hist.border = NULL, hist.density = NULL, hist.angle = 45, ...)
{
	## Here, we try to define all arguments that are specific to the histogram
	## (col, border, density and angle) with specific arguments to allow better
	## control of the appearance of the histograms independently from the other
	## panels
	par(new = TRUE)
	hist(x, breaks = breaks, col = hist.col, border = hist.border,
		density = hist.density, angle = hist.angle, axes = FALSE,
		xlab = "", ylab = "", main = "")
}

## QQ-plot agains a Normal distribution
panel.qqnorm <- function(x, pch = par("pch"), col = par("col"), bg = par("bg"),
cex = par("cex"), lwd = par("lwd"), qq.pch = pch, qq.col = col, qq.bg = bg,
qq.cex = cex, qqline.col = qq.col, qqline.lwd = lwd, ...)
{
	par(new = TRUE)
    ylim <- range(x, na.rm = TRUE)
	## Leave enough space for name of variables on top of the graph
	ylim[2] <- ylim[2] + (ylim[2] - ylim[1]) / 4
	qqnorm(x, axes = FALSE, xlab = "", ylab = "", main = "",
		ylim = ylim, col = qq.col, bg = qq.bg, pch = qq.pch, cex = qq.cex)
    qqline(x, col = qqline.col, lwd = qqline.lwd, ...)
}
