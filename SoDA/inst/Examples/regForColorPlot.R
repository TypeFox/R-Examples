## a slightly silly artificial regression problem to illustrate using
## color to code by size of residuals.

N = 150
set.seed(402)

x1 = rnorm(N, 15, 5)

x2 = rnorm(N, 15, 5)

y = .3 * x1 + 1.* (x1+x2)^1.5 + .1 * x2 + rnorm(N)

yreg = lm(y ~ x1+x2)

require(RColorBrewer)
require(lattice)

resids <- yreg$residuals

residGroups = cut(resids,  quantile(resids, seq(0,1,len=6)))

## using plotting character
pdf(file = "Examples/regResidPch.pdf", width =4, height = 4)

xyplot(x2~x1, groups = residGroups, pch = c("X", "=", "=", "=", "O"), col="black", cex=1.5)

cols = brewer.pal(5, "RdYlGn")

resCols <- cols[residGroups]

dev.off()

{
pdf(file = "Examples/regForColorPlot.pdf", width =
7.5, height =10)

par(mfcol=c(2,1), pty = "s")
myLwd = 1.5
myCex = 1.25

plot(x1, x2, col = resCols, main = "Color Brewer", pty = "s", lwd = myLwd, cex = myCex)

redGreen <- colorRampPalette(c("red", "blue","green"))

par(mar = c(7.1, 4.1, 2.1, 2.1))

plot(x1, x2, col = redGreen(5)[residGroups], main = "Color Ramp",
     xlab = "", sub = "Plate 1:  Coding residuals by color: \nGreen for large positive, red for large negative. " , pty = "s", lwd =myLwd, cex = myCex)

system("cp Examples/regForColorPlot.pdf Plate1.pdf; mv Plate1.pdf Chambers-color*")
system("cp Examples/regForColorPlot.R Chambers-color*")

dev.off()
}

## trellis.device("pdf", file = "Examples/regForColorLattice.pdf", width = 7.5, height = 10)

## ## pars <- trellis.par.get("plot.symbol")
## ## pars$col <- redGreen(5)
## ## trellis.par.set("plot.symbol", pars)

## xyplot(x2~ x1, cex=1.5, lwd=2,
##        groups = residGroups, col = redGreen(5), aspect = 1. , sub ="Plate 2:  Coding residuals by color\n with xyplot()" )

## # draw.key(simpleKey(levels(residGroups), col = redGreen(5)))
## dev.off()
