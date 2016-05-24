###
###
### Ramsey & Silverman (2002) Applied Functional Data Analysis (Springer)
###
### ch. 3.  Nondurable Goods Index
###
library(fda)
##
## sec. 3.1.  Introduction
##
# pp. 41-42, Figure 3.1.  US monthly nondurable good manufacturing

plot(nondurables)

##
## sec. 3.2.  Transformation and smoothing
##
# p. 43, Figure 3.2.  US nondurable goods production, log scale

plot(log10(nondurables))

nondur1929 <- window(nondurables, 1929, 1930)
max1929 <- max(nondur1929)
max1929t <- which(nondur1929==max1929)

library(zoo)
(Max1929 <- index(nondur1929)[max1929t])
# function 'index' is in 'zoo'
text(Max1929, log10(max1929)+.01, "Stock market crash", 0, srt=90)

restrMoney <- window(nondurables, 1928, 1940)
maxRestr <- max(restrMoney)
maxRestr.t <- which(restrMoney==maxRestr)
(MaxRestr <- mean(index(restrMoney)[maxRestr.t]))
text(MaxRestr, log10(maxRestr)+.01,
     "Restriction of money supply", 0, srt=90)

endVietnam <- (1975+4/12)
endVN <- min(window(nondurables, 1975, 1976))
text(endVietnam, log10(endVN)-0.01, "End of Vietnam War", 1,
     srt=90)

(nondurFit <- lm(log10(nondurables)~index(nondurables)) )
abline(nondurFit, lty=3)

# p. 44, Figure 3.3.  log nondurable goods index 1964 to 1967

#  Fit smooth per sec. 3.6.
goodsbasis <- create.bspline.basis(rangeval=c(1919,2000),
                                   nbasis=979, norder=8)

LfdobjNonDur     = int2Lfd(4);

#goodsfdPar = fdPar(goodsbasis, LfdobjNonDur, lambda=1e-11)
#lognondursmth = smooth.basis(durtime, coredata(lognondur), goodsfdPar);
logNondurSm <- smooth.basisPar(argvals=index(nondurables),
                y=log10(coredata(nondurables)), fdobj=goodsbasis,
                Lfdobj=LfdobjNonDur, lambda=1e-11)

#str(lognondursmth)

nondur1964.1967 <- window(nondurables, 1964, 1967)

plot(log10(nondur1964.1967), type="p", axes=FALSE, xlab="Year",
     ylab=expression(paste(log[10], " nondurable goods index")) )
axis(2)
axis(1, 1964:1967)
axis(1, seq(1964, 1967, by=0.5), labels=FALSE)

#durtimefine = linspace(1964,1967,101);
durtimefine <- seq(1964, 1967, length=181)

#fit = eval.fd(durtimefine, lognondursmth);
logNondurSm1964.67 = eval.fd(durtimefine, logNondurSm$fd);
lines(durtimefine, logNondurSm1964.67)
abline(v=1965:1966, lty=2)

##
## sec. 3.3.  Phase-plane plots
##

#p. 45-46.  Figure 3.4.  Phase-plane plot of sin(2*pi*t)

sin. <- expression(sin(2*pi*x))
D.sin <- D(sin., "x")
D2.sin <- D(D.sin, "x")

op <- par(pty="s")
# square plot region so we get a circle not an ellipse
with(data.frame(x=seq(0, 1, length=46)),
     plot(eval(D.sin), eval(D2.sin), type="l",
          xlim=c(-10, 10), ylim=c(-50, 50),
          xlab="Velocity", ylab="Acceleration") )
pi.2 <- (2*pi)
#lines(x=c(-pi.2,pi.2), y=c(0,0), lty=3)
abline(h=0, lty="longdash")
pi.2.2 <- pi.2^2
lines(x=c(0,0), y=c(-pi.2.2, pi.2.2), lty="longdash")

text(c(0,0), c(-47, 47), rep("no kinetic, max potential", 2))
text(c(-8.5,8.5), c(0,0), rep("max kinetic\nno potential", 2))
par(op)

##
## sec. 3.4.  Nondurable goods cycles
##
# Applied Functional Data Analysis, p. 56, sec. 3.6
# "Our final choice for lambda was 10^(-9.5)
lam9.5 <- 10^(-9.5)
goodsfdPar9.5 = fdPar(goodsbasis, LfdobjNonDur, lam9.5);
#lognondursmth = smooth_basis(durtime, lognondur, goodsfdPar);
lognondursm9.5 = smooth.basis(durtime, coredata(lognondur), goodsfdPar9.5);

# p. 47, Figure 3.5.  Nondurable phase-plane plot, 1964
##*** Need to add xlim and ylim to the following
## to match the plots in the book.

phaseplanePlot(1964, logNondurSm$fd)

# pp. 48-49, Figure 3.6.  Nondurable phase-plane plots, 1929-1931
phaseplanePlot(1929, logNondurSm$fd)
phaseplanePlot(1930, logNondurSm$fd)
phaseplanePlot(1931, logNondurSm$fd)

# pp. 48, 50, Figure 3.7.  Nondurable phase-plane plots, 1937-1938, 1943
phaseplanePlot(1937, logNondurSm$fd)
phaseplanePlot(1938, logNondurSm$fd)
phaseplanePlot(1939, logNondurSm$fd)

# pp. 48, 51, Figure 3.8.  Nondurable phase-plane plots, 1974-76
phaseplanePlot(1974, logNondurSm$fd)
phaseplanePlot(1975, logNondurSm$fd)
phaseplanePlot(1976, logNondurSm$fd)

# pp. 52-53, Figure 3.9.  Nondurable phase-plane plots, 1996-1998
phaseplanePlot(1996, logNondurSm$fd)
phaseplanePlot(1997, logNondurSm$fd)
phaseplanePlot(1998, logNondurSm$fd)

# p. 53, Figure 3.10.  Nondurable phase-plane plot, 1997, larger scale
phaseplanePlot(1997, logNondurSm$fd)

##
## sec. 3.5.  What have we seen?
##
# All discussion, no data analysis in this section

##
## sec. 3.6.  Smoothing data for phase-plane plots
##
# All discussion, no data analysis in this section

