filePrefix <- "G:/vg/vg0.2-2/R/"
source(paste(filePrefix,"vgChangePars.R",sep=""))
source(paste(filePrefix,"vgMode.R",sep=""))
source(paste(filePrefix,"ddvg.R",sep=""))
source(paste(filePrefix,"dvg.R",sep=""))
source(paste(filePrefix,"vgCalcRange.R",sep=""))
source(paste(filePrefix,"vgBreaks.R",sep=""))
source(paste(filePrefix,"pvg.R",sep=""))
source(paste(filePrefix,"rvg.R",sep=""))
source(paste(filePrefix,"ppvg.R",sep=""))
source(paste(filePrefix,"qvg.R",sep=""))
source(paste(filePrefix,"qqvg.R",sep=""))
source(paste(filePrefix,"vgCheckPars.R",sep=""))

filePrefix <- "G:/vg/vg0.2-2/unitTests/testFunction/graphical tests/"
setwd(paste(filePrefix,"graphOutputs",sep=""))
source(paste(filePrefix,"vgGraphsdpqr.R",sep=""))
source(paste(filePrefix,"vgGraphsfit.R",sep=""))
data(vgParam)
params = vgSmallShape
n <- 1000
method <- "Nelder-Mead"
hessian <- FALSE

# run dpqr group graphical tests
## graphs for dvg and ddvg
test.vgGraphdpqrddvg ()
## graphs for dvg, vgBreaks and vgCalcRange
test.vgGraphdpqrvgBreaks ()
## graph for pvg and qvg
test.vgGraphpdpqrvgqvg ()
## graph for pvg, qvg, ppvg and qqvg
test.vgGraphdpqrppvgqqvg ()
## graph for rvg
test.vgGraphdpqrrvg ()

# run fitting group graphical tests
test.vgGraphvgFitStart ()
test.vgGraphvgFit ()