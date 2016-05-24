magplot <-
function(x, y, log='', main='', side = 1:2, majorn = 5, minorn = 5, tcl = 0.5, ratio = 0.5, labels = TRUE, unlog = "Auto", mgp = c(2,0.5,0), mtline = 2, xlab = NULL, ylab = NULL, crunch = TRUE, logpretty = TRUE, prettybase = 10, hersh = FALSE, family = "sans", frame.plot=TRUE, usepar=FALSE, ...){
if(missing(y)) plot(x, axes=F, xlab='', ylab='', main=main, log=log, frame.plot=FALSE, ...)
else plot(x, y, axes=F, xlab='', ylab='', main=main, log=log, frame.plot=FALSE, ...)

if(side[1] !=FALSE){magaxis(side = side, majorn = majorn, minorn = minorn, tcl = tcl, ratio = ratio, labels = labels, unlog = unlog, mgp = mgp, mtline = mtline, xlab = xlab, ylab = ylab, crunch = crunch, logpretty = logpretty, prettybase = prettybase, hersh = hersh, family = family, frame.plot = frame.plot, usepar = usepar)}
}
