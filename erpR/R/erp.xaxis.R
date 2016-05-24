# this function is basically a wrapper for axis with msectopoints embedded and with some little additional features.
# for finer tuning use directly axis (see also par).

erp.xaxis <- function (length.erp = NULL, startmsec=-200, endmsec = 1200, x.tick=seq(-200, 1200, 200), x.labels = x.tick, x.pos = NA,  x.outer = FALSE, x.font = NA, x.lty = "solid", x.lwd = 1, x.lwd.ticks = 1, x.col = NULL, x.col.ticks = NULL, x.hadj = NA, x.padj = NA, x.tcl = -0.5, x.tick.both= FALSE, x.cex = 1)
{

if (is.null(length.erp)){
	stop("the length of the erp waveform need to be specified.\n The axis is drawn taking into account this length and the two parameters starmsec and endmsec", call.=F)
}

# tick and labels positions
ticksandlab=msectopoints(x.tick, length.erp, startmsec, endmsec)

axis(1, at= ticksandlab, labels=x.labels, pos= x.pos, tick = x.tick, outer = x.outer, font = x.font, lty = x.lty, lwd = x.lwd, lwd.ticks = x.lwd.ticks, col = x.col, col.ticks = x.col.ticks, hadj =, x.hadj, padj = x.padj, tcl = x.tcl, cex.axis = x.cex)

if (x.tick.both==TRUE){
axis(1, at= ticksandlab, labels= FALSE, pos= x.pos, tick = x.tick, outer = x.outer, font = x.font, lty = x.lty, lwd = x.lwd, lwd.ticks = x.lwd.ticks, col = x.col, col.ticks = x.col.ticks, tcl = - x.tcl) # NOTICE the minus in tcl: it is for plotting in the other direction
	
}


}