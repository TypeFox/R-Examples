## tsq <- scan(hh("datasets/tsq.dat"))
## tsq <- as.rts(tsq)
## if.R(s=attr(attr(tsq,"tspar"),"units") <- "days",
##      r={})

data(tsq)
tsq[96:100]   ## last 5 of n=100 observations

tsq.plot <- tsacfplots(tsq, main="tsq 1", lwd=1)
## trellis.device(file=hh("tser/figure/tsq1.ps"), postscript, horizontal=TRUE); strip.background0()
tsq.plot
## dev.off()
## export.eps(hh("tser/figure/tsq1.ps"))

tmp <- tsq.plot$acf.plots
acfs <- if.R(r={
               acfs <- t(sapply(tmp$panel.args, `[[`, "y"))
               dimnames(acfs) <- list(c("acf","pacf"), 0:(dim(acfs)[2]-1))
               acfs
             }, s=
             t(matrix(tmp$y, ncol=2,
                      dimnames=list(matrix(tmp$x, ncol=2)[,1], c("acf","pacf"))))
             )
round(acfs[,1:6], digits=3)

tsq.loop <- if.R(s=
               arma.loop(tsq, model=list(order=c(2,0,2)))
               ,r=
               arma.loop(tsq, order=c(2,0,2))
               )

tsq.loop

tsq.diag <-
  rearrange.diag.arma.loop(diag.arma.loop(tsq.loop, tsq))
tsq.diagplot <- tsdiagplot(armas=tsq.loop, ts.diag=tsq.diag, lwd=1)
## trellis.device(file=hh("tser/figure/tsq2.ps"), postscript, horizontal=TRUE); strip.background0()
tsq.diagplot
## dev.off()
## export.eps(hh("tser/figure/tsq2.ps"))
