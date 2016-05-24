## Nottingham Castle, mean monthly air temperature in degrees Fahrenheit
## January 1920 -- December 1939

library(HH)
library("MASS")  ## nottem is part of the MASS package

nottem.dataplot <- tsacfplots(nottem, lwd=1)
## trellis.device(file=hh("tser/figure/nottema.ps"), postscript, horizontal=TRUE); strip.background0()
nottem.dataplot
## dev.off()
## export.eps(hh("tser/figure/nottema.ps"))

nottem.diff.dataplot <- tsacfplots(diff(nottem), lwd=1)  ## not displayed in book
nottem.diff.dataplot

nottem.diff12.dataplot <- tsacfplots(diff(nottem, 12), lwd=1)
## trellis.device(file=hh("tser/figure/nottemb.ps"), postscript, horizontal=TRUE); strip.background0()
nottem.diff12.dataplot
## dev.off()
## export.eps(hh("tser/figure/nottemb.ps"))


nottem.loop <- if.R(s=
                    arma.loop(nottem,  list(list(order=c(2,0,2)),
                                       list(order=c(2,1,0), period=12)))
                    ,r=
                    arma.loop(nottem,  order=c(2,0,2),
                                       seasonal=list(order=c(2,1,0), period=12),
                              method="ML")
                    )
print(nottem.loop, digits=4)

nottem.diag <- rearrange.diag.arma.loop(diag.arma.loop(nottem.loop, nottem))
nottem.diagplot <- tsdiagplot(armas=nottem.loop, ts.diag=nottem.diag, lwd=1)
## trellis.device(file=hh("tser/figure/nottemc.ps"), postscript, horizontal=TRUE); strip.background0()
nottem.diagplot
## dev.off()
## export.eps(hh("tser/figure/nottemc.ps"))

nottem.loop[["1","0"]]

## if.R(s=detach("MASS"), r=detach("package:MASS"))
