## ------------------------------------------------------------------------
# test
library(Dark)
data(dark)
tmp<-dark
names(tmp)

## ----, fig.width=6, fig.height=6, fig.align='center'---------------------
par(las=1, bty='n',mfrow=c(1,1))
plot(tmp$time, tmp$thrs)

P<-Start(tmp, 2000)
MSC<-ModelSelect(tmp,P)
tmp<-BestFit(tmp,MSC, T)

## ----, fig.width=6, fig.height=6, fig.align='center'---------------------
par(las=1, bty='n')
tmp<-MultiStart(tmp,repeats = 25)
tmp<-BootDark(tmp,R = 150, T)


## ------------------------------------------------------------------------
tmp

