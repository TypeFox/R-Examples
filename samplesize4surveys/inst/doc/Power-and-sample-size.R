## ---- echo=FALSE, message=FALSE------------------------------------------
library(samplesize4surveys)

## ---- fig.show='hold', fig.width=7, fig.height=6, fig.retina=3-----------
ss4dpH(N = 1000, P1 = 0.5, P2 = 0.5, D=0.03, DEFF = 2, plot=TRUE)

## ---- message=FALSE, fig.show='hold', fig.width=7, fig.height=6, fig.retina=3----
b4dp(N = 1000, n = 873, P1 = 0.5, P2 = 0.5, D=0.03, DEFF = 2, plot=TRUE)

