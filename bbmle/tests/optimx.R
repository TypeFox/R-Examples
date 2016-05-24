library(bbmle)
old_opt <- options(digits=3)
require(optimx)
x <- 0:10
y <- c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
d <- data.frame(x,y)

## breaks, don't try this
## optimx(fn=Lfn,par=c(15,6),method="Rvmmin")

suppressWarnings(m1 <- mle2(minuslogl=y~dpois(lambda=ymax/(1+x/xhalf)),
     start=list(ymax=15,xhalf=6),data=d,
     optimizer="optimx",
           method=c("BFGS","Nelder-Mead","CG")))

## FIXME!! fails (although not with an error, because
##  errors are caught by profiling) due to npar now
## being restricted to >1 in optimx 2012.05.24 ...

suppressWarnings(head(as.data.frame(profile(m1))))
detach("package:optimx")
options(old_opt)
