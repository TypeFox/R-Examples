## from Eric Weese
library(bbmle)
f <- function(x=2,a=1) x^2 - a
f.g <- function(x,a) 2*x
f.g2 <- function(x,a) c(2*x,0)
options(digits=3)
mle2(f,fixed=list(a=1))
mle2(f,gr=f.g,fixed=list(a=1))
mle2(f,gr=f.g2,fixed=list(a=1))
