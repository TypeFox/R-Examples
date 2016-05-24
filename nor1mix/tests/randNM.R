library(nor1mix)

## From: "Jenifer Boshes" <boshes@mathpost.la.asu.edu>
## Date: Sun, 2 Dec 2007 20:16:28 -0700

## ................

## The use of sample(x, ...) inside rnorMix() typically was wrong when
## x was of length one

obj <- norMix(mu = c(1,0), sigma = c(5,1), w = c(.5,.5))

for(ss in round(runif(1000, min=0,max=1000))) {
    set.seed(ss)
    stopifnot(length(r <- rnorMix(1,obj)) == 1)
}

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
