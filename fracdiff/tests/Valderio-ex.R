library(fracdiff)

set.seed(1)
##  examples(fdSperio)
mem.long <- fracdiff.sim(1500, d = 0.3)
spm <- fdSperio(mem.long$series)
str(spm, digits=6)

set.seed(8)
##  examples(fdGPH)
mem.l2 <- fracdiff.sim(1024, d = 0.25)
fdGPH(mem.l2$series)

stopifnot(diffseries ( 1:20, d = 1) == c(-9.5, rep(1, 20-1)),
          diffseries(-10:10, d = 0) == -10:10)
set.seed(123)
## example(diffseries)
mem.l3 <- fracdiff.sim(80, d = 0.3)
mGPH <- fdGPH(mem.l3$series)
r <- diffseries(mem.l3$series, d = mGPH$d)
print(r, digits = 4)
print(acf(r)) # shouldn't show structure - ideally

cat("Time used: ", proc.time(),"\n") # for ``statistical reasons''
