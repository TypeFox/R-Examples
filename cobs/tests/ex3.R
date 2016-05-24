suppressMessages(library(cobs))

options(digits = 6)
## pdf("ex3.pdf")

data(women) # 15 obs.
attach(women)

## Interpolation! very easy problem,   BUT :
try( ## gives  "ifl = 5" !!!!!
cobw <- cobs(weight, height, knots = weight, nknots = length(weight))
)

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
