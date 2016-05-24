#Loading package
library(R0)

# Data taken from traced cases of H1N1 viruses.
data(H1N1.serial.interval)
est.GT(serial.interval=H1N1.serial.interval)

## Best fitting GT distribution is a gamma distribution with mean = 3.039437 and sd = 1.676551 .
## Discretized Generation Time distribution
## mean: 3.070303 , sd: 1.676531 
## [1] 0.0000000000 0.1621208802 0.2704857362 0.2358751176 0.1561845680 0.0888997193 0.0459909903 
## 0.0222778094 0.0102848887 0.0045773285 0.0019791984 0.0008360608 0.0003464431 0.0001412594


# The same result can be achieved with two vectors of dates of onset.
# Here we use the same data, but trick the function into thinking onset dates are all "0".
data(H1N1.serial.interval)
est.GT(infector.onset.dates=rep(0,length(H1N1.serial.interval)), 
       infectee.onset.dates=H1N1.serial.interval)
