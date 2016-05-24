#Loading package
library(R0)

## Outbreak during 1918 influenza pandemic in Germany)
data(Germany.1918)
mGT<-generation.time("gamma", c(3, 1.5))
estR0<-estimate.R(Germany.1918, mGT, begin=1, end=27, methods=c("EG", "ML", "TD", "AR", "SB"), 
                  pop.size=100000, nsim=100)

attributes(estR0)
## $names
## [1] "epid"      "GT"        "begin"     "end"       "estimates"
## 
## $class
## [1] "R0.sR"

## Estimates results are stored in the $estimates object
estR0
## Reproduction number estimate using  Exponential Growth  method.
## R :  1.525895[ 1.494984 , 1.557779 ]
## 
## Reproduction number estimate using  Maximum Likelihood  method.
## R :  1.383996[ 1.309545 , 1.461203 ]
## 
## Reproduction number estimate using  Attack Rate  method.
## R :  1.047392[ 1.046394 , 1.048393 ]
## 
## Reproduction number estimate using  Time-Dependent  method.
## 2.322239 2.272013 1.998474 1.843703 2.019297 1.867488 1.644993 1.553265 1.553317 1.601317 ...
## 
## Reproduction number estimate using  Sequential Bayesian  method.
## 0 0 2.22 0.66 1.2 1.84 1.43 1.63 1.34 1.52 ...


## If no date vector nor date of first observation are provided, results are the same
## except time values in $t are replaced by index
