##
## Calibration with badly-scaled initial weights (bug report by Takahiro Tsuchiya)
##
library(survey)
data <- data.frame(x=c(1,1,1,1,2,2,2,2,2,2), w=rep(10,10))
des <- svydesign(ids=~1, weights=~w, data=data)
des.c <- calibrate(des, ~factor(x), c(10000, 5000))
des.r <- calibrate(des, ~factor(x), c(10000, 5000), calfun='raking')
stopifnot(all.equal(svytotal(~factor(x), des.c), svytotal(~factor(x), des.r)))
