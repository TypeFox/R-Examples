## ----include=FALSE-------------------------------------------------------
require(knitr)
require(reval)
require(ggplot2)
require(dplyr)
opts_chunk$set(tidy=TRUE, message=FALSE)

## ----single--------------------------------------------------------------
# install.packages("rivr")
require(rivr)
myprofile = compute_profile(So = 0.001, n = 0.045, Q = 250, y0 = 2.5, 
  Cm = 1.486, g = 32.2, B = 100, SS = 0, stepdist = 50, totaldist = 3000)
head(myprofile)

## ----loop----------------------------------------------------------------
# loop through values of n
ns = seq(0.03, 0.06, by = 0.005)
results = vector("list", length=length(ns))
for(i in 1:length(ns)){
  results[[i]] = compute_profile(So = 0.001, n = ns[i], Q = 250, y0 = 2.5, 
    Cm = 1.486, g = 32.2, B = 100, SS = 0, stepdist = 50, totaldist = 3000)
  # add an identifier to the result
  results[[i]]["n"] = ns[i]
}
# combine outputs
combined = do.call(rbind.data.frame, results)

## ----reval-single--------------------------------------------------------
results = evalmany(compute_profile, n = seq(0.03, 0.06, by = 0.005),
  default.args = list(So = 0.001, Q = 250, y0 = 2.5, Cm = 1.486, g = 32.2, 
  B = 100, SS = 0, stepdist = 50, totaldist = 3000))

## ----plot-single, echo=FALSE, fig.width=10, fig.height=5, dpi=200--------
ggplot(results, aes(x = x, y = y, color = id)) + geom_line()

## ----reval-multi---------------------------------------------------------
results = evalmany(compute_profile, n = seq(0.03, 0.06, by = 0.005),
  So = seq(0.001, 0.0015, by = 0.00025), SS = seq(0, 6, by = 2),
  default.args = list(Q = 250, y0 = 2.5, Cm = 1.486, g = 32.2, B = 100, 
  stepdist = 50, totaldist = 3000), method = "permute", collate.id = "multi",
  clusters = 2, packages = "rivr")

## ----plot-curves, fig.width=10, dpi=200----------------------------------
require(ggplot2)
ggplot(results, aes(x = x, y = y, color = factor(n))) + geom_line() + 
  facet_grid(SS ~ So)  

## ----filter-data, results="hide"-----------------------------------------
require(dplyr)
filter(results, n == 0.045, So == 0.0015, SS == 2)

## ----first-five, echo=FALSE----------------------------------------------
kable(head(filter(results, n == 0.045, So == 0.0015, SS == 2)))

