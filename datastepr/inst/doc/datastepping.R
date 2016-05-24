## ---- message=FALSE------------------------------------------------------
library(dplyr)
library(datastepr)
library(magrittr)
library(knitr)

## ------------------------------------------------------------------------
step = dataStepClass$new()

## ------------------------------------------------------------------------
?dataStepClass()

## ------------------------------------------------------------------------
xFrame = data.frame(x = 0:9, group_id = 1:5)
kable(xFrame %>% arrange(group_id))

## ------------------------------------------------------------------------
yFrame = data.frame(y = c(-1, 1), group_id = 1)
kable(yFrame)

## ------------------------------------------------------------------------
stairs = function(...) {
  step$begin(environment())

  if (step$i == 1) step$set(yFrame, group_id)

  if (step$i > 1) lagx = x

  step$set(xFrame, group_id)
  
  if (step$i > 1) y = y + dydx*(x - lagx)
 
  dydx = x*y

  series_id = c(1, 2)
  
  step$output(list(
    result = y,
    type = "y",
    series_id = series_id))
  
  step$output(list(
    result = y^2,
    type = "y squared",
    series_id = series_id))
  
  step$end(stairs)
}

stairs()

## ------------------------------------------------------------------------
step$results %>%
  arrange(series_id, type) %>%
  kable

