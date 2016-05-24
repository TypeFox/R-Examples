## ----install, eval=FALSE-------------------------------------------------
#  install.packages("sweidnumbr")

## ----test, message=FALSE, warning=FALSE, eval=TRUE-----------------------
library(sweidnumbr)

## ----locale, eval=FALSE--------------------------------------------------
#  Sys.setlocale(locale="UTF-8")

## ----example1, message=FALSE, eval=TRUE----------------------------------
example_pin <- c("640823-3234", "6408233234", "19640823-3230")
example_pin <- as.pin(example_pin)
example_pin

## ----example2, message=FALSE, eval=TRUE----------------------------------
is.pin(example_pin)

## ----example3, message=FALSE, eval=TRUE----------------------------------
pin_ctrl(example_pin)

## ----example4, message=FALSE, eval=TRUE----------------------------------
pin_sex(example_pin)
pin_birthplace(example_pin)

## ----example5, message=FALSE, eval=TRUE----------------------------------
pin_age(example_pin)
pin_age(example_pin, date = "2000-01-01")

## ----example6, message=FALSE, eval=TRUE----------------------------------
format_pin(example_pin, "%Y-%m-%d-%N")
format_pin(example_pin, "%P")

## ----oin1, message=FALSE, eval=TRUE--------------------------------------
example_oin <- c("556000-4615", "232100-0156", "802002-4280")
example_oin <- as.oin(example_oin)
example_oin

## ----oin2, message=FALSE, eval=TRUE--------------------------------------
is.oin(example_oin)

## ----oin3, message=FALSE, eval=TRUE--------------------------------------
oin_ctrl(example_oin)

## ----oin4, message=FALSE, eval=TRUE--------------------------------------
oin_group(example_oin)

## ----citation, message=FALSE, eval=TRUE----------------------------------
citation("sweidnumbr")

## ----sessioninfo, message=FALSE, warning=FALSE---------------------------
sessionInfo()

