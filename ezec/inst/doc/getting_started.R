## ---- echo = FALSE, message = FALSE, warning = FALSE---------------------
knitr::opts_chunk$set(message = FALSE, warning = FALSE, tidy = FALSE)
knitr::opts_chunk$set(fig.align = "center", fig.show = 'asis', fig.height = 5,
                      fig.width = 5)
options(out.width = 100)

## ------------------------------------------------------------------------
packageVersion("drc")
citation("drc")

## ---- dummydata----------------------------------------------------------
library("ezec")
data("dummydata", library = "ezec")
head(dummydata) # the function head means "look at the top of the object"

## ---- ec_table-----------------------------------------------------------
library("ezec")
data(dummydata)
res <- EC_table(dummydata, form = response ~ dose)
print(res)

## ---- ec_table_par, fig.width = 7, fig.height = 4------------------------
par(mfrow = c(1, 2)) # set window to have 1 row and two columns
EC_table(dummydata, form = response ~ dose)
par(mfrow = c(1, 1)) # reset the window

## ------------------------------------------------------------------------
EC_table(dummydata, form = response ~ dose, plot = FALSE, result = "summary")

## ------------------------------------------------------------------------
EC_table(dummydata, form = response ~ dose, plot = FALSE, result = "model")

