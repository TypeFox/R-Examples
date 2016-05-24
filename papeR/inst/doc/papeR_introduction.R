## ----initialization, include=FALSE---------------------------------------
library("knitr")
options(knitr.table.format = 'markdown')
opts_chunk$set(message = FALSE, warning = FALSE, fig.path = "graphics/fig-",
               out.width = "480px")
set.seed(1234)

## ----install_CRAN, eval = FALSE------------------------------------------
#  install.packages("papeR")

## ----install_github, eval=FALSE------------------------------------------
#  install.packages("devtools")
#  library("devtools")
#  install_github("hofnerb/papeR")

## ----load_pkg------------------------------------------------------------
library("papeR")

## ------------------------------------------------------------------------
data(Orthodont, package = "nlme")
## keep the original data set for later use
Orthodont_orig <- Orthodont

## ------------------------------------------------------------------------
is.ldf(Orthodont)

## ------------------------------------------------------------------------
labels(Orthodont)

## ------------------------------------------------------------------------
labels(Orthodont) <- c("fissure distance (mm)", "age (years)", "Subject", "Sex")

## ------------------------------------------------------------------------
is.ldf(Orthodont)
class(Orthodont)

## ------------------------------------------------------------------------
labels(Orthodont)

## ------------------------------------------------------------------------
## set labels for distance and age
labels(Orthodont, which = c("distance", "age")) <- c("Fissure distance (mm)", "Age (years)")
## extract labels for age only
labels(Orthodont, which = "age")
## or for the first two variables (i.e., distance and age)
labels(Orthodont, which = 1:2)

## ------------------------------------------------------------------------
Orthodont2 <- convert.labels(Orthodont_orig)
class(Orthodont2)
labels(Orthodont2)

## ----plot_labeled_dataframe----------------------------------------------
par(mfrow = c(2, 2))
plot(Orthodont)

## ----grouped_plot--------------------------------------------------------
par(mfrow = c(1, 3))
plot(Orthodont, by = "Sex")

## ----with_x--------------------------------------------------------------
par(mfrow = c(1, 3))
plot(Orthodont, with = "distance")

## ----univariate_no_regressionline----------------------------------------
par(mfrow = c(1, 2))
plot(Orthodont, variables = -3, with = "distance", regression.line = FALSE)

## ------------------------------------------------------------------------
data(Orthodont, package = "nlme")
summarize(Orthodont, type = "numeric")

## ------------------------------------------------------------------------
summarize(Orthodont, type = "factor", variables = "Sex")

## ------------------------------------------------------------------------
summarize(Orthodont, type = "numeric", group = "Sex", test = FALSE)

## ------------------------------------------------------------------------
summarize(Orthodont, type = "numeric", group = "Sex")

## ---- echo = FALSE, results = 'asis'-------------------------------------
library("knitr")
kable(summarize(Orthodont, type = "numeric"))
kable(summarize(Orthodont, type = "factor", variables = "Sex", cumulative = TRUE))
kable(summarize(Orthodont, type = "numeric", group = "Sex"))

## ------------------------------------------------------------------------
linmod <- lm(distance ~ age + Sex, data = Orthodont)
## Extract pretty summary
(pretty_lm <- prettify(summary(linmod)))

## ---- results='hide'-----------------------------------------------------
xtable(pretty_lm)

## ---- results='asis'-----------------------------------------------------
kable(pretty_lm)

