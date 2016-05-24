## ----include = FALSE-----------------------------------------------------
library("papeR")
opts_chunk$set(message = FALSE, warning = FALSE, echo = TRUE, results = 'asis')
data(Orthodont, package = "nlme")

## ------------------------------------------------------------------------
library("papeR")
data(Orthodont, package = "nlme")
labels(Orthodont) <- c("fissure distance (mm)",
                       "age (years)", "Subject", "Sex")

## ------------------------------------------------------------------------
xtable(summarize(Orthodont, type = "numeric"))

## ------------------------------------------------------------------------
xtable(summarize(Orthodont, type = "numeric", group = "Sex"))

## ------------------------------------------------------------------------
xtable(summarize(Orthodont, type = "numeric", group = "Sex",
                 test = c("wilcox.test", "t.test")))

## ------------------------------------------------------------------------
xtable(summarize(Orthodont, type = "numeric", group = "Sex",
                 quantiles = FALSE))

## ------------------------------------------------------------------------
xtable(summarize(Orthodont, type = "factor", variables = "Sex"))

## ------------------------------------------------------------------------
print(xtable(summarize(Orthodont, type = "factor")),
      tabular.environment = "longtable")

## ------------------------------------------------------------------------
xtable(summarize(Orthodont, type = "factor", variables = "Sex",
                 cumulative = TRUE))

## ------------------------------------------------------------------------
Ortho_small <- subset(Orthodont,
                      Subject %in% c("M01", "M02", "F01", "F02"))
xtable(summarize(Ortho_small, type = "factor",
                 variables = "Subject", group = "Sex"))

## ------------------------------------------------------------------------
xtable(summarize(Ortho_small, type = "factor",
                 variables = "Subject", group = "Sex"),
       caption = "Example table for Fisher's exact test",
       label = "tab:Fisher")

