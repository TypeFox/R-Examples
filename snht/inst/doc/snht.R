## ----include=FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(
concordance=TRUE
)

## ----setup, include=FALSE, cache=FALSE, echo=FALSE----------------------------
library(knitr)
opts_chunk$set(fig.path = 'figure/', fig.align = 'center', fig.show = 'hold',
               warning = FALSE, message = FALSE, error = FALSE, tidy = FALSE,
               results = 'markup', eval = TRUE, echo = TRUE, cache = FALSE,
               fig.height = 4)
options(replace.assign = TRUE, width = 80)

## -----------------------------------------------------------------------------
set.seed(123)
baseData = rnorm(1000)
baseData[201:500] = baseData[201:500] + .4
baseData[501:600] = baseData[501:600] - .6

## ----echo=FALSE---------------------------------------------------------------
suppressWarnings(library(ggplot2))
p1 = qplot(1:1000, baseData) +
    geom_segment(aes(x = 0, xend = 200, y = 0, yend = 0),
                 color = "red", size = 1) +
    geom_segment(aes(x = 201, xend = 500, y = 0.4, yend = 0.4),
                 color = "red", size = 1) +
    geom_segment(aes(x = 501, xend = 600, y = -0.6, yend = -0.6),
                 color = "red", size = 1) +
    geom_segment(aes(x = 601, xend = 1000, y = 0, yend = 0),
                 color = "red", size = 1) +
    labs(x = "Time", y = "Data")
p1

## -----------------------------------------------------------------------------
library(snht)
snhtStatistic30 = snht(data = baseData, period = 30)
summary(snhtStatistic30)
snhtStatistic60 = snht(data = baseData, period = 60)
summary(snhtStatistic60)

## -----------------------------------------------------------------------------
plotSNHT(data = baseData, stat = snhtStatistic30, alpha = .05)

## -----------------------------------------------------------------------------
largestStatTime = which.max(snhtStatistic60$score)
snhtStatistic60[largestStatTime, ]

## -----------------------------------------------------------------------------
plotSNHT(data = baseData, stat = snhtStatistic30, alpha = 0.05)
plotSNHT(data = baseData, stat = snhtStatistic60, alpha = 0.05)

## -----------------------------------------------------------------------------
seasonalData = baseData + cos(1:200 * 2 * pi / 200)
seasonalData = seasonalData +
    rbinom(1000, p = .1, size = 1) * rnorm(1000, sd = 10)
qplot(1:1000, seasonalData) + labs(x = "Time", y = "Seasonal Data")

## -----------------------------------------------------------------------------
snhtStatistic = snht(data = seasonalData, period = 200)
plotSNHT(data = seasonalData, stat = snhtStatistic, alpha = 0.05)

## -----------------------------------------------------------------------------
snhtStatistic = snht(data = seasonalData, period = 200, robust = TRUE,
                     rmSeasonalPeriod = 200)
plotSNHT(data = seasonalData, stat = snhtStatistic, alpha = .05)

## -----------------------------------------------------------------------------
times = 1:60 + rnorm(60, sd = 3)
times = sort(times)
data = rnorm(60) + c(rep(0, 30), rep(1, 30))
snhtStatistic = snht(data = data, period = 5, time = times)
summary(snhtStatistic)

## -----------------------------------------------------------------------------
plotSNHT(data = data, stat = snhtStatistic, time = times, alpha = .05)

