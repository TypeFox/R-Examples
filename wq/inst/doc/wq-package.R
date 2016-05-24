## ----setup, echo=FALSE---------------------------------------------------
knitr::opts_chunk$set(fig.align="center", warning=FALSE)

## ------------------------------------------------------------------------
library(wq)

## ----eval=FALSE----------------------------------------------------------
#  sfbay <- read.csv("sfbay.csv", header = FALSE, as.is = TRUE,
#                    skip = 2)
#  names(sfbay) <- c('date', 'time', 'stn', 'depth', 'chl', 'dox',
#                    'spm', 'ext', 'sal', 'temp', 'nox', 'nhx')
#  sfbay <- subset(sfbay, stn %in% c(21, 24, 27, 30, 32, 36) &
#                  substring(date, 7, 10) %in% 1985:2004)

## ------------------------------------------------------------------------
head(sfbay)

## ------------------------------------------------------------------------
x <- sample(1:nrow(sfbay), 10)
sfbay[x, "dox"]
sfbay1 <- transform(sfbay,
                    dox = round(100 * dox/oxySol(sal, temp), 1))
sfbay1[x, "dox"]

## ------------------------------------------------------------------------
sfb <- wqData(sfbay, c(1, 3:4), 5:12, site.order = TRUE,
              type = "wide", time.format = "%m/%d/%Y")
head(sfb)

## ------------------------------------------------------------------------
summary(sfb)

## ---- fig.width=5, fig.height=3.1----------------------------------------
plot(sfb, vars = c('dox', 'temp'), num.col = 2)

## ------------------------------------------------------------------------
y <- tsMake(sfb, focus = "chl", layer = c(0, 5))
y[1:4, ]
tsp(y)

## ---- width=5, height=8--------------------------------------------------
plotTs(y[, 1:4], dot.size = 1.3, ylab = "Chlorophyll in San Francisco Bay",
      strip.labels = paste("Station", 21:24), ncol = 1, scales = "free_y")

## ------------------------------------------------------------------------
head(tsMake(sfb, focus = "chl", layer = c(0, 5), type = 'zoo'))

## ------------------------------------------------------------------------
chl27 <- sfbayChla[, 's27']
tsp(chl27)
chl27 <- round(chl27, 1)
head(ts2df(chl27))

## ------------------------------------------------------------------------
y <- window(sfbayChla, start = 2005,
            end = c(2009, 12))  # 5 years, 16 sites
round(mts2ts(y, seas = 2:4), 1)  # focus on Feb-Apr spring bloom

## ------------------------------------------------------------------------
chl27 <- sfbayChla[, "s27"]
chl27a <- interpTs(chl27, gap = 3)

## ----fig.height=3.7, fig.width=6-----------------------------------------
plot(chl27a, col = "red", lwd = .5, xlab = "")
lines(chl27, col = "blue", lwd = 1.5)

## ------------------------------------------------------------------------
mannKen(Nile)

## ---- fig.height=3.7, fig.width=6----------------------------------------
plot(Nile, ylab = "Flow", xlab = "")
abline(v=1898, col='blue')
pett(Nile)

## ------------------------------------------------------------------------
y <- ts.intersect(Nile, LakeHuron)
pett(y)

## ------------------------------------------------------------------------
y <- sfbayChla
y1 <- tsSub(y, seas = 2:4)  # focus on Feb-Apr spring bloom
y2 <- aggregate(y1, 1, mean, na.rm = FALSE)
signif(mannKen(y2), 3)

## ------------------------------------------------------------------------
chl27 <- sfbayChla[, "s27"]
seaKen(chl27)

## ------------------------------------------------------------------------
seaRoll(chl27, w = 10)

## ---- fig.height=8, fig.width=7------------------------------------------
x <- sfbayChla
seasonTrend(x, plot = TRUE, ncol = 2, scales = 'free_y')

## ------------------------------------------------------------------------
x <- sfbayChla[, 's27']
trendHomog(x)

## ------------------------------------------------------------------------
chl <- sfbayChla[, 1:12]  # first 12 stns have good data coverage
seaKen(mts2ts(chl, 2:4))  # regional trend in spring bloom

## ---- fig.height=3.7, fig.width=6----------------------------------------
chla1 <- aggregate(sfbayChla, 1, mean, na.rm = TRUE)
chla1 <- chla1[, 1:12]
eofNum(chla1)

## ------------------------------------------------------------------------
e1 <- eof(chla1, n = 1)
e1

## ---- fig.height=3.1, fig.width=5----------------------------------------
eofPlot(e1, type = "amp")

## ------------------------------------------------------------------------
chl27b <- interpTs(sfbayChla[, "s27"], gap = 3)
chl27b <- ts2df(chl27b, mon1 = 10, addYr = TRUE, omit = TRUE)
head(round(chl27b, 1))

## ---- fig.height=3.1, fig.width=5----------------------------------------
e2 <- eof(chl27b, n = 2, scale. = TRUE)
eofPlot(e2, type = "coef")

## ---- fig.height=6, fig.width=6------------------------------------------
chl27 <- sfbayChla[, "s27"]
d1 <- decompTs(chl27)
plot(d1, nc = 1, main = "Station 27 Chl-a decomposition")

## ---- fig.height=4.3, fig.width=7----------------------------------------
plotSeason(chl27, num.era = 3, same.plot = FALSE, ylab = 'Stn 27 Chl-a')

## ---- fig.height=4.3, fig.width=7----------------------------------------
plotSeason(chl27, num.era = 3, same.plot = TRUE, ylab = 'Stn 27 Chl-a')

## ---- fig.height=4.3, fig.width=7----------------------------------------
plotSeason(chl27, "by.month", ylab = 'Stn 27 Chl-a')

## ------------------------------------------------------------------------
chl27 <- sfbayChla[, 's27']
p1 <- phenoPhase(chl27)
head(p1)
p2 <- phenoPhase(chl27, c(1, 6))
head(p2)
p3 <- phenoAmp(chl27, c(1, 6))
head(p3)

## ------------------------------------------------------------------------
zchl <- tsMake(sfb, focus = "chl", layer = c(0, 5), type = 'zoo')
head(zchl)
zchl27 <- zchl[, 3]
head(phenoPhase(zchl27))
head(phenoPhase(zchl27, c(1, 6), out = 'doy'))
head(phenoPhase(zchl27, c(1, 6), out = 'julian'))

## ---- fig.height=3.7, fig.width=6----------------------------------------
chl <- aggregate(sfbayChla[, 1:6], 1, meanSub, 2:4, na.rm = TRUE)
plotTsAnom(chl, ylab = 'Chlorophyll-a',
           strip.labels = paste('Station', substring(colnames(chl), 2, 3)))

## ---- fig.height=4.3, fig.width=7----------------------------------------
chl27 <- sfbayChla[, "s27"]
plotTsTile(chl27)

