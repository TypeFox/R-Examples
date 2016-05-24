## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.show = "hold")
options(width = 80)

## -----------------------------------------------------------------------------
library(nestedRanksTest)
data(woodpecker_multiyear)
d.Qlob <- subset(woodpecker_multiyear, Species == "lobata")
d.Qagr <- subset(woodpecker_multiyear, Species == "agrifolia")
table(d.Qlob$Year, d.Qlob$Granary)
table(d.Qagr$Year, d.Qagr$Granary)

## -----------------------------------------------------------------------------
sum(d.Qlob$Distance > 500)
sum(d.Qagr$Distance > 500)

## ---- fig.width = 5, fig.height = 5-------------------------------------------
opa <- par(mfcol = c(2, 1), cex = 0.8, mar = c(3, 3, 0, 0), las = 1,
           mgp = c(2, 0.5, 0), tcl = -0.3)
stripchart(Distance ~ Granary, data = d.Qlob, subset = Year == "2002",
           xlab = "Distance (m)", ylab = "Granary", pch = 1, col = "blue",
           method = "jitter", jitter = 0.2, xlim = c(0, 500),
           frame.plot = FALSE)
stripchart(Distance ~ Granary, data = d.Qlob, subset = Year == "2004",
           pch = 1, col = "red", method = "jitter", jitter = 0.2, add = TRUE)
legend("bottomright", legend = c("2002", "2004"), pch = 1, bty = "n",
       col = c("blue", "red"), title = "Valley oak acorns")
stripchart(Distance ~ Granary, data = d.Qagr, subset = Year == "2006",
           xlab = "Distance (m)", ylab = "Granary", pch = 1, col = "blue",
           method = "jitter", jitter = 0.2, xlim = c(0, 500),
           frame.plot = FALSE)
stripchart(Distance ~ Granary, data = d.Qagr, subset = Year == "2007",
           pch = 1, col = "red", method = "jitter", jitter = 0.2, add = TRUE)
legend("bottomright", legend = c("2006", "2007"), pch = 1, bty = "n",
       col = c("blue", "red"), title = "Live oak acorns")
par(opa)

## -----------------------------------------------------------------------------
with(d.Qlob, unlist(lapply(lapply(split(Distance, Granary), unique), median)))
with(d.Qagr, unlist(lapply(lapply(split(Distance, Granary), unique), median)))

## -----------------------------------------------------------------------------
wilcox.test(Distance ~ Year, data = d.Qlob)
wilcox.test(Distance ~ Year, data = d.Qagr)

## -----------------------------------------------------------------------------
wilcox.p.value <- function(x) wilcox.test(Distance ~ Year, data = x, exact = FALSE)$p.value
round(unlist(lapply(split(d.Qlob, d.Qlob$Granary), wilcox.p.value)), 4)
round(unlist(lapply(split(d.Qagr, d.Qagr$Granary), wilcox.p.value)), 4)

## ---- eval = FALSE------------------------------------------------------------
#  result <- nestedRanksTest(Distance ~ Year | Granary, data = d.Qlob)
#  result <- nestedRanksTest(Distance ~ Year, groups = Granary, data = d.Qlob)
#  ## Because no data= for default interface, here we use with()
#  result <- with(d.Qlob, nestedRanksTest(Year, Distance, Granary))
#  result <- with(d.Qlob, nestedRanksTest(y = Distance, x = Year, groups = Granary))

## ---- eval = FALSE------------------------------------------------------------
#  ## Error: 'groups' missing
#  result <- with(d.Qlob, nestedRanksTest(Year, Distance))
#  ## Error: invalid group specification in formula (no "| grouping.variable")
#  result <- nestedRanksTest(Distance ~ Year, data = d.Qlob)
#  ## Error: groups are specified with '|' in formula or with groups= argument, but not both
#  result <- nestedRanksTest(Distance ~ Year | Granary, groups = Granary, data = d.Qlob)

## -----------------------------------------------------------------------------
result.Qlob <- nestedRanksTest(Distance ~ Year | Granary, data = d.Qlob, n.iter = 2000)
print(result.Qlob)
#
#
result.Qagr <- nestedRanksTest(Distance ~ Year | Granary, data = d.Qagr, n.iter = 2000)
print(result.Qagr)

## ---- fig.width = 6, fig.height = 3-------------------------------------------
plot(result.Qlob)

## ---- fig.width = 6, fig.height = 3-------------------------------------------
plot(result.Qagr, main = expression(italic("Quercus agrifolia")),
     col = "lightgreen", p.col = "darkgreen", p.lty = 1, p.lwd = 4,
     mgp = c(2.5, 0.5, 0), cex.axis = 0.8, las = 1, tcl = -0.3)

## -----------------------------------------------------------------------------
str(result.Qlob)

## ---- eval = FALSE------------------------------------------------------------
#  nestedRanksTest::nestedRanksTest_Z

## ---- eval = FALSE------------------------------------------------------------
#  nestedRanksTest::nestedRanksTest_weights

