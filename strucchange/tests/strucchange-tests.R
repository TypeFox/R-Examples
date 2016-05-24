library("strucchange")

## Nile data
data("Nile")

Nile.rcus <- efp(Nile ~ 1, type = "Rec-CUSUM")
Nile.ocus <- efp(Nile ~ 1, type = "OLS-CUSUM")
Nile.re <- efp(Nile ~ 1, type = "RE")
Nile.rmos <- efp(Nile ~ 1, type = "Rec-MOSUM", h = 0.5)
Nile.omos <- efp(Nile ~ 1, type = "OLS-MOSUM", h = 0.5)
Nile.me <- efp(Nile ~ 1, type = "ME", h = 0.5)
sctest(Nile.rcus)
sctest(Nile.ocus)
sctest(Nile.ocus, alt = TRUE)
sctest(Nile.re)
sctest(Nile.re, alt = TRUE)
sctest(Nile.re, fun = "range")
sctest(Nile.rmos)
sctest(Nile.omos)
sctest(Nile.omos, fun = "range")
sctest(Nile.me)
sctest(Nile.me, fun = "range")

Nile.fs <- Fstats(Nile ~ 1)
sctest(Nile.fs, type = "supF")
sctest(Nile.fs, type = "aveF")
sctest(Nile.fs, type = "expF")
breakpoints(Nile.fs)

## Seatbelt data
data("UKDriverDeaths")
seatbelt <- log10(UKDriverDeaths)
seatbelt <- cbind(seatbelt, lag(seatbelt, k = -1), lag(seatbelt, k = -12))
colnames(seatbelt) <- c("y", "ylag1", "ylag12")
seatbelt <- window(seatbelt, start = c(1970, 1), end = c(1984,12))

seat.rcus <- efp(y ~ ylag1 + ylag12, data = seatbelt, type = "Rec-CUSUM")
seat.ocus <- efp(y ~ ylag1 + ylag12, data = seatbelt, type = "OLS-CUSUM")
seat.re <- efp(y ~ ylag1 + ylag12, data = seatbelt, type = "RE")
seat.rmos <- efp(y ~ ylag1 + ylag12, data = seatbelt, type = "Rec-MOSUM", h = 0.5)
seat.omos <- efp(y ~ ylag1 + ylag12, data = seatbelt, type = "OLS-MOSUM", h = 0.5)
seat.me <- efp(y ~ ylag1 + ylag12, data = seatbelt, type = "ME", h = 0.5)
sctest(seat.rcus)
sctest(seat.ocus)
sctest(seat.ocus, alt = TRUE)
sctest(seat.re)
sctest(seat.re, alt = TRUE)
sctest(seat.re, fun = "range")
sctest(seat.rmos)
sctest(seat.omos)
sctest(seat.omos, fun = "range")
sctest(seat.me)
sctest(seat.me, fun = "range")

seat.fs <- Fstats(y ~ ylag1 + ylag12, data = seatbelt, from = 0.1)
sctest(seat.fs, type = "supF")
sctest(seat.fs, type = "aveF")
sctest(seat.fs, type = "expF")
breakpoints(seat.fs)

## German M1 data
data("GermanM1")
LTW.model <- dm ~ dy2 + dR + dR1 + dp + m1 + y1 + R1 + season
M1.model <- dm ~ dy2 + dR + dR1 + dp + ecm.res + season

M1.ocus <- efp(LTW.model, data = GermanM1, type = "OLS-CUSUM")
M1.re <- efp(LTW.model, data = GermanM1, type = "RE")
M1.fs <- Fstats(LTW.model, data = GermanM1, from = 0.1)
sctest(M1.ocus)
sctest(M1.ocus, alt = TRUE)
sctest(M1.re)
sctest(M1.re, fun = "range")
sctest(M1.re, alt = TRUE)
sctest(M1.fs, type = "supF")
sctest(M1.fs, type = "aveF")
sctest(M1.fs, type = "expF")

M1 <- historyM1
ols.efp <- efp(M1.model, type = "OLS-CUSUM", data = M1)
newborder <- function(k) 1.5778*k/118
ols.mefp <- mefp(ols.efp, period = 2)
ols.mefp2 <- mefp(ols.efp, border = newborder)
M1 <- GermanM1
ols.mon <- monitor(ols.mefp)
ols.mon2 <- monitor(ols.mefp2)
ols.mon
ols.mon2

## Grossarl data
data("Grossarl")
Grossarl.bp <- breakpoints(fraction ~ 1, data = Grossarl, h = 0.1)
summary(Grossarl.bp)
confint(Grossarl.bp)

