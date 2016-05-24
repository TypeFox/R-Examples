## Package and data
library("betareg")
data("GasolineYield", package = "betareg")

## Same results as in Table 1 in Kosmidis and Firth (2010, EJS,
## http://dx.doi.org/10.1214/10-EJS579)
gyML <- betareg(yield ~ batch + temp, data = GasolineYield, link.phi = "identity")
gyBC <- betareg(yield ~ batch + temp, data = GasolineYield, link.phi = "identity", type = "BC")
gyBR <- betareg(yield ~ batch + temp, data = GasolineYield, link.phi = "identity", type = "BR")

## Coefficients and standard errors
se <- function(obj, ...) sqrt(diag(vcov(obj, ...)))
sapply(list(gyML, gyBC, gyBR), coef)
sapply(list(gyML, gyBC, gyBR), se)

## Same results as in Table 3 in Kosmidis and Firth (2010, EJS,
## http://dx.doi.org/10.1214/10-EJS579). BR and BC estimates in the
## latter study were calculated in a different way than betareg
## computes them which provides some indication on the correctness of
## implementation
gyMLlog <- betareg(yield ~ batch + temp, data = GasolineYield, link.phi = "log")
gyBClog <- betareg(yield ~ batch + temp, data = GasolineYield, link.phi = "log", type = "BC")
gyBRlog <- betareg(yield ~ batch + temp, data = GasolineYield, link.phi = "log", type = "BR")
sapply(list(gyMLlog, gyBClog, gyBRlog), coef)
sapply(list(gyMLlog, gyBClog, gyBRlog), se)

## Fit with temp as a dispersion covariate
gy2ML <- betareg(yield ~ batch + temp | temp, data = GasolineYield)
gy2BC <- betareg(yield ~ batch + temp | temp, data = GasolineYield, type = "BC")
gy2BR <- betareg(yield ~ batch + temp | temp, data = GasolineYield, type = "BR")
sapply(list(gy2ML, gy2BC, gy2BR), coef)
sapply(list(gy2ML, gy2BC, gy2BR), se)

