###################
## PRELIMINARIES ##
###################

## packages
stopifnot(require("glogis") & require("fxregime"))

## convenience code
moments <- function(object, ...) UseMethod("moments")
moments.glogisfit <- function(object, ...) object$moments
moments.breakpoints.glogisfit <- function(object, breaks = NULL, ...)
  t(sapply(refit(object, breaks = breaks), "[[", "moments"))

## data
data("hicps", package = "glogis")
x <- na.omit(hicps[, "Portugal"])


##############
## MODELING ##
##############

## fit full-sample model
x_gf <- glogisfit(x)

## assess stability of full-sample model
x_efp <- gefp(x_gf, fit = NULL)

## estimate breakpoints (with minimal segment size of 2 years)
x_bp <- breakpoints(x_gf, h = 24)

## (minimal) number of breaks
x_nbreaks <- if(sctest(x_efp, functional = supLM(0.1))$p.value > 0.05) {
  0
} else {
  max(2, length(breakdates(x_bp)))
}

## refit segmented model
x_rf <- refit(x_bp, breaks = x_nbreaks)


#############
## RESULTS ##
#############

## parameter estimates and associated moments
coef(x_bp, breaks = x_nbreaks)
moments(x_bp, breaks = x_nbreaks)

## series and test
par(mfrow = c(1, 2))
plot(x, main = "Series with Fitted Mean", ylab = "hicps", xlab = "Time")
if(x_nbreaks > 0) lines(confint(x_bp, breaks = x_nbreaks))
lines(fitted(x_bp, breaks = x_nbreaks, type = "mean"), col = 4)
plot(x_efp, functional = supLM(0.1), main = "supLM test")
lines(x_bp, breaks = x_nbreaks)

## fitted segmented model
par(mfrow = c(1, x_nbreaks + 1))
for(i in 1:(x_nbreaks+1)) plot(x_rf[[i]], moments = TRUE,
  xlab = paste(names(x_rf)[i], "\n", "Goodness-of-fit p-value: ",
  format.pval(summary(x_rf[[i]])$chisq.test$p.value, digits = 4), sep = ""))
