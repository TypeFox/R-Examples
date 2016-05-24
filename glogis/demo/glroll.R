## package
stopifnot(require("glogis"))

## data
data("hicps", package = "glogis")
hicps <- hicps[, colnames(hicps) != "EU"]

## rolling mean from fitted generalized logistic distribution
glmean <- function(x) {
  if(any(is.na(x))) NA else glogisfit(x, hessian = FALSE)$moments[1]
}

## rolling GL mean with window width of one year for all countries
hicp_glmean <- rollapply(hicps, 12, glmean)

## visualization
gray_red <- rgb(c(0.2, 0.8), c(0.2, 0), c(0.2, 0), alpha = 0.3)
gray_red1 <- rgb(c(0.2, 0.8), c(0.2, 0), c(0.2, 0))
plot(hicps, plot.type = "single", lwd = 1.5, col = gray_red[1],
  xlab = "Time", ylab = "Monthly inflation rates")
for(i in 1:ncol(hicps)) lines(hicp_glmean[,i], lwd = 1.5, col = gray_red[2])
