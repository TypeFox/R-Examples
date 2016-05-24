ciTarLag <- function(y, x, model = c("tar", "mtar"), maxlag = 4, 
  thresh = 0, adjust = TRUE) 
{
  if(!is.ts(y) || !is.ts(x) ) stop("Please provide time series data.\n")
  if (!identical(tsp(y), tsp(x))) {
    stop("Properties of y and x are different.\n")
  }
  model <- match.arg(model)

  small <- c(start(y)[1], start(y)[2] + 1 + maxlag)
  lag <- totObs <- coinObs <- sse <- NULL
  aic <- bic <- LB4 <- LB8 <-LB12 <- NULL
  for (k in 0:maxlag) 
  {
    if (adjust) {
        curt <- ciTarFit(y = y, x = x, model = model, lag = k,
          thresh = thresh, small.win = small)
    } else {
        curt <- ciTarFit(y = y, x = x, model = model, lag = k, 
          thresh = thresh)
    }
    cur <- summary(curt)$dia
    lag[k + 1] <- k;             totObs[k + 1]<- cur[3, 2]
    coinObs[k + 1] <- cur[4, 2]; sse[k + 1] <- cur[5, 2]
    aic[k + 1] <- cur[6, 2];     bic[k + 1] <- cur[7, 2]
    LB4[k + 1] <- cur[12, 2];    LB8[k + 1] <- cur[13, 2]
    LB12[k + 1]<- cur[14, 2]
  } 
  path <- data.frame(lag, totObs, coinObs, sse, aic, bic, LB4, LB8, LB12)
  cr1 <- path[order(aic), ]; cr2 <- path[order(bic), ]
  BestLag.byAic <- cr1[1, 1]; BestAic <- cr1[1, 5] 
  BestLag.byBic <- cr2[1, 1]; BestBic <- cr2[1, 6]   
  
  Item <- c("model", "max lag", "threshold", "BestLag.byAic", 
    "BestLag.byBic", "Best AIC", "Best BIC")
  Value <- c(model, maxlag, thresh, BestLag.byAic, 
    BestLag.byBic, BestAic, BestBic)
  out <- data.frame(Item, Value)
  result <- listn(path, out)
  class(result) <- "ciTarLag"
  return(result)
}

print.ciTarLag <- function(x, ...) {print(x$out)}

plot.ciTarLag <- function(x, ...)
{
  xlabel <- paste("Lag value (model = ", x$out[1, 2], 
      ", threshold = ", x$out[3, 2], ")")
  par(mfrow = c(2, 1))
  plot(aic ~ lag, data = x$path, type = "l", col = "green", 
    xlab = xlabel, ylab = "AIC")
  plot(bic ~ lag, data = x$path, type = "l", col = "red", 
    xlab = xlabel, ylab = "BIC")
}