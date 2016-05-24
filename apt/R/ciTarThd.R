ciTarThd <- function(y, x, model = c('tar', 'mtar'), lag = 1, 
  th.range = 0.15, digits = 3)
{
  if(!is.ts(y) || !is.ts(x) ) stop("Please provide time series data.\n")
  if (!identical(tsp(y), tsp(x))) {
    stop("Properties of y and x are different.\n")
  }
  model <- match.arg(model)
  
  Q <- ciTarFit(y = y, x = x, model = model, lag = lag, thresh = 0)
  if(model == "tar")  { thresh.test <- sort(as.vector(Q$lz )) } 
  if(model == "mtar") { thresh.test <- sort(as.vector(Q$ldz)) }
  obs.tot <- NROW(y); obs.CI  <- nrow(Q$data.CI)
  a <- ceiling(obs.CI * th.range)
  b <- obs.CI - floor(obs.CI * th.range)
  
  path <- data.frame(matrix(0, ncol = 5, nrow = b - a + 1))
  colnames(path) <- paste("path", c("num", "thr", "sse", "aic", "bic"), 
    sep = '.')
  path$path.num <- 1:(b - a + 1)
  path$path.thr <- thresh.test[a:b]   
  for(i in a:b) {
    H <- ciTarFit(y = y, x = x, model = model, lag = lag, 
      thresh = thresh.test[i])
    path[(i-a+1), 3:5] <- c(H$sse, H$aic, H$bic) 
  }  
  sse.lowest <- min(path$path.sse)
  thresh.final <- path[which(sse.lowest == path$path.sse)[1], 2]

  Item <- c("lag", "thresh final", "thresh range", "sse.lowest", 
     "Total obs", "CI obs", "Lower obs", "Upper obs")
  Value <- round(c(lag, thresh.final, th.range, sse.lowest,
      obs.tot, obs.CI, a, b), digits = digits)
  basic <- data.frame(Item, Value) 
  colnames(basic)[2] <- model
  
  result <- listn(model, lag, th.range, th.final = thresh.final, 
      ssef = sse.lowest, obs.tot, obs.CI, basic, path)
  class(result) <- "ciTarThd"
  return(result)
}

print.ciTarThd <- function(x, ...) {print(x$basic)}

plot.ciTarThd <- function(x, ...)
{
  par(mfrow = c(2, 2))
  plot(path.sse ~ path.thr, data = x$path, type = "l", col = "blue", 
     ylab = "SSE", xlab = paste("threshold value for", x$model))
  plot(path.aic ~ path.thr, data = x$path, type = "l", col = "green",
     xlab = paste("threshold value for", x$model))
  plot(path.bic ~ path.thr, data = x$path, type = "l", col = "red",
     xlab = paste("threshold value for", x$model))
}