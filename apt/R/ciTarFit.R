ciTarFit <- function(y, x, model = c('tar','mtar'), lag = 1, thresh = 0,
  small.win = NULL)
{
  if(!is.ts(y) || !is.ts(x) ) stop("Please provide time series data.\n")
  if (!identical(tsp(y), tsp(x))) {
    stop("y and x have different properties.\n")
  }
  model <- match.arg(model)
  
  A <- start(y); B <- end(y); Q <- tsp(y)[3]
  name.y <- deparse(substitute(y)); name.x <- deparse(substitute(x))
  data.LR <- data.frame(cbind(y, x))
  colnames(data.LR) <- c(name.y, name.x)
  formula.LR <- as.formula(paste(name.y, "~", name.x, sep = ""))
  LR <- lm(formula = formula.LR, data = data.LR) 
  z <- ts(data = residuals(LR), start = A, end = B, frequency = Q)  
  lz <- lag(z, k = -1); dz <- diff(z); ldz <- lag(dz, k = -1)
  
  if(model == "tar")  {ind <- ifelse(test = lz >= thresh, yes = 1, no = 0)}
  if(model == "mtar") {ind <- ifelse(test =ldz >= thresh, yes = 1, no = 0)}
  pos <- lz * ind
  neg <- lz * (1 - ind)
  
  xx <- bsLag(h = dz, lag = lag, var.name = "diff.resid")
  if(tsp(xx)[1] >= tsp(pos)[1]) {sa <- start(xx)} else {sa <- start(pos)}       
  data.CI <- window(cbind(pos, neg, xx), start = sa, end = B, frequency=Q)
  colnames(data.CI) <- c("pos.resid.t_1", "neg.resid.t_1", colnames(xx)) 

  if (!is.null(small.win)){
    sw <- small.win[1] + (small.win[2] - 1) / Q
    if (sw <= time(y)[1]) {
      stop(paste("The value of small.win should be bigger than the start",
        "value of raw data to have a smaller window.", sep = "\n"))
    } else {    
      data.CI <- window(data.CI, start = small.win, end = B, frequency = Q)
    }
  }  
  
  CI <- lm(formula = diff.resid.t_0 ~ 0 + ., data = data.CI)
  f.phi <- car::linearHypothesis(CI, c("pos.resid.t_1 = 0", 
    "neg.resid.t_1 = 0")) 
  f.apt <- car::linearHypothesis(CI, "pos.resid.t_1 = neg.resid.t_1")
  sse <- deviance(CI)
  aic <- AIC(CI, k = 2)
  bic <- AIC(CI, k = log(nrow(data.CI)))    
  
  result <- listn(y, x, model, lag, thresh, data.LR, data.CI, z, lz, ldz, 
      LR, CI, f.phi, f.apt, sse, aic, bic)
  class(result) <- "ciTarFit"
  return(result)
} 

print.ciTarFit <- function(x, ...)
{
  cat("=== Long Run Regression\n") ; print(summary(x$LR))
  cat("=== Threshold Cointegration Regression\n"); print(summary(x$CI))
  cat("=== H1: No cointegration b/w two variables\n"); print(x$f.phi)
  cat("=== H2: Symmetric adjustment in the long run\n"); print(x$f.apt)    
} 