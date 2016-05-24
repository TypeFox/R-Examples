x0 <- function(object)
{
  if (class(object)[1] != "onls") stop("Input must be a result of 'onls'!")
  val <- as.vector(object$x0)
  if (!is.null(object$na.action)) 
    val <- napredict(object$na.action, val)
  lab <- "x0 values from orthogonal model"
  attr(val, "label") <- lab
  val  
}

y0 <- function(object)
{
  if (class(object)[1] != "onls") stop("Input must be a result of 'onls'!")
  val <- as.vector(object$y0)
  if (!is.null(object$na.action)) 
    val <- napredict(object$na.action, val)
  lab <- "y0 values from orthogonal model"
  attr(val, "label") <- lab
  val  
}

deviance_o <- function(object)
{
  val <- sum(object$resid_o^2)   
  lab <- "Deviance (RSS) of orthogonal residuals from orthogonal model"
  attr(val, "label") <- lab
  val  
}

residuals_o <- function(object)
{
  val <- object$resid_o 
  if (!is.null(object$na.action)) 
    val <- napredict(object$na.action, val)
  lab <- "Orthogonal residuals from orthogonal model"
  attr(val, "label") <- lab
  val  
}

logLik_o <- function(object) 
{
  res_o <- object$resid_o
  N <- length(res_o)
  if (is.null(w <- object$weights)) w <- rep_len(1, N)
  zw <- w == 0
  val_o <- -N * (log(2 * pi) + 1 - log(N) - sum(log(w + zw)) + 
                   log(sum(w * res_o^2)))/2
  val <- val_o
  attr(val, "df") <- 1L + length(coef(object))
  attr(val, "nobs") <- attr(val, "nall") <- sum(!zw)
  class(val) <- "logLik"
  val
  attr(val, "label") <- "Log-likelihood using orthogonal residuals from orthogonal model"
  class(val) <- c("logLik_o", "logLik") 
  val
}

print.logLik_o <- function(x, ...) 
{
  cat("'log Lik.' ", paste(format(c(x), digits = 7), collapse = ", "), 
      " (df=", format(attr(x, "df")), ")\n", 
      "\"", attr(x, "label"), "\"\n", sep = "")
  invisible(x)
}

check_o <- function(object, plot = TRUE)
{
  if (class(object)[1] != "onls") stop("Input must be a result of 'onls'!")
  X <- object$pred
  Y <- object$resp
  X0 <- x0(object)
  Y0 <- y0(object)
  MODEL <- object$model
  FORMULA <- object$formula[[3]]  
  nameX <- attr(X, "name")   
  nameY <- attr(Y, "name")   
  eps <- sqrt(.Machine$double.eps)
  
  ## slope of model at x0, i.e. first derivative
  UP <- X0 + eps
  MODEL[[nameX]] <- UP
  upEVAL <- eval(FORMULA, envir = MODEL)
  DOWN <- X0 - eps
  MODEL[[nameX]] <- DOWN
  downEVAL <- eval(FORMULA, envir = MODEL)
  modelSLOPE <- (upEVAL - downEVAL)/(UP - DOWN)
  
  ## slope of object using X, Y, X0, Y0
  dataSLOPE <- (Y - Y0)/(X - X0)
    
  ## angle between two slopes, must be 90 degrees for orthogonality  
  tanalpha <- abs((modelSLOPE - dataSLOPE)/(1 + modelSLOPE * dataSLOPE))
  alpha <- atan(tanalpha) * 360/(2*pi)  
  QUAL <- alpha > 89.95 & alpha < 90.05
  
  ## plot difference to 90 degrees
  if (plot) plot(as.numeric(alpha), type = "o",  col = ifelse(QUAL == TRUE, "black", "darkred"), 
                 pch = 16, cex = 1.5, ylab = "Alpha", xlab = "x")
  
  OUT <- data.frame(X, X0, Y, Y0, alpha, modelSLOPE, QUAL)
  colnames(OUT) <- c(nameX, "x0", nameY, "y0", "alpha", "df/dx", "Ortho")
  return(OUT)
}
