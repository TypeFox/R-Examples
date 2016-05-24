dfmix <- function(x = NULL, w = NULL, Theta = NULL, ...)
{
  digits <- getOption("digits"); options(digits = 15)
  
  if (is.null(x)) {
    stop(sQuote("x"), " must not be NULL!", call. = FALSE)
  }

  if ((!is.numeric(x)) && (!is.data.frame(x))) {
    stop(sQuote("x"), " numeric or data frame is requested!", call. = FALSE)
  }
  
  x <- as.matrix(x)
  
  d <- ncol(x)
  n <- nrow(x)  
  
  if (is.null(w)) {
    stop(sQuote("w"), " must not be NULL!", call. = FALSE)
  }  
  
  if (!is.numeric(w)) {
    stop(sQuote("w"), " numeric vector is requested!", call. = FALSE)
  }
  
  if (!all(w > 0.0)) {
    stop("all ", sQuote("w"), " must be greater than 0.0!", call. = FALSE)
  }
  
  if (is.null(Theta)) {
    stop(sQuote("Theta"), " must not be NULL!", call. = FALSE)
  }

  if (!is.list(Theta)) {
    stop(sQuote("Theta"), " list is requested!", call. = FALSE)
  }
  
  Names <- names(Theta)
  
  c <- length(w)

  pdf <- Theta[grep("pdf", Names)]

  theta1 <- Theta[grep("theta1", Names)]
  
  theta2 <- Theta[grep("theta2", Names)]

  f <- array(data = 0.0, dim = n, dimnames = NULL)

  for (i in 1:c) {
    fi <- rep(1.0, n)
    
    for (j in 1:d) {
      if (pdf[[i]][j] == .rebmix$pdf[1]) {
        fi <- fi * dnorm(as.numeric(x[, j]), mean = as.numeric(theta1[[i]][j]), sd = as.numeric(theta2[[i]][j]), ...)
      }
      else
      if (pdf[[i]][j] == .rebmix$pdf[2]) {
        fi <- fi * dlnorm(as.numeric(x[, j]), meanlog = as.numeric(theta1[[i]][j]), sdlog = as.numeric(theta2[[i]][j]), ...)
      }
      else
      if (pdf[[i]][j] == .rebmix$pdf[3]) {
        fi <- fi * dweibull(as.numeric(x[, j]), scale = as.numeric(theta1[[i]][j]), shape = as.numeric(theta2[[i]][j]), ...)
      }
      else
      if (pdf[[i]][j] == .rebmix$pdf[4]) {
        fi <- fi * dbinom(as.integer(x[, j]), size = as.integer(theta1[[i]][j]), prob = as.numeric(theta2[[i]][j]), ...)
      }
      else
      if (pdf[[i]][j] == .rebmix$pdf[5]) {
        fi <- fi * dpois(as.integer(x[, j]), lambda = as.numeric(theta1[[i]][j]), ...)
      }
      else
      if (pdf[[i]][j] == .rebmix$pdf[6]) {
        fi <- fi * ddirac(as.numeric(x[, j]), location = as.numeric(theta1[[i]][j]))
      }
      else
      if (pdf[[i]][j] == .rebmix$pdf[7]) {
        fi <- fi * dgamma(as.numeric(x[, j]), scale = as.numeric(theta1[[i]][j]), shape = as.numeric(theta2[[i]][j]), ...)
      }      
    }
    
    f <- f + as.numeric(w[i]) * fi
  }
  
  options(digits = digits)

  rm(list = ls()[!(ls() %in% c("f"))])
  
  invisible(f)
} ## dfmix
