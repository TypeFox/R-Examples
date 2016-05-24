#' @title Function to Give the Log Likelihood Function
#' @aliases qrLogLikFun
#' @description Function to compute the (negative) log likelihood for fitting cdfquantile models.
#' 
#' @param fd A string that specifies the parent distribution.
#' @param sd A string that specifies the sub-family distribution.
#'
#' @return The log-likelihood calculation function given a specified cdfquantile distribution.
#' 
#' @export
#' @examples
#' qrLogLikFun('t2','t2')
#' 
qrLogLikFun <- function(fd, sd) {
  # logit-XX
  if (fd == "logit") {
    
    if (sd == "logistic") {
      quantreg_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- ((mu - sigma * log(sigma)) - (-1 + sigma) * log(-1 + 1/y) - 
          2 * sigma * log(y + exp(mu/sigma) * (-1 + 1/y)^(1/sigma) * y))/sigma
        -sum(loglik, na.rm = TRUE)
      }
    }
    
    # logit-Cauchy
    if (sd == "cauchy") {
      quantreg_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- (mu + sigma * log(pi) - sigma * log(sigma) + cot(pi * y) - 
                     2 * sigma * log(1 + exp((mu + cot(pi * y))/sigma)) + 
                     2 * sigma * log(csc(pi * y)))/sigma
        -sum(loglik, na.rm = TRUE)
      }
    }
    
    # logit-t2
    if (sd == "t2") {
      quantreg_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        a1 <- a2 <- y
        
        for (i in 1:length(y)) {
          c1 <- ((y[i] >= 0) & (y[i] < 1/2))  #1st situation
          c2 <- (y[i] >= 1/2 & y[i] <= 1)  #2ed situation
          
          if (c1) {
            a1[i] <- -sqrt((1 - 2 * y[i])^2/(2 * (1 - y[i]) * y[i]))
          } else if (c2) {
            a1[i] <- sqrt((1 - 2 * y[i])^2/(2 * (1 - y[i]) * y[i]))
          } else {
            a1[i] <- 0
          }
          
          
          c3 <- (y[i] > 0 & y[i] < 1/2)  #3rd situation
          c4 <- (y[i] > 1/2 & y[i] < 1)  #4th situation
          c5 <- (y[i] == 1 | y[i] == 0)
          if (y[i] == 0.5) {
            a2[i] <- 2 * sqrt(2)
          } else if (c3) {
            a2[i] <- (1 - 2 * y[i])/sqrt(8 * (1 - 2 * y[i])^2 * ((1 - y[i]) * y[i])^3)
          } else if (c4) {
            a2[i] <- (2 * y[i] - 1)/sqrt(8 * (1 - 2 * y[i])^2 * ((1 - y[i]) * y[i])^3)
          } else if (c5) {
            a2[i] <- NA  # NO IDEA    
          } else {
            a2[i] <- 0
          }
        }
        
        lik <- (exp((mu + a1)/sigma) * a2)/((exp(mu/sigma) + exp(a1/sigma))^2 * sigma)
        loglik <- log(lik)
        -sum(loglik, na.rm = TRUE)
      }
    }
    
    if (sd == "burr8") {
      quantreg_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <-log(csc(pi * y)) + ((mu + sigma * log(pi) - sigma * log(sigma)) 
                  + log(tan((pi *   y)/2)) - 2 * sigma * log(exp(mu/sigma) + 
                  tan((pi * y)/2)^(1/sigma)))/sigma
        -sum(loglik, na.rm = TRUE)
      }
    }
  }
  
  # arcsinh-XX
  if (fd == "arcsinh") {
    
    # arcsinh-Cauchy
    if (sd == "cauchy") {
      quantreg_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- log(pi) + asinh((mu + cot(pi * y))/sigma) - 
          2 * log(1 + exp(asinh((mu + cot(pi *   y))/sigma))) - 
          (1/2) * log(mu^2 + sigma^2 + 2 * mu * cot(pi * y) + cot(pi * y)^2) + 
          2 * log(csc(pi * y))
        
        -sum(loglik, na.rm = TRUE)
      }
    }
    
    # arcsinh-t2
    if (sd == "t2") {
      quantreg_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        a1 <- a2 <- y
        
        for (i in 1:length(y)) {
          c1 <- ((y[i] >= 0) & (y[i] < 1/2))  #1st situation
          c2 <- (y[i] >= 1/2 & y[i] <= 1)  #2ed situation
          
          if (c1) {
            a1[i] <- -sqrt((1 - 2 * y[i])^2/(2 * (1 - y[i]) * y[i]))
          } else if (c2) {
            a1[i] <- sqrt((1 - 2 * y[i])^2/(2 * (1 - y[i]) * y[i]))
          } else {
            a1[i] <- 0
          }
          
          
          c3 <- (y[i] > 0 & y[i] < 1/2)  #3rd situation
          c4 <- (y[i] > 1/2 & y[i] < 1)  #4th situation
          c5 <- (y[i] == 1 | y[i] == 0)
          if (y[i] == 0.5) {
            a2[i] <- 2 * sqrt(2)
          } else if (c3) {
            a2[i] <- (1 - 2 * y[i])/sqrt(8 * (1 - 2 * y[i])^2 * ((1 - y[i]) * y[i])^3)
          } else if (c4) {
            a2[i] <- (2 * y[i] - 1)/sqrt(8 * (1 - 2 * y[i])^2 * ((1 - y[i]) * y[i])^3)
          } else if (c5) {
            a2[i] <- NA  # NO IDEA    
          } else {
            a2[i] <- 0
          }
        }
        
        lik <- (exp(asinh((-mu + a1)/sigma)) * a2)/((1 + exp(asinh((-mu + a1)/sigma)))^2 * sigma * 
                                                      sqrt(1 + (mu - a1)^2/sigma^2))
        loglik <- log(lik)
        -sum(loglik, na.rm = TRUE)
      }
    }
    
    # arcsinh-burr8
    if (sd == "burr8") {
      quantreg_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- log(pi) - log(sigma) + asinh((-mu + log(tan((pi * y)/2)))/sigma) - 
          2 * log(1 + exp(asinh((-mu + log(tan((pi * y)/2)))/sigma))) + 
          log(csc(pi * y)) - (1/2) * log(1 +(mu - log(tan((pi * y)/2)))^2/sigma^2)
        
        -sum(loglik, na.rm = TRUE)
      }
    }
    
  }
  
  # t2-XX
  if (fd == "t2") {
    # t2-Cauchy
    if (sd == "cauchy") {
      quantreg_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- log(pi) + log(sigma) - log(mu^2 + 2 * sigma^2 + 2 * mu * cot(pi * y)
                   + cot(pi * y)^2) - (1/2) * log(2 + (mu + cot(pi * y))^2/sigma^2) + 2 * log(csc(pi * y))
        -sum(loglik, na.rm = TRUE)
      }
    }
    
    
    # t2-t2
    if (sd == "t2") {
      quantreg_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- 2 * log(sigma) - (3/2) * log(1 + 2 * sqrt(2) * mu * sqrt(1 - y) * sqrt(y) + 2 * 
                   (-2 + mu^2 + 2 * sigma^2) * y - 
                  4 * mu * sqrt(2 - 2 * y) * y^(3/2) - 2 * (-2 + mu^2 + 2 * sigma^2) * y^2)
        -sum(loglik, na.rm = TRUE)
      }
    }
    
    # t2-burr7
    if (sd == "burr7") {
      quantreg_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- -log(2) + 2 * log(sigma) - (3/2) * log(mu^2 + 2 * sigma^2 + 
                   2 * mu * atanh(1 - 2 * y) + atanh(1 - 2 * y)^2) - log(1 - y) - log(y)
        
        -sum(loglik, na.rm = TRUE)
      }
    }
    
    # t2-burr8
    if (sd == "burr8") {
      quantreg_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- (log(pi) + 2 * log(sigma)) + log(csc(pi * y)) - 
          (3/2) * log(mu^2 + 2 * sigma^2 - 2 * mu * log(tan((pi * y)/2)) + log(tan((pi * y)/2))^2)
        -sum(loglik, na.rm = TRUE)
      }
      
    }
  }
  
  # burr8-XX
  if (fd == "burr8") {
    # burr8 Cauchy
    if (sd == "cauchy") {
      quantreg_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- (mu + sigma * log(2) - sigma * log(sigma) + cot(pi * y) - 
                     sigma * log(1 + exp((2 * (mu + cot(pi * y)))/sigma)) + 
                     2 * sigma * log(csc(pi * y)))/sigma
        
        -sum(loglik, na.rm = TRUE)
      }
    }
    
    # burr8 t2
    if (sd == "t2") {
      quantreg_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        a1 <- a2 <- y
        
        for (i in 1:length(y)) {
          c1 <- ((y[i] >= 0) & (y[i] < 1))  #1st situation
          c2 <- (y[i] == 0 & y[i] == 1)  #2ed situation
          
          if (c1) {
            a1[i] <- sech(((1 - 2 * y[i])/(sqrt(2) * sqrt((-(-1 + y[i])) * y[i])) + mu[i])/sigma[i])/(2 * sqrt(2) * pi * ((-(-1 + y[i])) * y[i])^(3/2) * sigma[i])
            # sech((mu[i] + (1 - 2*y[i])/(sqrt(2)*sqrt(1 - y[i])*sqrt(y[i])))/sigma[i])/
            # (2*sqrt(2)*pi*sigma[i]*(1 - y[i])^(3/2)*y[i]^(3/2))
          } else if (c2) {
            a1[i] <- NA
          } else {
            a1[i] <- 0
          }
        }
        
        loglik <- log(a1)
        -sum(loglik, na.rm = TRUE)
      }
    }
    
    # burr8 burr7
    if (sd == "burr7") {
      quantreg_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- -log(2 * pi * y * sigma - 2 * pi * y^2 * sigma) + 
          log(sech((mu + atanh(1 - 2 *  y))/sigma))
        -sum(loglik, na.rm = TRUE)
      }
    }
    
    # burr8 burr8
    if (sd == "burr8") {
      quantreg_loglik <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        mu <- hx
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        sigma <- exp(gz)
        loglik <- (mu + sigma * log(2) + sigma * log(csc(pi * y)) + log(tan(pi * y/2)) - sigma * 
                     log(exp((2 * mu)/sigma) * sigma + sigma * tan((pi * y)/2)^(2/sigma)))/sigma
        -sum(loglik, na.rm = TRUE)
      }
    }
  }
  
  if (fd == "km" | sd == "km") {
    quantreg_loglik <- function(h, y, x, z) {
      hx <- x %*% h[1:length(x[1, ])]
      a <- exp(hx)
      gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
      b <- exp(gz)
      loglik <- log(a) + log(b) + (a - 1) * log(y) + (b - 1) * log(1 - y^a)
      -sum(loglik, na.rm = TRUE)
    }
  }
  
  
  quantreg_loglik
} 
