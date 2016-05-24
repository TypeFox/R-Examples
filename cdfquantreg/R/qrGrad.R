#' @title Give the Gradient Function for CDF-Quantile Distribution Modles
#' @description Give the Gradient Function for CDF-Quantile Distribution Modles. 
#' @aliases qrGrad
#' @param fd A string that specifies the parent distribution.
#' @param sd A string that specifies the sub-family distribution.
#' 
#' @return grad The gradient function of parameter estimates, given a specified cdf-quantile distribution
#' 
#' @export
#' 
#' @examples
#' qrGrad('t2','t2')
qrGrad <- function(fd, sd) {
  # logit-XX
  if (fd == "logit") {
    
    if (sd == "logistic") {
      grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep(((1 - exp((hx)/exp(gz)) * (-1 + 1/y)^exp(-gz)))/(exp(gz) * 
          (1 + exp((hx)/exp(gz)) * (-1 + 1/y)^exp(-gz))), length(x[1, ])), 
          z * rep(((-(1/(1 + exp((hx)/exp(gz)) * (-1 + 1/y)^exp(-gz)))) * 
          (exp(gz) + exp((hx)/exp(gz) + gz) * (-1 + 1/y)^exp(-gz) + (1 - 
            exp((hx)/exp(gz)) * (-1 + 1/y)^exp(-gz)) * log(-1 + 1/y) + 
            hx - exp((hx)/exp(gz)) * (-1 + 1/y)^exp(-gz) * hx))/exp(gz), 
          length(z[1, ])))
        
        colSums(gd, na.rm = TRUE)
      }
    }
    
    # logit-Cauchy
    if (sd == "cauchy") {
      grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep(((-(-1 + exp((cot(pi * y) + hx)/exp(gz))))/(exp(gz) * 
          (1 + exp((cot(pi * y) + hx)/exp(gz))))), length(x[1, ])), z * rep(-(((exp(gz) + 
          exp((cot(pi * y) + hx + exp(gz) * gz)/exp(gz)) - (-1 + exp((cot(pi * 
          y) + hx)/exp(gz))) * cot(pi * y) + hx - exp((cot(pi * y) + hx)/exp(gz)) * 
          hx))/(exp(gz) * (1 + exp((cot(pi * y) + hx)/exp(gz))))), length(z[1, 
          ])))
        
        colSums(gd, na.rm = TRUE)
      }
    }
    
    # logit-t2
    if (sd == "t2") {
      grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        
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
          
        }
        
        
        gd <- cbind(x * rep(((exp(a1/exp(gz)) - exp((hx)/exp(gz))))/(exp(gz) * 
          (exp(a1/exp(gz)) + exp((hx)/exp(gz)))), length(x[1, ])), z * rep(-(((exp(a1/exp(gz) + 
          gz) + exp((hx)/exp(gz) + gz) + a1 * (-exp(a1/exp(gz)) + exp((hx)/exp(gz))) + 
          exp(a1/exp(gz)) * hx - exp((hx)/exp(gz)) * hx))/(exp(gz) * (exp(a1/exp(gz)) + 
          exp((hx)/exp(gz))))), length(z[1, ])))
        
        colSums(gd, na.rm = TRUE)
      }
    }
    
    if (sd == "burr8") {
      grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep(((1 - (2 * exp((hx)/exp(gz)))/(exp((hx)/exp(gz)) + 
          tan((pi * y)/2)^exp(-gz))))/exp(gz), length(x[1, ])), z * rep(-(((exp((hx)/exp(gz)) * 
          (exp(gz) + log(tan((pi * y)/2)) - hx) + tan((pi * y)/2)^exp(-gz) * 
          (exp(gz) - log(tan((pi * y)/2)) + hx)))/(exp(gz) * (exp((hx)/exp(gz)) + 
          tan((pi * y)/2)^exp(-gz)))), length(z[1, ])))
        
        colSums(gd, na.rm = TRUE)
      }
    }
  }
  
  # arcsinh-XX
  if (fd == "arcsinh") {
    
    # arcsinh-Cauchy
    if (sd == "cauchy") {
      grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        
        gd <- cbind(x * rep((-((cot(pi * y) + hx)/(exp(2 * gz) + cot(pi * 
          y)^2 + 2 * cot(pi * y) * hx + hx^2)) + 1/(exp(gz) * sqrt(1 + (cot(pi * 
          y) + hx)^2/exp(2 * gz))) - (2 * exp(asinh((cot(pi * y) + hx)/exp(gz)) - 
          gz))/((1 + exp(asinh((cot(pi * y) + hx)/exp(gz)))) * sqrt(1 + (cot(pi * 
          y) + hx)^2/exp(2 * gz)))), length(x[1, ])), z * rep(((-(exp(3 * 
          gz)/(exp(2 * gz) + cot(pi * y)^2 + 2 * cot(pi * y) * hx + hx^2)) - 
          (cot(pi * y) + hx)/sqrt(1 + (cot(pi * y) + hx)^2/exp(2 * gz)) + 
          (2 * exp(asinh((cot(pi * y) + hx)/exp(gz))) * (cot(pi * y) + hx))/((1 + 
          exp(asinh((cot(pi * y) + hx)/exp(gz)))) * sqrt(1 + (cot(pi * 
          y) + hx)^2/exp(2 * gz)))))/exp(gz), length(z[1, ])))
        colSums(gd, na.rm = TRUE)
      }
    }
    
    # arcsinh-t2
    if (sd == "t2") {
      grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        
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
        }
        
        
        gd <- cbind(x * rep((-(((-a1) * (1 + exp(asinh((a1 - hx)/exp(gz)))) + 
          (1 + exp(asinh((a1 - hx)/exp(gz)))) * hx - exp(gz) * (-1 + exp(asinh((a1 - 
          hx)/exp(gz)))) * sqrt((a1^2 + exp(2 * gz) - 2 * a1 * hx + hx^2)/exp(2 * 
          gz))))/((1 + exp(asinh((a1 - hx)/exp(gz)))) * (a1^2 + exp(2 * gz) - 
          2 * a1 * hx + hx^2))), length(x[1, ])), z * rep(-((exp(gz) * (exp(gz) + 
          exp(asinh((a1 - hx)/exp(gz)) + gz) + (a1 - hx) * sqrt((a1^2 + exp(2 * 
          gz) - 2 * a1 * hx + hx^2)/exp(2 * gz)) + exp(asinh((a1 - hx)/exp(gz))) * 
          (-a1 + hx) * sqrt((a1^2 + exp(2 * gz) - 2 * a1 * hx + hx^2)/exp(2 * 
          gz))))/((1 + exp(asinh((a1 - hx)/exp(gz)))) * (a1^2 + exp(2 * gz) - 
          2 * a1 * hx + hx^2))), length(z[1, ])))
        
        
        colSums(gd, na.rm = TRUE)
      }
    }
    
    # arcsinh-burr8
    if (sd == "burr8") {
      grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep(((-((-log(tan((pi * y)/2)) + hx)/(1 + (log(tan((pi * 
          y)/2)) - hx)^2/exp(2 * gz))) - exp(gz)/sqrt(1 + (log(tan((pi * 
          y)/2)) - hx)^2/exp(2 * gz)) + (2 * exp(asinh((log(tan((pi * y)/2)) - 
          hx)/exp(gz)) + gz))/((1 + exp(asinh((log(tan((pi * y)/2)) - hx)/exp(gz)))) * 
          sqrt(1 + (log(tan((pi * y)/2)) - hx)^2/exp(2 * gz)))))/exp(2 * 
          gz), length(x[1, ])), z * rep((-1 + (log(tan((pi * y)/2)) - hx)^2/(exp(2 * 
          gz) * (1 + (log(tan((pi * y)/2)) - hx)^2/exp(2 * gz))) + (2 * exp(asinh((log(tan((pi * 
          y)/2)) - hx)/exp(gz)) - gz) * (log(tan((pi * y)/2)) - hx))/((1 + 
          exp(asinh((log(tan((pi * y)/2)) - hx)/exp(gz)))) * sqrt(1 + (log(tan((pi * 
          y)/2)) - hx)^2/exp(2 * gz))) + (-log(tan((pi * y)/2)) + hx)/(exp(gz) * 
          sqrt(1 + (log(tan((pi * y)/2)) - hx)^2/exp(2 * gz)))), length(z[1, 
          ])))
        
        colSums(gd, na.rm = TRUE)
      }
    }
    
  }
  
  # t2-XX
  if (fd == "t2") {
    # t2-Cauchy
    if (sd == "cauchy") {
      grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep((-(3 * sqrt(y) * (-sqrt(2 - 2 * y) + 2 * sqrt(2 - 
          2 * y) * y - 2 * sqrt(y) * hx + 2 * y^(3/2) * hx))/(-1 + 4 * sqrt(2 - 
          2 * y) * y^(3/2) * hx - 2 * sqrt(2) * sqrt((-(-1 + y)) * y) * hx - 
          2 * y * (-2 + 2 * exp(2 * gz) + hx^2) + 2 * y^2 * (-2 + 2 * exp(2 * 
          gz) + hx^2))), length(x[1, ])), z * rep(2 - (12 * exp(2 * gz) * 
          (-1 + y) * y)/(-1 + 4 * sqrt(2 - 2 * y) * y^(3/2) * hx - 2 * sqrt(2) * 
          sqrt((-(-1 + y)) * y) * hx - 2 * y * (-2 + 2 * exp(2 * gz) + hx^2) + 
          2 * y^2 * (-2 + 2 * exp(2 * gz) + hx^2)), length(z[1, ])))
        
        colSums(gd, na.rm = TRUE)
      }
    }
    
    
    # t2-t2
    if (sd == "t2") {
      grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep((-(3 * sqrt(y) * (-sqrt(2 - 2 * y) + 2 * sqrt(2 - 
          2 * y) * y - 2 * sqrt(y) * hx + 2 * y^(3/2) * hx))/(-1 + 4 * sqrt(2 - 
          2 * y) * y^(3/2) * hx - 2 * sqrt(2) * sqrt((-(-1 + y)) * y) * hx - 
          2 * y * (-2 + 2 * exp(2 * gz) + hx^2) + 2 * y^2 * (-2 + 2 * exp(2 * 
          gz) + hx^2))), length(x[1, ])), z * rep(2 - (12 * exp(2 * gz) * 
          (-1 + y) * y)/(-1 + 4 * sqrt(2 - 2 * y) * y^(3/2) * hx - 2 * sqrt(2) * 
          sqrt((-(-1 + y)) * y) * hx - 2 * y * (-2 + 2 * exp(2 * gz) + hx^2) + 
          2 * y^2 * (-2 + 2 * exp(2 * gz) + hx^2)), length(z[1, ])))
        
        colSums(gd, na.rm = TRUE)
      }
    }
    
    # t2-burr7
    if (sd == "burr7") {
      grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep((-(3 * (atanh(1 - 2 * y) + hx))/(2 * exp(2 * 
          gz) + atanh(1 - 2 * y)^2 + 2 * atanh(1 - 2 * y) * hx + hx^2)), 
          length(x[1, ])), z * rep(2 * (1 - (3 * exp(2 * gz))/(2 * exp(2 * 
          gz) + atanh(1 - 2 * y)^2 + 2 * atanh(1 - 2 * y) * hx + hx^2)), 
          length(z[1, ])))
        
        colSums(gd, na.rm = TRUE)
      }
    }
    
    # t2-burr8
    if (sd == "burr8") {
      grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep((3 * (log(tan((pi * y)/2)) - hx))/(2 * exp(2 * 
          gz) + log(tan((pi * y)/2))^2 - 2 * hx * log(tan((pi * y)/2)) + 
          hx^2), length(x[1, ])), z * rep(2 * (1 - (3 * exp(2 * gz))/(2 * 
          exp(2 * gz) + log(tan((pi * y)/2))^2 - 2 * hx * log(tan((pi * y)/2)) + 
          hx^2)), length(z[1, ])))
        
        colSums(gd, na.rm = TRUE)
      }
      
    }
  }
  
  # burr8-XX
  if (fd == "burr8") {
    # burr8 Cauchy
    if (sd == "cauchy") {
      grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep(((-1) * ((-1 + exp((2 * (cot(pi * y) + hx))/exp(gz))))/(exp(gz) * 
          (1 + exp((2 * (cot(pi * y) + hx))/exp(gz))))), length(x[1, ])), 
          z * rep(-(((exp(gz) + exp((2 * cot(pi * y) + 2 * hx + exp(gz) * 
          gz)/exp(gz)) - (-1 + exp((2 * (cot(pi * y) + hx))/exp(gz))) * 
          cot(pi * y) + hx - exp((2 * (cot(pi * y) + hx))/exp(gz)) * hx))/(exp(gz) * 
          (1 + exp((2 * (cot(pi * y) + hx))/exp(gz))))), length(z[1, ])))
        
        colSums(gd, na.rm = TRUE)
      }
    }
    
    # burr8 t2
    if (sd == "t2") {
      grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        
        for (i in 1:length(y)) {
          c1 <- ((y[i] <= 0) | (y[i] >= 1))  #1st situation
          if (c1) {
          gd <- NA
          }
        }
        gd <- cbind(x * rep((-exp(-gz)) * tanh(((1 - 2 * y)/(sqrt(2) * sqrt((-(-1 + 
          y)) * y)) + hx)/exp(gz)), length(x[1, ])), z * rep(-(((2 * exp(gz) * 
          sqrt((-(-1 + y)) * y) + tanh(((1 - 2 * y)/(sqrt(2) * sqrt((-(-1 + 
          y)) * y)) + hx)/exp(gz)) * (-sqrt(2) + 2 * sqrt(2) * y - 2 * sqrt((-(-1 + 
          y)) * y) * hx)))/(exp(gz) * (2 * sqrt((-(-1 + y)) * y)))), length(z[1, 
          ])))
        
        
        colSums(gd, na.rm = TRUE)
      }
    }
    
    # burr8 burr7
    if (sd == "burr7") {
      grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep((-exp(-gz)) * tanh((atanh(1 - 2 * y) + hx)/exp(gz)), 
          length(x[1, ])), z * rep((-exp(-gz)) * (exp(gz) - atanh(1 - 2 * 
          y) * tanh((atanh(1 - 2 * y) + hx)/exp(gz)) - hx * tanh((atanh(1 - 
          2 * y) + hx)/exp(gz))), length(z[1, ])))
        colSums(gd, na.rm = TRUE)
      }
    }
    
    # burr8 burr8
    if (sd == "burr8") {
      grad <- function(h, y, x, z) {
        hx <- x %*% h[1:length(x[1, ])]
        gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
        gd <- cbind(x * rep((-1) * (((exp((2 * hx)/exp(gz)) - tan((pi * y)/2)^(2/exp(gz))))/(exp(gz) * 
          (exp((2 * hx)/exp(gz)) + tan((pi * y)/2)^(2/exp(gz))))), length(x[1, 
          ])), z * rep(-(((exp((2 * hx)/exp(gz)) * (exp(gz) + log(tan(pi * 
          y/2)) - hx) + tan((pi * y)/2)^(2/exp(gz)) * (exp(gz) - 2 * log(tan((pi * 
          y)/2)) + log(tan(pi * y/2)) + hx)))/(exp(gz) * (exp((2 * hx)/exp(gz)) + 
          tan((pi * y)/2)^(2/exp(gz))))), length(z[1, ])))
        colSums(gd, na.rm = TRUE)
      }
    }
  }
  
  if (fd == "km" | sd == "km") {
    grad <- function(h, y, x, z) {
      hx <- x %*% h[1:length(x[1, ])]
      gz <- z %*% h[length(x[1, ]) + 1:length(z[1, ])]
      gd <- cbind(x * rep(((-1 + y^exp(hx) + exp(hx) * (-1 + exp(gz) * y^exp(hx)) * 
        log(y)))/(-1 + y^exp(hx)), length(x[1, ])), z * rep((1 + exp(gz) * 
        log(1 - y^exp(hx))), length(z[1, ])))
      colSums(gd, na.rm = TRUE)
    }
  }
  
  
  grad
} 
