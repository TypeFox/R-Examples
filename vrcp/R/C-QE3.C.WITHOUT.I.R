# Continuous QE3
# Without initial values
# Common variance

llsearch.QE3.CC.WITHOUT.I <- function(x, y, n, jlo, jhi){
  
  loglik <- array(NA,n)
  s2_list <- array(NA,n)
  
  for (j in jlo:jhi){ 
    
    est_x <- x[1:j]
    est_y <- y[1:j]
    est_fun <- lm(est_y ~ est_x + I(est_x^2))
    est_a0 <- est_fun$coe[1] 
    est_a1 <- est_fun$coe[2]
    est_a2 <- est_fun$coe[3]
    ypred <- est_a0 + est_a1 * x[j] + est_a2 * x[j]^2
    
    seg2 <- data.frame(x[(j+1):n],y[(j+1):n])
    names(seg2)=c("est_x2","est_y2")
    
    positive_y2_value <- length(y[y[(j+1):n] > ypred])
    negative_y2_value <- length(y[y[(j+1):n] < ypred])
    
    # determine the direction of the 2nd segment
    if (positive_y2_value >= negative_y2_value){
      est_y2 <- seg2$est_y2[seg2$est_y2 > ypred]
      est_x2 <- seg2$est_x2[seg2$est_y2 > ypred]
      est_log_y2 <- log(est_y2-ypred)
    }
    else{
      est_y2 <- seg2$est_y2[seg2$est_y2 < ypred]
      est_x2 <- seg2$est_x2[seg2$est_y2 < ypred]
      est_log_y2 <- log(-est_y2+ypred)
    }
    
    est_x2_shift <- est_x2 - x[j]
    est_log_fun <- lm(est_log_y2 ~ est_x2_shift)
    
    est_b2 <- est_log_fun$coe[2]
    if (positive_y2_value < negative_y2_value){
      est_b1 <- -exp(est_log_fun$coe[1])
    }
    else{
      est_b1 <- exp(est_log_fun$coe[1]) # for shifted model 
    }
    est_b0 <- ypred - est_b1 # continuous # for shifted model 
    
    s2 <- (sum((y[1:j] - est_a0 - est_a1 * x[1:j] - est_a2 * x[1:j]^2)^2) +
             sum((y[(j+1):n] - est_b0 - est_b1 * exp(est_b2 * (x[(j+1):n]-x[j])))^2))/n
    q1 <- n * log(sqrt(2 * pi))
    q2 <- 0.5 * n  * (1 + log(s2))
    
    s2_list[j]<- s2
    loglik[j] <- -(q1+q2)
    est_j <- which.max(loglik) # find the change point
  }
  
  est_x <- x[1:est_j]
  est_y <- y[1:est_j]
  est_fun <- lm(est_y ~ est_x + I(est_x^2))
  est_a0 <- est_fun$coe[1] # initial a0
  est_a1 <- est_fun$coe[2] # initial a1
  est_a2 <- est_fun$coe[3] # initial a2
  ypred <- est_a0 + est_a1 * x[est_j] + est_a2 * x[est_j]^2
  
  seg2_with_j <- data.frame(x[(est_j+1):n],y[(est_j+1):n])
  names(seg2_with_j)=c("est_x2_with_j","est_y2_with_j")
  
  positive_y2_value <- length(y[y[(est_j+1):n] > ypred])
  negative_y2_value <- length(y[y[(est_j+1):n] < ypred])
  
  if (positive_y2_value >= negative_y2_value){
    est_y2_with_j <- seg2_with_j$est_y2_with_j[seg2_with_j$est_y2_with_j > ypred]
    est_x2_with_j <- seg2_with_j$est_x2_with_j[seg2_with_j$est_y2_with_j > ypred]
    est_log_y2_with_j <- log(est_y2_with_j-ypred) 
  }
  if (positive_y2_value < negative_y2_value){
    est_y2_with_j <- seg2_with_j$est_y2_with_j[seg2_with_j$est_y2_with_j < ypred]
    est_x2_with_j <- seg2_with_j$est_x2_with_j[seg2_with_j$est_y2_with_j < ypred]
    est_log_y2_with_j <- log(-est_y2_with_j+ypred) 
  } 
  
  est_x2_shift_with_j <- est_x2_with_j - x[est_j]
  est_log_fun_with_j <- lm(est_log_y2_with_j ~ est_x2_shift_with_j)
  
  est_b2 <- est_log_fun_with_j$coe[2] # initial b2
  
  if (positive_y2_value >= negative_y2_value){
    est_b1 <- exp(est_log_fun_with_j$coe[1])
  }
  if (positive_y2_value < negative_y2_value){
    est_b1 <- -exp(est_log_fun_with_j$coe[1])
  }
  
  est_b0 <- ypred - est_b1 * (exp(est_b2*(x[est_j]-x[est_j]))) ## continuous
  
  list(a0=est_a0,a1=est_a1,a2 = est_a2, b0=est_b0,b1=est_b1,b2=est_b2,
       maxloglik = loglik[est_j],
       sigma2=s2_list[est_j],xj=x[est_j])
}
