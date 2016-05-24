L_hNV_rho <- function(p, y = y, X = X, sc = sc, w = w, sigmau2_sar = sigmau2_sar) {

    b <- p[4:length(p)]
    rho  <- p[1]
    
    sigmau2_dmu <- p[2]
    sigmav2 <- p[3]
    sigmau2_sar <- sigmau2_sar

    
    DIM_w <- dim(w)[1]
    I_m <- diag(DIM_w)
    epsilon <- (I_m - rho *w)%*%(y - X%*%b)
    sc <- sc
    
    N <- length(y)
    
    ret = (N * log(sqrt(2) / sqrt(pi)) + N * log(1 / (sqrt(sigmau2_dmu + sigmau2_sar + sigmav2)))
            + sum(log(pnorm(-sc*(epsilon * (sqrt(sigmau2_dmu + sigmau2_sar) / sqrt(sigmav2))) / 
                              (sqrt(sigmau2_dmu + sigmau2_sar + sigmav2)))))
            - 1 / (2 * (sigmau2_dmu + sigmau2_sar + sigmav2)) * sum(epsilon^2))
    
    names(ret) <- "Log-Lik normal/half-normal distribution"
    return(ret)
  }
