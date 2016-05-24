# The code is based on the Modified Transformation-Based Kernel Estimator (MTK) 
# for copula densities, which is proposed in Wen and Wu (2015).
# The code is written by Kuangyu Wen, Email: kweneco@gmail.com
# Affiliation: International School of Economics and Management, Capital 
# University of Economics and Business, Beijing, China

eval_tt <- function(uev, udata, bw) {
    # This function returns copula density estimates
    # given the bandwidth matrix H = h^2 * (1, rho \\
    # rho, 1) and the additional tuning parameters
    # theta = (theta_1, theta_2).  udata is pseudo-data
    # calculated from the rescaled rank of the orginal
    # data; udata should be a matrix of dimension N*2.
    # The first arguments of the evaluation points: u
    # The second arguments of the evaluation points: v
    # In this function, u and v should have the same
    # length.
    u <- uev[, 1]
    v <- uev[, 2]
    h <- bw[1]
    rho <- bw[2]
    theta <- bw[3:4]
    
    n <- dim(udata)[1]
    
    Si <- qnorm(udata[, 1])
    Ti <- qnorm(udata[, 2])
    s <- qnorm(u)
    t <- qnorm(v)
    
    delta <- sqrt(h^4 * (1 - rho^2) * (4 * theta[1]^2 -  theta[2]^2) + 
                      2 * h^2 * (2 * theta[1] + rho * theta[2]) + 1)
    eta <- mean(exp(-((4 * h^2 * theta[1]^2 - h^2 * 
                           theta[2]^2 + 2 * theta[1]) * (Si^2 + Ti^2) + 
                          (2 * rho * h^2 * theta[2]^2 - 8 * rho * h^2 * 
                             theta[1]^2 + 2 * theta[2]) * Si * Ti)/2/delta^2))/delta
    
    arg1 <- outer(s, Si, function(p, q) (p - q)/h)
    arg2 <- outer(t, Ti, function(p, q) (p - q)/h)
    temp <- exp(-(arg1^2 + arg2^2 - 2 * rho * arg1 * 
                      arg2)/2/(1 - rho^2))/2/pi/sqrt(1 - rho^2)
    if (is.null(dim(temp))) 
        temp <- mean(temp) else temp <- rowMeans(temp)
    
    temp <- temp/h^2/eta
    temp * exp(-theta[1] * (s^2 + t^2) - theta[2] * 
                   s * t)/dnorm(s)/dnorm(t)
}
