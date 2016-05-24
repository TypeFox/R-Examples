"epi.cluster2size" <- function(nbar, R = NA, n, mean, sigma2.x, sigma2.y = NA, sigma2.xy = NA, epsilon.r, method = "mean", conf.level = 0.95){
    N. <- 1 - ((1 - conf.level) / 2)
    z <- qnorm(N., mean = 0, sd = 1)
    
    if (method == "total") {
        if (length(n) != 2) 
           stop("Error: n must be of length 2")
        if (length(mean) != 2) 
           stop("Error: mean must be of length 2")
        if (length(sigma2.x) != 2) 
           stop("Error: sigma2.x must be of length 2")
        
        sigma2.1x <- sigma2.x[1]; sigma2.2x <- sigma2.x[2]
        M <- n[1]; N <- n[2]
        Xbar <- mean[1]; xbar <- mean[2]
        
        # Equation 10.6, page 292 Levy and Lemeshow:
        numerator <- (sigma2.1x / Xbar^2) * (M / (M - 1)) + (1 / nbar) * (sigma2.2x / xbar^2) * ((N - nbar) / (N - 1))
        denominator <- (epsilon.r^2 / z^2) + (sigma2.1x / (Xbar^2 * (M - 1)))
        rval <- round(numerator/denominator, digits = 0)
        }
    
    if (method == "mean") {
        if (length(n) != 2) 
           stop("Error: n must be of length 2")
        if (length(mean) != 2) 
           stop("Error: mean must be of length 2")
        if (length(sigma2.x) != 2) 
           stop("Error: sigma2.x must be of length 2")
           
        sigma2.1x <- sigma2.x[1]; sigma2.2x <- sigma2.x[2]
        M <- n[1]; N <- n[2]
        Xbar <- mean[1]; xbar <- mean[2]

        # Equation 10.6, page 292 Levy and Lemeshow:
        numerator <- ((sigma2.1x / Xbar^2) * (M / (M - 1))) + ((1 / nbar) * (sigma2.2x / xbar^2) * ((N - nbar) / (N - 1)))
        denominator <- (epsilon.r^2 / z^2) + (sigma2.1x / (Xbar^2 * (M - 1)))
        rval <- round(numerator/denominator, digits = 0)
        }
       
    if (method == "proportion") {
      if (length(n) != 2) 
        stop("Error: n must be of length 2")
      if (length(mean) != 2) 
        stop("Error: mean must be of length 2")
      if (length(sigma2.x) != 2) 
        stop("Error: sigma2.x must be of length 2")
      if (length(sigma2.y) != 2) 
        stop("Error: sigma2.y must be of length 2")
      if (length(sigma2.xy) != 2) 
        stop("Error: sigma2.xy must be of length 2")
      
      sigma2.1x <- sigma2.x[1]; sigma2.2x <- sigma2.x[2]
      sigma2.1y <- sigma2.y[1]; sigma2.2y <- sigma2.y[2]
      sigma.1xy <- sigma2.xy[1]; sigma.2xy <- sigma2.xy[2]
      M <- n[1]; N <- n[2]
      Xbar <- mean[1]; xbar <- mean[2]
      
      sigmasq.r1 <- sigma2.1y + (R^2 * sigma2.1x) - (2 * R * sigma.1xy)
      sigmasq.r2 <- sigma2.2y + (R^2 * sigma2.2x) - (2 * R * sigma.2xy)
      
      # Equation 10.7, page 292 Levy and Lemeshow:
      numerator <- ((sigmasq.r1 / Xbar^2) * (M / (M - 1))) + ((1 / nbar) * (sigmasq.r2 / xbar^2) * ((N - nbar) / (N - 1)))
      denominator <- (epsilon.r^2 / z^2) + (sigmasq.r1 / (Xbar^2 * (M - 1)))
      rval <- round(numerator/denominator, digits = 0)
       }
    return(rval)
}
