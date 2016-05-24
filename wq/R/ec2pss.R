ec2pss <-
function (ec, t, p = 0) {
	# Define conductivity ratio
  R <- ec/42.914
    
  # Estimate temperature correction (valid for -2 < t < 35)
  c <- c(0.6766097, 0.0200564, 0.0001104259, -6.9698e-07, 1.0031e-09)
  rt <- c[1] + c[2] * t + c[3] * t^2 + c[4] * t^3 + c[5] * t^4
    
	# Estimate pressure correction (validity range varies with t and S)
  d <- c(0.03426, 0.0004464, 0.4215, -0.003107)
  e <- c(2.07e-05, -6.37e-10, 3.989e-15)
  Rp <- 1 + p * (e[1] + e[2] * p + e[3] * p^2)/(1 + d[1] * t + 
    d[2] * t^2 + (d[3] + d[4] * t) * R)
       
	# Estimate salinity (valid for 2 < S < 42 and -2 < t < 35).       
  Rt <- R/(Rp * rt)
  a <- c(0.008, -0.1692, 25.3851, 14.0941, -7.0261, 2.7081)
  b <- c(5e-04, -0.0056, -0.0066, -0.0375, 0.0636, -0.0144)
  ft <- (t - 15)/(1 + 0.0162 * (t - 15))
  S <- a[1] + a[2] * Rt^0.5 + a[3] * Rt + a[4] * Rt^1.5 + a[5] * 
    Rt^2 + a[6] * Rt^2.5 + ft * (b[1] + b[2] * Rt^0.5 + b[3] * 
    Rt + b[4] * Rt^1.5 + b[5] * Rt^2 + b[6] * Rt^2.5)
        
	# Estimate salinity correction for S < 2
  x <- 400 * Rt
  y <- 100 * Rt
  ifelse(S >= 2, S, S - a[1]/(1 + 1.5 * x + x^2) - b[1] * ft/(1 + 
    y^0.5 + y + y^1.5))
}