qnorm.app3 <- function(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE) {
  ## 26.2.22 in Abramowitz and Segun, Dover, 1968, |error| < 0.003
  a0 <- 2.30753e0
  a1 <- 0.27061e0
  b1 <- 0.99229e0
  b2 <- 0.04481e0
  vorz <- rep(1,length(p))
  if (log.p) p <- exp(p)
  vorz <- ifelse((p < 0.5e0) == lower.tail, vorz, -vorz)
  p    <- ifelse(p < 0.5e0, p, 1 - p)
  t    <- ifelse(p <= 0, NaN, sqrt(-2.0e0*log(p)))
  ifelse(p <= 0, NaN, mean - (t - (a1*t +a0)/((b2*t + b1)*t + 1))*sd*vorz)
}

qnorm.app4 <- function(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE) {
  ## 26.2.23 in Abramowitz and Segun, Dover, 1968, |error| < 0.00045
  c0 <- 2.515517e0
  c1 <- 0.802853e0
  c2 <- 0.010328e0
  d1 <- 1.432788e0
  d2 <- 0.189269e0
  d3 <- 0.001308e0
  vorz <- rep(1,length(p))
  if (log.p) p <- exp(p)
  vorz <- ifelse((p < 0.5e0) == lower.tail, vorz, -vorz)
  p    <- ifelse(p < 0.5e0, p, 1 - p)
  t    <- ifelse(p <= 0, NaN, sqrt(-2.0e0*log(p)))
  ifelse(p <= 0, NaN, mean - (t - ((c2*t + c1)*t +c0)/(((d3*t + d2)*t + d1)*t + 1))*sd*vorz)
}

qnorm.ap16 <- function(p, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE) {
  # algorithm AS241  Applied Statistics (1988) Vol. 37, No. 3
  half   <- 0.5e0
  split1 <- 0.425e0
  split2 <- 5.0e0
  const1 <- 0.180625e0
  const2 <- 1.6e0
#   coefficients for p close to 1/2
  a0 <- 3.3871328727963666080e0
  a1 <- 1.3314166789178437745e2
  a2 <- 1.9715909503065514427e3
  a3 <- 1.3731693765509461125e4
  a4 <- 4.5921953931549871457e4
  a5 <- 6.7265770927008700853e4
  a6 <- 3.3430575583588128105e4
  a7 <- 2.5090809287301226727e3
  b1 <- 4.2313330701600911252e1
  b2 <- 6.8718700749205790830e2
  b3 <- 5.3941960214247511077e3
  b4 <- 2.1213794301586595867e4
  b5 <- 3.9307895800092710619e4
  b6 <- 2.8729085735721942674e4
  b7 <- 5.2264952788528545610e3
  # a0+a1+a2+a3+a4+a5+a6+a7+b1+b2+b3+b4+b5+b6+b7
# hash sum ab = 55.8831928806149014439
  #             55.883192880614892
  #                             014448 OK
 
# coefficients for p neither close to 1/2 nor 0 nor 1
  c0 <- 1.42343711074968357734e0
  c1 <- 4.63033784615654529590e0
  c2 <- 5.76949722146069140550e0
  c3 <- 3.64784832476320460504e0
  c4 <- 1.27045825245236838258e0
  c5 <- 2.41780725177450611770e-1
  c6 <- 2.27238449892691845833e-2
  c7 <- 7.74545014278341407640e-4
  d1 <- 2.05319162663775882187e0
  d2 <- 1.67638483018380384940e0
  d3 <- 6.89767334985100004550e-1
  d4 <- 1.48103976427480074590e-1
  d5 <- 1.51986665636164571966e-2
  d6 <- 5.47593808499534494600e-4
  d7 <- 1.05075007164441684324e-9
  # c0+c1+c2+c3+c4+c5+c6+c7+d1+d2+d3+d4+d5+d6+d7
# hash sum cd = 49.33206503301610289036
  #             49.33206503301611
  #                            10289036
 
# coefficients for p near 0 or 1
  e0 <- 6.65790464350110377720e0
  e1 <- 5.46378491116411436990e0
  e2 <- 1.78482653991729133580e0
  e3 <- 2.96560571828504891230e-1
  e4 <- 2.65321895265761230930e-2
  e5 <- 1.24266094738807843860e-3
  e6 <- 2.71155556874348757815e-5
  e7 <- 2.01033439929228813265e-7
  f1 <- 5.99832206555887937690e-1
  f2 <- 1.46929880922835805310e-1
  f3 <- 1.48753612908506148525e-2
  f4 <- 7.86869131145613259100e-4
  f5 <- 1.84631831751005468180e-5
  f6 <- 1.42151175831644588870e-7
  f7 <- 2.04426310338993978564e-15
  # e0+e1+e2+e3+e4+e5+e6+e7+f1+f2+f3+f4+f5+f6+f7
# hash sum cd = 47.52583317549289671629
  #             47.625833175493895
  #                            89671629
    
  vorz <- rep(1,length(p))
  if (log.p) p <- exp(p)
  Q <-  p - half
  vorz <- ifelse((Q < 0) == lower.tail, vorz, -vorz)
  p    <- ifelse(Q < 0, p, 1 - p)
  Q    <- abs(Q)

  t <- ifelse(Q <= split1, 0, sqrt(-log(p)))
  r <- ifelse(Q <= split1, const1 - Q*Q, ifelse(t <= split2, t - const2, t - split2))

  res <- ifelse (p <= 0,NaN,
           ifelse (Q <= split1,
             Q*(((((((a7*r + a6)*r + a5)*r + a4)*r + a3)*r + a2)*r + a1)*r + a0)/(((((((b7*r + b6)*r + b5)*r + b4)*r + b3)*r + b2)*r + b1)*r + 1),
             ifelse (t <= split2,
               (((((((c7*r + c6)*r + c5)*r + c4)*r + c3)*r + c2)*r + c1)*r + c0)/(((((((d7*r + d6)*r + d5)*r + d4)*r + d3)*r + d2)*r + d1)*r + 1),
               (((((((e7*r + e6)*r + e5)*r + e4)*r + e3)*r + e2)*r + e1)*r + e0)/(((((((f7*r + f6)*r + f5)*r + f4)*r + f3)*r + f2)*r + f1)*r + 1)
             )))
  mean - res*sd*vorz
}
