pistonSimulation <-
function(m = 60, s = 0.02, v0 = 0.01, 
  k = 5000, p0 = 110000, t = 296, t0 = 360, each = 50, seed = NA, check = TRUE){
  if(check){
    if(any(m < 30) || any(m > 60))
      stop("Piston weight m is out of range, [30, 60] kg")
    
    if(any(s < 0.005) || any(s > 0.02))
      stop("Piston surface area s is out of range, [0.005, 0.020] m^2")
    
    if(any(v0 < 0.002) || any(v0 > 0.01))
      stop("Value of initial gas volume v0 is out of range, [0.002, 0.010] m^3")
    
    if(any(k < 1000) || any(k > 5000))
      stop("Value of spring coefficient k is out of range, [1000, 5000] N/m")
    
    if(any(p0 < 90000) || any(p0 > 110000))
      stop("Value of atmospheric pressure p0 is out of range, [90000, 110000] N/m^2") 
    
    if(any(t < 290) || any(t > 296))
      stop("Value of ambient temperature t is out of range, [290, 296] K")
    
    if(any(t0 < 340) || any(t0 > 360))
      stop("Value of filling gas temperature t0 is out of range, [340, 360] K")
  }
  
  Lm <- length(m)
  if(any(c(length(s) != Lm, length(v0) != Lm, length(k) != Lm, 
         length(p0) != Lm, length(t) != Lm, length(t0) != Lm)))
    stop("m, s, v0, k, p0, t, t0 should be of same length")
  rm(Lm)
  
  m <-  rep(m,  each = each)
  s <-  rep(s,  each = each)
  v0 <- rep(v0, each = each)
  k <-  rep(k,  each = each)
  p0 <- rep(p0, each = each)
  t <- rep(t, each = each)
  t0 <- rep(t0, each = each)
  n <- length(m)
  
  X <- matrix(NA, 14, n)
  Sd <- c(0.1, 0.01, 0.00025, 50, 0.01, 0.13, 0.13)
  if(!is.na(seed))
    set.seed(seed)
  # .runifsum(n, min -3, max = 3, k = 6) * sqrt(2) * Sd[j]
  # is quite close to rnorm(n, 0, Sd[j])
  # for(col in 1:n){
  #   for(row in 1:7){
  #     X[row, col] <- .runifsum(1, min = -3, max = 3, k = 6) * sqrt(2) * Sd[row]      
  #   }
  # }
  for(row in 1:7){
    X[row, ] <- rnorm(n, mean=0, sd=Sd[row])
  }

  X[8,] <- m
  X[9,] <- s
  X[10,] <- v0
  X[11,] <- k
  X[12,] <- p0
  X[13,] <- t
  X[14,] <- t0

  res <- as.data.frame(t(apply(X, 2, .tCycle)))
  class(res) <- c(class(res), "mistatSimulation", "pistonSimulation")
  return(res)
}


# a random number from uniform sum distribution [0, k]
.runifs <- function(k){
  return(sum(runif(k)))
}
# n random numbers from uniform sum distribution [min, max], k is parameter n
.runifsum <- function(n, min = 0, max = 1, k = 6){
  if(max <= min || n < 1 || k < 1)
    stop("wrong parameters")
  
  x <- replicate(n, .runifs(k = k))
  return(((x/k)*(max-min))+min)
}

.tCycle <- function(x){
  Ms <- x[8] + x[1] # a value and its error
  Ss <- x[9] + x[2]
  if(Ss < 0)
    Ss <- 0.00001
  V0s <- x[10] + x[3]
  if(V0s < 0)
    V0s <- 0.001
  Ks <- x[11] + x[4]
  P0s <- x[12] + x[5]
  if(P0s < 0)
    P0s <- 0.001
  Tms <- x[13] + x[6]
  T0s <- x[14] + x[7]
  Mg <- Ks*0.2 #X0
  A <- (P0s*Ss) + (2*Mg) - (Ks*V0s/Ss)
  V <- Ss*(sqrt(((A)^2) + (4*Ks*P0s*V0s*Tms/T0s))-A)/(2*Ks)
  res <- 2*pi * sqrt(Ms / (Ks + ((Ss ^ 2) * P0s * V0s * Tms) / (T0s * V * V)))
  res <- c(x[8:14], res)
  names(res) <- c("m", "s", "v0", "k", "p0", "t", "t0", "seconds") 
  return(res)
}
