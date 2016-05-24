a <- function(t) {
  a <- get("aaa")(t)
}

b <- function(t) {
  b <- get("bbb")(t)
}
S <- function(t) {
  S <- get("SSS")(t)
}

Sp <- function(t) {
  Sp <- get("SSSp")(t)
}

cc <- function(t) {
  s2 <- get("sigma2")
    cc <- s2 + 0.0*t
}

a1 <- function(x,t) {
  a1 <- a(t)*x + b(t)
}

a2 <- function(t) {
  a2 <- cc(t)
}


mdt <- function(t,y,tau) {
  mean.vector <- get("mp")
  v.vector <- get("vp")
  t0 <- get("t0")
  i <- floor((t-t0)/deltat+1.001)
  j <- floor((tau-t0)/deltat+1.001)
  mdt <- mean.vector[i] + v.vector[i] * (y - mean.vector[j])/v.vector[j]
}


vdt <- function(t,tau) {
  u.vector <- get("up")
  v.vector <- get("vp")
  t0 <- get("t0")
  i <- floor((t-t0)/deltat+1.001)
  j <- floor((tau-t0)/deltat+1.001)
  vdt <- v.vector[i] *(u.vector[i] * v.vector[j] - u.vector[j]*v.vector[i])/v.vector[j]
}

fdt <- function(x,t,y,tau) {
  temp1 <- (x-mdt(t,y,tau))^2
  temp2 <- 2.0*vdt(t,tau)
  temp3 <- sqrt(pi*temp2)
  fdt <- exp(-temp1/temp2)/temp3
}

psi <- function(t,y,tau) {
  temp1 <- Sp(t) - a1(S(t),t)
  temp2 <- (S(t) - mdt(t,y,tau))/vdt(t,tau)
  f1 <- temp1 - a2(t) * temp2
  f2 <- fdt(S(t),t,y,tau)
  psi <- f1*f2
}
