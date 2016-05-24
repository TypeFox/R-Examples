# These functions are
# Copyright (C) 1998-2015 T.W. Yee, University of Auckland.
# All rights reserved.





DeLury <- function(catch, effort,
                   type = c("DeLury", "Leslie"),
                   ricker = FALSE) {
  type <- match.arg(type, c("DeLury", "Leslie"))[1]
  if (!is.logical(ricker))
    stop("bad input for argument 'ricker'")
  if ((LLL <- Lcatch <- length(catch)) != (Leffort <- length(effort)))
    stop("length(catch) != length(effort)")

  CPUE <- catch / effort
  if (type == "DeLury") {
    Et <- cumsum(effort) - ifelse(ricker, 0.5, 1) * effort
    logCPUE <- log(CPUE)
    lmfit <- lm(logCPUE ~ Et, x = TRUE)
    myq <- catchabilityCoefficient <- -coef(lmfit)[2]
    N0 <- exp(coef(lmfit)["(Intercept)"]) / myq
  } else {
    Kt <- cumsum(catch) - ifelse(ricker, 0.5, 1) * catch
    lmfit <- lm(CPUE ~ Kt, x = TRUE)
    myq <- catchabilityCoefficient <- -coef(lmfit)[2]
    N0 <- coef(lmfit)["(Intercept)"] / myq
  }

  rlist <-
  list(catch = catch,
       effort = effort,
       type = type,
       N0 = N0,
       CPUE = CPUE,
       lmfit = lmfit)
  if (type == "DeLury") {
    rlist$E <- Et
  } else {
    rlist$K <- Kt
  }
  rlist
}






wffc.P1     <- function(length, c1 = 100, min.eligible = 0.18, ppm = 2000)
  ifelse(length >= min.eligible, c1 + (ppm/100) *
         ceiling(  signif(100 * length, digits = 8)  ), 0)


wffc.P1star <- function(length, c1 = 100, min.eligible = 0.18, ppm = 2000)
  ifelse(length >= min.eligible, c1 + ppm * length, 0)
















wffc.P2star <- function(length, c1 = 100, min.eligible = 0.18, ppm = 2000,
                        c.quad = 12700)
  wffc.P1star(length, c1 = c1, min.eligible = min.eligible, ppm = ppm) +
  ifelse(length > min.eligible, c.quad * (length - min.eligible)^2, 0)





wffc.P2     <- function(length, c1 = 100, min.eligible = 0.18, ppm = 2000,
                        c.quad = 12700)
  wffc.P2star(ifelse(length > min.eligible,
                     ceiling(100 * length) / 100,
                     length),
              c1 = c1,
              min.eligible = min.eligible,
              ppm = ppm,
              c.quad = c.quad)









wffc.P3star <-
  function(length, c1 = 100, min.eligible = 0.18, ppm = 2000) {
  kay <- floor(length / min.eligible)
  ans <- ifelse(kay >= 1, c1, length * 0)  # Handles NAs

  ans <- ans +
         ifelse(kay >= 1, ppm * min.eligible, 0) +
         ifelse(kay >= 1, ppm * min.eligible * kay*(kay-1)/2, 0) +
         ifelse(kay >= 1, ppm * (length - kay * min.eligible) * kay, 0)


  ans
}



wffc.P3 <- function(length, c1 = 100, min.eligible = 0.18, ppm = 2000) {


  wffc.P3star(ifelse(length > min.eligible,
                     ceiling(100 * length) / 100,
                     length),
              c1 = c1,
              min.eligible = min.eligible,
              ppm = ppm)
}






wffc.P4star <-
  function(length, c1 = 100, min.eligible = 0.18, ppm = 2000) {

  kay <- floor(length / (min.eligible / 2))
  km1 <- kay - 1

  ans <- ifelse(length >= min.eligible, c1, length * 0)  # Handles NAs

  ans <- ans +
       ifelse(km1 >= 1, ppm * min.eligible, 0) +
       ifelse(km1 >= 1, ppm * (min.eligible/2) * km1*(km1-1)/2, 0) +
       ifelse(km1 >= 1, ppm * (length - (km1+1) * min.eligible/2) * km1, 0)
}



wffc.P4     <-
  function(length, c1 = 100, min.eligible = 0.18, ppm = 2000) {


  wffc.P4star(ifelse(length > min.eligible,
                     ceiling(100 * length) / 100,
                     length),
              c1 = c1,
              min.eligible = min.eligible,
              ppm = ppm)
}













