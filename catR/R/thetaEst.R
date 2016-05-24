thetaEst<-function (it, x, model = NULL, D = 1, method = "BM", priorDist = "norm", 
    priorPar = c(0, 1), range = c(-4, 4), parInt = c(-4, 4, 33), 
    constantPatt = NULL, current.th = 0, bRange = c(-2, 2)) 
{
    constantPattern <- function(t) ifelse(sum(t) == 0 | sum(t) == 
        length(t), TRUE, FALSE)
    METHOD <- NULL
    if (!constantPattern(x) | !is.null(model) | is.null(constantPatt)) 
        METHOD <- method
    else {
        if (sum(constantPatt == c("BM", "EAP", "WL")) == 1) 
            METHOD <- constantPatt
        else {
            if (sum(x) == 0) 
                res <- switch(constantPatt, fixed4 = current.th - 
                  0.4, fixed7 = current.th - 0.7, var = 0.5 * 
                  (current.th + bRange[1]))
            else res <- switch(constantPatt, fixed4 = current.th + 
                0.4, fixed7 = current.th + 0.7, var = 0.5 * (current.th + 
                bRange[2]))
        }
    }
ind<-which(!is.na(x))
it<-it[ind,]
x<-x[ind]
    if (!is.null(METHOD)) {
        if (METHOD == "EAP") 
            res <- eapEst(it, x, model = model, D = D, priorDist = priorDist, 
                priorPar = priorPar, lower = parInt[1], upper = parInt[2], 
                nqp = parInt[3])
        else {
            if (is.null(model)) {
                r0 <- function(th, it, D = 1, method = "BM", 
                  priorDist = "norm", priorPar = c(0, 1)) {
                  if (method == "BM") 
                    res <- switch(priorDist, norm = (priorPar[1] - 
                      th)/priorPar[2]^2, unif = 0, Jeffreys = sum(Ii(th, 
                      it, D = D)$dIi)/(2 * sum(Ii(th, it, D = D)$Ii)))
                  else res <- switch(method, ML = 0, WL = sum(Ji(th, 
                    it, D = D)$Ji)/(2 * sum(Ii(th, it, D = D)$Ii)))
                  return(res)
                }
                r <- function(th, it, x, D = 1) {
                  pr <- Pi(th, it, D = D)
                  P <- pr$Pi
                  Q <- 1 - P
                  dP <- pr$dPi
                  res <- sum(dP * (x - P)/(P * Q))
                  return(res)
                }
                T <- function(th, it, x, D = 1, method = "BM", 
                  priorDist = "norm", priorPar = c(0, 1)) {
                  r0(th, it, D = D, method = method, priorDist = priorDist, 
                    priorPar = priorPar) + r(th, it, x, D = D)
                }
                if (METHOD == "BM" & priorDist == "unif") 
                  f <- function(th) T(th, it, x, D = D, method = "ML")
                else f <- function(th) T(th, it, x, D = D, method = METHOD, 
                  priorDist = priorDist, priorPar = priorPar)
                if (METHOD == "BM" & priorDist == "unif") 
                  RANGE <- priorPar
                else RANGE <- range
                if (f(RANGE[1]) < 0 & f(RANGE[2]) < 0) 
                  res <- RANGE[1]
                else {
                  if (f(RANGE[1]) > 0 & f(RANGE[2]) > 0) 
                    res <- RANGE[2]
                  else res <- uniroot(f, RANGE)$root
                }
            }
            else {
                met <- switch(METHOD, ML = 1, BM = 2, WL = 3, 
                  EAP = 4)
                pd <- switch(priorDist, norm = 1, unif = 2, Jeffreys = 3)
                if (met == 2 & pd == 2) 
                  RANGE <- c(max(c(priorPar[1], range[1])), min(c(priorPar[2], 
                    range[2])))
                else RANGE <- range
                dl <- function(th) {
                  p <- Pi(th, it, model = model, D = D)
                  pr <- p$Pi
                  dpr <- p$dPi
                  res <- 0
                  for (i in 1:length(x)) res <- res + dpr[i, 
                    x[i] + 1]/pr[i, x[i] + 1]
                  return(res)
                }
                f0 <- function(th) {
                  if (met == 2) 
                    res <- switch(pd, `1` = (priorPar[1] - th)/priorPar[2]^2, 
                      `2` = 0, `3` = sum(Ii(th, it, model = model, 
                        D = D)$dIi)/(2 * sum(Ii(th, it, model = model, 
                        D = D)$Ii)))
                  else res <- switch(met, `1` = 0, `2` = 0, `3` = sum(Ji(th, 
                    it, model = model, D = D)$Ji)/(2 * sum(Ii(th, 
                    it, model = model, D = D)$Ii)))
                  return(res)
                }
                optF <- function(th) dl(th) + f0(th)
                if (optF(RANGE[1]) < 0 & optF(RANGE[2]) < 0) 
                  res <- RANGE[1]
                else {
                  if (optF(RANGE[1]) > 0 & optF(RANGE[2]) > 0) 
                    res <- RANGE[2]
                  else res <- uniroot(optF, RANGE)$root
                }
            }
        }
    }
    return(res)
}

