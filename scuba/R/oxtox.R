#
# 	oxtox.R
#
#	$Revision: 1.5 $	$Date: 2013/12/10 12:58:54 $
#
###################################################################
#  
# Oxygen toxicity and nitrox limits
#

"ead" <-
function(depth, g) {
  if(!is.gas(g))
    g <- nitrox(g)
  if(!is.nitrox(g))
    stop("Not defined for trimix")
  fO2 <- g$fO2
  (10 + depth) * (1- fO2)/0.79 - 10
}

"END" <-
function(depth, g) {
  stopifnot(is.gas(g))
  fHe <- g$fHe
  pmax(0, (10 + depth) * (1- fHe) - 10)
}

"maxmix" <-
function(depth, ppO2max = 1.4) {
  fo <- pmin(1, ppO2max/(1 + depth/10))
  fo <- signif(fo, 2)
  nitrox(fo)
}

"mod" <-
function(g, ppO2max=1.4) {
  if(!is.gas(g))
    g <- nitrox(g)
  fO2 <- g$fO2
  10 * (ppO2max/fO2 - 1)
}

"eadtable" <-
function(g, ppO2max=1.4) {
  if(!is.gas(g))
    g <- nitrox(g)
  if(!is.nitrox(g))
    stop("Not yet implemented for trimix")
  fO2 <- g$fO2
  depth <- 0:40
  EAD <- round(ead(depth, fO2), 1)
  EAD <- ifelse(EAD >= 0, as.character(EAD), " ")
  EAD <- ifelse(depth > mod(fO2, ppO2max), "Warning", EAD)
  data.frame(depth=depth,EAD=EAD)
}

ppO2 <- function(d) {
  times <- times.dive(d)
  depths <- depths.dive(d)
  fO2 <- d$data$fO2
  pp <- fO2 * (depths/10 + 1)
  data.frame(time=times, ppO2=pp)
}

"oxtox" <- 
function(d, progressive=FALSE, warn=TRUE)
{
  times <- times.dive(d)
  depths <- depths.dive(d)
  n <- length(depths)
  fO2 <- d$data$fO2

  # fO2[i] applies to interval (times[i], times[i+1])
  pO2start <- fO2 * (depths/10 + 1)
  pO2end   <- c(fO2[-n] * (depths[-1]/10 + 1), NA)
  pO2high  <- pmax(pO2start, pO2end, na.rm=TRUE)
  maxpO2 <- max(pO2high)
  if(warn) {
    if(maxpO2 > 1.6)
      warning("O2 partial pressure exceeded 1.6")
    else if(maxpO2 > 1.5)
      warning("O2 partial pressure exceeded 1.5")
    else if(maxpO2 > 1.4)
      warning("O2 partial pressure exceeded 1.4")
  }
  
  # determine which intervals contain toxicity contributions
  toxicstart <- (pO2start > 0.5)
  toxicend   <- (pO2end   > 0.5)
  toxic <- (pO2high > 0.5)
  if(!any(toxic)) {
    if(!progressive) return(0) else return(rep(0, n))
  }

  # calculations for each interval
  durations <- diff(times)
  flat <- (diff(depths) == 0)
  powerit <- function(x) { exp(0.83 * log(x)) }
  Powerit <- function(x) { exp(1.83 * log(x)) }

  # vector of toxicity contributions
  dotu <- rep(0, n)
  
  # compute toxicity contributions for each toxic interval
  for(i in seq(n-1)[toxic[-n]]) {
    dura <- durations[i]
    p0 <- pO2start[i]
    p1 <- pO2end[i]
    if(flat[i])
      # integrate the constant (2 *(p02-0.5))^0.83
      dotu[i+1] <- dura * powerit(2 * (p0 - 0.5))
    else if(toxicstart[i] && toxicend[i]) {
      # entire interval is in the toxic zone
      # integrate (a + bx)^0.83)
      a <- 2 * p0 - 1
      b <- 2 * (p1-p0)/dura
      dotu[i+1] <- (1/(1.83 * b)) * (Powerit(a + b * dura) - Powerit(a))
    } else {
      # interval is only partly toxic
      # compute toxic subinterval
      t0 <- times[i]
      t1 <- times[i+1]
      tx <- t0 + (t1-t0) * (0.5-p0)/(p1-p0)
      # restrict to subinterval
      if(p0 < p1) {
        # [tx, t1]
        t0 <- tx
        p0 <- 0.5
        dura <- t1-tx
      } else {
        # [t0, tx]
        t1 <- tx
        p1 <- 0.5
        dura <- tx-t0
      }
      # integrate (a + bx)^0.83)
      a <- 2 * p0 - 1
      b <- 2 * (p1-p0)/dura
      dotu[i+1] <- (1/(1.83 * b)) * (Powerit(a + b * dura) - Powerit(a))
    }
  }
  if(progressive)
    return(cumsum(dotu))
  else
    return(sum(dotu))
}

