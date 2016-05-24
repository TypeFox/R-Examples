#
# theoretical optimal solution for double dive
#
#  $Revision: 1.6 $   $Date: 2011/08/10 02:13:42 $
#
# Originally written for Baddeley & Bassom paper
#
bestdoubledive <- function(d1, d2, surface=30, verbose=FALSE,
                           model="D", asdive=TRUE) {
  if(is.character(model))
    model <- pickmodel(model)
  else
    stopifnot(inherits(model, "hm"))
  m <- model
  # case (a1): first dive omitted
  t2 <- ndl(d2, model=m)
  out <- data.frame(case="a1", t1=0, t2=t2, phi=d2*t2,
                    stringsAsFactors=FALSE)
  # case (b2): stationary point
  Q0 <- 0.79
  Q1 <- 0.79 * d1/10
  Q2 <- 0.79 * d2/10
  halft <- param(m, "N2", "HalfT")
  M0    <- param(m, "N2", "M0")
  ratecoefs <- log(2)/halft
  Q1i <- Q1 * exp(- surface * ratecoefs)
  rhs1 <- Q1i * (d1-d2)/(d1 * (Q1i - Q2))
  t1 <- t2 <- rep(NA, length(ratecoefs))
  sane <- is.finite(rhs1) & (rhs1 > 0)
  t1[sane] <- log(rhs1[sane])/ratecoefs[sane]
  rhs2 <- - (Q1i-Q2) * d2/((d1-d2) * (M0 - Q0 - Q2))
  sane <- is.finite(rhs2) & (rhs2 > 0)
  t2[sane] <- log(rhs2[sane])/ratecoefs[sane]
  ok <- !is.na(t1) & !is.na(t2) & t1 >= 0 & t2 >= 0
  if(any(ok)) {
    if(verbose)
      cat(paste("Case (b2):\nCandidates for stationary point: tissue",
                paste(which(ok), collapse=", "), "\n"))
    for(j in seq(along=t1)[ok]) {
      T1 <- t1[j]
      T2 <- t2[j]
      enddive1 <- Q0 + Q1 - Q1 * exp(-ratecoefs*T1)
      enddive2 <- Q0 + Q2 + (Q1i - Q2) * exp(-ratecoefs * T2) -
                     Q1 * exp(-ratecoefs * (T1+T2+surface))
      valid1 <- (enddive1 <= M0)
      valid2 <- (enddive2 <= M0)
      ok[j] <- ok[j] && all(valid1[-j]) && all(valid2[-j])
    }
  }
  if(any(ok)) {
    if(verbose)
      cat(paste("Tissues with stationary point in feasible interval:",
                paste(which(ok), collapse=", "), "\n"))
    t1 <- t1[ok]
    t2 <- t2[ok]
    phi <- d1*t1 + d2*t2
    best <- which.max(phi)
    out <- rbind(out,
                 data.frame(case="b2", t1=t1[best], t2=t2[best], phi=max(phi),
                            stringsAsFactors=FALSE))
  }
  # case (c2): two controlling tissues on second dive
  # ..................................................
  # Determine which tissues can be controlling
  ntissues <- nrow(as.data.frame(m))
  excess <- M0 - Q0
  crosspoint <- - log(1-ifelse(excess < Q1i, excess/Q1i, NA))/ratecoefs
  mint1 <- maxt1 <- numeric(ntissues)
  ifA <- (excess < Q2) & (excess >= Q1i)
  ifB <- (excess < Q2) & (excess < Q1i)
  ifC <- (excess > Q2) & (excess >= Q1i)
  ifD <- (excess > Q2) & (excess < Q1i)
  # under condition A, tissue constraint eventually bites (for some t2)
  # no matter the value of t1
  mint1[ifA] <- 0
  maxt1[ifA] <- Inf
  # under condition B, tissue constraint eventually bites (for some t2)
  # if t1 is sufficiently small. 
  mint1[ifB] <- 0
  maxt1[ifB] <- crosspoint[ifB]
  # no solutions under condition C
  mint1[ifC] <- maxt1[ifC] <- NA
  # under condition D, tissue constraint eventually bites (for some t2)
  # if t1 is sufficiently large
  mint1[ifD] <- crosspoint[ifD]
  maxt1[ifD] <- Inf
  #
  if(verbose) {
    cat("Case (c2):\n")
    cat("Range of relevant values of t1 for each tissue:\n")
    df <- data.frame(tissue=1:ntissues,
                     type=ifelse(ifA, "a",
                       ifelse(ifB, "B",
                              ifelse(ifC, "C", "D"))),
                     t1min=mint1,
                     t1max=maxt1)
    print(df[complete.cases(df),])
    cat("\n")
  }
    
  # constrain t1 to NDL
  ndl1 <- ndl(d1, model=m)
  maxt1 <- pmin(maxt1, ndl1)
  if(any(missed <- !is.na(mint1) & !is.na(maxt1) & (mint1 > maxt1)))
    mint1[missed] <- maxt1[missed] <- NA
  #
  relevant <- seq(ntissues)[!is.na(mint1) & !is.na(maxt1)]
  if(length(relevant) > 1) {
    if(verbose) {
      cat(paste("Restricting to [0,T] where T = NDL(d1) = ",
                round(ndl1,4), "\n"))
      cat("Search domain of t1 values for each relevant tissue:\n")
      df2 <- data.frame(tissue=relevant,
                        t1min=mint1[relevant],
                        t1max=maxt1[relevant])
      print(df2)
      cat("\n")
    }
    objective <- function(x, M0i, M0j, Q0, Q1i, Q1j, Q2, ki, kj) {
      -(1/ki) * log((M0i-Q0-Q2)/(Q1i * (1 - exp(-x*ki)) - Q2)) + 
       (1/kj) * log((M0j-Q0-Q2)/(Q1j * (1 - exp(-x*kj)) - Q2))
    }
    ff <- function(x, M0i, Q0, Q1i, Q2, ki) {
      -(1/ki) * log((M0i-Q0-Q2)/(Q1i * (1 - exp(-x*ki)) - Q2))
    }
    # loop through pairs of tissues
    for(i in relevant) {
      jj <- relevant[relevant > i]
      if(length(jj) != 0) {
        for(j in jj) {
          # determine whether the domains intersect
          lo <- max(mint1[c(i,j)])
          hi <- min(maxt1[c(i,j)])
          if(lo < hi) {
            # creep inward to avoid NA's
#            dd <- (hi-lo)/500
#            lo <- lo+dd
#            hi <- hi-dd
            # test whether graphs cross
            olo <- objective(lo,
                             M0i=M0[i],
                             M0j=M0[j],
                             Q0=Q0,
                             Q1i=Q1i[i],
                             Q1j=Q1i[j],
                             Q2=Q2,
                             ki=ratecoefs[i],
                             kj=ratecoefs[j])
            ohi <- objective(hi,
                             M0i=M0[i],
                             M0j=M0[j],
                             Q0=Q0,
                             Q1i=Q1i[i],
                             Q1j=Q1i[j],
                             Q2=Q2,
                             ki=ratecoefs[i],
                             kj=ratecoefs[j])
            if(ohi * olo < 0) {
              # OK, find solution
              if(verbose)
                cat(paste("Crossing between tissues", i, "and", j, "\n"))
              roo <- uniroot(objective, c(lo, hi),
                             M0i=M0[i],
                             M0j=M0[j],
                             Q0=Q0,
                             Q1i=Q1i[i],
                             Q1j=Q1i[j],
                             Q2=Q2,
                             ki=ratecoefs[i],
                             kj=ratecoefs[j])
              t1 <- roo$root
              t2 <- ff(t1,
                       M0i=M0[i],
                       Q0=Q0,
                       Q1i=Q1i[i],
                       Q2=Q2,
                       ki=ratecoefs[i])
              if(verbose)
                cat(paste("\tt1 =", round(t1, 4), ", t2 =", round(t2,4), "\n"))
              # check that proposed solution satisfies no-deco constraints
              if(t1 > ndl(d1, model=m))
                stop("Internal error: t1 exceeds ndl in case (c2)")
              firstdive <- dive(ascent(time=0.001), descent(time=0.001),
                                c(d1, t1), c(0, surface))
              t2max <- ndl(d2, model=m, prevstate=haldane(firstdive, model=m))
              #
              if(t2 <= t2max) {
                if(verbose) cat("\tfeasible\n")
                # OK, found a feasible optimum
                out <- rbind(out,
                             data.frame(case="c2",
                                        t1=t1, t2=t2, phi=d1*t1+d2*t2,
                                        stringsAsFactors=FALSE))
              } else {
                if(verbose)
                  cat(paste("\tinfeasible: t2 =", round(t2,3),
                            ">", round(t2max, 3), " = NDL(d2, x(t1+s))\n"))
              }
            } 
          }
        }
      }
    }
  }
  if(verbose) {
    cat("\n")
    nrel <- length(relevant)
    nexam <- nrel * (nrel-1)/2
    cat(paste(nexam, "possible constraint crossings examined.\n"))
    nc2 <- sum(out$case == "c2")
    if(nc2 == 0) {
      cat("No constraint crossings are feasible optima.\n")
    } else if(nc2 == 1) {
      cat("One constraint crossing is a feasible optimum.\n")
    } else 
      cat(paste(nc2, "constraint crossings are feasible optima.\n"))
    cat("End case (c2).\n\n")
  }
  # case (c3): NDL on both dives
  t1 <- ndl(d1, model=m)
  firstdive <- dive(ascent(time=0.001), descent(time=0.001),
                      c(d1, t1), c(0, surface))
  t2 <- ndl(d2, model=m, prevstate=haldane(firstdive, model=m))
  out <- rbind(out,
               data.frame(case="c3", t1=t1, t2=t2, phi=d1*t1+d2*t2,
                          stringsAsFactors=FALSE))
  if(verbose) cat("Done.\n\n") else {
    # Find THE optimal row
    ibest <- which.max(out$phi)
    out <- out[ibest, , drop=FALSE]
    if(asdive) 
      out <- dive(c(d1, out$t1), c(0, surface), c(d2, out$t2))
  }
  #
  return(out)
}

