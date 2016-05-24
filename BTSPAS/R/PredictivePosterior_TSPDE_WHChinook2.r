# 2010-03-28 CJS initial creation of function for the second Wild-Hatchery Chinook problem of Eric Logan

PredictivePosterior.TSPDE.WHCH2 <- function (time, n1, m2, 
       u2.A.YoY, u2.N.YoY, u2.A.1, u2.N.1,  clip.frac.H.YoY, clip.frac.H.1, p, 
       U.W.YoY, U.H.YoY, U.W.1, U.H.1, hatch.after.YoY) {
#  Generate Predictive Posterior Plot (Bayesian p-value) given the data
#  for a TimeStratified Petersen with Diagonal Elements and error
#    n1, m2, u2.A.YoY, u2.N.YoY, u2.A.1, u2.N.1  = vectors of input data
#    p, U.W.YoY, U.H.YoY, U.W.1, U.H.1          = matrix of values (rows=number of posterior samples, columns=strata)
#               These are returned from the call to OpenBugs/ WinBugs
#
#cat("Call to PredictivePosterior for Wild vs Hatchery and YoY vs Age1 \n")
#browser()
discrep <- matrix(0, nrow=0, ncol=16)
select.m2 <- !is.na(m2)
select.u2.A.YoY <- !is.na(u2.A.YoY)
select.u2.N.YoY <- !is.na(u2.N.YoY)
select.u2.A.1   <- !is.na(u2.A.1)
select.u2.N.1   <- !is.na(u2.N.1)
for(i in 1:nrow(p)){
   # generate sample data
   gen.m2 <- rbinom(ncol(p), n1, p[i,])
   gen.u2.A.YoY <- rbinom(ncol(p), U.H.YoY[i,], p[i,]*clip.frac.H.YoY)  # only hatchery fish can generate adipose clipped fish
   gen.u2.N.YoY <- rbinom(ncol(p), U.W.YoY[i,], p[i,]) +
                   rbinom(ncol(p), U.H.YoY[i,], p[i,]*(1-clip.frac.H.YoY)) # wild and hatchery fish generate non-clipped fish
   gen.u2.A.1   <- rbinom(ncol(p), U.H.1  [i,], p[i,]*clip.frac.H.1)   # only hatchery fish can generate adipose clipped fish
   gen.u2.N.1   <- rbinom(ncol(p), U.W.1  [i,], p[i,]) +
                   rbinom(ncol(p), U.H.1  [i,], p[i,]*(1-clip.frac.H.1)) # wild and hatchery fish generate non-clipped fish
   # compute a discrepancy measure
   # Observed vs expected values for recaptures of marked fish
     temp <- sqrt(m2) - sqrt(n1*p[i,])
     d1.m2.o <- sum( temp[select.m2]^2, na.rm=TRUE)
     temp <- sqrt(gen.m2) - sqrt(n1*p[i,])
     d1.m2.s <- sum( temp[select.m2]^2, na.rm=TRUE)

   # Observed vs expected values for captures of unmarked but clipped fish
     temp <- sqrt(u2.A.YoY) - sqrt(U.H.YoY[i,]*p[i,]*clip.frac.H.YoY) # YoY fish. Recall clipped YoY only come after hatch.after
     temp <- temp[time>hatch.after.YoY]
     d1.u2.A.YoY.o <- sum( temp[select.u2.A.YoY]^2, na.rm=TRUE)
     temp <- sqrt(gen.u2.A.YoY) - sqrt(U.H.YoY[i,]*p[i,]*clip.frac.H.YoY)
     temp <- temp[time>hatch.after.YoY]
     d1.u2.A.YoY.s <- sum( temp[select.u2.A.YoY]^2, na.rm=TRUE)

     temp <- sqrt(u2.A.1  ) - sqrt(U.H.1[i,]*p[i,]*clip.frac.H.1)   # age 1 fish
     d1.u2.A.1.o <- sum( temp[select.u2.A.1]^2, na.rm=TRUE)
     temp <- sqrt(gen.u2.A.1) - sqrt(U.H.1[i,]*p[i,]*clip.frac.H.1)
     d1.u2.A.1.s <- sum( temp[select.u2.A.1]^2, na.rm=TRUE)

   # Observed vs expected values for captures of unmarked fish with NO adipose clips (a mixture of wild and hatchery fish)
     temp <- sqrt(u2.N.YoY) - sqrt(U.W.YoY[i,]*p[i,]+U.H.YoY[i]*p[i,]*(1-clip.frac.H.YoY))
     d1.u2.N.YoY.o <- sum( temp[select.u2.N.YoY]^2, na.rm=TRUE)
     temp <- sqrt(gen.u2.N.YoY) - sqrt(U.W.YoY[i,]*p[i,]+U.H.YoY[i]*p[i,]*(1-clip.frac.H.YoY))
     d1.u2.N.YoY.s <- sum( temp[select.u2.N.YoY]^2, na.rm=TRUE)

     temp <- sqrt(u2.N.1) - sqrt(U.W.1[i,]*p[i,]+U.H.1[i]*p[i,]*(1-clip.frac.H.1))
     d1.u2.N.1.o <- sum( temp[select.u2.N.1]^2, na.rm=TRUE)
     temp <- sqrt(gen.u2.N.1) - sqrt(U.W.1[i,]*p[i,]+U.H.1[i]*p[i,]*(1-clip.frac.H.1))
     d1.u2.N.1.s <- sum( temp[select.u2.N.1]^2, na.rm=TRUE)

   # combined discrepancy measures
     d1.YoY.o <- d1.u2.A.YoY.o + d1.u2.N.YoY.o  # observed data total discrepancy  for YoY
     d1.YoY.s <- d1.u2.A.YoY.s + d1.u2.N.YoY.s  # simulated data total discrepancy for YoY
     d1.1.o   <- d1.u2.A.1.o   + d1.u2.N.1.o    # observed data total discrepancy  for Age1
     d1.1.s   <- d1.u2.A.1.s   + d1.u2.N.1.s    # simulated data total discrepancy for Age1
     d1.o     <- d1.m2.o + d1.YoY.o + d1.1.o    # observed data total discrepancy  all data
     d1.s     <- d1.m2.s + d1.1.s   + d1.1.s    # simulated data total discrepancy all data
   # update the array
     discrep <- rbind(discrep, 
              c(d1.m2.o, d1.m2.s,     
                d1.u2.A.YoY.o, d1.u2.A.YoY.s, 
                d1.u2.N.YoY.o, d1.u2.N.YoY.s, 
                d1.u2.A.1.o,   d1.u2.A.1.s, 
                d1.u2.N.1.o,   d1.u2.N.1.s, 
                d1.YoY.o,      d1.YoY.s,
                d1.1.o,        d1.1.s,
                d1.o   , d1.s         
                )) 
}
#browser()
discrep
}

