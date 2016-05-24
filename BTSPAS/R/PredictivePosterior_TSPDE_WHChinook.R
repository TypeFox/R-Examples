PredictivePosterior.TSPDE.WHCH <- function (time, n1, m2, u2.A, u2.N, clip.frac.H, p, U.W, U.H, hatch.after) {
#  Generate Predictive Posterior Plot (Bayesian p-value) given the data
#  for a TimeStratified Petersen with Diagonal Elements and error
#    n1, m2, u2.A, u2.N  = vectors of input data
#    p, U.W. U.H         = matrix of values (rows=number of posterior samples, columns=strata)
#               These are returned from the call to OpenBugs/ WinBugs
#
#cat("Call to PredictivePosterior\n")
#browser()
discrep <- matrix(0, nrow=0, ncol=8)
select.m2 <- !is.na(m2)
select.u2.A <- !is.na(u2.A)
select.u2.N <- !is.na(u2.N)
for(i in 1:nrow(p)){
   # generate sample data
   gen.m2 <- rbinom(ncol(p), n1, p[i,])
   gen.u2.A <- rbinom(ncol(p), U.H[i,], p[i,]*clip.frac.H)  # only hatchery fish can generate adipose clipped fish
   gen.u2.N <- rbinom(ncol(p), U.W[i,], p[i,]) +
               rbinom(ncol(p), U.H[i,], p[i,]*(1-clip.frac.H)) # wild and hatchery fish generate non-clipped fish
   # compute a discrepancy measure
   # Observed vs expected values for recaptures of marked fish
     temp <- sqrt(m2) - sqrt(n1*p[i,])
     d1.m2.o <- sum( temp[select.m2]^2, na.rm=TRUE)
     temp <- sqrt(gen.m2) - sqrt(n1*p[i,])
     d1.m2.s <- sum( temp[select.m2]^2, na.rm=TRUE)

   # Observed vs expected values for captures of unmarked but clipped fish
     temp <- sqrt(u2.A) - sqrt(U.H[i,]*p[i,]*clip.frac.H)
     temp <- temp[time>hatch.after]
     d1.u2.A.o <- sum( temp[select.u2.A]^2, na.rm=TRUE)
     temp <- sqrt(gen.u2.A) - sqrt(U.H[i,]*p[i,]*clip.frac.H)
     temp <- temp[time>hatch.after]
     d1.u2.A.s <- sum( temp[select.u2.A]^2, na.rm=TRUE)

   # Observed vs expected values for captures of unmarked fish with NO adipose clips (a mixture of wild and hatchery fish)
     temp <- sqrt(u2.N) - sqrt(U.W[i,]*p[i,]+U.H[i]*p[i,]*(1-clip.frac.H))
     d1.u2.N.o <- sum( temp[select.u2.N]^2, na.rm=TRUE)
     temp <- sqrt(gen.u2.N) - sqrt(U.W[i,]*p[i,]+U.H[i]*p[i,]*(1-clip.frac.H))
     d1.u2.N.s <- sum( temp[select.u2.N]^2, na.rm=TRUE)

   # combined discrepancy measures
     d1.o <- d1.m2.o + d1.u2.A.o + d1.u2.N.o  # observed data total discrepancy
     d1.s <- d1.m2.s + d1.u2.A.s + d1.u2.N.s  # simulated data total discrepancy
   # update the array
     discrep <- rbind(discrep, 
              c(d1.m2.o, d1.m2.s,     
                d1.u2.A.o, d1.u2.A.s, 
                d1.u2.N.o, d1.u2.N.s, 
                d1.o   , d1.s         
                )) 
}
#browser()
discrep
}

