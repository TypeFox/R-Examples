PredictivePosterior.TSPDE.WHSteel <- function (time, n1, m2, u2.W.YoY, u2.W.1, u2.H.1, p, U.W.YoY, U.W.1, U.H.1, hatch.after) {
#  Generate Predictive Posterior Plot (Bayesian p-value) given the data
#  for a TimeStratified Petersen with Diagonal Elements and error
#    n1, m2, u2.*   = vectors of input data
#    p, U.*         = matrix of values (rows=number of posterior samples, columns=strata)
#               These are returned from the call to OpenBugs/ WinBugs
#
#cat("Call to PredictivePosterior\n")
#browser()
discrep <- matrix(0, nrow=0, ncol=10)
select.m2       <- !is.na(m2)
select.u2.W.YoY <- !is.na(u2.W.YoY)
select.u2.W.1   <- !is.na(u2.W.1)
select.u2.H.1   <- !is.na(u2.H.1)

for(i in 1:nrow(p)){
   # generate sample data
   gen.m2       <- rbinom(ncol(p), n1,          p[i,])
   gen.u2.W.YoY <- rbinom(ncol(p), U.W.YoY[i,], p[i,])  
   gen.u2.W.1   <- rbinom(ncol(p), U.W.1  [i,], p[i,])  
   gen.u2.H.1   <- rbinom(ncol(p), U.H.1  [i,], p[i,])  
 
   # compute a discrepancy measure
   # Observed vs expected values for recaptures of marked fish
     temp <- sqrt(m2) - sqrt(n1*p[i,])
     d1.m2.o <- sum( temp[select.m2]^2, na.rm=TRUE)
     temp <- sqrt(gen.m2) - sqrt(n1*p[i,])
     d1.m2.s <- sum( temp[select.m2]^2, na.rm=TRUE)

   # Observed vs expected values for observed data
     temp          <- sqrt(u2.W.YoY) - sqrt(U.W.YoY[i,]*p[i,])
     d1.u2.W.YoY.o <- sum( temp[select.u2.W.YoY]^2, na.rm=TRUE)
     temp          <- sqrt(u2.W.1)    - sqrt(U.W.1[i,]*p[i,])
     d1.u2.W.1.o   <- sum( temp[select.u2.W.1]^2, na.rm=TRUE)
     temp          <- sqrt(u2.H.1)    - sqrt(U.H.1[i,]*p[i,])
     temp          <- temp[time>hatch.after]  # don't include terms before the hatchery fish are released
     d1.u2.H.1.o   <- sum( temp[select.u2.H.1]^2, na.rm=TRUE)

   # Observed vs expected values for simulated data
     temp          <- sqrt(gen.u2.W.YoY) - sqrt(U.W.YoY[i,]*p[i,])
     d1.u2.W.YoY.s <- sum( temp[select.u2.W.YoY]^2, na.rm=TRUE)
     temp          <- sqrt(gen.u2.W.1)    - sqrt(U.W.1[i,]*p[i,])
     d1.u2.W.1.s   <- sum( temp[select.u2.W.1]^2, na.rm=TRUE)
     temp          <- sqrt(gen.u2.H.1)    - sqrt(U.H.1[i,]*p[i,])
     temp          <- temp[time>hatch.after]  # don't include terms before the hatchery fish are released
     d1.u2.H.1.s   <- sum( temp[select.u2.H.1]^2, na.rm=TRUE)


   # combined discrepancy measures
     d1.o <- d1.m2.o + d1.u2.W.YoY.o + d1.u2.W.1.o + d1.u2.H.1.o  # observed data total discrepancy
     d1.s <- d1.m2.s + d1.u2.W.YoY.s + d1.u2.W.1.s + d1.u2.H.1.s # simulated data total discrepancy
   # update the array
     discrep <- rbind(discrep, 
              c(d1.m2.o,       d1.m2.s,      
                d1.u2.W.YoY.o, d1.u2.W.YoY.s, 
                d1.u2.W.1.o,   d1.u2.W.1.s,    
                d1.u2.H.1.o,   d1.u2.H.1.s,    
                d1.o   , d1.s
                )) 
}
#browser()
discrep
}
