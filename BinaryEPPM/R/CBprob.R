CBprob <-
function(parameter,nt) {
   p   <- parameter[1] 
   q   <- 1 - p 
   rho <- parameter[2] 
#  m should equal nt + 1 
   m  <- nt + 1
   probability <- rep(0,m)
# limits on rho given in Kupper & Haseman Biometrics (1978) 
   if (nt>1) {
      if (p<0.5) { llrho <- - 2*p/(nt*(nt-1)*(1-p)) 
                 } else { llrho <- - 2*(1-p)/(nt*(nt-1)*p) }
      wkv <- c(0:nt)
      wkv <- (wkv - (nt - 1)*p - 0.5)**2 
      gamma0 <- min(wkv) 
      ulrho <- 2*p*(1-p)/((nt-1)*p*(1-p) + 0.25 - gamma0) 
      if (is.finite(ulrho)==FALSE) { ulrho <- 0 }
              } else {
         llrho <- -1 
         ulrho <- +1 } 

  if ((rho>=llrho) & (rho<=ulrho)) {

     if (round(rho,digits=14)==0) {
        probability <- dbinom(c(0:nt),nt,p,log=FALSE)
                                       } else {
        np  <- nt*p 
        np2 <- np*p 
        twopm1  <- 2*p - 1 
        twop2q2 <- 2*(p*q)**2 
        for ( i in 1:m ) {
           im1 <- i - 1 
           wks <- 1 + rho*( (im1-np)**2 + im1*twopm1 - np2 )/twop2q2 
           if ((wks>=0) & (p>0) & (q>0)) { 
              probability[i] <- wks*exp(im1*log(p) + (nt-im1)*log(q))
                        }} # end of for loop
#  Calculating factorials 
        vfac <- choose(c(rep(nt,m)),c(0:nt))
        probability <- probability*vfac 
#  Checking for valid probability distribution
#  Total probability can be slightly greater than 1 therefore rounding
       total.prob <- round(sum(probability),digits=10)
       valid.prob <- sum(((probability>=0) & (probability<=1)))
       if (((total.prob>1) | (valid.prob<m) | (is.na(total.prob)==TRUE) | 
          (is.na(valid.prob)==TRUE))) { probability <- rep(NA,m) } 
                                                } # end of if (rho==0)
          } else { probability <- rep(NA,m) } # end of llrho, ulrho check
   output <- list(probability=probability,parameter=parameter,
                  limit=matrix(c(llrho,ulrho),ncol=2))
   return(output)
                                    }
