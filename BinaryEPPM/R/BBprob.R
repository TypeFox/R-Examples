BBprob <-
function(parameter,nt) {
   p     <- parameter[1]
   theta <- parameter[2]
   m <- nt + 1
   probability <- rep(0,m)
# limits as in Smith & Ridout A Remark on Algorithm AS 189 
# max(n) <- nt here   
   if (nt>1) { lltheta <- -p/(nt-1)
      if (p>0.5) { lltheta <- -lltheta - 1/(nt-1) }
             } else { lltheta <- -1 } # end if nt

   if (theta>=lltheta) {

       if (round(theta,digits=14)==0) {
           probability <- dbinom(c(0:nt),nt,p,log=FALSE)
                                       } else {
# Calculating probabilities using Smith AS189 formula
# product of (1 + r*theta) r=0 to r=n-1; constant
# as a single sample of grouped binary data
      ntm1  <- nt - 1
      if (ntm1>0) { denom <- sum(log(c(rep(1,ntm1))+c(1:ntm1)*c(rep(theta,ntm1))))
          } else { denom <- 0 }
      for ( i in 0:nt ) { 
# product of (p + r*theta) r=0 to r=x-1
         xm1 <- i - 1
         if (xm1>-1) { numer_one <- sum(log(c(rep(p,i))+c(0:xm1)*c(rep(theta,i))))
                     } else { numer_one <- 0 }
# product of (1 - p + r*theta) r=0 to r=n-x-1
         nmx   <- nt - i
         nmxm1 <- nmx - 1
         if ( nmxm1>-1) { numer_two <- sum(log(c(rep((1-p),nmx))+c(0:nmxm1)*c(rep(theta,nmx))))
                        } else { numer_two <- 0 } 
         probability[i+1] <- exp(numer_one + numer_two - denom)
                    } # end of for loop
#  Calculating factorials  
   vfac <- choose(c(rep(nt,m)),c(0:nt))
   probability <- probability*vfac 
#  Checking for valid probability distribution
#  Total probability can be slightly greater than 1 therefore rounding
   total.prob <- round(sum(probability),digits=10)
   valid.prob <- sum(((probability>=0) & (probability<=1)))
   if (((total.prob>1) | (valid.prob<m) | (is.na(total.prob)==TRUE) | 
       (is.na(valid.prob)==TRUE))) { probability <- rep(NA,m) } 
                                             } # end of if (theta==0)
       } else { probability <- rep(NA,m) } # end of lltheta check
   output <- list(probability=probability,parameter=parameter,
                  limit=c(lltheta))
   return(output)
                                    }
