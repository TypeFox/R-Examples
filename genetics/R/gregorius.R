# $Id: gregorius.R 114 2003-05-22 17:25:23Z warnesgr $
#
# Code contributed by David Duffy <davidD@qumr.edu.au>.
#
# Gregorius, H.-R. 1980. The probability of losing an allele when
# diploid genotypes are sampled.  Biometrics 36, 643-652.
#
# Formula from "Corollary 2" and "Corollary 3" of that paper
#
# N is the number of genotypes sampled,
# freq=frequency of least common allele to be detected by the study,
# missprob=the probability of missing at least one allele
#
# tol=smallest term in series to be accumulated
#
gregorius <- function(freq, N, missprob, tol=1.0e-10, maxN=1e4, maxiter=100,
                      showiter=FALSE)
{
   

  find.alpha <- function(N, freq, tol) #, showiter=FALSE)
    {
      n<- floor(1/freq)
      i<-1
      sgn<- -1
      term<-1.0
      res<-0.0
      while(abs(term)>tol && i<n) {
        sgn<- (-1) ^ (i+1)
        term<- exp( lchoose(n-1,i) +
                   log(exp(N*log(1-i*freq))+
                       exp(i+N*log(freq)+(N-1)*log(n-i))))
        res<-res+sgn*term
        i<-i+1

#        if(showiter)
#          {
#            cat("i=",i,"\n")
#            cat("sgn=",sgn,"\n")
#            cat("term=",term,"\n")
#            cat("res=",res,"\n")
#          }
      }

      max(min(res,1),0)
    }


  retval <- list()
  retval$call <- match.call()
    
  
  if(!missing(N) && missing(missprob) )
    {
      retval$method <- "Compute missprob given N and freq"
      retval$freq <- freq
      retval$N <- N
      retval$missprob <- find.alpha(N=N,freq=freq,tol=tol) 
    }
  else if(missing(N) && !missing(missprob) )
    {
      retval$method <- "Determine minimal N given missprob and freq"
      retval$freq <- freq
      val <- binsearch( function(N) find.alpha(N=N, freq=freq, tol=tol),
                       range=c(1, maxN), target=missprob, showiter=showiter,
                       maxiter=maxiter )
      if(length(val$where)==2)
        {
          retval$N <- val$where[2]
          retval$missprob <- val$value[2]
        }
      else
        {
          retval$N <- val$where[1]
          retval$missprob <- val$value[1] 
       }
    }
  else
    stop("Exactly two of N, freq, and missprob must be specified")

  return(retval)
}

