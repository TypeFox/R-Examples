"pmoms" <-
function(x) {
  n  <- length(x)

  if(n < 4) {
    warning("Four or more (hopefully much more) data values are needed")
    return(NULL)
  }

  MU <- mean(x)
  SD <- sd(x) # returns with the proper n - 1
              # division for bias correction of variance
  my.term  <- exp(lgamma((n-1)/2) - lgamma(n/2)) # this construct--avoid overflow
  SD.prime <- my.term / sqrt(2)
  SD.prime <- SD.prime * sqrt(sum((x - MU)^2))

  hosking.M2 <- sum((x-MU)^2) / n # HW(1997, eq. 2.20)
  M2 <- sqrt(hosking.M2) # theoretical definition of standard deviation

  M3 <- sum((x-MU)^3) / n
  classic.skew  <- M3 / M2^3 # theoretical definition of skew
  skew <- n^2*M3 / ((n-1)*(n-2)*SD^3) # correct for sample size bias

  hosking.M4   <- M4 <- sum((x-MU)^4) / n # HW(1997, eq. 2.20)
  classic.kurt <- M4 / M2^4 # theoretical definition of kurtosis

  # Hosking and Wallis (1997, eq. 2.23) have a slight rearrangement of
  # algebra M2 and M4 are defined in HW(1997, eq. 2.20)
  tmp <- n^2 / ((n-2)*(n-3)) # double checked 5/8/2013
  hosking.k <- (n+1)/(n-1) * hosking.M4 - 3*hosking.M2^2 # reimplemented, 5/8/2013
  hosking.k <- hosking.k * tmp
  hosking.k <- hosking.k/SD^4 + 3 # HW(1997, eq. 2.24), confirmed SD needed 5/8/2013
  #cat(c("DEBUG: Hosking k",hosking.k,"\n"))
  # Hosking and Wallis (1997) definition of kurtosis is repeated correctly in
  # Asquith (2011, ISBN 978-1463508418

  # S.L. Dingman (1992) "Physical Hydrology," Appendix C.
  #dingman.k <- n^3*M4 / ((n-1)*(n-2)*(n-3)*SD^4) # correct for sample size bias
  #cat(c("DEBUG: Dingman k",dingman.k,"\n")) # DINGMAN APPEARS WRONG

  # Additions 5/8/2013 with discussion with Sergio Gomez
  #Rimoldini.M2 <- hosking.M2; Rimoldini.M4 <- M4
  #Rimoldini.k <- n^2*(n+1)/((n-1)*(n-2)*(n-3)) * Rimoldini.M4
  #Rimoldini.k <- Rimoldini.k - 3*n^2/((n-2)*(n-3)) * Rimoldini.M2^2
  #Rimoldini.k <- Rimoldini.k/SD^4 # thanks Sergio Gomez
  #Rimoldini.k <- Rimoldini.k + 3 # swap away from excess
  #cat(c("DEBUG: Rimoldini k",Rimoldini.k,"\n")) # this appears to be in error
  # when one looks at all the DEBUG outputs.

  # Addition 5/8/2013 with discussion with Sergio Gomez
  #Wiki.k <- (n+1)*n/((n-1)*(n-2)*(n-3))*(n*M4)/SD^4
  #Wiki.k <- Wiki.k - 3 * (n-1)^2/((n-2)*(n-3))
  # Note, http://en.wikipedia.org/wiki/Kurtosis#Estimators_of_population_kurtosis
  # says "k2 is is the unbiased estimate of the second cumulant (identical to
  # the unbiased estimate of the sample variance)" However, it seems then that
  # re-squaring k2 as k^2_2 as wiki notation show is incorrect. The words
  # "sample variance"
  #Wiki.k <- Wiki.k + 3 # swap away from excess
  #cat(c("DEBUG: Wiki k",Wiki.k,"\n"))

  # Addition 5/8/2013 with discussion with Sergio Gomez
  # http://www.ats.ucla.edu/stat/mult_pkg/faq/general/kurtosis.htm
  #idre.k <- n*(n+1)/((n-1)*(n-2)*(n-3))*(n*M4/SD^4)
  #cat(c("DEBUG: idre k",idre.k,"\n")) # this appears to be in error
  # when one looks at all the DEBUG outputs.  # IDRE APPEARS WRONG

  # Addition 5/8/2013 with discussion with Sergio Gomez
  # http://www.tc3.edu/instruct/sbrown/stat/shape.htm#KurtosisCompute
  # D. N. Joanes and C. A. Gill., 1998, Comparing Measures of Sample Skewness
  # and Kurtosis. The Statistician 47(1):183-189.
  #Brown.a4 <- hosking.M4 / (hosking.M2)^2
  #Brown.g2 <- Brown.a4 - 3 # excess kurtosis
  #Brown.k  <- (n-1)/((n-2)*(n-3)) * ((n+1)*Brown.g2 + 6)
  #Brown.k  <- Brown.k + 3 # swap away from excess
  #cat(c("DEBUG: Brown k",Brown.k,"\n"))

  # Addition 5/8/2013 with discussion with Sergio Gomez
  # Let us look at the Joanes and Gill (1998) definition
  # tmp <- n^2/((n-1)*(n-2)*(n-3))
  #JG.k <- tmp*((n+1)*hosking.M4 - 3*(n-1)*hosking.M2^2)
  #JG.k <- JG.k/(n*hosking.M2/(n-1))^2
  #JG.k <- JG.k + 3 # swap away from excess
  #cat(c("DEBUG: JG k",JG.k,"\n"))

  z <- list(moments      = c(MU, SD,     M3,       M4),
            ratios       = c(NA, SD/MU,skew,hosking.k),
            sd           = SD,
            umvu.sd      = SD.prime,
            skew         = skew,
            kurt         = hosking.k,
            excesskurt   = hosking.k - 3,
            classic.sd   = M2,
            classic.skew = classic.skew,
            classic.kurt = classic.kurt,
            classic.excesskurt = classic.kurt - 3,
            message = c("The 'classic' values should be slightly different (the bias).",
                        "sd is not truly unbiased; although, surprizingly sd^2 ",
                        "(variance) is."),
            source = "pmoms")
  return(z)
}
