"tau34sq.normtest" <-
function(x, alpha=0.05, pvalue.only=FALSE, getlist=TRUE,
         useHoskingZt4=TRUE, verbose=FALSE, digits=4) {
   n   <- length(x)
   lmr <- lmoms(x); T3 <- lmr$ratios[3]; T4 <- lmr$ratios[4]
   Zt3 <- T3 * (1/sqrt(exp(log(0.1866) - log(n)) + exp(log(0.8) - 2*log(n))))
   Zt4 <- exp(log(0.0883) - log(n))

   if(useHoskingZt4) Zt4 <- Zt4 + exp(log(0.68) - 2*log(n)) + exp(log(4.9) - 3*log(n))
   Zt4 <- (T4 - 0.1226) * (1/sqrt(Zt4))

   t34sq <- Zt3^2 + Zt4^2

   chiSQ     <- pchisq(t34sq, df=2)
   pvalue    <- 1 - chiSQ

   t34sqrd  <- round(t34sq,  digits=digits)
   chiSQrd  <- round(chiSQ,  digits=digits)*100
   pvaluerd <- round(pvalue, digits=digits)
   Zt3rd    <- round(Zt3,    digits=digits)
   Zt4rd    <- round(Zt4,    digits=digits)
   issig    <- ifelse(pvalue <= alpha, TRUE, FALSE)
   txt      <- ifelse(issig, "<= alpha: reject Ho, conclude 'non-normal data'.",
                             "> alpha: accept Ho, conclude 'normal data'.")
   if(verbose) {
      message("\n        *** Harri-Coble Tau34-squared Test for Normality ***\n")
        message("  --A normality test using sample L-skew (T3) and L-kurtosis (T4)--\n")
      if(! useHoskingZt4) {
         message("Warning: Hosking's personal approximation for Z-score(L-kurtosis) is *NOT* being used.\n")
      }
      message("  Z-score(T3) = ", Zt3rd, " and Z-score(T4) = ", Zt4rd)
      message("    (Mapping of T3 and T4 to standard normal variates)\n")
      message("  Tau34-squared test statistic, T34sq = Z(T3)^2 + Z(T4)^2")
      message("    (A squared Euclidean distance that is Chi-Squared(df=2) distributed.)")
      message("  T34sq = ", t34sqrd ," |--> ChiSquared(T34sq, df=2) = pchisq(T34sq, df=2)\n")
      message("  ChiSquared(", t34sqrd, ", 2 degrees freedom) = ", chiSQrd, " percentile\n", appendLF = FALSE)
      message("    p-value = ",pvaluerd," [complement of the pchisq(T34sq, df=2, lower.tail=TRUE)]")
      message("    The p-value is ",txt,"\n")
      message("  Reference:")
      message("   Harri, A. and Coble, K.H., 2011, Normality testing---\n",
              "     Two new tests using L-moments: J. Appl. Stat. 38(7), 1369-1379.\n")
   }
   if(pvalue.only) return(pvalue)
   if(getlist) {
      z <- list(SampleTau3=T3, SampleTau4=T4,
                Ztau3=Zt3, Ztau4=Zt4, Tau34sq=t34sq,
                ChiSq.2df=chiSQ, pvalue=pvalue, isSig=issig,
                source="tau34sq.normtest")
      return(as.data.frame(z, row.names="Harri-Coble test"))
   }
   return(ifelse(issig, TRUE, FALSE))
}

