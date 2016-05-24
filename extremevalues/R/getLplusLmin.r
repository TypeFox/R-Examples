# compute limits.
# mvdl 10.12.2009

getLplusLmin <- function(sigmaE, alpha)
{
   Lplus <- Inf
   Lmin <- -Inf
   if ( !is.na(alpha[2]) )
      Lplus <- sqrt(2)*sigmaE*invErf(1-2*alpha[2])
   if ( !is.na(alpha[1]) )
      Lmin  <- sqrt(2)*sigmaE*invErf(2*alpha[1]-1)
   return(list(Lplus=Lplus,Lmin=Lmin));
}

