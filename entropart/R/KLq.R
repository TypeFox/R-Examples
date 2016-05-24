KLq <-
function(Ps, Pexp, q = 1, CheckArguments = TRUE) 
{
  if (CheckArguments)
    CheckentropartArguments()
  
  dataKLq <- Ps^q*lnq(Ps/Pexp, q)
  # Ignore if Pexp==0 (then P==0 too), and set 0log0=0
  dataKLq[Pexp == 0 | Ps == 0] <- 0
  return (sum(dataKLq)) # equals sum(P* (P^(q-1)-Pexp^(q-1))/(q-1), na.rm = TRUE)
}
