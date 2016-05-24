distance <- function(P,Q)
{
  if (length(P)!=length(Q))
  {
    print("erreur")
    return(0)
  } else {
    tailcop=length(P)
    KL=sum(P*log(P/Q),na.rm=TRUE)
    L2=sum((P-Q)^2)*tailcop
    L1=sum(abs(P-Q))
    APE=sum(abs(1-P/Q))/tailcop
    return(cbind(KL,L2,L1,APE))
  }
}
