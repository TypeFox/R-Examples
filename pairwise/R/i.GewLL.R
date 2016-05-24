### ------------------------- GewLL ------------------------------
GewLL <- function(np,sb,fw,m,  k=length(m),mscs=sum(m-1)){ # Verbess., f2 und w(eighted)LL
  # argumente:
  # np vector der länge mscs+1 (anzahl der rohwertgruppen) Anzahl der Personen mit Scoresumme r
  # sb matrix mit  bedingten Schwierigkeitsparametern (sb) aus tau: t(apply(U$tau,1,function(x){cumsum(x)}))
  # fw vector der länge mscs+1 (anzahl der rohwertgruppen) mit fähigkeitswerten fw
  # ergaenzt um : 
  # mscs integer maximale rohwertsumme 
  # m vector der länge k mit anzahl der kategorien (ehemals k)
  # k integer anzahl der items (ehemals m)
  verb <- f2 <- rep(0.,mscs+1); 
  wLL = 0.
  for(r0 in 0:mscs){ r = r0+1; L0=L1=f3=0.
                     for(i in 1:k){ N0=1.; N1=N2=N3=0.
                                    for(x in 1:(m[i]-1)){ e = exp(x*fw[r]-sb[i,x])
                                                          N0 = N0+e; N1 = N1+x*e
                                                          N2 = N2+x^2*e; N3 = N3+x^3*e } # Ende x
                                    f2[r] = f2[r]-N2/N0+(N1/N0)^2
                                    L0 = L0+log(N0); L1 = L1+N1/N0
                                    f3 = f3+N3/N0-3*N2*N1/N0^2+2*(N1/N0)^3 } # Ende i
                     wLL = wLL+np[r]*( r0*fw[r]-L0+log(-f2[r])/2 )
                     f1 = r0-L1-f3/f2[r]/2
                     verb[r] = f1/f2[r] } # Ende r0
  erg<-list(wLL=wLL, verb=verb, f2=f2) 
  return(erg)
}